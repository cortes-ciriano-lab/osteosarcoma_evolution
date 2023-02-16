#! /bin/bash
#Jose Espejo Valle-Inclan 2022

#Required: Argparser
# ##https://github.com/nhoffman/argparse-bash
# wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
# chmod +x argparse.bash

ARGPARSE=/nfs/research/icortes/SOFTWARE/argparse_bash/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-n','--normal-bam', required=True, help='Path to normal/germline BAM')
parser.add_argument('-o', '--output-dir', required=True, help='Output folder')
parser.add_argument('-s', '--sample-id', help='Sample ID to append to jobs')
parser.add_argument('-r','--reference', default="/nfs/research/icortes/DATA/hg38_noALT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz", help='Reference genome (must not have HLA ALT contigs) [/nfs/research/icortes/DATA/hg38_noALT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz]')
parser.add_argument('-rv','--ref-version', default="V38", help='Reference genome version (V37, V38, HG19)[V38]')
parser.add_argument('-l','--lilac', default="/nfs/research/icortes/SOFTWARE/hmftools/lilac_v1.2.jar", help='Path to Lilac jar [/nfs/research/icortes/SOFTWARE/hmftools/lilac_v1.2.jar]')
parser.add_argument('-lr','--lilac-resource', default="/nfs/research/icortes/SOFTWARE/hmftools/HMFTools-Resources/Lilac/", help='Path to Lilac resource dir [/nfs/research/icortes/SOFTWARE/hmftools/HMFTools-Resources/Lilac/]')
parser.add_argument('-ha','--hla-alts', default="/nfs/research/icortes/SOFTWARE/hmftools/HMFTools-Resources/Lilac/38/hla_alt.bed", help='Path to HLA alts BED file [/nfs/research/icortes/SOFTWARE/hmftools/HMFTools-Resources/Lilac/38/hla_alt.bed]')
parser.add_argument('-hl','--hla-loci', default="/nfs/research/icortes/SOFTWARE/hmftools/HMFTools-Resources/Lilac/38/hla.38.bed", help='Path to HLA loci BED file [/nfs/research/icortes/SOFTWARE/hmftools/HMFTools-Resources/Lilac/38/hla.38.bed]')
EOF

#Check required parameters

if [ -z ${NORMAL_BAM} ]; then
  echo "Missing -n|--normal-bam parameter"
  exit
fi
if [ -z ${OUTPUT_DIR} ]; then
  echo "Missing -o|--output-dir parameter"
  exit
fi


#1. PREPARE
prepare_function() {

  #Check files and directories
  echo "Checking files and directories"

  if [ ! -f ${NORMAL_BAM} ]; then
    echo "ERROR: ${NORMAL_BAM} does not exist"
    exit
  elif [ ! -f ${NORMAL_BAM}.bai ] && [ ! -f ${NORMAL_BAM/.bam/.bai} ]; then
    echo "ERROR: ${NORMAL_BAM} is not indexed"
    exit
  else
    NORMAL_BAM=$(realpath ${NORMAL_BAM})
  fi
  if [ ! -f ${LILAC} ]; then echo "ERROR: ${LILAC} does not exist"; exit; fi
  if [ ! -d ${LILAC_RESOURCE} ]; then echo "ERROR: ${LILAC_RESOURCE} does not exist"; exit; fi
  if [ ! -f ${REFERENCE} ]; then echo "ERROR: ${REFERENCE} does not exist"; exit; fi
  if [ ! -f ${HLA_ALTS} ]; then echo "ERROR: ${HLA_ALTS} does not exist"; exit; fi
  if [ ! -f ${HLA_LOCI} ]; then echo "ERROR: ${HLA_LOCI} does not exist"; exit; fi

  for tool in java samtools bwa
	do
		if ! command -v ${tool} &> /dev/null; then
	    	echo "${tool} not found, did you load the Conda environment?"
	    	exit
	    else echo "${tool}: OK"
		fi
	done

  OUTPUT_DIR=$(realpath ${OUTPUT_DIR})
  TMPDIR=${OUTPUT_DIR}/tmp
  LOGDIR=${OUTPUT_DIR}/log
  mkdir -p ${OUTPUT_DIR} ${LOGDIR} ${TMPDIR}
  SAMPLE=$(samtools view -H ${NORMAL_BAM} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq | head -1)
  echo "Sample for this run is ${SAMPLE}"

  if [ -z ${SAMPLE_ID} ]; then
    RAND="${SAMPLE}_$(cat /dev/urandom | tr -cd 'a-zA-Z0-9' | head -c 5)"
  else
    RAND="${SAMPLE_ID}_$(cat /dev/urandom | tr -cd 'a-zA-Z0-9' | head -c 5)"
  fi
  echo "Identifier string for this run is ${RAND}"

}

remap_normal_function() {
  REMAP_HLA_N_DIR=${TMPDIR}/normal
  REMAP_HLA_N_TMP_DIR=${REMAP_HLA_N_DIR}/tmp
  REMAP_HLA_N_OUT=${REMAP_HLA_N_DIR}/${SAMPLE}.normal.hla.sorted.bam
  REMAP_HLA_N_DONE=${REMAP_HLA_N_OUT}.done
  REMAP_HLA_N_JOB="REMAP_HLA_N_${RAND}"
  REMAP_HLA_N_FLAG=true
  if [ -e ${REMAP_HLA_N_DONE} ]; then
		echo "REMAP_HLA_N output exists at ${REMAP_HLA_N_OUT}"
		REMAP_HLA_N_FLAG=false
	else
		echo "Creating REMAP_HLA_N output directory in ${REMAP_HLA_N_DIR}"
    mkdir -p ${REMAP_HLA_N_DIR} ${REMAP_HLA_N_TMP_DIR}
    REMAP_HLA_N_CMD="samtools view -huL ${HLA_ALTS} ${NORMAL_BAM} | \
samtools collate -Oun128 -@ 2 - \
${REMAP_HLA_N_TMP_DIR}/${SAMPLE}_normal.collate.tmp | \
samtools fastq -OT RG,BG -@ 2 - | \
bwa mem -p -t 2 -C ${REFERENCE} - | \
samtools sort -@ 2 -m 500m \
-T ${REMAP_HLA_N_TMP_DIR} \
-o ${REMAP_HLA_N_OUT} - \
&& samtools index ${REMAP_HLA_N_OUT} \
&& touch ${REMAP_HLA_N_DONE}
"
    echo "REMAP_HLA_N command is:"
    echo "${REMAP_HLA_N_CMD}" | tee "${LOGDIR}/${REMAP_HLA_N_JOB}.cmd"
    REMAP_HLA_N_BSUB=(bsub
			-M 10G
      -n 8
			-o "${LOGDIR}/${REMAP_HLA_N_JOB}.o"
			-e "${LOGDIR}/${REMAP_HLA_N_JOB}.e"
			-J "${REMAP_HLA_N_JOB}"
			"${REMAP_HLA_N_CMD}")
    echo "${REMAP_HLA_N_BSUB[@]}" > "${LOGDIR}/${REMAP_HLA_N_JOB}.bsub"
    "${REMAP_HLA_N_BSUB[@]}"
  fi
}

merge_normal_function() {
  MERGE_HLA_N_OUT=${REMAP_HLA_N_OUT/.bam/.merged.bam}
  MERGE_HLA_N_DONE=${MERGE_HLA_N_OUT}.done
  MERGE_HLA_N_JOB="MERGE_HLA_N_${RAND}"
  MERGE_HLA_N_FLAG=true
  if [ -e ${MERGE_HLA_N_DONE} ]; then
		echo "MERGE_HLA_N output exists at ${MERGE_HLA_N_OUT}"
		MERGE_HLA_N_FLAG=false
	else
    MERGE_HLA_N_CMD="samtools merge -@ 8 -h ${NORMAL_BAM} \
-L ${HLA_LOCI} ${MERGE_HLA_N_OUT} ${NORMAL_BAM} ${REMAP_HLA_N_OUT} \
&& samtools index ${MERGE_HLA_N_OUT} \
&& touch ${MERGE_HLA_N_DONE}
"
    echo "MERGE_HLA_N command is:"
    echo "${MERGE_HLA_N_CMD}" | tee "${LOGDIR}/${MERGE_HLA_N_JOB}.cmd"
    MERGE_HLA_N_BSUB=(bsub
			-M 4G
      -n 8
			-o "${LOGDIR}/${MERGE_HLA_N_JOB}.o"
			-e "${LOGDIR}/${MERGE_HLA_N_JOB}.e"
			-J "${MERGE_HLA_N_JOB}"
    )
    if [ "${REMAP_HLA_N_FLAG}" = "true" ]; then
			MERGE_HLA_N_BSUB+=("-w done(${REMAP_HLA_N_JOB})")
		fi
		MERGE_HLA_N_BSUB+=("${MERGE_HLA_N_CMD}")
    echo "${MERGE_HLA_N_BSUB[@]}" > "${LOGDIR}/${MERGE_HLA_N_JOB}.bsub"
    "${MERGE_HLA_N_BSUB[@]}"
  fi
}

lilac_function(){
  LILAC_DONE=${OUTPUT_DIR}/LILAC.done
  LILAC_JOB="LILAC_${RAND}"
  LILAC_FLAG=true
  if [ -e ${LILAC_DONE} ]; then
		echo "LILAC output exists at ${OUTPUT_DIR}"
		LILAC_FLAG=false
	else
    LILAC_CMD="java -jar ${LILAC} \
-sample ${SAMPLE} \
-ref_genome ${REFERENCE} \
-ref_genome_version ${REF_VERSION} \
-resource_dir ${LILAC_RESOURCE} \
-reference_bam ${MERGE_HLA_N_OUT} \
-output_dir ${OUTPUT_DIR} && \
touch ${LILAC_DONE}"
    echo "LILAC command is:"
    echo "${LILAC_CMD}" | tee "${LOGDIR}/${LILAC_JOB}.cmd"
    LILAC_BSUB=(bsub
      -M 10G
      -o "${LOGDIR}/${LILAC_JOB}.o"
      -e "${LOGDIR}/${LILAC_JOB}.e"
      -J "${LILAC_JOB}"
    )
    LILAC_DEPENDENCY=()
    if [ "${MERGE_HLA_T_FLAG}" = "true" ]; then
      LILAC_DEPENDENCY+=("done(${MERGE_HLA_T_JOB})")
    fi
    if [ "${MERGE_HLA_N_FLAG}" = "true" ]; then
      LILAC_DEPENDENCY+=("done(${MERGE_HLA_N_JOB})")
    fi
    if [ ! ${#LILAC_DEPENDENCY[@]} -eq 0 ]; then
  			LILAC_BSUB+=("-w $(echo ${LILAC_DEPENDENCY[@]} | sed 's: : \&\& :g')")
  	fi
		LILAC_BSUB+=("${LILAC_CMD}")
    echo "${LILAC_BSUB[@]}" > "${LOGDIR}/${LILAC_JOB}.bsub"
		"${LILAC_BSUB[@]}"
  fi
}

echo "###Running preparations"
prepare_function
echo "###Preparations finished"
echo ""
echo "###Preparing REMAP_HLA_N"
remap_normal_function
echo "###"
echo ""
echo "###Preparing MERGE_HLA_N"
merge_normal_function
echo "###"
echo ""
echo "###Preparing LILAC"
lilac_function
echo "###"
echo ""
