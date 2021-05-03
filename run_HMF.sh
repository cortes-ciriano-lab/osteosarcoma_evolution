#!/bin/bash

#Jose Espejo Valle-Inclan 2021


#TO DO#
#Add kraken2db
#TEST & TEST & TEST
#######


installDir=$(realpath $(dirname ${BASH_SOURCE[0]}))

#Defaults
iniFile="${installDir}/hmf.ini"
reference="/hps/research1/icortes/DATA/hg38/Homo_sapiens_assembly38.fasta"
condaLoadedFlag="false"
condaName=hmf
condaSrc=/hps/research1/icortes/jespejo/conda-envs/miniconda3/



#Took the arg parsing from gridss.sh (https://github.com/PapenfussLab/gridss)
USAGE_MESSAGE="Usage: runHMF.sh [options] -t <tumor.bam> -n <normal.bam> -o <outputDir> -i <iniFile>

Required parameters:
	-t/--tumorBam: path to tumor BAM file.
	-n/--normalBam: path to normal BAM file.
	-o/--outputDir: path to output directory (will be created).
	

Optional parameters:
	-i/--iniFile: path to ini file [${iniFile}]
	-h/--help: show this usage help.
	-r/--reference: reference genome to use [${reference}].
	--condaLoaded: flag; your env is already pre-loaded, don't load another one [${condaLoadedFlag}].
	--condaName: Name of your conda env [${condaName}]
	--condaSrc Path to conda install directory, needed to properly source [${condaSrc}]:
	"
usage () {
	echo "${USAGE_MESSAGE}"
	exit 1
}

OPTIONS="hr:t:n:o:i:"
LONGOPTS="help,reference:,tumorBam:,normalBam:,outputDir:,iniFile:,condaLoaded,condaName:,condaSrc:"
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	usage
fi
eval set -- "$PARSED"
while true; do
    case "$1" in
    	-h|--help)
		    usage
    		shift # past argument
   			;;
        -r|--reference)
            reference="$2"
            shift 2
            ;;
        -t|--tumorBam)
            tumorBam="$2"
            shift 2
            ;;
        -n|--normalBam)
            normalBam="$2"
            shift 2
            ;;
        -o|--outputDir)
            outputDir="$2"
            shift 2
            ;;
        -i|--iniFile)
            iniFile="$2"
            shift 2
            ;;
        --condaLoaded)
            condaLoadedFlag="true"
            shift
            ;;
        --condaName)
            condaName="$2"
            shift 2
            ;;
        --condaSrc)
            condaSrc="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            usage
            ;;
    esac
done

#Check require parameters
if [ -z ${tumorBam} ]; then
	echo "Missing -t|--tumorBam parameter"
 	usage
fi
if [ -z ${normalBam} ]; then
	echo "Missing -n|--normalBam parameter"
 	usage
fi
if [ -z ${outputDir} ]; then
	echo "Missing -o|--outputDir parameter"
 	usage
fi



prepare_function() {

	#First steps
	echo "Install dir is ${installDir}"
	if [ -z ${iniFile}]
	echo "Sourcing ini file at ${iniFile}"
	source ${iniFile}
	RAND=$(cat /dev/urandom | tr -cd 'a-zA-Z0-9' | head -c 10)
	echo "Random string for this run is ${RAND}"

	#ACTIVATE CONDA
	if [ "$condaLoadedFlag" = "false" ]; then
		echo "Sourcing ${condaSrc} and then loading ${condaName}"
		source "${condaSrc}/etc/profile.d/conda.sh"
		conda activate "${condaName}"
	else echo "Not activating any conda env"
	fi

	# Test Programs
	echo "Testing if all tools in path:"
	for tool in SAGE gridss gripss AMBER COBALT PURPLE linx samtools bcftools circos snpEff
	do
		if ! command -v ${tool} &> /dev/null; then
	    	echo "${tool} not found, did you load the Conda environment?"
	    	exit
	    else echo "${tool}: OK"
		fi
	done

	#Prepare outDir
	outputDir=$(realpath ${outputDir})
	echo "Creating output directory: ${outputDir}"
	mkdir -p ${outputDir}
	if [ ! -d ${outputDir} ]; then
	    echo "${outputDir} does not exist"
	    exit
	fi
	logDir=${outputDir}/log
	echo "Creating log directory: ${logDir}"
	mkdir -p ${logDir}
	if [ ! -d ${logDir} ]; then
	    echo "${logDir} does not exist"
	    exit
	fi

	#Prepare version logs
	versionLog=${logDir}/version.log
	versionCmd="conda list --explicit > ${versionLog}"
	echo "Submitting version log job to have package versions on ${versionLog}"
	bsub -M 1G -o ${logDir}/versionLog.o -e ${logDir}/versionLog.e "${versionCmd}"

	#Get sample names
	echo "Getting sample names"
	tumorBam=$(realpath ${tumorBam})
	normalBam=$(realpath ${normalBam})
	tumorSample=$(samtools view -H ${tumorBam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
	normalSample=$(samtools view -H ${normalBam} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
	echo "Tumor sample: $tumorSample
Normal sample: $normalSample" | tee ${logDir}/samples.log
}


#####SAGE
sage_function() {
	echo "Testing if needed files exist"
	for file in ${tumorBam} ${normalBam} ${reference} ${hmfActionableCodingPanel} ${hmfSomaticHotspots} ${hmfHighConfidence}; do
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/sage"
			exit
		fi
		echo "${file} exists"
	done

	sageDir=${outputDir}/sage
	sageOutput=${sageDir}/${tumorSample}_${normalSample}.sage.vcf.gz
	sageDone=${sageOutput}.done
	sageJob="SAGE_${RAND}"
	if [ -e ${sageDone} ]; then
		echo "SAGE output exists at ${sageOutput}"
	else
		echo "Creating SAGE output directory in ${sageDir}"
		mkdir -p ${sageDir}
		
		sageCmd="SAGE \
-Xms${sageXms} \
-Xmx${sageXmx} \
-tumor ${tumorSample} \
-tumor_bam ${tumorBam} \
-reference ${normalSample} \
-reference_bam ${normalBam} \
-out ${sageOutput} \
-ref_genome ${reference} \
-threads ${sageThreads} \
-assembly ${sageAssembly} \
-hotspots ${hmfSomaticHotspots} \
-panel_bed ${hmfActionableCodingPanel} \
-high_confidence_bed ${hmfHighConfidence} \
${sageExtraParameters} \
&& touch ${sageDone}"
		echo "SAGE command is:"
		echo "${sageCmd}" | tee "${logDir}/${sageJob}.cmd"
		bsub -M "${sageMem}" -n "${sageThreads}" -o "${logDir}/${sageJob}.o" -e "${logDir}/${sageJob}.e" -J "${sageJob}" "${sageCmd}"
	fi
}


####SAGE-POSTPROCESSING
sage_postprocessing_function(){
	echo "Testing if needed files exist"
	for file in ${hmfGermlinePon}; do
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/sage"
			exit
		fi
		echo "${file} exists"
	done
	sageAnnotatedOutput=${sageOutput/.vcf.gz/.annotated.vcf.gz}
	sageFilteredOutput=${sageOutput/.vcf.gz/.filtered.vcf.gz}
	sageFilteredDone=${sageFilteredOutput}.done
	sageFilterJob="SAGE_FILTER_${RAND}"

	if [ -e ${sageFilteredDone} ]; then
		echo "SAGE filtered output exists at ${sageFilteredOutput}"
	else
		echo "SAGE filtering:"
		sageFilterCmd="bcftools annotate \
-a ${hmfGermlinePon} \
-c PON_COUNT,PON_MAX \
-O z \
-o $sageAnnotatedOutput \
${sageOutput} \
&& bcftools filter \
-e 'PON_COUNT!=\".\" && INFO/TIER=\"HOTSPOT\" && PON_MAX>=5 && PON_COUNT >= 5' -s PON -m+ \
${sageAnnotatedOutput} -O u | \
bcftools filter \
-e 'PON_COUNT!=\".\" && INFO/TIER=\"PANEL\" && PON_MAX>=5 && PON_COUNT >= 2' -s PON -m+ -O u | \
bcftools filter \
-e 'PON_COUNT!=\".\" && INFO/TIER!=\"HOTSPOT\" && INFO/TIER!=\"PANEL\" && PON_COUNT >= 2' -s PON -m+ -O z -o ${sageFilteredOutput} \
&& touch ${sageFilteredDone}"
		echo "SAGE filtering command is:"
		echo "${sageFilterCmd}" | tee "${logDir}/${sageFilterJob}.cmd"
		bsub -w "done(${sageJob})" -M "${sageFilterMem}" -o "${logDir}/${sageFilterJob}.o" -e "${logDir}/${sageFilterJob}.e" \
		-J "${sageFilterJob}" "${sageFilterCmd}"
	fi

	sageSnpeffOutput=${sageFilteredOutput/.vcf.gz/.snpeff.vcf}
	sageSnpeffDone=${sageSnpeffOutput}.done
	sageSnpeffJob="SAGE_SNPEFF_${RAND}"
	if [ -e ${sageSnpeffDone} ]; then
		echo "SAGE snpEff annotated output exists at ${sageSnpeffOutput}"
	else
		echo "SnpEff annotation:"
		sageSnpeffCmd="snpEff -i vcf -o vcf ${sageSnpeffGenomeBuild} \
${sageFilteredOutput} ${sageSnpeffExtraParameters} \
> ${sageSnpeffOutput}; \
gzip ${sageSnpeffOutput} \
&& touch ${sageSnpeffDone}"
		echo "SAGE SnpEff command is:"
		echo "${sageSnpeffCmd}" | tee "${logDir}/${sageSnpeffJob}.cmd"
		bsub -w "done(${sageFilterJob})" -M "${sageSnpeffMem}" \
		-o "${logDir}/${sageSnpeffJob}.o" -e "${logDir}/${sageSnpeffJob}.e" \
		-J "${sageSnpeffJob}" "${sageSnpeffCmd}"
	fi

	sageFinalOutput=${sageSnpeffOutput}.gz

}


#####AMBER
amber_function() {
	echo "Testing if needed files exist"
	for file in ${tumorBam} ${normalBam} ${hmfHetPonLoci}; do		
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. 
			Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/amber"
			exit
		fi
		echo "${file} exists"
	done
	
	amberDir=${outputDir}/amber
	amberDone=${amberDir}/AMBER.done
	amberJob="AMBER_${RAND}"

	if [ -e ${amberDone} ]; then
		echo "AMBER output exists at ${amberDir}"
	else
		echo "Creating AMBER output directory in ${amberDir}"
		mkdir -p ${amberDir}
		
		amberCmd="AMBER \
-Xms${amberXms} \
-Xmx${amberXmx} \
-tumor ${tumorSample} \
-tumor_bam ${tumorBam} \
-reference ${normalSample} \
-reference_bam ${normalBam} \
-output_dir ${amberDir} \
-threads ${amberThreads} \
-loci ${hmfHetPonLoci} \
${amberExtraParameters} \
&& touch ${amberDone}"
		echo "AMBER command is:"
		echo "${amberCmd}" | tee "${logDir}/${amberJob}.cmd"
		bsub -M "${amberMem}" -n "${amberThreads}" -o "${logDir}/${amberJob}.o" -e "${logDir}/${amberJob}.e" -J "${amberJob}" "${amberCmd}"
	fi
}


#####COBALT
cobalt_function() {
	echo "Testing if needed files exist"
	for file in ${tumorBam} ${normalBam} ${cobaltLoci}; do		
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. 
			Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/cobalt"
			exit
		fi
		echo "${file} exists"
	done
	
	cobaltDir=${outputDir}/cobalt
	cobaltDone=${cobaltDir}/COBALT.done
	cobaltJob="COBALT_${RAND}"

	if [ -e ${cobaltDone} ]; then
		echo "COBALT output exists at ${cobaltDir}"
	else
		echo "Creating COBALT output directory in ${cobaltDir}"
		mkdir -p ${cobaltDir}
		
		cobaltCmd="COBALT \
-Xms${cobaltXms} \
-Xmx${cobaltXmx} \
-tumor ${tumorSample} \
-tumor_bam ${tumorBam} \
-reference ${normalSample} \
-reference_bam ${normalBam} \
-output_dir ${cobaltDir} \
-threads ${cobaltThreads} \
-gc_profile ${hmfGcProfile} \
${cobaltExtraParameters} \
&& touch ${cobaltDone}"
		echo "COBALT command is:"
		echo "${cobaltCmd}" | tee "${logDir}/${cobaltJob}.cmd"
		bsub -M "${cobaltMem}" -n "${cobaltThreads}" -o "${logDir}/${cobaltJob}.o" -e "${logDir}/${cobaltJob}.e" -J "${cobaltJob}" "${cobaltCmd}"
	fi
}

#####GRIDSS
gridss_function() {
	
	echo "Testing if needed files exist"
	for file in ${tumorBam} ${normalBam} ${reference}; do		
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. 
			Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/cobalt"
			exit
		fi
		echo "${file} exists"
	done
		
	gridssDir=${outputDir}/gridss
	gridssTmpDir=${gridssDir}/gridssTmp
	gridssOutput=${gridssDir}/${tumorSample}_${normalSample}.gridss.vcf.gz
	gridssAssembly=${gridssDir}/${tumorSample}_${normalSample}.assembly.gridss.bam
	gridssDone=${gridssOutput}.done
	gridssJob="GRIDSS_${RAND}"

	if [ -e ${gridssDone} ]; then
		echo "GRIDSS output exists at ${gridssOutput}"
	else
		echo "Creating GRIDSS output directory in ${gridssDir}"
		mkdir -p ${gridssDir}
		mkdir -p ${gridssTmpDir}
		
		gridssCmd="gridss \
--jvmheap ${gridssJvmHeap} \
--otherjvmheap ${gridssOtherJvmHeap} \
--reference ${reference} \
--output ${gridssOutput} \
--assembly ${gridssAssembly} \
--threads ${gridssThreads} \
--workingdir ${gridssTmpDir} \
--labels ${tumorSample},${normalSample} \
${gridssExtraParameters} \
$tumorBam \
$normalBam \
&& touch ${gridssDone}"
		echo "GRIDSS command is:"
		echo "${gridssCmd}" | tee "${logDir}/${gridssJob}.cmd"
		bsub -M "${gridssMem}" -n "${gridssThreads}" -o "${logDir}/${gridssJob}.o" -e "${logDir}/${gridssJob}.e" -J "${gridssJob}" "${gridssCmd}"
	fi
}

####GRIDSS-POSTPROCESSING
gridss_postprocessing_function() {
	#Find GRIDSS directory
	#Taken from the conda gridss.sh script, to find jar
	gridssSource="$(which gridss)"
	while [ -h "$gridssSource" ]; do # resolve $gridssSource until the file is no longer a symlink
	    gridssDir="$( cd -P "$( dirname "$gridssSource" )" && pwd )"
	    gridssSource="$(readlink "$gridssSource")"
	    [[ $gridssSource != /* ]] && gridssSource="$gridssDir/$gridssSource" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	gridssDir="$( cd -P "$( dirname "$gridssSource" )" && pwd )"
	gridssJar="${gridssDir}/gridss.jar"
	repeatMaskerScript="${gridssDir}/gridss_annotate_vcf_repeatmasker.sh"
	krakenScript="${gridssDir}/gridss_annotate_vcf_kraken2.sh"

	#Repeatmasker

	gridssRmOutput=${gridssOutput/.vcf.gz/.repeatmasker.vcf.gz}
	gridssRmDone=${gridssRmOutput}.done
	gridssRmJob="GRIDSS_RM_${RAND}"
	if [ -e ${gridssRmDone} ]; then
		echo "GRIDSS-RM output exists at ${gridssRmOutput}"
	else
		gridssRmCmd="bash ${repeatMaskerScript} \
-j ${gridssJar} \
-o ${gridssRmOutput} \
-w ${gridssTmpDir} \
-t ${gridssRmThreads} ${gridssRmExtraParameters} ${gridssOutput} \
&& touch ${gridssRmDone}"

		echo "GRIDSS-REPEATMASKER command is:"
		echo "${gridssRmCmd}" | tee "${logDir}/${gridssRmJob}.cmd"

		bsub -w "done(${gridssJob})" -M "${gridssRmMem}" \
		-n "${gridssRmThreads}" -o "${logDir}/${gridssRmJob}.o" -e "${logDir}/${gridssRmJob}.e" \
		-J "${gridssRmJob}" "${gridssRmCmd}"
	fi

	#Kraken
	gridssKrOutput=${gridssRmOutput/.vcf.gz/.kraken2.vcf.gz}
	gridssKrDone=${gridssKrOutput}.done
	gridssKrJob="GRIDSS_KR_${RAND}"
	if [ -e ${gridssKrDone} ]; then
		echo "GRIDSS-KR output exists at ${gridssKrOutput}"
	else
		gridssKrCmd="bash ${krakenScript} \
--kraken2db ${kraken2db} \
-o ${gridssKrOutput} \
-j ${gridssJar} \
-w ${gridssTmpDir} \
-t ${gridssKrThreads} \
${gridssKrExtraParameters} \
${gridssRmOutput} \
&& touch ${gridssKrDone}"
		echo "GRIDSS-KRAKEN command is:"
		echo "${gridssKrCmd}" | tee "${logDir}/${gridssKrJob}.cmd"
		bsub -w "done(${gridssRmJob})" -M "${gridssKrMem}" \
		-n "${gridssKrThreads}" -o "${logDir}/${gridssKrJob}.o" -e "${logDir}/${gridssKrJob}.e" \
		-J "${gridssKrJob}" "${gridssKrCmd}"
	fi
}

####GRIPSS
gripss_function() {
	echo "Testing if needed files exist"
	for file in ${reference} ${hmfBreakendPon} ${hmfBreakpointPon} ${hmfBreakpointHotspot}; do		
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. 
			Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/gripss"
			exit
		fi
		echo "${file} exists"
	done
	
	gripssDir=${outputDir}/gripss
	gripssOutput=${gripssDir}/${tumorSample}_${normalSample}.somatic.vcf.gz
	gripssFilteredOutput=${gripssOutput/.vcf.gz/.filtered.vcf.gz}
	gripssDone=${gripssFilteredOutput}.done
	gripssJob="GRIPSS_${RAND}"

	#Find GRIPSS jar
	#Taken from the conda gridss.sh script, to find jar
	gripssSource="$(which gripss)"
	while [ -h "$gripssSource" ]; do # resolve $gripssSource until the file is no longer a symlink
	    gripssDir="$( cd -P "$( dirname "$gripssSource" )" && pwd )"
	    gripssSource="$(readlink "$gripssSource")"
	    [[ $gripssSource != /* ]] && gripssSource="$gripssDir/$gripssSource" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	gripssDir="$( cd -P "$( dirname "$gripssSource" )" && pwd )"
	gripssJar="${gripssDir}/gripss.jar"

	if [ -e ${gripssFilteredDone} ]; then
		echo "GRIPSS output exists at ${gripssOutput}"
	else
		echo "Creating GRIPSS output directory in ${gripssDir}"
		mkdir -p ${gripssDir}
		
		gripssCmd="gripss \
-Xms${gripssXms} \
-Xmx${gripssXmx} \
-tumor ${tumorSample} \
-reference ${normalSample} \
-input_vcf ${gridssKrOutput} \
-output_vcf ${gripssOutput} \
-ref_genome ${reference} \
-breakend_pon ${hmfBreakendPon} \
-breakpoint_pon ${hmfBreakpointPon} \
-breakpoint_hotspot ${hmfBreakpointHotspot} \
&& java \
-Xms${gripssXms} \
-Xmx${gripssXmx} \
-cp ${gripssJar} com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
-input_vcf ${gripssOutput} \
-output_vcf ${gripssFilteredOutput} \
&& touch ${gripssDone}"
		echo "GRIPSS command is:"
		echo "${gripssCmd}" | tee "${logDir}/${gripssJob}.cmd"
		bsub -w "done(${gridssKrJob})" \
		-M "${gripssMem}" -o "${logDir}/${gripssJob}.o" -e "${logDir}/${gripssJob}.e" \
		-J "${gripssJob}" "${gripssCmd}"
	fi
}

####PURPLE
purple_function() {
	echo "Testing if needed files exist"
	for file in ${reference} ${hmfGcProfile} ${hmfSomaticHotspots} ${hmfGermlineHotspots} ${hmfDriverGenePanel}; do		
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. 
			Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/purple"
			exit
		fi
		echo "${file} exists"
	done
	purpleDir=${outputDir}/purple
	purpleDone=${purpleDir}/PURPLE.done
	purpleJob="PURPLE_${RAND}"

	if [ -e ${purpleDone} ]; then
		echo "PURPLE output exists at ${purpleDir}"
	else
		echo "Creating PURPLE output directory in ${purpleDir}"
		mkdir -p ${purpleDir}
		
		purpleCmd="PURPLE \
-Xms${purpleXms} \
-Xmx${purpleXmx} \
-tumor ${tumorSample} \
-reference ${normalSample} \
-output_dir ${purpleDir} \
-amber ${amberDir} \
-cobal ${cobaltDir} \
-gc_profile ${hmfGcProfile} \
-ref_genome ${reference} \
-somatic_vcf ${sageFinalOutput} \
-structural_vcf ${gripssFilteredOutput} \
-sv_recovery_vcf ${gripssOutput} \
-driver_catalog \
-somatic_hotspots ${hmfSomaticHotspots} \
-driver_gene_panel ${hmfDriverGenePanel} \
${purpleExtraParameters} \
&& touch ${purpleDone}"
		
		echo "PURPLE command is:"
		echo "${purpleCmd}" | tee "${logDir}/${purpleJob}.cmd"
		bsub -w "done(${sageSnpeffJob}) && done(${amberJob}) && done(${cobaltJob}) && done(${gripssJob})" \
		-n "${purpleThreads}" -M "${purpleMem}" -o "${logDir}/${purpleJob}.o" -e "${logDir}/${purpleJob}.e" \
		-J "${purpleJob}" "${purpleCmd}"
	fi
}

####LINX
linx_function() {
	echo "Testing if needed files exist"
	for file in ${reference} ${hmfFragileSites} ${hmfLineElements} ${hmfReplicationOrigins} ${hmfDriverGenePanel} \
	${hmfViralHosts} ${hmfKnownFusions}; do		
		if [ ! -f ${file} ]; then
			bfile=$(basename ${file})
			echo "File ${bfile} not found, looked at path ${file} does not exist. 
			Please provide correct path or check https://github.com/hartwigmedical/hmftools/tree/master/gripss"
			exit
		fi
		echo "${file} exists"
	done
	if [ ! -d ${hmfGeneTranscripts} ]; then
		echo "Directory with ensembl cache does not exist, check https://github.com/hartwigmedical/hmftools/tree/master/linx"
		exit
	fi
	echo "${hmfGeneTranscripts} also exists"

	purpleOutputVcf=${purpleDir}/${tumorSample}.purple.sv.vcf.gz
	linxDir=${outputDir}/linx
	linxDone=${linxDir}/LINX.done
	linxJob="LINX_${RAND}"
	if [ -e ${linxDone} ]; then
		echo "LINX output exists at ${linxDir}"
	else
		echo "Creating LINX output directory in ${linxDir}"
		mkdir -p ${linxDir}
		linxCmd="linx \
-Xms${linxXms} \
-Xmx${linxXmx} \
-sample ${tumorSample} \
-sv_vcf ${purpleOutputVcf} \
-purple_dir ${purpleDir} \
-output_dir ${linxDir} \
-ref_genome_version ${linxGenomeVersion} \
-fragile_site_file ${hmfFragileSites} \
-line_element_file ${hmfLineElements} \
-replication_origins_file ${hmfReplicationOrigins} \
-viral_hosts_file ${hmfViralHosts} \
-gene_transcripts_dir ${hmfGeneTranscripts} \
-check_fusions \
-known_fusion_file ${hmfKnownFusions} \
-check_drivers \
-log_debug \
-write_vis_data \
${linxExtraParameters} \
&& touch ${linxDone}"
		echo "LINX command is:"
		echo "${linxCmd}" | tee "${logDir}/${linxJob}.cmd"
		bsub -w "done(${purpleJob})" \
		-M "${linxMem}" -o "${logDir}/${linxJob}.o" -e "${logDir}/${linxJob}.e" \
		-J "${linxJob}" "${linxCmd}"
	fi

	linxVizPlotDir=${linxDir}/plot
	linxVizDataDir=${linxDir}/data
	linxVizDone=${linxDir}/LINX-VIZ.done
	linxVizJob="LINX_VIZ_${RAND}"

	if [ -e ${linxVizDone} ]; then
		echo "LINX-VIZ output exists at ${linxDir}"
	else
		echo "Creating LINX-VIZ output directory in ${linxDir}"
		mkdir -p ${linxVizPlotDir}
		mkdir -p ${linxVizDataDir}

		#Find LINX jar
		#Taken from the conda gridss.sh script, to find jar
		linxSource="$(which linx)"
		while [ -h "$linxSource" ]; do # resolve $linxSource until the file is no longer a symlink
		    linxDir="$( cd -P "$( dirname "$linxSource" )" && pwd )"
		    linxSource="$(readlink "$linxSource")"
		    [[ $linxSource != /* ]] && linxSource="$linxDir/$linxSource" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
		done
		linxDir="$( cd -P "$( dirname "$linxSource" )" && pwd )"
		linxJar="${linxDir}/sv-linx.jar"
		linxVizCmd="java \
-Xms${linxVizXms} \
-Xmx${linxVizXmx} \
-cp ${linxJar} com.hartwig.hmftools.linx.visualiser.SvVisualiser \
-sample ${tumorSample} \
-gene_transcripts_dir ${hmfGeneTranscripts} \
-plot_out ${linxVizPlotDir} \
-data_out ${linxVizDataDir} \
-segment ${linxDir}/${tumorSample}.linx.vis_segments.tsv \
-link ${linxDir}/${tumorSample}.linx.vis_sv_data.tsv \
-exon ${linxDir}/${tumorSample}.linx.vis_gene_exon.tsv \
-cna ${linxDir}/${tumorSample}.linx.vis_copy_number.tsv \
-protein_domain ${linxDir}/${tumorSample}.linx.vis_protein_domain.tsv \
-fusion ${linxDir}/${tumorSample}.linx.vis_fusion.tsv \
-threads ${linxVizThreads} \
${linxVizExtraParameters} \
&& touch ${linxVizDone}"
		echo "LINX-VIZ command is:"

		echo "${linxVizCmd}" | tee "${logDir}/${linxVizJob}.cmd"
		bsub -w "done(${linxJob})" -n "${linxVizThreads}" \
		-M "${linxVizMem}" -o "${logDir}/${linxVizJob}.o" -e "${logDir}/${linxVizJob}.e" \
		-J "${linxVizJob}" "${linxVizCmd}"
	fi
}

echo "###Running preparations"
prepare_function
echo "###Preparations finished"
echo ""
echo "###Preparing SAGE"
sage_function
echo "###"
echo ""
echo "###Preparing SAGE-POSTPROCESSING"
sage_postprocessing_function
echo "###"
echo ""
echo "###Preparing AMBER"
amber_function
echo "###"
echo ""
echo "###Preparing COBALT"
cobalt_function
echo "###"
echo ""
echo "###Preparing GRIDSS"
gridss_function
echo "###"
echo ""
echo "###Preparing GRIDSS postprocessing"
gridss_postprocessing_function
echo "###"
echo ""
echo "###Preparing GRIPSS"
gripss_function
echo "###"
echo ""
echo "###Preparing PURPLE"
purple_function
echo "###"
echo ""
echo "###Preparing LINX"
linx_function
echo "###"
echo ""






