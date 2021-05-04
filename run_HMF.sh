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
	versionLog="${logDir}/version_${RAND}.log"
	versionCmd="conda list --explicit > ${versionLog}"
	versionJob="version_${RAND}"
	echo "Submitting version log job to have package versions on ${versionLog}"
	bsub -M 1G -J "${versionJob}" -o "${logDir}/${versionJob}.o" -e "${logDir}/${versionJob}.e" "\"${versionCmd}\""

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
	sageFlag=true
	if [ -e ${sageDone} ]; then
		echo "SAGE output exists at ${sageOutput}"
		sageFlag=false
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
		sageBsub=(bsub
			-M "${sageMem}"
			-n "${sageThreads}"
			-o "${logDir}/${sageJob}.o"
			-e "${logDir}/${sageJob}.e"
			-J "${sageJob}"
			"${sageCmd}")
	"${sageBsub[@]}"
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
	sageFilterFlag=true
	if [ -e ${sageFilteredDone} ]; then
		echo "SAGE filtered output exists at ${sageFilteredOutput}"
		sageFilterFlag=false
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

		sageFilterBsub=(bsub
		-M "${sageFilterMem}"
		-o "${logDir}/${sageFilterJob}.o"
		-e "${logDir}/${sageFilterJob}.e"
		-J "${sageFilterJob}"
		)
		if [ "${sageFlag}" = "true" ]; then
			sageFilterBsub+=("-w done(${sageJob})")
		fi
		sageFilterBsub+=("${sageFilterCmd}")
		"${sageFilterBsub[@]}"
		
	fi

	sageSnpeffOutput=${sageFilteredOutput/.vcf.gz/.snpeff.vcf}
	sageSnpeffDone=${sageSnpeffOutput}.gz.done
	sageSnpeffJob="SAGE_SNPEFF_${RAND}"
	sageSnpeffFlag=true
	if [ -e ${sageSnpeffDone} ]; then
		echo "SAGE snpEff annotated output exists at ${sageSnpeffOutput}.gz"
		sageSnpeffFlag=false
	else
		echo "SnpEff annotation:"
		sageSnpeffCmd="snpEff \
-Xms${sageSnpeffXms} \
-Xmx${sageSnpeffXmx} \
-i vcf \
-o vcf \
${sageSnpeffGenomeBuild} \
${sageFilteredOutput} \
${sageSnpeffExtraParameters} \
> ${sageSnpeffOutput} \
&& gzip ${sageSnpeffOutput} \
&& touch ${sageSnpeffDone}"
		echo "SAGE SnpEff command is:"
		echo "${sageSnpeffCmd}" | tee "${logDir}/${sageSnpeffJob}.cmd"
		sageSnpeffBsub=(bsub
			-M "${sageSnpeffMem}"
			-o "${logDir}/${sageSnpeffJob}.o"
			-e "${logDir}/${sageSnpeffJob}.e"
			-J "${sageSnpeffJob}"
			)
		if [ "${sageFilterFlag}" = "true" ]; then
			sageSnpeffBsub+=("-w done(${sageFilterJob})")
		fi
		sageSnpeffBsub+=("${sageSnpeffCmd}")
		"${sageSnpeffBsub[@]}"
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
	amberFlag=true

	if [ -e ${amberDone} ]; then
		echo "AMBER output exists at ${amberDir}"
		amberFlag=false
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
		amberBsub=(bsub
			-M "${amberMem}"
			-n "${amberThreads}"
			-o "${logDir}/${amberJob}.o"
			-e "${logDir}/${amberJob}.e"
			-J "${amberJob}"
			"${amberCmd}"
		)
		"${amberBsub[@]}"
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
	cobaltFlag=true

	if [ -e ${cobaltDone} ]; then
		echo "COBALT output exists at ${cobaltDir}"
		cobaltFlag=false
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
		cobaltBsub=(bsub
			-M "${cobaltMem}"
			-n "${cobaltThreads}"
			-o "${logDir}/${cobaltJob}.o"
			-e "${logDir}/${cobaltJob}.e"
			-J "${cobaltJob}"
			"${cobaltCmd}")
		"${cobaltBsub[@]}"
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
	gridssFlag=true

	if [ -e ${gridssDone} ]; then
		echo "GRIDSS output exists at ${gridssOutput}"
		gridssFlag=false
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
		gridssBsub=(bsub
			-M "${gridssMem}"
			-n "${gridssThreads}"
			-o "${logDir}/${gridssJob}.o"
			-e "${logDir}/${gridssJob}.e"
			-J "${gridssJob}"
			"${gridssCmd}")
		"${gridssBsub[@]}"
	fi
}

####GRIDSS-POSTPROCESSING
gridss_postprocessing_function() {
	#Find GRIDSS directory
	#Taken from the conda gridss.sh script, to find jar
	gridssSource="$(which gridss)"
	while [ -h "${gridssSource}" ]; do # resolve $gridssSource until the file is no longer a symlink
	    gridssSrcDir="$( cd -P "$( dirname "${gridssSource}" )" && pwd )"
	    gridssSource="$(readlink "${gridssSource}")"
	    [[ ${gridssSource} != /* ]] && gridssSource="${gridssSrcDir}/${gridssSource}" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	gridssSrcDir="$( cd -P "$( dirname "${gridssSource}" )" && pwd )"
	gridssJar="${gridssSrcDir}/gridss.jar"
	repeatMaskerScript="${gridssSrcDir}/gridss_annotate_vcf_repeatmasker.sh"
	krakenScript="${gridssSrcDir}/gridss_annotate_vcf_kraken2.sh"

	#Repeatmasker

	gridssRmOutput=${gridssOutput/.vcf.gz/.repeatmasker.vcf.gz}
	gridssRmDone=${gridssRmOutput}.done
	gridssRmJob="GRIDSS_RM_${RAND}"
	gridssRmFlag=true
	if [ -e ${gridssRmDone} ]; then
		echo "GRIDSS-RM output exists at ${gridssRmOutput}"
		gridssRmFlag=false
	else
		gridssRmCmd="bash ${repeatMaskerScript} \
-j ${gridssJar} \
-o ${gridssRmOutput} \
-w ${gridssTmpDir} \
-t ${gridssRmThreads} ${gridssRmExtraParameters} ${gridssOutput} \
&& touch ${gridssRmDone}"

		echo "GRIDSS-REPEATMASKER command is:"
		echo "${gridssRmCmd}" | tee "${logDir}/${gridssRmJob}.cmd"
		gridssRmBsub=(bsub
			-M "${gridssRmMem}"
			-n "${gridssRmThreads}"
			-o "${logDir}/${gridssRmJob}.o"
			-e "${logDir}/${gridssRmJob}.e"
			-J "${gridssRmJob}"
			)
		if [ "${gridssFlag}" = "true" ]; then
			gridssRmBsub+=("-w done(${gridssJob})")
		fi
		gridssRmBsub+=("${gridssRmCmd}")
		"${gridssRmBsub[@]}"
	fi

	#Kraken
	gridssKrOutput=${gridssRmOutput/.vcf.gz/.kraken2.vcf.gz}
	gridssKrDone=${gridssKrOutput}.done
	gridssKrJob="GRIDSS_KR_${RAND}"
	gridssKrFlag=true
	if [ -e ${gridssKrDone} ]; then
		echo "GRIDSS-KR output exists at ${gridssKrOutput}"
		gridssKrFlag=false
	else
		gridssKrCmd="sh ${krakenScript} \
--kraken2db ${kraken2db} \
-o ${gridssKrOutput} \
-j ${gridssJar} \
-t ${gridssKrThreads} \
${gridssKrExtraParameters} \
${gridssRmOutput} \
&& touch ${gridssKrDone}"
		echo "GRIDSS-KRAKEN command is:"
		echo "${gridssKrCmd}" | tee "${logDir}/${gridssKrJob}.cmd"
		gridssKrBsub=(bsub
			-M "${gridssKrMem}"
			-n "${gridssKrThreads}"
			-o "${logDir}/${gridssKrJob}.o"
			-e "${logDir}/${gridssKrJob}.e"
			-J "${gridssKrJob}"
			)
		if [ "${gridssRmFlag}" = "true" ]; then
			gridssKrBsub+=("-w done(${gridssRmJob})")
		fi
		gridssKrBsub+=("${gridssKrCmd}")
		"${gridssKrBsub[@]}"
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
	gripssFlag=true
	#Find GRIPSS jar
	#Taken from the conda gridss.sh script, to find jar
	gripssSource="$(which gripss)"
	while [ -h "${gripssSource}" ]; do # resolve $gripssSource until the file is no longer a symlink
	    gripssSrcDir="$( cd -P "$( dirname "${gripssSource}" )" && pwd )"
	    gripssSource="$(readlink "${gripssSource}")"
	    [[ ${gripssSource} != /* ]] && gripssSource="${gripssSrcDir}/${gripssSource}" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	gripssSrcDir="$( cd -P "$( dirname "${gripssSource}" )" && pwd )"
	gripssJar="${gripssSrcDir}/gripss.jar"

	if [ -e ${gripssDone} ]; then
		echo "GRIPSS output exists at ${gripssFilteredOutput}"
		gripssFlag=false
	else
		echo "Creating GRIPSS output directory in ${gripssDir}"
		mkdir -p "${gripssDir}"
		
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
		gripssBsub=(bsub
			-M "${gripssMem}"
			-o "${logDir}/${gripssJob}.o"
			-e "${logDir}/${gripssJob}.e"
			-J "${gripssJob}"
			)
		if [ "${gridssKrFlag}" = "true" ]; then
			gripssBsub+=("-w done(${gridssKrJob})")
		fi
		gripssBsub+=("${gripssCmd}")
		"${gripssBsub[@]}"
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
	purpleFlag=true

	if [ -e ${purpleDone} ]; then
		echo "PURPLE output exists at ${purpleDir}"
		purpleFlag=false
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
-cobalt ${cobaltDir} \
-gc_profile ${hmfGcProfile} \
-ref_genome ${reference} \
-somatic_vcf ${sageFinalOutput} \
-structural_vcf ${gripssFilteredOutput} \
-sv_recovery_vcf ${gripssOutput} \
-driver_catalog \
-somatic_hotspots ${hmfSomaticHotspots} \
-driver_gene_panel ${hmfDriverGenePanel} \
-circos circos \
${purpleExtraParameters} \
&& touch ${purpleDone}"		
		echo "PURPLE command is:"
		echo "${purpleCmd}" | tee "${logDir}/${purpleJob}.cmd"
		purpleBsub=(bsub
		-n "${purpleThreads}"
		-M "${purpleMem}"
		-o "${logDir}/${purpleJob}.o"
		-e "${logDir}/${purpleJob}.e"
		-J "${purpleJob}"
		)

		purpleDependency=()
		if [ "${sageSnpeffFlag}" = "true" ]; then
			purpleDependency+=("done(${sageSnpeffJob})")
		fi
		if [ "${amberFlag}" = "true" ]; then
			purpleDependency+=("done(${amberJob})")
		fi
		if [ "${cobaltFlag}" = "true" ]; then
			purpleDependency+=("done(${cobaltJob})")
		fi
		if [ "${gripssFlag}" = "true" ]; then
			purpleDependency+=("done(${gripssJob})")
		fi
		if [ ! ${#purpleDependency[@]} -eq 0 ]; then
			purpleBsub+=("-w $(echo ${purpleDependency[@]} | sed 's: : \&\& :g')")
		fi
		purpleBsub+=("${purpleCmd}")
		"${purpleBsub[@]}"
	fi
}

####LINX
linx_function() {
	echo "Testing if needed files exist"
	for file in ${reference} ${hmfFragileSites} ${hmfLineElements} ${hmfReplicationOrigins} ${hmfDriverGenePanel} \
	${hmfViralHosts} ${hmfKnownFusionData}; do		
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
	linxFlag=true
	if [ -e ${linxDone} ]; then
		echo "LINX output exists at ${linxDir}"
		linxFlag=false
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
-known_fusion_file ${hmfKnownFusionData} \
-check_drivers \
-log_debug \
-write_vis_data \
${linxExtraParameters} \
&& touch ${linxDone}"
		echo "LINX command is:"
		echo "${linxCmd}" | tee "${logDir}/${linxJob}.cmd"
		linxBsub=(bsub
		-M "${linxMem}"
		-o "${logDir}/${linxJob}.o"
		-e "${logDir}/${linxJob}.e"
		-J "${linxJob}"
		)
		if [ "${purpleFlag}" = "true" ]; then
			linxBsub+=("-w done(${purpleJob})")
		fi
		linxBsub+=("${linxCmd}")
		"${linxBsub[@]}"
	fi

	linxVizPlotDir=${linxDir}/plot
	linxVizDataDir=${linxDir}/data
	linxVizDone=${linxDir}/LINX-VIZ.done
	linxVizJob="LINX_VIZ_${RAND}"
	linxVizFlag=true

	if [ -e ${linxVizDone} ]; then
		echo "LINX-VIZ output exists at ${linxDir}"
		linxVizFlag=false
	else
		echo "Creating LINX-VIZ output directory in ${linxDir}"
		mkdir -p ${linxVizPlotDir}
		mkdir -p ${linxVizDataDir}

		#Find LINX jar
		#Taken from the conda gridss.sh script, to find jar
		linxSource="$(which linx)"
		while [ -h "${linxSource}" ]; do # resolve $linxSource until the file is no longer a symlink
		    linxSrcDir="$( cd -P "$( dirname "${linxSource}" )" && pwd )"
		    linxSource="$(readlink "${linxSource}")"
		    [[ $linxSource != /* ]] && linxSource="${linxSrcDir}/${linxSource}" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
		done
		linxSrcDir="$( cd -P "$( dirname "${linxSource}" )" && pwd )"
		linxJar="${linxSrcDir}/sv-linx.jar"
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
-circos circos
${linxVizExtraParameters} \
&& touch ${linxVizDone}"
		echo "LINX-VIZ command is:"

		echo "${linxVizCmd}" | tee "${logDir}/${linxVizJob}.cmd"
		linxVizBsub=(bsub
			-n "${linxVizThreads}"
			-M "${linxVizMem}"
			-o "${logDir}/${linxVizJob}.o"
			-e "${logDir}/${linxVizJob}.e"
			-J "${linxVizJob}"
			)
		if [ "${linxFlag}" = "true" ]; then
			linxVizBsub+=("-w done(${linxJob})")
		fi
		linxVizBsub+=("${linxVizCmd}")
		"${linxVizBsub[@]}"
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






