#!/bin/bash

#Jose Espejo Valle-Inclan 2021

#Clean up script for HMF pipeline

#Defaults
# mode="soft"
forceFlag="false"

USAGE_MESSAGE="Usage: clean_HMF.sh [options] -o <outputDir>
Cleans up the HMF run directory, leaving just the input files needed to regenerate \
PURPLE and LINX output.

Required parameters:
	-o/--outputDir: path to output directory of the HMF run

Optional parameters:
	-h/--help: show this usage help.
	-f/--force: Ignore done files sanity check and force removal.
	"
  # -m/--mode: Clean up mode. Either soft, hard or extreme [soft].
  #   soft: Removes intermediate files but leaves the output of each step
  #   hard: Also removes output files that are unlikely to be useful, leaving only
  #       the somatic output of PURPLE and LINX. Leaves the done files
  #   extreme: Remove even the done files. Leaves only the purple and linx folders
usage () {
	echo "${USAGE_MESSAGE}"
	exit 1
}

OPTIONS="hfo:"
LONGOPTS="help,force,outputDir:"
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
    # -m|--mode)
    #   mode="$2"
    #   shift 2
    #   ;;
    -o|--outputDir)
      outputDir="$2"
      shift 2
      ;;
    -f|--force)
      forceFlag="true"
      shift
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

#Check required parameters
if [ -z ${outputDir} ]; then
	echo "Missing -o|--outputDir parameter"
 	usage
fi

echo "Cleaning up HMF run for directory: ${outputDir}"
# echo "Mode is ${mode}"

if compgen -G ${outputDir}/HMF_*.done > /dev/null; then
    echo "Final done file exists in ${outputDir}"
  else
    if [ ${forceFlag} == "false" ]; then
      echo "No done file in ${outputDir}, did you run finish? Otherwise, use --force to delete"
    else
      echo "No done file in ${outputDir} but --force is present, forcing deletion"
    fi
fi

# Soft mode
echo "Looking for SAGE intermediate files"
fileList=(${outputDir}/sage/*.sage.vcf.gz
${outputDir}/sage/*.sage.vcf.gz.tbi
${outputDir}/sage/*.sage.filtered.vcf.gz
${outputDir}/sage/*.sage.annotated.vcf.gz
)

echo "Looking for GRIDSS intermediate files"
fileList+=(${outputDir}/gridss/*.assembly.gridss.bam
${outputDir}/gridss/gridssTmp/
${outputDir}/gridss/*.gridss.vcf.gz
${outputDir}/gridss/*.gridss.vcf.gz.tbi
${outputDir}/gridss/*.gridss.repeatmasker.vcf.gz
${outputDir}/gridss/*.gridss.repeatmasker.vcf.gz.tbi
)
echo "Deleting files: ${fileList[@]}"
rm -r "${fileList[@]}"
