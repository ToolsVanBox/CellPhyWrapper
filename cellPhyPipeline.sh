#!/bin/bash

#SBATCH --time=8:0:0
#SBATCH -c 8


# define default parameter settings ---------------
ERROR="true"
MODEL="GT16+FO"
INPUTTYPE="GL"
THREADS=8
OUTPUTDIR=$(pwd)
SEED=2
NBOOTSTRAP=100
CELLPHY=/hpc/pmc_vanboxtel/tools/external/cellphy/cellphy.sh
SUPPORT_MAP=/hpc/pmc_vanboxtel/tools/external/cellphy/script/support-map.R
OUTGROUP="NONE"
PTATODIR="NONE"

#SOURCE=$(dirname ${BASH_SOURCE[0]} )
#echo "${SLURM_JOB_ID}"
if [ -n $SLURM_JOB_ID ];  then
        # check the original location through scontrol and $SLURM_JOB_ID
        SOURCE=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | awk '{print $1}' |  sed 's:/[^/]*$::')
        #echo "Sbatch start"
        #echo "${SCRIPT_PATH}"
#else
#        # otherwise: started with bash. Get the real location.
#        SCRIPT_PATH=$(realpath $0)
#        echo "Bash start"
#        echo "${SCRIPT_PATH}"
fi
#echo "${SCRIPT_PATH}"


# Load module since R is not standard available on the HPC
module load R/4.3.0

# user message ---------------
usage_msg="Usage: Cellphy Wrapper

Required command line arguments:
        --input: input multi-sample VCF
	    --outgr: the full name of the sample that should be used to root the tree (e.g. a control sample, like MSCs)
Optional parameters:
        --ptato: directory with final PTATO vcfs folders                            
        --error: whether an error/FP rate should be assumed. Turn off when having only clonally-expanded cells ($ERROR)
        --model: the model to use. "GT16+FO" is the most comprehensive model. "GT10+FO" can be selected when no phasing info is present ($MODEL)
        --inputtype: the type of input data. Should be genotype-likelihood ("GL") for best results. If not present genotype ("GT") can be selected ($INPUTTYPE)
        --threads: number of threads to use ($THREADS)
        --outputdir: output directory, if not current ($OUTPUTDIR)
        --seed: seed to start from ($SEED)
        --nbootstrap: number of bootstraps ($NBOOTSTRAP)
        --prefix: a prefix for naming output names. If no prefix is provided, the vcf name will be taken
        --percentage: a percentage of cells that can be missed in a specific branch
        --plottingonly: Whether only the plotting part needs to be rerun (TRUE/FALSE)
        --help: print this message
"
usage() {
        echo "$usage_msg" 1>&2
        exit $EX_USAGE
}
if [[ $# -eq 0 ]] ; then
    usage
    exit 0
fi

echo $@
POSITIONAL=()

# capture user input ---------------
while [[ $# -gt 0 ]]; do
        case "$1" in
                --input)
                        INPUT="$2"
                        shift 2
                        ;;
	        --outgr)
	            	OUTGROUP="$2"
                        shift 2
                        ;;
                --ptato)
                        PTATODIR="$2"
                        shift 2
                        ;;
                --error)
                        ERROR="$2"
                        shift 2
                        ;;
                --model)
                        MODEL="$2"
                        shift 1
                        ;;
                --inputtype)
                        INPUTTYPE="$2"
                        shift 2
                        ;;
                ---threads)
                        THREADS="$2"
                        shift 2
                        ;;
                --outputdir)
                        OUTPUTDIR=$(realpath "$2")
                        shift 2
                        ;;
                --seed)
                        SEED="$2"
                        shift 2
                        ;;
                --nbootstrap)
                        NBOOTSTRAP="$2"
                        shift 2
                        ;;
                --prefix)
                        PREFIX="$2"
                        shift 2
                        ;;
                --percentage)
                        PERCENT="$2"
                        shift 2
                        ;;
                --plottingonly)
                        PLOTONLY="$2"
                        shift 2
                        ;;
                --help)
                        usage
                        exit
                        ;;
                *)
                        echo "Command line parsing error ($1)"
                        echo "$@"
                        exit 3
                        ;;
        esac
done

# check if all files and directories exist ---------------
if [[ ! -f "$INPUT" ]] ; then
        echo "File $INPUT not found. Specify using the command line argument --input"
        exit
fi
if [[ ! -d "$OUTPUTDIR" ]] ; then
	mkdir -p $OUTPUTDIR
fi
if [[ $ERROR == "true" ]]; then
        MODEL=${MODEL}"+E"
fi
if [[ $OUTGROUP == "NONE" ]]; then
	echo "please define an outgroup using --outgr"
	#exit
fi

# Create prefix if no prefix is given
if [ -z "$PREFIX" ] ; then 
        PREFIX=${INPUT/%.vcf/}
        PREFIX=$(basename ${PREFIX})
fi
# Create percentage if no percentage is given
if [ -z "$PERCENT" ] ; then 
        PERCENT=0.4
fi
# Create percentage if no percentage is given
if [ -z "$PLOTONLY" ] ; then 
        PLOTONLY=FALSE
fi

# process the input variables ---------------
# build the command to build the tree: add "prob-msa off" when cellphy should be run on GT instead of the default PL
TREE_COMM="bash ${CELLPHY} RAXML --msa ${INPUT} --msa-format VCF --model ${MODEL} --seed ${SEED} --threads ${THREADS} --prefix ${PREFIX}.Tree"
if [[ $INPUTTYPE != "GL" ]]; then
        TREE_COMM=${TREE_COMM}" --prob-msa off"
fi


cd $OUTPUTDIR

if [ $PLOTONLY == FALSE ]; then
        # run cellphy ---------------
        # find the best tree
        eval "$TREE_COMM"
        # do bootstrapping to determine confidence
        bash ${CELLPHY} RAXML --bootstrap --msa ${INPUT} --model ${MODEL} --seed ${SEED} --threads ${THREADS} --bs-trees ${NBOOTSTRAP} --prefix ${PREFIX}.Boot
        # summarize bootstrap results
        bash ${CELLPHY} RAXML --support -tree ${PREFIX}.Tree.raxml.bestTree --bs-trees ${PREFIX}.Boot.raxml.bootstraps --prefix ${PREFIX}.Support --threads ${THREADS} --redo
        # map mutations to the tree
        bash ${CELLPHY} RAXML --mutmap --msa ${INPUT} --msa-format VCF -model ${PREFIX}.Tree.raxml.bestModel -tree ${PREFIX}.Tree.raxml.bestTree --prefix ${PREFIX}.Support --threads ${THREADS} --opt-branches off

        # Default cellphy tree  
        Rscript ${SUPPORT_MAP} ${PREFIX}.Support.raxml.support ${OUTGROUP}
fi



# make the tree in R ---------------
Rscript --vanilla ${SOURCE}/cellPhyPlotTree.R $OUTPUTDIR $INPUT $PTATODIR $OUTGROUP $PERCENT


# Create upset plot ---------------
Rscript --vanilla ${SOURCE}/Upsetplot.R $INPUT $OUTPUTDIR ${PREFIX} $OUTGROUP





