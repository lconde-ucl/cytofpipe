#!/bin/bash -l


set -o pipefail


configfile="${CYTOFPIPE_HOME}/aux_files/clustering_config.txt"
gatingtemplate="${CYTOFPIPE_HOME}/aux_files/gating_template_transform.csv"
inputfiles="";
outputfiles="";
markersfile="";
ref="";
conditions="";
transform="-";
merge="-";
downsample="-";
asinh="-";
displayAll="-";
medians="-";


command=$1
arguments=$@

case "$command" in
  # Parse options to the install sub command
  --clustering)
    mode=$1; shift  # Remove '--clustering' from the argument list

    # Process mode options
    while getopts ":i:o:m:-:" opt; do
      case ${opt} in
	-)
            case "${OPTARG}" in
                displayAll)
                    displayAll="yes";
                    ;;
                flow)
                    transform="autoLgcl";
                    ;;
                cytof)
                    transform="arcsinh";
                    ;;
                all)
                    merge="all";
                    ;;
                downsample)
                    downsample="${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                config)
                    configfile="${PWD}/${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                *)
                    if [ "$OPTERR" = 1 ] && [ "${optspec:0:1}" != ":" ]; then
                        echo "Unknown option --${OPTARG}" >&2
                    fi
                    ;;
            esac;;
        i )
	  inputfiles=$OPTARG
          ;;
        o )
	  outputfiles=$OPTARG
          ;;
        m )
	  markersfile=$OPTARG
          ;;
        \? )
          echo "Invalid Option: -$OPTARG" 1>&2
          exit 1
          ;;
        : )
          echo "Invalid Option: -$OPTARG requires an argument" 1>&2
          exit 1
          ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
  --scaffold)
    mode=$1; shift  # Remove '--scaffold' from the argument list

    # Process mode options
    while getopts ":i:o:m:-:" opt; do
      case ${opt} in
	-)
            case "${OPTARG}" in
                flow)
                    asinh="150";
                    ;;
                cytof)
                    asinh="5";
                    ;;
                all)
                    merge="all";
                    ;;
                downsample)
                    downsample="${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                ref)
                    ref="${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                *)
                    if [ "$OPTERR" = 1 ] && [ "${optspec:0:1}" != ":" ]; then
                        echo "Unknown option --${OPTARG}" >&2
                    fi
                    ;;
            esac;;
        i )
	  inputfiles=$OPTARG
          ;;
        o )
	  outputfiles=$OPTARG
          ;;
        m )
	  markersfile=$OPTARG
          ;;
        \? )
          echo "Invalid Option: -$OPTARG" 1>&2
          exit 1
          ;;
        : )
          echo "Invalid Option: -$OPTARG requires an argument" 1>&2
          exit 1
          ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
  --citrus)
    mode=$1; shift  # Remove '--citrus' from the argument list

    # Process mode options
    while getopts ":i:o:m:-:" opt; do
      case ${opt} in
	-)
            case "${OPTARG}" in
                flow)
                    asinh="150";
                    ;;
                cytof)
                    asinh="5";
                    ;;
                all)
                    merge="all";
                    ;;
                downsample)
                    downsample="${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                cond)
                    conditions="${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                medians)
                    medians="${!OPTIND}";
		    OPTIND=$(( $OPTIND + 1 ))
                    ;;
                *)
                    if [ "$OPTERR" = 1 ] && [ "${optspec:0:1}" != ":" ]; then
                        echo "Unknown option --${OPTARG}" >&2
                    fi
                    ;;
            esac;;
        i )
	  inputfiles=$OPTARG
          ;;
        o )
	  outputfiles=$OPTARG
          ;;
        m )
	  markersfile=$OPTARG
          ;;
        \? )
          echo "Invalid Option: -$OPTARG" 1>&2
          exit 1
          ;;
        : )
          echo "Invalid Option: -$OPTARG requires an argument" 1>&2
          exit 1
          ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
esac


## Required by function die
trap "exit 1" TERM
export TOP_PID=$$
 

function die {
  echo "******************************************************"
  echo "["`date`"] ERROR while running $1"
  echo "******************************************************"
  kill -s TERM $TOP_PID
}
 

function run_clustering {

	JOB=${RAND_ID} R CMD BATCH --vanilla --no-save ${CYTOFPIPE_HOME}/code/cytofpipe_clustering.R  ${PWD}/${outputfiles}/log_R.txt
	
	for i in $PWD/$outputfiles/*_mean_*; do
	  if [ -f "$i" ]; 
		then 
		        rm ${PWD}/${outputfiles}/*cluster_mean_heatmap.*
		        rm ${PWD}/${outputfiles}/*cluster_mean_data.*
			break;
	 fi
	done
	
	for i in $PWD/$outputfiles/*Rphenograph*; do
	  if [ -f "$i" ]; 
		then 
			mkdir -p $PWD/$outputfiles/Rphenograph
			mv $PWD/$outputfiles/*Rphenograph*\.* $PWD/$outputfiles/Rphenograph/.
			break;
	 fi
	done
	
	for i in $PWD/$outputfiles/*FlowSOM*; do
	  if [ -f "$i" ]; 
		then 
			mkdir -p $PWD/$outputfiles/FlowSOM
			mv $PWD/$outputfiles/*FlowSOM*\.* $PWD/$outputfiles/FlowSOM/.
			break;
	 fi
	done
	
	for i in $PWD/$outputfiles/*DensVM*; do
	  if [ -f "$i" ]; 
		then 
			mkdir -p $PWD/$outputfiles/DensVM
			mv $PWD/$outputfiles/*DensVM*\.* $PWD/$outputfiles/DensVM/.
			break;
	 fi
	done
	
	for i in $PWD/$outputfiles/*ClusterX*; do
	  if [ -f "$i" ]; 
		then 		
			mkdir -p $PWD/$outputfiles/ClusterX
			mv $PWD/$outputfiles/*ClusterX*\.* $PWD/$outputfiles/ClusterX/.
			break;
	 fi
	done
	
	for i in $PWD/$outputfiles/gating*; do
	   if [ -f "$i" ]; 
	 	then 
	 		mkdir -p $PWD/$outputfiles/Gating
	 		mv $PWD/$outputfiles/gating* $PWD/$outputfiles/Gating/.
	 		break;
	  fi
	done

	for i in $PWD/$outputfiles/*_sample_level_plot*; do
	   if [ -f "$i" ]; 
	 	then 
	 		mkdir -p $PWD/$outputfiles/Marker_level_plots_by_sample
	 		mv $PWD/$outputfiles/*_sample_level_plot* $PWD/$outputfiles/Marker_level_plots_by_sample/.
	 		break;
	  fi
	done	
	
	cp ${CYTOFPIPE_HOME}/code/summary_clustering.Rmd ${PWD}/${outputfiles}/.
	R --vanilla  -e "rmarkdown::render('${PWD}/${outputfiles}/summary_clustering.Rmd',params=list(rscript='${CYTOFPIPE_HOME}/code/cytofpipe_clustering.R',rdata='${PWD}/${outputfiles}/cytofpipe.RData',inputparams='${FILE}'))"
	rm -rf ${PWD}/${outputfiles}/summary_clustering.Rmd
	
	for i in $PWD/Rplots*; do
	   if [ -f "$i" ]; 
	 	then 
	 		rm -rf $PWD/Rplots*
	 		break;
	  fi
	done
	
	if [ -e "$PWD/${RAND_ID}.txt" ]
	  then
	      rm $PWD/${RAND_ID}.txt
	fi
}

function run_scaffold {

	JOB=${RAND_ID} R CMD BATCH --vanilla --no-save ${CYTOFPIPE_HOME}/code/cytofpipe_scaffold.R  ${PWD}/${outputfiles}/log_R.txt
	
	if [ -d "$PWD/$inputfiles/downsampled_*" ]; 
	   then 
	 	mv $PWD/$inputfiles/downsampled_* $PWD/$outputfiles/.
	fi

	for i in $PWD/$outputfiles/downsampled_*/*clustered*; do
	   if [ -f "$i" ]; 
	 	then 
	 		mkdir -p $PWD/$outputfiles/clustering
	 		mv $PWD/$outputfiles/downsampled_*/*clustered* $PWD/$outputfiles/clustering/.
	 		break;
	  fi
	done

	cp ${CYTOFPIPE_HOME}/code/summary_scaffold.Rmd ${PWD}/${outputfiles}/.
	R --vanilla  -e "rmarkdown::render('${PWD}/${outputfiles}/summary_scaffold.Rmd',params=list(rscript='${CYTOFPIPE_HOME}/code/cytofpipe_scaffold.R',rdata='${PWD}/${outputfiles}/cytofpipe.scaffold',inputparams='${FILE}'))"
	rm -rf ${PWD}/${outputfiles}/summary_scaffold.Rmd
	

	if [ -e "$PWD/${RAND_ID}.txt" ]
	  then
	      rm $PWD/${RAND_ID}.txt
	fi
}

function run_citrus {

	JOB=${RAND_ID} R CMD BATCH --vanilla --no-save ${CYTOFPIPE_HOME}/code/cytofpipe_citrus.R  ${PWD}/${outputfiles}/log_R.txt
	
	cp ${CYTOFPIPE_HOME}/code/summary_citrus.Rmd ${PWD}/${outputfiles}/.
	R --vanilla  -e "rmarkdown::render('${PWD}/${outputfiles}/summary_citrus.Rmd',params=list(rscript='${CYTOFPIPE_HOME}/code/cytofpipe_citrus.R',inputparams='${FILE}'))"
	rm -rf ${PWD}/${outputfiles}/summary_citrus.Rmd


	if [ -e "$PWD/${RAND_ID}.txt" ]
	  then
	      rm $PWD/${RAND_ID}.txt
	fi
}



echo "======================================================"
echo " HOST: "`hostname`
echo " TIME: "`date`
echo " ARGS: ${arguments[*]}"
echo "======================================================"


FILE=${PWD}/${RAND_ID}.txt

if [ $command == "--clustering" ]
then
	
/bin/cat <<EOM >$FILE
[ paramsclustering ]
INPUTFILE = ${PWD}/${inputfiles}
OUTPUTFILE = ${PWD}/${outputfiles}
MARKERSFILE = ${PWD}/${markersfile}
CONFIGFILE = $configfile
GATINGFILE = $gatingtemplate
TRANSFORM = $transform
MERGE = $merge
DOWNSAMPLE = $downsample
DISPLAY_ALL = $displayAll
ARGS = ${arguments[*]}
EOM

	mkdir -p $PWD/$outputfiles	
	run_clustering;

elif [ $command == "--scaffold" ]
then

/bin/cat <<EOM >$FILE
[ paramsscaffold ]
INPUTFILE = ${PWD}/${inputfiles}
OUTPUTFILE = ${PWD}/${outputfiles}
MARKERSFILE = ${PWD}/${markersfile}
REF = ${ref}
ASINH = $asinh
MERGE = $merge
DOWNSAMPLE = $downsample
ARGS = ${arguments[*]}
EOM

	mkdir -p $PWD/$outputfiles	
	run_scaffold;

elif [ $command == "--citrus" ]
then

	
/bin/cat <<EOM >$FILE
[ paramscitrus ]
INPUTFILE = ${PWD}/${inputfiles}
OUTPUTFILE = ${PWD}/${outputfiles}
MARKERSFILE = ${PWD}/${markersfile}
CONDITIONS = ${PWD}/${conditions}
ASINH = $asinh
MERGE = $merge
DOWNSAMPLE = $downsample
MEDIANS = ${PWD}/${medians}
ARGS = ${arguments[*]}
EOM

	mkdir -p $PWD/${outputfiles}	
	run_citrus;

else
	echo "This shoud not be happening: command $command";
	die $*;
fi

echo
echo "======================================================"
echo "["`date`"] Done."
echo "======================================================"

