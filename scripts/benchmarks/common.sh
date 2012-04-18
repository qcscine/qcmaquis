#!/bin/bash

## locations ##
SCRIPT="$BENCHMARK_SCRIPTS_DIR/common.sh"
BENCHMARK_DIR=$ROOT_DIR/$TARGET/$BUILD_NAME/regression/performance
SCRATCH_DIR="$BENCHMARK_DIR/$BENCHMARK.`date +\"%m.%d.%H.%M\"`"
OUTPUT_LOG="output.log"
REPRODUCER="reproducer.sh"

## execution ##
COMMAND="( mpiexec -np \$NP \$EXECUTABLE \$NT \$WL | grep GFlops | tr -d 'GFlops ' ) 2> /dev/null"
[[ -n "$EXECUTABLE" ]] || EXECUTABLE="./$BENCHMARK.out"

## Handle SIGINT (^C) ##
## Terminate all runs ##
control_c(){ echo; echo; echo -n " "; kill $$; }
trap control_c SIGINT

## prepare the environment ##
pushd . &> /dev/null
mkdir $SCRATCH_DIR
touch $SCRATCH_DIR/$OUTPUT_LOG
touch $SCRATCH_DIR/$REPRODUCER
cp $BENCHMARK_DIR/$EXECUTABLE $SCRATCH_DIR
cd $SCRATCH_DIR

## generating reproducer ##
head -n 1 $SCRIPT > $REPRODUCER
echo "## auto-environment ##" >> $REPRODUCER
echo "## `uname -a` ##" >> $REPRODUCER
set >> $REPRODUCER
echo "## end of auto-environment ##" >> $REPRODUCER
tail +2 $SCRIPT >> $REPRODUCER

## main execution ## 
echo | tee $OUTPUT_LOG
echo " `date +\"%m.%d, %H:%M:%S\"`" | tee -a $OUTPUT_LOG

## loop over process counts ##
NP_STEP_BOUND=$NP_LOWER_BOUND
for NP in `seq $NP_LOWER_BOUND $NP_S_INCREMENT $NP_UPPER_BOUND`
do
    [[ $NP < $NP_STEP_BOUND ]] && continue
    NP_STEP_BOUND="$(( $NP*$NP_L_INCREMENT ))"
    echo | tee -a $OUTPUT_LOG
    echo " ------------------------------------------------------------------------------------------ " | tee -a $OUTPUT_LOG
    echo " $NAME: processes: $NP                                                                      " | tee -a $OUTPUT_LOG
    echo -n " ----" | tee -a $OUTPUT_LOG

    ## loop over thread counts ##
    NT_STEP_BOUND=$NT_LOWER_BOUND
    for NT in `seq $NT_LOWER_BOUND $NT_S_INCREMENT $NT_UPPER_BOUND`
    do
        [[ $NT < $NT_STEP_BOUND ]] && continue
        NT_STEP_BOUND="$(( $NT*$NT_L_INCREMENT ))"
        echo      "-------------------------------------------------------------------------------------- " | tee -a $OUTPUT_LOG
        echo "     $NAME: threads: $NT                                                                    " | tee -a $OUTPUT_LOG
        echo "     -------------------------------------------------------------------------------------- " | tee -a $OUTPUT_LOG

        ## loop over workloads ##
        WL_STEP_BOUND=$WL_LOWER_BOUND
        for WL in `seq $WL_LOWER_BOUND $WL_S_INCREMENT $WL_UPPER_BOUND`
        do
            [[ $WL < $WL_STEP_BOUND ]] && continue
            WL_STEP_BOUND="$(( $WL*$WL_L_INCREMENT ))"
        
            echo -n "     $NAME (${NP}x${NT} on ${WL}): " | tee -a $OUTPUT_LOG
            eval $COMMAND | tee -a $OUTPUT_LOG
        done
        echo -n "     " | tee -a $OUTPUT_LOG
    done
done

## exiting / cleaning ##
echo | tee -a $OUTPUT_LOG
echo " `date +\"%m.%d, %H:%M:%S\"`" | tee -a $OUTPUT_LOG
echo | tee -a $OUTPUT_LOG
popd &> /dev/null
