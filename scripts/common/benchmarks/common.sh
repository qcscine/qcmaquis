#!/bin/bash -l

## locations ##
SCRIPT="$BENCHMARK_SCRIPTS_DIR/common.sh"
BENCHMARK_DIR=$ROOT_DIR/$TARGET/$BUILD_NAME/regression/performance
SCRATCH_DIR="$BENCHMARK_DIR/$BENCHMARK.`date +\"%m.%d.%H.%M\"`"
OUTPUT_LOG="output.log"
PART_PLOT_OUTPUT_LOG="temp.data"
FULL_PLOT_OUTPUT_LOG="plot.data"
FULL_PLOT_OUTPUT_LOG_TMP="plot.data.tmp"
REPRODUCER="reproducer.sh"

## execution ##
COMMAND="( mpiexec -np \$NP \$EXECUTABLE \$NT \$WL | grep GFlops | tr -d 'GFlops ' ) 2> /dev/null"
[[ -n "$EXECUTABLE" ]] || EXECUTABLE="./$BENCHMARK.out"

## Handle SIGINT (^C) ##
## Terminate all runs ##
control_c(){ echo; echo; echo -n " "; kill $$; }
trap control_c SIGINT


## generating partial chart ##
## $1: title, $2: input, $3: output
part_plot(){
gnuplot << EOF
set terminal postscript eps color enhanced
set output "$3.eps"
set xlabel "size"
set ylabel "GFlops"
set title "$1"
plot "$2" using 1:2 notitle w l
EOF
epstopdf $3.eps; rm $3.eps
}
## generating full chart ##
## $1: title, $2: input, $3: legend, $4: output
full_plot(){
gnuplot << EOF
set terminal postscript eps color 
set output "$4.eps"
set xlabel "size"
set ylabel "GFlops"
set title "$1"

#red family
set style line 1  lt 1 lw 3 pt 1 lc rgb "#FF0000" #red
set style line 2  lt 1 lw 3 pt 1 lc rgb "#FF1493" #deep pink
set style line 3  lt 1 lw 3 pt 1 lc rgb "#FF6347" #tomato
set style line 4  lt 1 lw 3 pt 1 lc rgb "#FFC8CB" #pink
                  
#blue family      
set style line 5  lt 1 lw 3 pt 1 lc rgb "#0000FF" #blue
set style line 6  lt 1 lw 3 pt 1 lc rgb "#1E90FF" #dodger blue
set style line 7  lt 1 lw 3 pt 1 lc rgb "#4682B4" #steelblue
set style line 8  lt 1 lw 3 pt 1 lc rgb "#E6E6FA" #lavender

#green family
set style line 9  lt 1 lw 3 pt 1 lc rgb "#008000" #green
set style line 10 lt 1 lw 3 pt 1 lc rgb "#3CB371" #medium seagreen
set style line 11 lt 1 lw 3 pt 1 lc rgb "#7FFFd4" #Aquamarine
set style line 12 lt 1 lw 3 pt 1 lc rgb "#7FFF00" #chartreuse:w


`N=2; IFS=';' read -ra LABELS <<< "$3"; for i in "${LABELS[@]}"; do echo "set out \"$4.eps\";"; [[ $N == 2 ]] && echo -n "plot" || echo -n "replot"; echo -n " \"$2\" using 1:$N t \"$i\" w l ls $N;"; ((N++)); done`
set output
EOF
epstopdf $4.eps; rm $4.eps
mv $2 $4.data
}

## $1: line, $2: value ##
plot_axis_append(){
    [[ `cat $FULL_PLOT_OUTPUT_LOG | wc -l` < $1 ]] && echo "$2 " >> $FULL_PLOT_OUTPUT_LOG
    echo -n "$2 " >> $PART_PLOT_OUTPUT_LOG
}
## $1: line, $2: value ##
plot_data_append(){
    echo $2 >> $PART_PLOT_OUTPUT_LOG
    awk -vn="$1" -vsp=" $2" 'NR==n{$0=$0 sp}1' $FULL_PLOT_OUTPUT_LOG > $FULL_PLOT_OUTPUT_LOG_TMP
    mv $FULL_PLOT_OUTPUT_LOG_TMP $FULL_PLOT_OUTPUT_LOG
}


## prepare the environment ##
pushd . &> /dev/null
mkdir $SCRATCH_DIR
touch $SCRATCH_DIR/$OUTPUT_LOG
touch $SCRATCH_DIR/$REPRODUCER
cp $BENCHMARK_DIR/$EXECUTABLE $SCRATCH_DIR
cd $SCRATCH_DIR
touch $FULL_PLOT_OUTPUT_LOG

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
    LEGEND=""
    for NT in `seq $NT_LOWER_BOUND $NT_S_INCREMENT $NT_UPPER_BOUND`
    do
        [[ $NT < $NT_STEP_BOUND ]] && continue
        NT_STEP_BOUND="$(( $NT*$NT_L_INCREMENT ))"
        echo      "-------------------------------------------------------------------------------------- " | tee -a $OUTPUT_LOG
        echo "     $NAME: threads: $NT                                                                    " | tee -a $OUTPUT_LOG
        echo "     -------------------------------------------------------------------------------------- " | tee -a $OUTPUT_LOG

        echo > $PART_PLOT_OUTPUT_LOG
        ## loop over workloads ##
        WL_STEP_BOUND=$WL_LOWER_BOUND
        PLOT_OUTPUT_LINE=1
        for WL in `seq $WL_LOWER_BOUND $WL_S_INCREMENT $WL_UPPER_BOUND`
        do
            [[ $WL < $WL_STEP_BOUND ]] && continue
            WL_STEP_BOUND="$(( $WL*$WL_L_INCREMENT ))"
            VALUE="`eval $COMMAND`"
            echo "     $NAME (${NP}x${NT} on ${WL}): $VALUE" | tee -a $OUTPUT_LOG

            plot_axis_append $PLOT_OUTPUT_LINE $WL
            plot_data_append $PLOT_OUTPUT_LINE $VALUE
            ((PLOT_OUTPUT_LINE++))
        done
        echo -n "     " | tee -a $OUTPUT_LOG
        part_plot "$TITLE on $NP processes x $NT threads" $PART_PLOT_OUTPUT_LOG "plot.${NP}x${NT}"
        LEGEND="${LEGEND}thread count: ${NT};"
    done
    full_plot "$TITLE on $NP proc" $FULL_PLOT_OUTPUT_LOG "$LEGEND" "plot.${NP}"
done

## exiting / cleaning ##
echo | tee -a $OUTPUT_LOG
echo " `date +\"%m.%d, %H:%M:%S\"`" | tee -a $OUTPUT_LOG
echo | tee -a $OUTPUT_LOG
rm $PART_PLOT_OUTPUT_LOG

popd &> /dev/null
