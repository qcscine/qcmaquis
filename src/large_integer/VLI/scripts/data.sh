# !/bin/bash
#path
PATH_FILE=${PWD}/../build/Testing/Temporary
PATH_PLOT=${PWD}
#create files
touch $PATH_PLOT/out.data
#command to get the good value GS or TE
COMMAND="(sed -n "\${k}p" out.data.inter)"

#select the data from the log files
for i in 2 3 4
do
    for j in 11 12 13 14 15 
    do
      cat $PATH_FILE/LastTest.log | grep 'InnerProduct_VLI_SIZE_'$i'_POLY_ORDER_'$j | tee -a "out.data"
    done
done
     cut -d' ' -f2- $PWD/out.data | tee -a out.data.inter 

k=1

for i in 2 3 4
do
    for j in 11 12 13 14 15 
    do
        VALUE=`eval $COMMAND`
        k=$(($k+1))
        echo $i $j $VALUE | tee -a data.gnuplot.$i
    done
done

rm out.data
rm out.data.inter 
