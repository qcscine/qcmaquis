#!/bin/bash -l
## externals ##
## HOST
## ROOT_DIR
## CONFIG_DIR 
## CONFIG_NAME (if specific)

## targets ##
AMBIENT=ambient
DMRG=dmrg
TARGETS="ambient dmrg"

RUN_PRESETS="1n.micro:             ~/maquis2013/benchmarks/dmrg_gs/1n/micro/parms              ~/maquis2013/benchmarks/dmrg_gs/1n/micro/model
             1n.short:             ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L6/parms        ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L6/model
             1n.spinless.L8.6000:  ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L8/parms.6000   ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L8/model
             1n.spinless.L8.8000:  ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L8/parms.8000   ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L8/model
             1n.spinless.L8.10000: ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L8/parms.10000  ~/maquis2013/benchmarks/dmrg_gs/1n/spinless.L8/model
             1n.spinfull.L8.6000:  ~/maquis2013/benchmarks/dmrg_gs/1n/spinfull.L8/parms.6000   ~/maquis2013/benchmarks/dmrg_gs/1n/spinfull.L8/model
             1n.spinfull.L8.8000:  ~/maquis2013/benchmarks/dmrg_gs/1n/spinfull.L8/parms.8000   ~/maquis2013/benchmarks/dmrg_gs/1n/spinfull.L8/model
             1n.spinfull.L8.10000: ~/maquis2013/benchmarks/dmrg_gs/1n/spinfull.L8/parms.10000  ~/maquis2013/benchmarks/dmrg_gs/1n/spinfull.L8/model
             fermiwideladder: ~/maquis2013/benchmarks/dmrg_gs/fermiwideladder/parms  ~/maquis2013/benchmarks/dmrg_gs/fermiwideladder/model"

## settings ##

BUILD_NAME=${PREFIX}_${COMPILER}_${MPI_WRAPPER}
[[ $SHARED_FS ]] && BUILD_NAME=${HOST}/${BUILD_NAME}
BENCHMARK_SCRIPTS_DIR=${ROOT_DIR}/scripts/benchmarks
SELF=$BUILD_NAME

add_ambient(){
    local defines_common=`get_defines COMMON`
    local defines_target=`get_defines AMBIENT`
    eval cmake $defines_common $defines_target ${ROOT_DIR}/ambient
}

add_dmrg(){
    local defines_common=`get_defines COMMON`
    local defines_target=`get_defines DMRG`
    eval cmake $defines_common $defines_target ${ROOT_DIR}/dmrg
}

add_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: configuring"
    echo " ------------------------------------------------------------------------------------------ "
    [[ $SHARED_FS ]] && mkdir ${ROOT_DIR}/${!target}/${HOST} &> /dev/null
    pushd . &> /dev/null; 
    mkdir ${ROOT_DIR}/${!target}/${BUILD_NAME} &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_NAME}
    eval add_${1}
    popd &> /dev/null
    set_state ${1} configure
}

remove_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " $self: cleaning "
    set_state ${1} void
    rm -rf ${ROOT_DIR}/${!target}/${BUILD_NAME}
}

build_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: building"
    echo " ------------------------------------------------------------------------------------------ "
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_NAME}
    make -j
    popd &> /dev/null
    set_state ${1} build
}

test_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: testing"
    echo " ------------------------------------------------------------------------------------------ "
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_NAME}
    make test
    popd &> /dev/null
}

run_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: testing"
    echo " ------------------------------------------------------------------------------------------ "
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_NAME}
    make test
    popd &> /dev/null
}

dash_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: testing (dashboard)"
    echo " ------------------------------------------------------------------------------------------ "
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_NAME}
    ctest -S Dashboards/site.cmake
    popd &> /dev/null
}

use_dashboards(){
    mkdir Dashboards &> /dev/null
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo "set(PREDEFINED_CTEST_SITE \"${HOST}\")"                                          >  ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BUILD_NAME \"${PREFIX}_${COMPILER}_${MPI_WRAPPER}\")"       >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_SOURCE_DIRECTORY \"${ROOT_DIR}/${!target}\")"               >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BINARY_DIRECTORY \"${ROOT_DIR}/${!target}/${BUILD_NAME}\")" >> ./Dashboards/site.cmake
    cat ${ROOT_DIR}/scripts/common/ctest/site.cmake                                        >> ./Dashboards/site.cmake
    cp  ${ROOT_DIR}/scripts/common/ctest/cmake_common.cmake                                   ./Dashboards/
}

## auxiliary functions ##

die(){
    echo " Error: ${1}" 2>&1;
    kill -SIGINT $$
}

get_defines(){
    local input=`eval echo \\${!\${1}_*}`
    local i; for i in $input; do
        local d=`echo $i | sed "s/${1}_/D/"`
        echo -n " -$d=\"${!i}\""
    done
}

lookup(){
    local command=${1}
    local target=${command%% *}
    if [ ! -f $target ]
    then
        target=`basename $target`
        local i; for i in $TARGETS; do
            if [ -d ${ROOT_DIR}/${i}/${BUILD_NAME} ]; then
                local result=`find ${ROOT_DIR}/${i}/${BUILD_NAME} -name $target -type f -executable -print`
                [[ -n $result ]] && value=$result
            fi
        done
    else
        value=`readlink -f $target`
    fi
    echo "$value"
}

check_state(){ # 1: target
    if [ -n "${1}" ]
    then
        [[ -f "benchmarks/${1}" ]] && return
        target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
        value="`eval echo \$\`eval \"echo STATE_${target}\"\``"
        [[ -n "$value" ]] || die "unknown target $1"
    fi
}

get_state(){ # 1: target
    local value=$STATE
    local target="${1}"
    if [ -n "${1}" ]
    then
        target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
        value="`eval echo \$\`eval \"echo STATE_${target}\"\``"
        [[ -n "$value" ]] || die "unknown target"
    fi
    echo "$value"
}

set_state(){ # 1: target (optional, def: all) # 2: state
    if [ -n "${2}" ]
    then
        local state=`get_state ${1}` # check if target exists
        local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
        eval "$(echo STATE_${target})=${2}"
    else
        STATE="${1}"
    fi
    write_states
}

write_states(){
    [[ -z "${CONFIG_NAME}" ]] && export CONFIG_NAME=`basename $0`
    ORIG="${CONFIG_DIR}/${CONFIG_NAME}"
    MOD="${CONFIG_DIR}/${CONFIG_NAME}.mod"
    cp $ORIG $MOD # keeping permissions
    echo "#!/bin/bash -l"                   >  $MOD
    echo "STATE=\"$STATE\""                 >> $MOD
    echo "STATE_AMBIENT=\"$STATE_AMBIENT\"" >> $MOD
    echo "STATE_DMRG=\"$STATE_DMRG\""       >> $MOD
    tail -n +5 $ORIG                        >> $MOD
    mv $MOD $ORIG
}

## target wrappers ##

function clean(){
   local state=`get_state ${1}`
   if [ -n "${1}" ] 
   then
      remove_target ${1} 
   else
      echo " $SELF ($state): cleaning build tree             "
      echo " -------------------------------------------------------------------------------------------------------------------<< "
      remove_target ambient
      remove_target dmrg
      set_state void
      echo " $SELF (void): build cleaning finished ...       "
   fi
}

function configure(){
    clean ${1} # cleaning every configuration
    local state=`get_state ${1}`
    if [ -n "${1}" ]
    then
        add_target ${1}
    else
        echo " $SELF ($state): configuring build tree        "
        add_target ambient
        add_target dmrg
        set_state configure
        echo " $SELF (configure): build configure is done    "
    fi
}

function build(){
    local state=`get_state ${1}`
    [[ "$state" == "void" ]] && configure ${1}
    if [ -n "${1}" ]
    then
        build_target ${1}
    else
        echo " $SELF ($state): building source tree          "
        build_target ambient
        build_target dmrg
        set_state build
        echo " $SELF (build): build has finished             "
    fi
}

function test(){
    local state=`get_state ${1}`
    [[ "$state" != "build" ]] && build ${1}
    if [ -n "${1}" ]
    then
        test_target ${1}
    else
        echo " $SELF ($state): running all tests             "
        test_target ambient
        test_target dmrg
        echo " $SELF (build): testing has finished           "
    fi
}

function expend(){
    echo "short"
}

function run(){  
    [[ ! -n "${1}" ]] && die "please specify the binary      "
    args=${*:2}
    if [ -z ${3} ]; then
        echo $RUN_PRESETS | grep "${2}:" &> /dev/null
        [[ $? -eq 0 ]] && args=`(read line; echo $line;) < <(echo "${RUN_PRESETS#*${2}: }")`
    fi

    local target=`lookup ${1}`
    [[ ! -n $target ]] && die "couldn't find the binary      "

    local socketlist="0 1"
    local proclist="0,1,2,3,4,5,6,7,8,9,10,11"
    local rank_var="OMPI_COMM_WORLD_NODE_RANK"

    [[ ! -z $VALGRIND         ]] && VALGRIND="valgrind --error-limit=no"
    [[ ! -z $CILK_NUM_THREADS ]] && CILK_NUM_THREADS="CILK_NWORKERS=$CILK_NUM_THREADS"
    [[ ! -z $MPI_NUM_PROCS    ]] && HWLOC="hwloc-bind socket:\${socketlist[\$$rank_var]}" # --mempolicy firsttouch 
    if [ ! -z $OMP_NUM_THREADS  ]; then
        [[ $OMP_NUM_THREADS -eq 1  ]] && proclist="0 6 1 7 2 8 3 9 4 10 5 11"
        [[ $OMP_NUM_THREADS -eq 2  ]] && proclist="0,6 1,7 2,3 8,9 4,5 10,11"
        [[ $OMP_NUM_THREADS -eq 3  ]] && proclist="0,1,2 6,7,8 3,4,5 9,10,11"
        [[ $OMP_NUM_THREADS -eq 4  ]] && proclist="0,1,2,3 8,9,10,11 4,5,6,7"
        [[ $OMP_NUM_THREADS -eq 6  ]] && proclist="0,1,2,3,4,5 6,7,8,9,10,11"
        OMP_NUM_THREADS="KMP_AFFINITY=verbose,proclist=[\${proclist[\$$rank_var]}],explicit OMP_NUM_THREADS=$OMP_NUM_THREADS"
    fi
    
    local command="$OMP_NUM_THREADS 
                   $CILK_NUM_THREADS
                   $HWLOC
                   $VALGRIND 
                   $target $args"; 

    if [ ! -z "$MPI_NUM_PROCS" ]; then
        command="export proclist=($proclist); 
                 export socketlist=($socketlist); 
                 command=\"$command\"; 
                 echo \$command;
                 echo
                 eval \$command"
        
        echo "#!/bin/bash
        $command; rm -f bootstrap.sh" &> bootstrap.sh; chmod +x bootstrap.sh;
        command="mpiexec -np $MPI_NUM_PROCS bootstrap.sh"
    fi

    rm -f *.h5*
    echo $command
    echo
    eval $command
    rm -f *.h5*
}

function dash(){
    local state=`get_state ${1}`
    [[ "$state" != "build" ]] && build ${1}
    if [ -n "${1}" ]
    then
        dash_target ${1}
    else
        echo " $SELF ($state): running all tests (dashboard) "
        dash_target ambient
        dash_target dmrg
        echo " $SELF (build): dashboard refresh has finished "
    fi
}

function benchmark(){
    [[ -n "${1}" ]] || die "please supply the name of the benchmark"
    source $BENCHMARK_SCRIPTS_DIR/${1}
    local state=`get_state ${TARGET}`
    [[ "$state" != "build" ]] && build ${TARGET}
    pushd . &> /dev/null
    source $BENCHMARK_SCRIPTS_DIR/common.sh
    popd &> /dev/null
}

function execute(){
    echo
    if [ "$1" == "" ]
    then
        echo " $SELF: set the environment"
    else
        action=`echo $1 | sed "s/dashboard/dash/"` &&
        [[ "$action" != "clean" ]] && [[ "$action" != "configure" ]] && 
        [[ "$action" != "build" ]] && [[ "$action" != "test"      ]] &&
        [[ "$action" != "dash"  ]] && [[ "$action" != "benchmark" ]] &&
        [[ "$action" != "run"   ]] && 
        echo "  Usage: ./config {clean, configure, build, test, run, dashboard, benchmark} [targets]" &&
        echo "  Note: in order to set the environment use \`source ./config\`" && echo && exit

        if [ "$action" == "run" ]; then
            run ${*:2}
        else
            local i; for i in ${*:2}""; do
                check_state ${i} # safe check
                eval $action ${i}
            done
        fi
    fi
    echo
}
