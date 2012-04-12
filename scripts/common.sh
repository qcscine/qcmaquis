#!/bin/bash

## targets ##

AMBIENT=ambient
DMRG=dmrg
TYPES=types

## settings ##

BUILD_DIR=${SITE}_${COMPILER}_${MPI_WRAPPER}
SELF=${BUILD_DIR}

add_ambient(){
    local defines_common=`get_defines MAQUIS_COMMON`
    local defines_target=`get_defines MAQUIS_AMBIENT`
    eval cmake $defines_common $defines_target ..
}

add_types(){
    local defines_common=`get_defines MAQUIS_COMMON`
    local defines_target=`get_defines MAQUIS_TYPES`
    use_dashboards types
    eval cmake $defines_common $defines_target ..
}

add_dmrg(){
    local defines_common=`get_defines MAQUIS_COMMON`
    local defines_target=`get_defines MAQUIS_DMRG`
    eval cmake $defines_common $defines_target ..
}

add_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: configuring"
    echo " ------------------------------------------------------------------------------------------ "
    mkdir ${ROOT_DIR}/${!target}/${BUILD_DIR}
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_DIR}
    eval add_${1}
    popd &> /dev/null
    set_state ${1} configure
}

remove_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " $self: cleaning "
    set_state ${1} void
    rm -rf ${ROOT_DIR}/${!target}/${BUILD_DIR}
}

build_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: building"
    echo " ------------------------------------------------------------------------------------------ "
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_DIR}
    make -j
    popd &> /dev/null
    set_state ${1} build
}

run_target(){
    local self="${SELF}/${1} (`get_state ${1}`)"
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo " ------------------------------------------------------------------------------------------ "
    echo " $self: testing"
    echo " ------------------------------------------------------------------------------------------ "
    pushd . &> /dev/null
    cd ${ROOT_DIR}/${!target}/${BUILD_DIR}
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
    cd ${ROOT_DIR}/${!target}/${BUILD_DIR}
    ctest -S Dashboards/site.cmake
    popd &> /dev/null
}

use_dashboards(){
    mkdir Dashboards
    local target="`echo ${1} | tr '[:lower:]' '[:upper:]'`"
    echo "set(PREDEFINED_CTEST_SITE \"${SITE}\")"                                         >  ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BUILD_NAME \"${COMPILER}_${MPI_WRAPPER}\")"                >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_SOURCE_DIRECTORY \"${ROOT_DIR}/${!target}\")"              >> ./Dashboards/site.cmake
    echo "set(PREDEFINED_CTEST_BINARY_DIRECTORY \"${ROOT_DIR}/${!target}/${BUILD_DIR}\")" >> ./Dashboards/site.cmake
    cat ../Dashboards/site.cmake                                                          >> ./Dashboards/site.cmake
    cp ../Dashboards/cmake_common.cmake ./Dashboards/
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

check_state(){ # 1: target
    if [ -n "${1}" ]
    then
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
    ORIG="`basename $0`"
    MOD=".`basename $0`.mod"
    cp $ORIG $MOD # keeping permissions
    echo "#!/bin/bash"                      >  $MOD
    echo "STATE=\"$STATE\""                 >> $MOD
    echo "STATE_AMBIENT=\"$STATE_AMBIENT\"" >> $MOD
    echo "STATE_TYPES=\"$STATE_TYPES\""     >> $MOD
    echo "STATE_DMRG=\"$STATE_DMRG\""       >> $MOD
    tail -n +6 $ORIG                        >> $MOD
    mv $MOD $ORIG
}

## target wrappers ##

clean(){
   local state=`get_state ${1}`
   if [ -n "${1}" ] 
   then
      remove_target ${1} 
   else
      echo " $SELF ($state): cleaning build tree             "
      echo " -------------------------------------------------------------------------------------------------------------------<< "
      remove_target ambient
      remove_target types
      remove_target dmrg
      set_state void
      echo " $SELF (void): build cleaning finished ...       "
   fi
}

configure(){
    clean ${1} # cleaning every configuration
    local state=`get_state ${1}`
    if [ -n "${1}" ]
    then
        add_target ${1}
    else
        echo " $SELF ($state): configuring build tree        "
        add_target ambient
        add_target types
        add_target dmrg
        set_state configure
        echo " $SELF (configure): build configure is done    "
    fi
}

build(){
    local state=`get_state ${1}`
    [[ "$state" == "void" ]] && configure ${1}
    if [ -n "${1}" ]
    then
        build_target ${1}
    else
        echo " $SELF ($state): building source tree          "
        build_target ambient
        build_target types
        build_target dmrg
        set_state build
        echo " $SELF (build): build has finished             "
    fi
}

run(){
    local state=`get_state ${1}`
    [[ "$state" != "build" ]] && build ${1}
    if [ -n "${1}" ]
    then
        run_target ${1}
    else
        echo " $SELF ($state): running all tests             "
        run_target ambient
        run_target types
        run_target dmrg
        echo " $SELF (build): testing has finished           "
    fi
}

dash(){
    local state=`get_state ${1}`
    [[ "$state" != "build" ]] && build ${1}
    if [ -n "${1}" ]
    then
        dash_target ${1}
    else
        echo " $SELF ($state): running all tests (dashboard) "
        dash_target ambient
        dash_target types
        dash_target dmrg
        echo " $SELF (build): dashboard refresh has finished "
    fi
}

execute(){
    echo
    if [ "$0" == "-bash" ]
    then
        echo " $SELF: set the environment"
    else
        action=`echo $1 | sed "s/dashboard/dash/" | sed "s/test/run/"` &&
        [[ "$action" != "clean" ]] && [[ "$action" != "configure" ]]   && 
        [[ "$action" != "build" ]] && [[ "$action" != "run"       ]]   &&
        [[ "$action" != "dash"  ]] && 
        echo "  Usage: ./config {clean, configure, build, test, dashboard} [targets]"  &&
        echo "  Note: in order to set the environment use \`source ./config\`" && echo && exit

        local i; for i in ${*:2}""; do
            check_state ${i} # safe check
            eval $action ${i}
        done
    fi
    echo
}
