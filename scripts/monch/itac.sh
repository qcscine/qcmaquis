# INTEL_LICENSE_FILE:
if [ -n "/opt/intel/licenses:/root/intel/licenses" ]; then
   if [ -n "${INTEL_LICENSE_FILE}" ]; then
       INTEL_LICENSE_FILE="/opt/intel/licenses:/root/intel/licenses:${INTEL_LICENSE_FILE}"
   else
       INTEL_LICENSE_FILE="/opt/intel/licenses:/root/intel/licenses"
   fi
   export INTEL_LICENSE_FILE
fi

# VT_MPI:
VT_MPI=$1
if [ -z "${VT_MPI}" ]; then
   VT_MPI=impi4
fi

# VT_ROOT:
VT_ROOT="/apps/monch/intel//itac/8.1.3.037"
export VT_ROOT

# VT_ARCH:
VT_ARCH=intel64
export VT_ARCH

# Temporary internal variable VT_ROOT_ARCH:
if [ -z "${VT_ARCH}" ]; then
   VT_ROOT_ARCH="${VT_ROOT}"
else
   VT_ROOT_ARCH="${VT_ROOT}/${VT_ARCH}"
fi

# Check:
if [ ! -d "${VT_ROOT_ARCH}/itac" ]; then
   # if only Intel(R) Trace Analyzer is installed:
   if [ -d "${VT_ROOT_ARCH}/bin" ]; then
        # PATH:
       if [ -n "${PATH}" ]; then
           PATH="${VT_ROOT_ARCH}/bin:${PATH}"
       else
           PATH="${VT_ROOT_ARCH}/bin"
       fi
       export PATH
       unset VT_ROOT VT_ARCH VT_ROOT_ARCH
       return 0
   fi
   echo "*** itacvars.sh: invalid VT_ROOT '${VT_ROOT}'" >&2
   unset VT_ROOT VT_ARCH VT_ROOT_ARCH
   return 1
fi

# VT_LIB_DIR, VT_SLIB_DIR, VT_ADD_LIBS:
case "${VT_MPI}" in
   mpich2|impi3|impi4)
          VT_MPI="impi4"
          VT_LIB_DIR="${VT_ROOT_ARCH}/itac/lib_impi4"
          VT_SLIB_DIR="${VT_ROOT_ARCH}/itac/slib_impi4"
          VT_MIC_SLIB_DIR="${VT_ROOT}/mic/itac/slib_impi4"
          VT_ADD_LIBS="-ldwarf -lelf -lvtunwind -lnsl -lm -ldl -lpthread"
          ;;
   *) echo "*** itacvars.sh: '${VT_MPI}' is not a supported MPI" >&2
      unset VT_ROOT VT_ARCH VT_ROOT_ARCH VT_MPI VT_LIB_DIR VT_SLIB_DIR VT_ADD_LIBS
      return 1
esac
export VT_MPI VT_LIB_DIR VT_SLIB_DIR VT_ADD_LIBS

# Check:
if [ ! -d "${VT_LIB_DIR}" ] || [ ! -d "${VT_SLIB_DIR}" ]; then # || [ ! -d "${VT_MIC_SLIB_DIR}" ]; then
   echo "*** itacvars.sh: '${VT_MPI}' is not installed" >&2
   unset VT_ROOT VT_ARCH VT_ROOT_ARCH VT_MPI VT_LIB_DIR VT_SLIB_DIR VT_ADD_LIBS
   return 1
fi

# LD_LIBRARY_PATH:
if [ -n "${LD_LIBRARY_PATH}" ]; then
   LD_LIBRARY_PATH="${VT_SLIB_DIR}:${LD_LIBRARY_PATH}"
else
   LD_LIBRARY_PATH="${VT_SLIB_DIR}"
fi
#LD_LIBRARY_PATH="${VT_MIC_SLIB_DIR}:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH

# PATH:
if [ -n "${PATH}" ]; then
   PATH="${VT_ROOT_ARCH}/bin:${PATH}"
else
   PATH="${VT_ROOT_ARCH}/bin"
fi
export PATH

# MANPATH:
if [ -n "${MANPATH}" ]; then
   MANPATH="${VT_ROOT}/man:${MANPATH}"
else
   MANPATH="${VT_ROOT}/man":$(manpath)
fi
export MANPATH

# CLASSPATH:
if [ -n "${CLASSPATH}" ]; then
   CLASSPATH="${VT_LIB_DIR}:${CLASSPATH}"
else
   CLASSPATH="${VT_LIB_DIR}"
fi
export CLASSPATH

# Unset temporary variables
unset VT_ROOT_ARCH
