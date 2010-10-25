# global constants for all modules
# weak checks on values (i.e.no check on reasonable size)

# number of good (additive) QUANTUM_NUMBERS

QUANTUM_NUMBERS = 1

if not isinstance(QUANTUM_NUMBERS, int):
    raise ConstantError("QUANTUM NUMBERS not integer")
if QUANTUM_NUMBERS < 0:
    raise ConstantError("QUANTUM NUMBERS must be non-negative integer")

# maximum allowed number of bosons per site (positive integer)
# local physical dimension is MAX_BOSONS_PER_SITE+1

MAX_BOSONS_PER_SITE = 4

if not isinstance(MAX_BOSONS_PER_SITE,int):
    raise ConstantError("MAX_BOSONS_PER_SITE not integer")
if MAX_BOSONS_PER_SITE <= 0:
    raise ConstantError("MAX_BOSONS_PER_SITE not positive")

# length of spin: allowed values are 0.5,1.0,1.5,2.0,...

SPIN_LENGTH = 0.5

if not isinstance(SPIN_LENGTH, float):
    raise ConstantError("SPIN_LENGTH not float")
if (SPIN_LENGTH <= 0.0) or (int(2 * SPIN_LENGTH) != 2 * SPIN_LENGTH) :
    raise ConstantError("SPIN_LENGTH must be 0.5,1.0,1.5,...")

# lower cutoff for retained singular values to ensure stability
# default value in NumPy.pinv (Moore-Penrose pseudoinverse)
# is 1.0e-15
 
SV_TOL = 1.0e-08

if not isinstance(SV_TOL, float):
    raise ConstantError("SV_TOL not float")
if (SV_TOL <= 0.0) or (SV_TOL >= 1.0):
    raise ConstantError("SV_TOL outside meaningful range")

# tolerance for norm deviations

NORM_TOL = 1.0e-08

if not isinstance(NORM_TOL, float):
    raise ConstantError("NORM_TOL not float")
if (NORM_TOL <= 0.0) or (NORM_TOL >= 1.0):
    raise ConstantError("NORM_TOL outside meaningful range")

# maximal bond dimension to catch runaways

BOND_MAX = 1024

if not isinstance(BOND_MAX, int):
    raise ConstantError("BOND_MAX not integer")
if BOND_MAX<=0:
    raise ConstantError("BOND_MAX not positive")



    
