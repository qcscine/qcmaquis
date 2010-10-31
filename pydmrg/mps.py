from numpy import *
from constants import *
from myexceptions import *
from numpy.linalg import norm, qr
from numpy.random import rand
from numpy import tensordot, compress
from numpy.dual import svd
   
# ######################################################################################    

# class for local MPS tensors M^{s}_{a,b}

class MPSTensor:
    """
    local MPS tensor
    ----------------
    
    structure:    
    self.tensor: complex 3-index tensor; self.tensor[s,a,b]=M^{s}_{a,b} Eq. (75)
    self.normalization: normalization of tensor; 'U', 'L', 'R' 
        'U': unnormalized or status unknown
        'L': left-normalized \sum_s M^{s\dagger}M^{s}=I -> A in formulae; Eqs. (38), (70)
        'R': right-normalized \sum_s M^{s}M^{s\dagger}=I -> B in formulae; Eqs. (45), (71)
    
    arguments of constructor:
    site_dim: dimension of s, int>0
    row_dim:  dimension of a, int>0
    col_dim:  dimension of b, int>0
    normalization: 'U' (default), 'L','R': normalization status of M^{s}_{a,b}
        only force 'L', 'R' manually if entirely certain it is true 
    init_random=False (default), True: initializes self.tensor to 0.0 or random complex numbers
    
    attributes:
    self.site_dim(): dimension of local physical space
    self.row_dim():  row dimension of MPS matrices
    self.col_dim():  column dimension of MPS matrices
        return int
    self.isobccompatible(pos): if pos=='L', True if row_dim==1; if pos=='R', True if col_dim==1;
        first and last matrices must be row and column vectors
        return bool
    
    self.isleftnormalized(test=False): MPS matrices left-normalized? 
        updates self.normalization ONLY if explicit calculation yields 'L' 
    self.isrightnormalized(test=False): MPS matrices right-normalized? 
        updates self.normalization ONLY if explicit calculation yields 'R' 
    self.isnormalized(test=False): MPS matrices either left- or right-normalized? 
        updates self.normalization ONLY if explicit calculation carried out 
        test==False: status read off from self.normalization
        test==True:  status verified by explicit calculation
        return bool
    
    self.normalize_left(decomposition_method='qr',multiplied=True,truncation=0.0,bond_dim=BOND_MAX): 
        left-normalize tensor (one step of Sec. 4.4.1 or 4.5.1)
        decomposition_method='qr' (default), 'svd': QR is faster, but does not provide truncation information
        multiplied=True: if SVD is used, returns SV+ multiplied (default), or S, V+
        truncation=0.0: if SVD is used, discards all singular values <= truncation (e.g. all zeros) (Sec. 4.5.1)
        bond_dim=BOND_MAX: if SVD is used, limits maximal column dimension to bond_dim
        the more stringent of the last two criteria dominates
        returns matrix (matrices) to be multiplied to the tensor to the right
    self.normalize_right(decomposition_method='qr',multiplied=True,truncation=0.0,bond_dim=BOND_MAX): 
        right-normalize tensor (one step of Sec. 4.4.2 or 4.5.1)
        as self.normalize_left, with following (straightforward) changes:
        multiplied=True: if SVD is used, returns US multiplied (default), or U, S
        returns matrix (matrices) to be multiplied to the tensor to the left
    self.multiply_from_left(N): multiply N * M^{s}; N returned from self.normalize_left; Eq. (136)
        N is 2D array
    self.multiply_from_right(N): multiply M^{s} * N; N returned from self.normalize_right; Eq. (137)
        N is 2D array
    self.multiply_by_scalar(scalar): multiply by scalar
    self.scalar_overlap(set2): returns scalar overlap \sum_{s,a,b} (self)^{s*}_{a,b} (set2)^{s}_{a,b})
    self.scalar_norm(): returns scalar norm \sqrt (\sum_{s,a,b} M^{s*}_{a,b} M^{s}_{a,b})
    self.copy(): 'true' copy by value
    self.mps_tensor_ismatched(set2): True if column dimension of self and row dimension of set2 match
    self.mps_tensor_physical_ismatched(set2): True if physical dimensions of self, set2 match
    self.mps_tensor_bra_left_ismatched(mat): True if bra of overlap matrix and left tensor dim. match
    self.mps_tensor_ket_left_ismatched(mat): True if ket of overlap matrix and left tensor dim. match
    self.mps_tensor_bra_right_ismatched(mat): True if bra of overlap matrix and right tensor dim. match
    self.mps_tensor_ket_right_ismatched(mat): True if ket of overlap matrix and right tensor dim. match
    self.mps_tensor_dummy(): returns (1x1x1) with entry 1.0
    
    """
    
# constructor
    
    def __init__(self, site_dim, row_dim, col_dim, init_random=False, normalization='U'):
        
# error handling in arguments
        
        if (not isinstance(site_dim, int)) or site_dim < 1:
            raise MPSTensorError("Site dimension must be positive integer.")
        if (not isinstance(row_dim, int)) or row_dim < 1:
            raise MPSTensorError("Row dimension must be positive integer.")
        if (not isinstance(col_dim, int)) or col_dim < 1:
            raise MPSTensorError("Column dimension must be positive integer")
        if not isinstance(init_random, bool):
            raise MPSTensorError("init_random must be Boolean")            
        if not normalization in ['U', 'L', 'R']:
            raise MPSTensorError("Normalization must be 'U', 'L', 'R'")

        if init_random == False:
            self.tensor = zeros((site_dim, row_dim, col_dim),complex)
        else:
            self.tensor = (rand(site_dim, row_dim, col_dim) + 1.0j *\
                           rand(site_dim, row_dim, col_dim)) / sqrt(row_dim * col_dim)
        self.normalization = normalization
        return
        

# tensor dimension queries

    def site_dim(self):
        return size(self.tensor, 0)


    def row_dim(self):
        return size(self.tensor, 1)


    def col_dim(self):
        return size(self.tensor, 2)
    
    
    def isobccompatible(self, pos):
        if (pos == 'L' and self.row_dim() == 1) or\
           (pos == 'R' and self.col_dim() == 1):
            return True
        else:
            return False        


# normalization queries

# check (test) for left normalization; set if confirmed by explicit test

    def isleftnormalized(self, test=False):
        if test == False:
            if self.normalization == 'L':
                return True
            else:
                return False
        else:
            if norm(identity(self.col_dim(), complex) - \
                    tensordot(self.tensor.conj(), self.tensor, axes = ([1, 0], [1, 0]))) < NORM_TOL:
                self.normalization = 'L'
                return True
            else:
                return False


# check (test) for right normalization; set if confirmed by explicit test

    def isrightnormalized(self, test=False):
        if test == False:
            if self.normalization == 'R':
                return True
            else:
                return False
        else:       
            if norm(identity(self.row_dim(), complex) - \
                    tensordot(self.tensor, self.tensor.conj(), axes = ([2, 0], [2, 0]))) < NORM_TOL:
                self.normalization= 'R' 
                return True
            else:
                return False


# check (test) for normalization, if explicit tests fail in confirming normalization, update to 'U'

    def isnormalized(self, test=False):
        if self.isleftnormalized(test) or self.isrightnormalized(test):
            return True
        else:
            self.normalization = 'U'
            return False


# normalization methods

    def normalize_left(self, decomposition_method='qr', 
                       multiplied=True, truncation=0.0, bond_dim=BOND_MAX):

        if not decomposition_method in ['qr', 'svd']: 
            raise DecompositionError("Decomposition method must be either 'qr' or 'svd'.")
        if (not isinstance(truncation, float)):
            raise SingularValueTruncationError("Truncation limit not valid.")
        if (not isinstance(bond_dim, int)) or bond_dim < 1:
            raise SingularValueTruncationError("Dimensional limit must be positive integer.")

        if decomposition_method == 'qr': # here with qr, faster, but less information 
            if multiplied == False:
                raise DecompositionError("QR decomposition returns only one matrix")
            if truncation != 0.0:
                raise DecompositionError("QR cannot truncate non-trivially")
            q, r = qr(self.tensor.reshape(self.site_dim() * self.row_dim(), self.col_dim())) # Eq. (28), Sec. 4.4.1
            min_dim = min(self.site_dim() * self.row_dim(), self.col_dim()) # thin QR, Eq. (29) - typo!!
            self.tensor = q.reshape(self.site_dim(), self.row_dim(), min_dim) # reshape back
            self.normalization = 'L'
            return r
        else: # now with svd, slower, full information 
            u, s, vd = svd(self.tensor.reshape(self.site_dim() * self.row_dim(), self.col_dim()),0)
            min_dim = min(self.site_dim() * self.row_dim(), self.col_dim())
            self.tensor = u.reshape(self.site_dim(), self.row_dim(), min_dim)
            self.normalization = 'L'
            # restrict to the (at most bond_dim) SV greater than truncation
            # truncation=0.0: discard SVs that are strictly zero (default)
            # truncation<0.0: discard none
            # take into account that sum of SV^2 does not add to 1.0
            norm_factor = norm(s, 0)
            compress_list = [item > truncation * norm_factor for item in list(s)]
            if len(compress_list) > bond_dim: 
                compress_list[bond_dim : ] = [] # relies on ordering of SVs by svd
            if len(compress_list) < size(s, 0):
                s = compress(compress_list, s, 0)
                vd = compress(compress_list, vd, 0)
            if multiplied == True:
                for k in range(size(vd, 0)):
                    vd[k, : ] *= s[k]
                return vd
            else:
                return s, vd
                

    def normalize_right(self, decomposition_method='qr',
                        multiplied=True, truncation=0.0, bond_dim=BOND_MAX):
        if not decomposition_method in ['qr', 'svd']: 
            raise DecompositionError("Decomposition method must be either 'qr' or 'svd'.")
        if (not isinstance(truncation, float)) or truncation < 0.0:
            raise SingularValueTruncationError("Truncation limit not in reasonable range.")
        if (not isinstance(bond_dim, int)) or bond_dim < 1:
            raise SingularValueTruncationError("Dimensional limit must be positive integer.")
        if decomposition_method == 'qr': # here with qr, faster, but less information
            if multiplied == False:
                raise DecompositionError("QR returns only one matrix.")
            if truncation != 0.0:
                raise DecompositionError("QR cannot truncate nontrivially.")
            q, r = qr(self.tensor.\
                      transpose(0, 2, 1).conj().\
                      reshape(self.site_dim() * self.col_dim(), self.row_dim())) # Sec. 4.4.2
            min_dim = min(self.row_dim(), self.site_dim() * self.col_dim())
            self.tensor = q.conj().reshape(self.site_dim(), self.col_dim(), min_dim).transpose(0, 2, 1)
            self.normalization = 'R'
            return r.transpose(1, 0).conj()
        else: # now with svd, slower, full information 
            u, s, vd = svd(self.tensor.\
                           transpose(1, 0, 2).\
                           reshape(self.row_dim(), self.site_dim() * self.col_dim()),0)
            min_dim = min(self.row_dim(), self.site_dim() * self.col_dim())
            self.tensor = vd.reshape(min_dim, self.site_dim(), self.col_dim()).transpose(1, 0, 2)
            self.normalization = 'R'
            # restrict to the (at most bond_dim) SV greater than truncation
            # truncation=0.0: discard SVs that are strictly zero (default)
            # truncation<0.0: discard none
            norm_factor = norm(s, 0)
            compress_list = [item > truncation * norm_factor for item in list(s)]
            if len(compress_list) > bond_dim:
                compress_list[bond_dim : ] = []
            if len(compress_list) < size(s, 0):
                s = compress(compress_list, s, 0)
                u = compress(compress_list, u, 1)
            if multiplied == True:
                for k in range(size(u, 1)):
                    u[ : , k] *= s[k]
                return u
            else:
                return u,s


# multiply matrices to a set of M^{s} from left or right
# attention: tensordot orders remaining dimensions from 'left to right', may force transpose
# new normalization status unknown
            
    def multiply_from_left(self, mat):
        if size(mat, 1) != self.row_dim():
            raise DimensionMismatchError("column and row dimensions do not match")
        self.tensor = tensordot(mat, self.tensor, axes = (1, 1)).transpose(1, 0, 2)
        self.normalization = 'U'
        return
    

    def multiply_from_right(self, mat):
        if size(mat, 0) != self.col_dim():
            raise DimensionMismatchError("column and row dimensions do not match")
        self.tensor = tensordot(self.tensor, mat, axes=(2, 0))
        self.normalization = 'U'
        return
                

    def multiply_by_scalar(self, scalar):
        self.tensor *= scalar
        self.normalization = 'U'
        return
    

    def scalar_overlap(self, set2): 
        if self.site_dim() != set2.site_dim() or \
           self.row_dim() != set2.row_dim() or \
           self.col_dim() != set2.col_dim():
            raise DimensionMismatchError("column and row dimensions do not match")
        return tensordot(self.tensor.conj(), set2.tensor, ([2,1,0],[2,1,0])) # shorter notation around? 
    

    def scalar_norm(self):
        return abs(sqrt(MPSTensor.scalar_overlap(self, self)))
        

# copy method

    def copy(self):
        _acopy_ = MPSTensor(self.site_dim(), self.row_dim(), self.col_dim(),\
                            self.init_random, self.normalization)
        _acopy_.tensor = self.tensor.copy()
        return _acopy_


# consistency checks for 2 MPSTensors
# mps_tensor_ismatched checks in bond space
# mps_tensor_physical_ismatched checks in physical space
    
    def mps_tensor_ismatched(self, set2):

# do the dimensions of adjoining local tensors match for legs to be contracted?
#    
# arguments:
# self, set2: two MPSTensors, self to left of set 2
#    
# returns bool
    
        if self.col_dim() == set2.row_dim():
            return True
        else:
            return False
    
    
    def mps_tensor_physical_ismatched(self, set2):
        
# do the physical dimensions of two local tensors match for contraction?
# 
# arguments:
# self, set2: two MPSTensors, self to left of set 2
#
# returns bool
     
        if self.site_dim() == set2.site_dim():
            return True
        else:
            return False

        
    def mps_tensor_bra_left_ismatched(self, lc):
    
        # do the bra dimension of contracted c and left dimension of tensor match?
        # lc: left contracted matrix
        # returns bool
        
        if size(lc, 0) == self.row_dim():
            return True
        else:
            return False
    

    def mps_tensor_ket_left_ismatched(self, lc):

        # do the ket dimension of contracted c and left dimension of tensor match?
        # lc: left contracted matrix
        # returns bool

        if size(lc, 1) == self.row_dim():
            return True
        else:
            return False
    

    def mps_tensor_bra_right_ismatched(self, rc):

        # do the bra dimension of contracted c and right dimension of tensor match?
        # lc: left contracted matrix
        # returns bool
        
        if size(rc, 0) == self.col_dim():
            return True
        else:
            return False
    

    def mps_tensor_ket_right_ismatched(self, rc):

        # do the ket dimension of contracted c and right dimension of tensor match?
        # lc: left contracted matrix
        # returns bool
        
        if size(rc, 1) == self.col_dim():
            return True
        else:
            return False


# last, but not least: return a dummy local MPS tensor.         
        
    def mps_tensor_dummy():
        _dummy_ = MPSTensor((1, 1, 1))
        _dummy_.tensor[0, 0, 0] = 1.0
        return _dummy_
        
        
# #################################################################################

class LocalTransferOperator:
    """
    local transfer operator constructed from bra- and -ket LocalMPSSet's
    --------------------------------------------------------------------
    
    structure:
    self[lb, lk, rb, rk] = \sum_{s} M1^{s*}_{lb,rb} * M2^{s}_{lk,rk}
        4 dimensional array of complex
        for eigenvalue, eigenstate analysis it must be flattened
    
    arguments for constructor:
    bra_local, ket_local: two MPSTensors, corresponding to M1, M2 in formula
    
    attributes:
    self.left_dim(): left dimension
    self.right_dim(): right dimension
    self.left_bra_dim(): left bra dimension
    self.left_ket_dim(): left ket dimension
    self.right_bra_dim(): right bra dimension
    self.right_ket_dim(): right ket dimension
    self.flatten(): returns matrix mat[(left_bra,left_ket),(right_bra,right_ket)] flatted from self TODO
    self.unflatten_as(mat): returns unflattened mat with dimensions of self TODO
    
    self.copy(): deep copy
    
    no others implemented, as computationally inefficient object
    is however conceptually very important, so new methods might come 
    
    """
    
    def __init__(self, bra, ket):

        if not local_mps_physical_ismatched(bra, ket):
            raise MatrixProductStateDimensionError("Physical dimensions of MPS tensors do not match.")
        
        self.tensor=tensordot(bra.tensor.conj(), ket.tensor, (0, 0)).transpose(0, 2, 1, 3)
        return
    

    def left_dim(self):
        return size(self.tensor, 0) * size(self.tensor, 1)

    
    def right_dim(self):
        return size(self.tensor, 2) * size(self.tensor, 3)
    

    def left_bra_dim(self):
        return size(self.tensor, 0)

    
    def left_ket_dim(self):
        return size(self.tensor, 1)

    
    def right_bra_dim(self):
        return size(self.tensor, 2)
    
    
    def right_ket_dim(self):
        return size(self.tensor, 3)

 
# flatten local transfer operator into matrix, grouping left and right indices

    def flatten():
        flat_matrix = self.tensor.reshape(self.left_bra_dim() * self.left_ket_dim(),
                                          self.right_bra_dim() * self.right_ket_dim())
        return flat_matrix
        
        
# unflatten matrix into local transfer operator, aligning structure of self

    def unflattenas(flat_matrix):
        if self.left_bra_dim() * self.left_ket_dim() != size(flat_matrix, 0) or \
           self.right_bra_dim() * self.right_ket_dim() != size(flat_matrix, 1):
            raise DimensionMismatchError("Dimensions of matrix to be unflattened and model transfer operator do not match.")
        new_lto = LocalTransferOperator(mps_tensor_dummy, mps_tensor_dummy)
        new_lto.tensor = flat_matrix.reshape(self.left_bra_dim(), self.left_ket_dim(), \
                                             self.right_bra_dim(), self.right_ket_dim())
        return new_lto
    

    def copy(self):
        
        _ltocopy_= LocalTransferOperator(MPSTensor((1, 1, 1)), MPSTensor((1, 1, 1))) # dummy
        _ltocopy_.tensor = self.tensor.copy() # deep copy
        return _ltocopy_
    

# #################################################################################

# outside the class MPSTensor: steps in overlap calculations

    
# single steps of carrying zip one step in overlap forward from left or from right
# potentially including a local operator O_i
#
# IMPORTANT: for notational consistency, contracted objects growing from RIGHT
# contain the COMPLEX CONJUGATE in deviation from review notation

def overlap_left_step(bra_tensor, ket_tensor, c, local_operator=[]):
    """
    overlap matrix c carried one site from left
    -------------------------------------------

    arguments:
    bra_tensor: MPSTensor from bra - will be complex conjugated
    ket_tensor: MPSTensor from ket
    c: left contracted overlap matrix to be carried forward
    local_operator: 2d-array of local operator, nothing done if [] (default) 
    
    output:
    d:  left contracted overlap matrix carried one site forward
    see Eqs. (98), (99)
    """
    
    # error handling
    
    if not bra_tensor.mps_tensor_bra_left_ismatched(c):
        raise MatrixProductStateDimensionError("Dimensions of overlap matrix and bra set do not match.")
    if not ket_tensor.mps_tensor_ket_left_ismatched(c):
        raise MatrixProductStateDimensionError("Dimensions of overlap matrix and ket set do not match.")
    if not bra_tensor.mps_tensor_physical_ismatched(ket_tensor):    
        raise MatrixProductStateDimensionError("Local physical dimensions of the states do not match.")
    if local_operator != []:
        if size(local_operator, 0) != bra_tensor.site_dim() or \
           size(local_operator, 1) != ket_tensor.site_dim():
            raise MatrixProductStateDimensionError("Local physical dimensions of the states and operator do not match.")
           
    if local_operator == []: # no sandwiched local operator
        c1 = tensordot(c, ket_tensor.tensor, (1, 1))
        d = tensordot(bra_tensor.tensor.conj(), c1, ([1, 0], [0, 1]))
    else: # sandwiched local operator
        c1 = tensordot(c, ket_tensor.tensor, (1, 1))
        c2 = tensordot(local_operator, c1, (1, 1))
        d = tensordot(bra_tensor.tensor.conj(), c2, ([1, 0], [1, 0]))
    return d

                 
def overlap_right_step(bra_tensor, ket_tensor, d, local_operator=[]):
    """
    overlap matrix d carried one site from right
    --------------------------------------------

    arguments:
    bra_tensor: MPSTensor from bra - will be complex conjugated
    ket_tensor: MPSTensor from ket
    d: right contracted overlap matrix to be carried forward
    local_operator: 2d-array of local operator, nothing done if [] (default)
    
    output:
    c: right contracted overlap matrix carried one site forward
    see Eqs. (102), (103)
    
    IMPORTANT: for bra/ket consistency, transpose of (102),(103), and COMPLEX CONJUGATE
    to be taken into account when using result!
    """

    # error handling
    
    if not bra_tensor.mps_tensor_bra_right_ismatched(d):
        raise MatrixProductStateDimensionError("Dimensions of overlap matrix and bra set do not match.")
    if not ket_tensor.mps_tensor_ket_right_ismatched(d):
        raise MatrixProductStateDimensionError("Dimensions of overlap matrix and ket set do not match.")
    if not bra_tensor.mps_tensor_physical_ismatched(ket_tensor):    
        raise MatrixProductStateDimensionError("Local physical dimensions of the states do not match.")
    if local_operator != []:
        if size(local_operator, 0) != bra_tensor.site_dim() or \
           size(local_operator, 1) != ket_tensor.site_dim():
            raise MatrixProductStateDimensionError("Local physical dimensions of the states and operator do not match.")
           
    if local_operator == []:
        d1 = tensordot(d, ket_tensor.tensor.conj(), (1, 2))
        c = tensordot(bra_tensor.tensor, d1, ([2, 0], [0, 1]))
    else:
        d1 = tensordot(d, ket_tensor.tensor.conj(), (1, 2))
        d2 = tensordot(local_operator.conj(), d1, (1, 1))
        c = tensordot(bra_tensor.tensor, d2, ([2, 0], [1, 0]))        
    return c

   
# #################################################################################

class MPS:
    """
    matrix product state
    --------------------
    
    structure: see Eqs. (75), (77) -- essentially a list of tensors; dummy for site 0   
    self.prefactor: complex scalar for scaling normalized AAAAA or BBBBB structure
    self.m: list of (length+1) MPSTensors describing M^{s}_{a,b}, where length is the number of sites
        m[0] will be (1 x 1 x 1) tensor with entry 1.0, please do not change
        
    attributes:
    self.length():    length of system described by MPS
    self.site_dim(i): dimension of physical state space at site i
    self.row_dim(i):  row dimension of matrices at site i
    self.col_dim(i):  column dimension of matrices at site i
    self.site_dims(): list of dimensions of physical state spaces
    self.row_dims():  list of row dimensions of matrices
    self.col_dims():  list of column dimensions of matrices
    
    self.ismixedcanonical(center,test): checks AAAAAXBBBBB form, where X is on site (center).
        center: integer site, to the right of L-normalized matrices and left of R-normalized matrices
        choosing center=0, or length+1, defaults to all R-normalized, all L-normalized
    self.isleftcanonical(test):   checks L-normalization of MPS, Eqs. (37), (38)
    self.isrightcanonical(test):  checks R-normalization of MPS, Eqs. (44), (45)
    self.ismixednormalized(center,test): as self.ismixedcanonical, and checks norm of X.
    self.isleftnormalized(test):  checks L-normalization of MPS, Eqs. (37), (38), and self.prefactor==1
    self.isrightnormalized(test): checks R-normalization of MPS, Eqs. (44), (45), and self.prefactor==1
        test=False (default), True: determines whether real check is calculated
    self.mps_ismatched(mps2): determines whether the physical dimensions of 2 MPS match
        
    self.copy(): deep copy
    
    self.canonize_left_step(i,test): if site i is not L-normalized: 
        make it so, return remaining matrix
    self.canonize_right_step(i,test): if site i is not R-normalized: 
        make it so, return remaining matrix
    self.canonize(center,test):  L-normalize before site 'center' and R-normalize after it
    self.canonize_left(test):    L-normalize all sites
    self.canonize_right(test):   R-normalize all sites

    self.normalize(center,test): like self.canonize, but normalize site 'center' matrices 
        as a normal wave function (because the L- and R-matrices span ONB)
    self.normalize_left(test):   like self.canonize_left, but set self.prefactor=1.0
    self.normalize_right(test):  like self.canonize_right, but set self.prefactor=1.0
    
    self.norm2(test): squared norm with (un)tested accelerations if local normalizations hold
    self.overlap(ket,withprefactors=True,direction='L'): overlap of self as bra with ket, including
        prefactors (or not), zipping from Left or Right

    self.compress_svd_left_site(i,sv_min,bond_dim,test): 
        checks if state in form AAAAAXBBB, X at site i, takes it to AAAAAAXBB, truncates 
    self.compress_svd_right_site(i,sv_min,bond_dim,test):
        checks if state in form AAAAAXBBB, X at site i, takes it to AAAAXBBBB, truncates 
    self.compress_svd_left(sv_min,bond_dim,test): 
        starting from left on R-canonical state (is checked), makes L-canonical, truncates      
    self.compress_svd_right(sv_min,bond_dim,test): 
        starting from right on L-canonical state (is checked), makes R-canonical, truncates
        sv_min=0.0: only singular values larger are considered
        bond_dim=BOND_MAX: at most dimension bond_dim is allowed
        test=False: test=True if explicit test of matrix normalization desired
    """    
    
# initialization

    def __init__(self):

        
# self.m[0].tensor is special dummy

        self.m = []
        self.m.append(MPSTensor(1,1,1))
        self.prefactor = 1.0 
        return
        

# size queries


    def length(self):
        return len(self.m) - 1
    

    def site_dim(self, i):
        return self.m[i].site_dim()
    

    def row_dim(self, i):
        return self.m[i].row_dim()
    

    def col_dim(self, i):
        return self.m[i].col_dim()
    

    def site_dims(self):
        return [item.site_dim() for item in self.m]
        

    def row_dims(self):
        return [item.row_dim() for item in self.m]


    def col_dims(self):
        return [item.col_dim() for item in self.m]
    

    def ismixedcanonical(self, center, test=False):
        for i in range(1, min(center, self.length() + 1)):
            if not self.m[i].isleftnormalized(test):
                return False
        for i in range(max(1, center + 1), self.length() + 1):
            if not self.m[i].isrightnormalized(test):
                return False
        return True


    def isleftcanonical(self, test=False):
        return self.ismixedcanonical(self.length() + 1, test)
    

    def isrightcanonical(self, test=False):
        return self.ismixedcanonical(0, test)
    

    def isnormalized(self, center, test=False):
        if not self.ismixedcanonical(center, test):
            return False
        else:
            if not center in range(1, self.length() + 1):
                return abs(abs(self.prefactor) - 1.0) < NORM_TOL
            else:
                return abs(self.m[center].scalar_norm() - 1.0) < NORM_TOL 


    def isleftnormalized(self, test=False):
        return self.isnormalized(self.length() + 1, test)


    def isrightnormalized(self, test=False):
        return self.isnormalized(0, test)


    def mps_ismatched(self, mps2):
        if self.length() != mps2.length():
            return False
        else:
            for i in (1, self.length() + 1):
                if not mps_tensor_physical_ismatched(self.m[i], mps2.m[i]):
                    return False
            return True
        
        
# copy method 

    def copy(self):
        _psicopy_ = MPS()
        _psicopy_.prefactor = self.prefactor
        for i in range(1, len(self.m)):
            _psicopy_.m.append(self.m[i])
            _psicopy_.m[i] = self.m[i].copy()
        return _psicopy_
    

# transform to left-normalization at site i (unless already done)
# AAAAAX???->AAAAAA???
# changes also content of site i+1, destroying normalization
# if i+1==length+1 then adjust prefactor instead

    def canonize_left_step(self, i, test=False):
        if self.m[i].isleftnormalized(test) == False:
            r = self.m[i].normalize_left()
            if i != self.length():
                self.m[i+1].multiply_from_left(r)
                self.m[i+1].normalization = 'U'
            else:
                self.prefactor *= r[0,0]
        return
                

# transform to right-normalization at site i (unless already done)
# ???XBBBBB->???BBBBBB
# changes also content of site i-1, destroying normalization
# if i==1 then adjust prefactor

    def canonize_right_step(self, i, test=False):
        if self.m[i].isrightnormalized(test) == False:
            r = self.m[i].normalize_right()
            if i != 1:
                self.m[i-1].multiply_from_right(r)
                self.m[i-1].normalization = 'U'
            else:
                self.prefactor *= r[0,0]
        return
        

# transform to mixed canonical representation, AAAAAAXBBBB where X is site 'center'
    
    def canonize(self, center, test=False):
        for i in range(1, min(center, self.length() + 1)):
            self.canonize_left_step(i, test)
        for i in range(self.length(), max(0, center), -1):
            self.canonize_right_step(i, test)
        return


# transform to canonical representation with L-normalization, Sec. 4.4.1
# in the end, leftover factor multiplied into self.prefactor
    
    def canonize_left(self, test=False):
        self.canonize(self.length() + 1, test)
        return

    
# transform to canonical representation with R-normalization, Sec. 4.4.2
# in the end, leftover factor multiplied into self.prefactor

    def canonize_right(self, test=False):
        self.canonize(0, test)
        return
         

# normalize from left or right is like canonize, 
# but with the prefactor forced to be 1.0 (if L- or R-canonical) or center matrices normalized
# then the state has norm 1.0


    def normalize(self, center, test=False):
        if not self.ismixedcanonical(center, test):
            self.canonize(center, test)
        if not center in range(1, self.length()+1): # left or right-canonical
            if self.prefactor != 0.0:
                self.prefactor = 1.0
                return
            else:
                raise NormalizationError("State not normalizable.")
        else: # mixed-canonical
            norm_factor = self.m[center].scalar_norm()
            if norm_factor != 0.0 and self.prefactor != 0.0:
                norm_factor = 1.0 / norm_factor
                self.m[center].multiply_by_scalar(norm_factor)
                self.prefactor = 1.0
                return
            else:
                raise NormalizationError("State not normalizable.")
            

    def normalize_left(self, test=False):
        self.normalize(self.length()+1, test)
        return
    

    def normalize_right(self, test=False):
        self.normalize(0, test)
        return


# calculates squared norm of a state
# if string of L-normalized matrices at the left -> trivial contraction from left
# if string of R-normalized matrices at the right -> trivial contraction from right
    
    def norm2(self, test=False):
        
        l_end = 1
        r_end = self.length()
        while l_end <= self.length() and self.m[l_end].isleftnormalized(test): 
            l_end += 1
        while r_end >= 1 and self.m[r_end].isrightnormalized(test):
            r_end -= 1
        if l_end == self.length() + 1 or r_end == 0:
            return self.prefactor.conj() * self.prefactor
        else:
            c = identity(self.row_dim(l_end), complex)
            for i in range(l_end, r_end + 1):
                c = overlap_left_step(self.m[i], self.m[i], c)
            return trace(c)


# calculates overlap <bra|ket>
# arguments:
# bra: bra MPS
# ket: ket MPS
# withprefactors=True: take *.prefactors into account
# direction='L': zip from 'L'eft or 'R'ight
# output:
# <bra|ket> 
# if bra=ket, ket.norm2() is faster, analyzes simplifications of ON-structure
        
    def overlap(self,ket,withprefactors=True,direction='L'):
        
        # error handling
        
        if self.length() != ket.length():
            raise MatrixProductStateDimensionError("Lengths of the states do not match.")
        
        c = ones((1, 1), complex)
        if direction == 'L':
            for i in range(1, ket.length() + 1):
                c = overlap_left_step(self.m[i], ket.m[i], c).copy() # is this copy() necessary? check!
        else:
            for i in range(ket.length(), 0, -1):
                c = overlap_right_step(self.m[i], ket.m[i], c).copy() # is this copy() necessary? check!
            c[0, 0] = conj(c[0, 0]) # complex conjugate because of construction
        if withprefactors == True:
            c[0, 0] *= conj(self.prefactor) * ket.prefactor
        return c[0, 0]


# compression by svd; Sec. 4.5.1
# on site X, if AAAAAXBBB form: L-normalize and compress X and change site to right: AAAAAA(RB)BB

    def compress_svd_left_site(self, i, sv_min=0.0, bond_dim=BOND_MAX, test=False):
    
        if not self.ismixedcanonical(i, test):
            raise CanonicalFormMismatchError("State is not in right canonical form for SVD compression.")
        
        r = self.m[i].normalize_left('svd', True, sv_min, bond_dim)        
        if i != self.length():
            self.m[i+1].multiply_from_left(r)
        else:
            self.prefactor *= r[0,0]
        return


# on site X, if AAAAAXBBB form: R-normalize and compress X and change site to left: AAAA(AR)BBBB

    def compress_svd_right_site(self, i, sv_min=0.0, bond_dim=BOND_MAX, test=False):
    
        if not self.ismixedcanonical(i, test):
            raise CanonicalFormMismatchError("State is not in right canonical form for SVD compression.")
        
        r = self.m[i].normalize_right('svd', True, sv_min, bond_dim)        
        if i != 1:
            self.m[i-1].multiply_from_right(r)
        else:
            self.prefactor *= r[0,0]
        return
    

    def compress_svd_left(self, sv_min=0.0, bond_dim=BOND_MAX, test=False):
        
        for i in range(1, self.length()+1):
            self.compress_svd_left_site(i, sv_min, bond_dim, test)


    def compress_svd_right(self, sv_min=0.0, bond_dim=BOND_MAX, test=False):
        
        for i in range(self.length(), 0, -1):
            self.compress_svd_right_site(i, sv_min, bond_dim, test)
            
            
    def compress_iterative_left(self, psi_guess, sv_min=0.0, bond_max=BOND_MAX, test=False):

# check consistency of self and psi_guess
        
        if self.length() != psi_guess.length():
            raise MatrixProductStateDimensionError("length of state to be compressed and initial guess do not match")
        if self.site_dims() != psi_guess.site_dims():
            raise MatrixProductStateDimensionError("local physical dimensions of state to be compressed and initial guess do not match") 
        if not (psi_guess.isrightcanonical(test) or psi_guess.ismixedcanonical(1,test)):
            psi_guess.canonize_right()
            
        list_lc = []
        list_rc = []
        
        rc = array((1, 1),complex)
        rc[0, 0] = 1.0
        lc = array((1, 1),complex)
        lc[0, 0] = 1.0
        list_rc.append(rc)
        list_lc.append(lc)
        
        convergence = False
        
# create list of right-contractions        

        for i in range(self.length(), 1, -1):
            rc = overlap_right_step(psi_guess.m[i], self.m[i], rc)
            list_rc.append(rc)
        
# now dynamical sweeping forth and back

        while convergence == False:
        
# from left to right ...

            for i in range(1, self.length()):
                rc = list_rc.pop()
                psi_guess.m[i].tensor = tensordot(tensordot(list_lc[-1], self.m[i].tensor, (1, 1)), rc.conj(), (2, 1)).copy()
                psi_guess.m[i].normalization = 'U'
                psi_guess.canonize_left_step(i, False)
                list_lc.append(overlap_left_step(psi_guess.m[i], self.m[i], list_lc[-1]))
                            
# from right to left ...

            for i in range(self.length(), 1, -1):
                lc = list_lc.pop()
                psi_guess.m[i].tensor = tensordot(tensordot(lc, self.m[i].tensor, (1, 1)), list_rc[-1], (2, 1)).copy()
                psi_guess.m[i].normalization = 'U'
                psi_guess.canonize_right_step(i, False)
                list_lc.append(overlap_right_step(psi_guess.m[i], self.m[i], list_rc[-1]))

        return psi_guess

    def compress_iterative_right(self, psi_guess, sv_min=0.0, bond_max=BOND_MAX, test=False):

# check consistency of self and psi_guess
        
        if self.length() != psi_guess.length():
            raise MatrixProductStateDimensionError("length of state to be compressed and initial guess do not match")
        if self.site_dims() != psi_guess.site_dims():
            raise MatrixProductStateDimensionError("local physical dimensions of state to be compressed and initial guess do not match") 
        if not (psi_guess.isleftcanonical(test) or psi_guess.ismixedcanonical(psi_guess.length(),test)):
            psi_guess.canonize_right()
            
        list_lc = []
        list_rc = []
        
        rc = array((1, 1),complex)
        rc[0, 0] = 1.0
        lc = array((1, 1),complex)
        lc[0, 0] = 1.0
        list_rc.append(rc)
        list_lc.append(lc)
        
        convergence = False
        
# create list of left-contractions        

        for i in range(1, self.length()):
            lc = overlap_left_step(psi_guess.m[i], self.m[i], lc)
            list_lc.append(lc)
        
# now dynamical sweeping forth and back

        while convergence == False:
        
# from right to left ...

            for i in range(self.length(), 1, -1):
                lc = list_lc.pop()
                psi_guess.m[i].tensor = tensordot(tensordot(lc, self.m[i].tensor, (1, 1)), list_rc[-1], (2, 1)).copy()
                psi_guess.m[i].normalization = 'U'
                psi_guess.canonize_right_step(i, False)
                list_lc.append(overlap_right_step(psi_guess.m[i], self.m[i], list_rc[-1]))

# from left to right ...

            for i in range(1, self.length()):
                rc = list_rc.pop()
                psi_guess.m[i].tensor = tensordot(tensordot(list_lc[-1], self.m[i].tensor, (1, 1)), rc.conj(), (2, 1)).copy()
                psi_guess.m[i].normalization = 'U'
                psi_guess.canonize_left_step(i, False)
                list_lc.append(overlap_left_step(psi_guess.m[i], self.m[i], list_lc[-1]))
                            
        return psi_guess

# ############

# matching checks for LeftContracted, RightContracted, LocalMPSSets, in bra and ket.


# ################ current leftovers ##################

class LeftContracted:
    """
    left contraction part of two MPS
    --------------------------------
    
    structure:
    self.tensor: 2 dimensional array of complex
    
    arguments: 
    none at the moment, as the initialization is just a dummy matrix (1x1) = 1.0
    
    attributes:
    self.bra_dim(): dimension of bra link
    self.ket_dim(): dimension of ket link
    self.copy(): deep copy
    """
    
    def __init__(self):
        self.tensor=array((1, 1), complex)
        return
    

    def bra_dim(self):
        return size(self.tensor, 0)
    

    def ket_dim(self):
        return size(self.tensor, 1)
    
    
    def copy(self):        
        _lccopy_=LeftContracted()
        _lccopy_.tensor=self.tensor.copy()
        return _lccopy_
    
    
class RightContracted:
    """
    right contraction part of two MPS
    ---------------------------------
    
    structure:
    self.tensor: 2 dimensional array of complex
    
    arguments: 
    none at the moment, as the initialization is just a dummy matrix (1x1) = 1.0
    
    attributes:
    self.bra_dim(): dimension of bra link
    self.ket_dim(): dimension of ket link
    self.copy(): deep copy
    """
    
    def __init__(self):
        self.tensor=array((1, 1), complex)
        return
    

    def bra_dim(self):
        return size(self.tensor, 0)
    

    def ket_dim(self):
        return size(self.tensor, 1)    
    

    def copy(self):
        _rccopy_=RightContracted()
        _rccopy_.tensor=self.tensor.copy()
        return _rccopy_
    


        
        
