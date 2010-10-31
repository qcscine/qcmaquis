from numpy import *
from constants import *
from myexceptions import *
from numpy.linalg import norm,qr
from numpy.random import rand
from numpy import tensordot,compress
from numpy.dual import svd

# ######################################################################################    

# class for local MPO tensors W^{s,t}_{a,b}

class MPOTensor:
    """
    local MPO tensor
    ----------------
    
    structure:
    self.tensor: complex 4-index tensor, self.tensor[s,t,a,b]=W^{s,t}_{a,b} Eq. (177)
    self.normalization: normalization of tensor; 'U', 'L', 'R'
        'U': unnormalized or status unknown
        'L': left-normalized \sum_{s,t} M^{s,t\dagger}M^{s,t}=I by analogy to MPS
        'R': right-normalized \sum_{s,t} M^{s,t}M^{s,t\dagger}=I by analogy to MPS
    
    arguments of constructor:
    site_dim: dimension of s, int>0
    row_dim:  dimension of a, int>0
    col_dim:  dimension of b, int>0
    normalization: 'U' (default), 'L','R': normalization status of W^{s,t}_{a,b}
        only force 'L', 'R' manually if entirely certain it is true 
    init_random=False (default), True: initializes self.tensor to 0.0 or random complex numbers

    attributes:
    
    self.site_bra_dim(): dimension of local physical bra space
    self.site_ket_dim(): dimension of local physical ket space    
    self.row_dim(): row dimension of MPO matrices
    self.col_dim(): column dimension of MPO matrices
    self.isobccompatible(pos): if pos=='L', True if row_dim==1; if pos=='R', True if col_dim==1;
        first and last matrices must be row and column vectors
        return bool
    self.isleftnormalized(test=False): MPO matrices left-normalized? 
    self.isrightnormalized(test=False): MPO matrices right-normalized? 
    self.isnormalized(test=False): MPO matrices either left- or right-normalized? 
        if test==False only from self.normalization; if test==True numerical verification
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
    self.multiply_from_left(N): multiply N * M^{s,t}; N returned from self.normalize_left; Eq. (136)
        N is 2D array
    self.multiply_from_right(N): multiply M^{s,t} * N; N returned from self.normalize_right; Eq. (137)
        N is 2D array
    self.multiply_by_scalar(scalar): multiply by scalar
    self.copy(): deep copy
    self.mpo_tensor_ismatched(set2): True if column dimension of self and row dimension of set2 match
    self.mpo_tensor_physical_ismatched(set2): True if physical dimensions of self, set2 match
    self.mpo_tensor_dummy(): returns (1x1x1x1) with entry 1.0
    """
    
# initialization:
    
    def __init__(self, site_bra_dim, site_ket_dim, row_dim, col_dim, \
                 init_random=False, normalization='U'):

# error handling
        
        if (not isinstance(site_bra_dim, int)) or site_bra_dim < 1:
            raise MPOTensorError("Physical bra dimension must be positive integer.")
        if (not isinstance(site_ket_dim,int)) or site_ket_dim < 1:
            raise MPOTensorError("Physical ket dimension must be positive integer.")
        if (not isinstance(row_dim, int)) or row_dim < 1:
            raise MPOTensorError("Row dimension must be positive integer.")
        if (not isinstance(col_dim, int)) or col_dim < 1:
            raise MPOTensorError("Column dimension must be positive integer.")
        if not isinstance(init_random, bool):
            raise MPOTensorError("Variable init_random must be boolean.")            
        if not normalization in ['U', 'L', 'R']:
            raise MPOTensorError("Normalization must be 'U', 'L', 'R'.")

        if init_random == False:
            self.tensor = zeros((site_bra_dim, site_ket_dim, row_dim, col_dim), complex)
        else:
            self.tensor = (rand(site_bra_dim, site_ket_dim, row_dim, col_dim) + 1.0j * \
                           rand(site_bra_dim, site_ket_dim, row_dim, col_dim)) / sqrt(row_dim * col_dim)
        self.normalization = normalization
        return
    
        
# simple queries

    def site_bra_dim(self):
        return size(self.tensor, 0)


    def site_ket_dim(self):
        return size(self.tensor, 1)


    def row_dim(self):
        return size(self.tensor, 2)


    def col_dim(self):
        return size(self.tensor, 3)
    

    def isobccompatible(self, pos):
        if (pos == 'L' and self.row_dim() == 1) or \
           (pos == 'R' and self.col_dim() == 1):
            return True
        else:
            return False        


# advanced queries

# if normalization is checked explicitly, the result updates self.normalization

    def isleftnormalized(self, test=False):
        if test == False:
            if self.normalization == 'L':
                return True
            else:
                return False
        else:
            if norm(identity(self.col_dim(), complex) -\
                    tensordot(self.tensor.conj(), self.tensor, axes=([2, 1, 0], [2, 1, 0]))) < NORM_TOL:
                self.normalization = 'L'
                return True
            else:
                return False


    def isrightnormalized(self,test=False):
        if test == False:
            if self.normalization == 'R':
                return True
            else:
                return False
        else:       
            if norm(identity(self.row_dim(), complex) - \
                    tensordot(self.tensor, self.tensor.conj(), axes=([3, 1, 0], [3, 1, 0]))) < NORM_TOL:
                self.normalization = 'R'
                return True
            else:
                return False


    def isnormalized(self, test=False):
        if self.isleftnormalized(test) or self.isrightnormalized(test):
            return True
        else:
            self.normalization = 'U'
            return False


# normalization methods

    def normalize_left(self, decomposition_method='qr', multiplied=True, 
                       truncation=0.0, bond_dim=BOND_MAX):

        if not decomposition_method in ['qr', 'svd']: 
            raise DecompositionError("Decomposition method must be either 'qr' or 'svd'.")
        if not isinstance(truncation, float):
            raise SingularValueTruncationError("Truncation limit not valid.")
        if (not isinstance(bond_dim, int)) or bond_dim < 1:
            raise SingularValueTruncationError("Dimensional limit not in reasonable range.")

        if decomposition_method == 'qr': # here with qr, faster, but less information 
            if multiplied == False:
                raise DecompositionError("QR returns only one matrix.")
            if truncation != 0.0:
                raise DecompositionError("QR cannot truncate nontrivially.")
            q, r = qr(self.tensor.reshape(self.site_bra_dim() * self.site_ket_dim() * self.row_dim(),
                                          self.col_dim())) 
            min_dim = min(self.site_bra_dim() * self.site_ket_dim() * self.row_dim(), self.col_dim()) 
            self.tensor = q.reshape(self.site_bra_dim(), self.site_ket_dim(), self.row_dim(), min_dim) # reshape back
            self.normalization = 'L'
            return r
        else: # now with svd, slower, full information 
            u, s, vd = svd(self.tensor.reshape(self.site_bra_dim() * self.site_ket_dim() * 
                                               self.row_dim(), self.col_dim()),0)
            min_dim = min(self.site_bra_dim() * self.site_ket_dim() * self.row_dim(), self.col_dim())
            self.tensor = u.reshape(self.site_bra_dim(), self.site_in_dim(), self.row_dim(), min_dim)
            self.normalization = 'L'
            # restrict to the (at most bond_dim) SV greater than truncation
            # truncation=0.0: discard SVs that are strictly zero (default)
            # truncation<0.0: discard none
            norm_factor = norm(s, 0)
            compress_list = [item > truncation * norm_factor for item in list(s)]
            if len(compress_list) > bond_dim:
                compress_list[bond_dim : ] = []
            if len(compress_list) < size(s, 0):
                s = compress(compress_list, s, 0)
                vd = compress(compress_list, vd, 0)
            if multiplied == True:
                for k in range(size(vd, 0)):
                    vd[k, : ] *= s[k]
                return vd
            else:
                return s, vd
                

    def normalize_right(self, decomposition_method='qr', multiplied=True, 
                        truncation=0.0, bond_dim=BOND_MAX):
        if not decomposition_method in ['qr', 'svd']: 
            raise DecompositionError("Decomposition method must be either 'qr' or 'svd'.")
        if (not isinstance(truncation,float)) or truncation < 0.0:
            raise SingularValueTruncationError("Truncation limit not in reasonable range.")
        if (not isinstance(bond_dim, int)) or bond_dim < 1:
            raise SingularValueTruncationError("Dimensional limit not in reasonable range.")
        if decomposition_method == 'qr': # here with qr, faster, but less information
            if multiplied == False:
                raise DecompositionError("QR returns only one matrix.")
            if truncation != 0.0:
                raise DecompositionError("QR cannot truncate nontrivially.")
            q, r = qr(self.tensor.\
                   transpose(0, 1, 3, 2).conj().\
                   reshape(self.site_bra_dim() * self.site_ket_dim() * self.col_dim(), self.row_dim())) 
            min_dim = min(self.row_dim(), self.site_bra_dim() * self.site_ket_dim() * self.col_dim())
            self.tensor = q.conj().reshape(self.site_bra_dim(), self.site_ket_dim(), self.col_dim(), 
                                           min_dim).transpose(0, 1, 3, 2)
            self.normalization = 'R'
            return r.transpose(1, 0).conj()
        else: # now with svd, slower, full information 
            u, s, vd = svd(self.tensor.\
                           transpose(2, 0, 1, 3).\
                           reshape(self.row_dim(), 
                                   self.site_bra_dim() * self.site_ket_dim() * self.col_dim()), 0)
            min_dim = min(self.row_dim(), self.site_bra_dim() * self.site_ket_dim() * self.col_dim())
            self.tensor = vd.reshape(min_dim, self.site_bra_dim(), 
                                     self.site_ket_dim(), self.col_dim()).transpose(1, 2, 0, 3)
            self.normalization = 'R'
            # restrict to the (at most bond_dim) SV greater than truncation
            # truncation=0.0: discard SVs that are strictly zero (default)
            # truncation<0.0: discard none
            norm_factor = norm(s, 0)
            compress_list = [item > truncation * norm_factor for item in list(s)]
            if len(compress_list) > bond_dim:
                compress_list[bond_dim :] = []
            if len(compress_list) < size(s, 0):
                s = compress(compress_list, s, 0)
                u = compress(compress_list, u, 1)
            if multiplied == True:
                for k in range(size(u, 1)):
                    u[:,k] *= s[k]
                return u
            else:
                return u,s


# multiply matrices to a set of M^{s,t} from left or right
# attention: tensordot orders remaining dimensions from 'left to right', may force transpose
# new normalization status unknown
            

    def multiply_from_left(self, mat):
        if size(mat, 1) != self.row_dim():
            raise DimensionMismatchError("Column and row dimensions do not match.")
        self.tensor = tensordot(mat, self.tensor, axes=([1], [2])).transpose(1,2,0,3)
        self.normalization = 'U'
        return
    

    def multiply_from_right(self, mat):
        if size(mat, 0) != self.col_dim():
            raise DimensionMismatchError("Column and row dimensions do not match.")
        self.tensor = tensordot(self.tensor, mat, axes=([3],[0]))
        self.normalization = 'U'
        return
                

    def multiply_by_scalar(self, scalar):
        self.tensor *= scalar
        self.normalization = 'U'
        return
    

# copy method

    def copy(self):
        _acopy_ = MPOTensor(self.site_bra_dim(), self.site_ket_dim(),\
                            self.row_dim(), self.col_dim(),\
                            self.init_random, self.normalization)
        _acopy_.tensor = self.tensor.copy()
        return _acopy_


# consistency checks for 2 MPOTensors
# mpo_tensor_ismatched checks in bond space
# mpo_tensor_physical_ismatched checks in physical space
    
    def mpo_tensor_ismatched(self, set2):

# do the dimensions of adjoining local tensors match for legs to be contracted?
#    
# arguments:
# self, set2: two MPOTensors, self to left of set 2
#    
# returns bool
    
        if self.col_dim() == set2.row_dim():
            return True
        else:
            return False
    
    
    def mpo_tensor_physical_ismatched(self, set2):
        
# do the physical dimensions of two local tensors match for contraction?
# 
# arguments:
# self, set2: two MPOTensors, self to left of set 2
#
# returns bool
     
        if self.site_dim() == set2.site_dim():
            return True
        else:
            return False

        
# last, but not least: return a dummy local MPO tensor.         
        
    def mpo_tensor_dummy():
        _dummy_ = MPOTensor((1, 1, 1, 1))
        _dummy_.tensor[0, 0, 0, 0] = 1.0
        return _dummy_
        

# ###############################################################################

# IMPORTANT: for notational consistency, contracted objects growing from RIGHT
# contain the COMPLEX CONJUGATE in deviation from review notation

def overlap_mpo_left_step(bra_tensor, ket_tensor, mpo_tensor, c):
    """
    overlap with sandwiched MPO: 3-tensor c carried one site from left
    ------------------------------------------------------------------

    arguments:
    bra_tensor: MPSTensor from bra - will be complex conjugated
    ket_tensor: MPSTensor from ket
    mpo_tensor: MPOTensor from MPO
    c: overlap 3-tensor to be carried forward
    
    output:
    c3: overlap 3-tensor carried one site forward
    see Eqs. (197)
    """
    
    # error handling
    
    if bra_tensor.site_dim() != mpo_tensor.site_bra_dim():
        raise MatrixProductStateDimensionError("Local physical dimensions of bra and operator do not match.")
    if ket_tensor.site_dim() != mpo_tensor.site_ket_dim():
        raise MatrixProductStateDimensionError("Local physical dimensions of ket and operator do not match.")
    if size(c, 0) != bra_tensor.row_dim() or\
       size(c, 1) != mpo_tensor.row_dim() or\
       size(c, 2) != ket_tensor.row_dim():
        raise MatrixProductStateDimensionError("Dimensions of overlap matrix and local matrix sets do not match.")
           
    c1 = tensordot(c, ket_tensor.tensor, (2, 1))
    c2 = tensordot(mpo_tensor.tensor,c1, ([2, 1], [1, 2]))
    c3 = tensordot(bra_tensor.tensor.conj(), c2, ([1, 0], [2, 0]))
    
    return c3


def overlap_mpo_right_step(bra_tensor, ket_tensor, mpo_tensor, c):
    """
    overlap with sandwiched MPO: 3-tensor c carried one site from right
    -------------------------------------------------------------------

    arguments:
    bra_tensor: MPSTensor from bra - will be complex conjugated
    ket_tensor: MPSTensor from ket
    mpo_tensor: MPOTensor from MPO
    c: overlap 3-tensor to be carried forward
    
    output:
    c3: overlap 3-tensor carried one site forward
    see Eqs. (197) adapted to growth from right and the notational extension
    """
    
    # error handling
    
    if bra_tensor.site_dim() != mpo_tensor.site_bra_dim():
        raise MatrixProductStateDimensionError("Local physical dimensions of bra and operator do not match.")
    if ket_tensor.site_dim() != mpo_tensor.site_ket_dim():
        raise MatrixProductStateDimensionError("Local physical dimensions of ket and operator do not match.")
    if size(c, 0) != bra_tensor.col_dim() or\
       size(c, 1) != mpo_tensor.col_dim() or\
       size(c, 2) != ket_tensor.col_dim():
        raise MatrixProductStateDimensionError("Dimensions of overlap matrix and local matrix sets do not match.")
           
    c1 = tensordot(c,ket_tensor.tensor.conj(), (2, 2))
    c2 = tensordot(mpo_tensor.tensor.conj(), c1, ([3, 1], [1, 2]))
    c3 = tensordot(bra_tensor.tensor, c2, ([2, 0],[2, 0]))
    
    return c3

                 
# ########################################################################## 

class MPO:
    """
    matrix product operator
    -----------------------
    
    structure: Eq. (177), Sec. 5
    
    * self.prefactor: complex scalar for scaling normalized AAAAA or BBBBB structure
    * self.w: list of (length+1) MPOTensors describing W^{s,t}_{a,b}, where length is the number of sites
    
    attributes:
    self.length(): length of system described by MPS
    self.site_bra_dim(i): dimension of bra physical state space at site i
    self.site_ket_dim(i): dimension of ket physical state space at site i
    self.row_dim(i): row dimension of matrices at site i
    self.col_dim(i): column dimension of matrices at site i
    self.site_bra_dims(): list of dimensions of bra physical state spaces
    self.site_ket_dims(): list of dimensions of ket physical state spaces
    self.row_dims(): list of row dimensions of matrices
    self.col_dims(): list of column dimensions of matrices
    self.ismixedcanonical(center,test): checks AAAAAXBBBBB form, where X is on site (center).
    self.isleftcanonical(test): checks L-normalization of MPO in analogy to MPS
    self.isrightcanonical(test): checks R-normalization of MPO in analogy to MPS
    self.ismixednormalized(center,test): as ismixedcanonical, and checks norm of X.
    self.isleftnormalized(test): checks L-normalization of MPO, and self.prefactor==1
    self.isrightnormalized(test): checks R-normalization of MPO, and self.prefactor==1
        test=True entails real verification
    self.copy(): deep copy
    self.canonize_left_step(i,test): if site i is not L-normalized (test!?), 
        make it so, return remaining matrix
    self.canonize_right_step(i,test): if site i is not R-normalized (test!?), 
        make it so, return remaining matrix
    self.canonize(center,test): L-normalize before site (center) and R-normalize after it
    self.canonize_left(test): L-normalize all sites
    self.canonize_right(test): R-normalize all sites
    self.normalize(center,test): L-normalize up to site (center-1), R-normalize from (center+1) to end
        and normalize center-matrices as a normal wave function (because the L- and R- span ONB)
    self.normalize_left(test): like canonize_left, but set self.prefactor=1.0
    self.normalize_right(test): like canonize_right, but set self.prefactor=1.0
    self.compress_svd_left_site(i,sv_min,bond_dim,test): 
        checks if operator in form AAAAAXBBB, X center, takes it to AAAAAAXBB, truncating as indicated below
    self.compress_svd_right_site(i,sv_min,bond_dim,test):
        checks if operator in form AAAAAXBBB, X center, takes it to AAAAXBBBB, truncating as indicated below
    self.compress_svd_left(sv_min,bond_dim,test): 
        starting from left on R-canonical state (!), makes L-canonical and truncates      
    self.compress_svd_right(sv_min,bond_dim,test): 
        starting from right on L-canonical state (!), makes R-canonical and truncates
        bond_dim_max=BOND_MAX,sv_min=0.0,test=False: truncates down to the at most bond_dim_max singular values
        that are larger than sv_min. test=True if full test of matrix normalization desired
    """    
    
# initialization

    def __init__(self):

        self.prefactor = 1.0
        self.w = []
        self.w.append(MPOTensor(1, 1, 1, 1))
        
            
# simple queries

    def length(self):
        return len(self.w) - 1
    
    
    def site_bra_dim(self, i):
        return size(self.w[i].tensor, 0)

    
    def site_ket_dim(self, i):
        return size(self.w[i].tensor, 1)
    
  
    def row_dim(self, i):
        return size(self.w[i].tensor, 2)
    
 
    def col_dim(self, i):
        return size(self.w[i].tensor, 3)
    

    def site_bra_dims(self):
        return [size(item.tensor, 0) for item in self.w]

    
    def site_ket_dims(self):
        return [size(item.tensor, 1) for item in self.w]
        

    def row_dims(self):
        return [size(item.tensor, 2) for item in self.w]


    def col_dims(self):
        return [size(item.tensor, 3) for item in self.w]
    

    def ismixedcanonical(self, center, test=False):
        for i in range(1, min(center, self.length() + 1)):
            if not self.w[i].isleftnormalized(test):
                return False
        for i in range(max(1, center+1), self.length() + 1):
            if not self.w[i].isrightnormalized(test):
                return False
        return True


    def isleftcanonical(self, test=False):
        return self.ismixedcanonical(self.length() + 1,test)
    

    def isrightcanonical(self, test=False):
        return self.ismixedcanonical(0, test)
    

    def isnormalized(self, center, test=False):
        if not self.ismixedcanonical(center, test):
            return False
        else:
            if not center in range(1, self.length() + 1):
                return abs(abs(self.prefactor) - 1.0) < NORM_TOL
            else:
                return abs(norm(self.w[center].tensor) - 1.0) < NORM_TOL 

    def isleftnormalized(self, test=False):
        return self.isnormalized(self.length() + 1, test)


    def isrightnormalized(self, test=False):
        return self.isnormalized(0, test)


# copy method

    def copy(self):
        _psicopy_ = MPO()
        _psicopy_.prefactor = self.prefactor
        for i in range(1, len(self.w)):
            _psicopy_.w.append(self.w[i])
            _psicopy_.w[i] = self.w[i].copy()
        return _psicopy_

# transform to left-normalization at site i (unless already done)
# AAAAAX???->AAAAAA???
# changes also content of site i+1, destroying normalization
# if i+1==length+1 then adjust prefactor


    def canonize_left_step(self, i, test=False):
        if self.w[i].isleftnormalized(test) == False:
            r = self.w[i].normalize_left()
            if i != self.length():
                self.w[i+1].multiply_from_left(r)
                self.w[i+1].normalization = 'U'
            else:
                self.prefactor *= r[0, 0]
        return
                
# transform to right-normalization at site i (unless already done)
# ???XBBBBB->???BBBBBB
# changes also content of site i-1, destroying normalization
# if i==1 then adjust prefactor


    def canonize_right_step(self, i, test=False):
        if self.w[i].isrightnormalized(test)==False:
            r = self.w[i].normalize_right()
            if i != 1:
                self.w[i-1].multiply_from_right(r)
                self.w[i-1].normalization = 'U'
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
        self.canonize(self.length() + 1,test)
        return


# transform to canonical representation with R-normalization, Sec. 4.4.2
# in the end, leftover factor multiplied into self.prefactor

    def canonize_right(self, test=False):
        self.canonize(0, test)
        return
         

# normalize from left or right is like canonize, 
# but with the prefactor forced to be 1 (if L- or R-canonical) or center matrices normalized
# then the state has norm 1

    def normalize(self, center, test=False):
        if not self.ismixedcanonical(center, test):
            self.canonize(center, test)
        if not center in range(1, self.length() + 1): # left or right-canonical
            if self.prefactor != 0.0:
                self.prefactor = 1.0
                return
            else:
                raise NormalizationError("State not normalizable.")
        else: # mixed-canonical
            norm_factor = norm(self.w[center].tensor)
            if norm_factor != 0.0 and self.prefactor != 0.0:
                norm_factor = 1.0 / norm_factor
                self.w[center].tensor *= normfactor
                self.prefactor = 1.0
                return
            else:
                raise NormalizationError("State not normalizable.")
            

    def normalize_left(self, test=False):
        self.normalize(self.length() + 1, test)
        return
    

    def normalize_right(self, test=False):
        self.normalize(0, test)
        return

    
# compression by svd; Sec. 4.5.1
# on site X, if AAAAAXBBB form: L-normalize and compress X and change site to right: AAAAAA(RB)BB

    def compress_svd_left_site(self, i, sv_min=0.0, bond_dim=BOND_MAX, test=False):
    
        if not self.ismixedcanonical(i, test):
            raise CanonicalFormMismatchError("Operator is not in right canonical form for SVD compression.")
        
        r = self.w[i].normalize_left('svd', True, sv_min, bond_dim)        
        if i != self.length():
            self.w[i+1].multiply_from_left(r)
        else:
            self.prefactor *= r[0,0]
        return


# on site X, if AAAAAXBBB form: R-normalize and compress X and change site to left: AAAA(AR)BBBB

    def compress_svd_right_site(self, i, sv_min=0.0, bond_dim=BOND_MAX, test=False):
    
        if not self.ismixedcanonical(i, test):
            raise CanonicalFormMismatchError("Operator is not in right canonical form for SVD compression")
        
        r=self.w[i].normalize_right('svd', True, sv_min, bond_dim)        
        if i != 1:
            self.w[i-1].multiply_from_right(r)
        else:
            self.prefactor *= r[0,0]
        return
    

    def compress_svd_left(self, sv_min=0.0, bond_dim=BOND_MAX, test=False):
        
        for i in range(1,self.length()+1):
            self.compress_svd_left_site(i, sv_min, bond_dim, test)



    def compress_svd_right(self, sv_min=0.0, bond_dim=BOND_MAX, test=False):
        
        for i in range(self.length(),0,-1):
            self.compress_svd_right_site(i, sv_min, bond_dim, test)
            




