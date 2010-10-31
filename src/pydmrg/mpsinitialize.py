from numpy import *
from constants import *
from myexceptions import *
from mps import *

# This file contains general routines to initialize matrix product states.
# Three cases are covered:
#
#     1) mps_from_dimensions(site_dims, row_dims, col_dims, init_random=False)
#     ------------------------------------------------------------------------ 
#     Dimensions of the local physical state spaces and the bond state spaces
#     are available, corresponding to the row and column dimensions of the local matrices.
#     Initialization is to 0.0, unless init_random=True.
#
#     2) mps_from_site_dimensions(site_dims, bond_dim, init_random=False)
#     -------------------------------------------------------------------
#     Dimensions of the local physical state spaces are available and used to construct 
#     corresponding MPS matrices of growing dimensions, bounded by a maximum bond dimension.
#     Initialization is to 0.0, unless init_random=True.
#
#     3) mps_from_tensor_list(tensor_list, prefactor=1.0)
#     -----------------------------------------------
#     A list of MPS tensors is combined to form an MPS. A prefactor to the MPS can be given.
#
#     IMPORTANT: All lists (site_dims, row_dims, col_dims, tensor_list) must contain a 
#     dummy element list[0], to get in line with the standard physical notation of labelling
#     the first site of a chain by 1, and the output by other routines. This avoids counterintuitive
#     notation elsewhere.
#
#     Moreover, a number of further generation routines for specific useful quantum states are
#     provided.

def mps_from_dimensions(site_dims, row_dims, col_dims, init_random=False):
    """
    generates random- or zero-initialized MPS from dimensional lists
    
    arguments:
    site_dims: list of local physical dimensions (int), implies length, with dummy entry on 0
    row_dims:  list of local matrix row dimensions (int)
    col_dims:  list of local matrix column dimensions (int)
    init_random: False (default), True: initialize with zeros (default) or random complex numbers
    
    output:
    returns MPS class object
    """
    
    # consistency and error checks
    
    if not isinstance(site_dims, list):
        raise MatrixProductStateDimensionError("Site dimensions do not form a list.")
    if not isinstance(row_dims, list):
        raise MatrixProductStateDimensionError("Row dimensions do not form a list.")
    if not isinstance(col_dims, list):
        raise MatrixProductStateDimensionError("Column dimensions do not form a list.")
    if len(site_dims) != len(row_dims):
        raise MatrixProductStateDimensionError("Lengths of site and row dimension lists do not match.")
    if len(site_dims) != len(col_dims):
        raise MatrixProductStateDimensionError("Lengths of site and column dimension lists do not match.")
    for i in range(1, len(site_dims)):
        if not isinstance(site_dims[i], int):
            raise MatrixProductStateDimensionError("Site dimension is not an integer.")
        if not isinstance(row_dims[i], int):
            raise MatrixProductStateDimensionError("Row dimension is not an integer.")
        if not isinstance(col_dims[i], int):
            raise MatrixProductStateDimensionError("Column dimension is not an integer.")
    for i in range(1, len(site_dims) - 1):
        if col_dims[i] != row_dims[i+1]:
            raise MatrixProductStateDimensionError("Dimensional mismatch of row and column dimensions.")
    if row_dims[1] != 1:
            raise MatrixProductStateDimensionError("Left boundary not consistent with OBC.")
    if col_dims[-1] != 1:
            raise MatrixProductStateDimensionError("Right boundary not consistent with OBC.")
        
    # actual construction of MPS
        
    psi = MPS()
    for i in range(1, len(site_dims)):
        psi.m.append(MPSTensor(site_dims[i], row_dims[i], col_dims[i], init_random, 'U'))
    return psi



# generate the dimensional structure for matrices in MPS


def mps_matrix_sizes(site_dims, bond_dim):
    """
    generates matrix sizes for MPS
    
    arguments:
    site_dims: list of local physical dimensions from sites 0 (dummy), 1 through length,
        to tie in with normal physical notation
    bond_dim: cut on matrix dimensions
    
    outcome:
    returns two length consistent lists:
    list of row dimensions of matrices (1 on dummy site),
    list of column dimensions of matrices (1 on dummy site)     
    """
    
    # consistency and error checks 
    
    if (not isinstance(bond_dim, int)) or bond_dim < 1:
        raise MatrixProductStateDimensionError("Bond dimension is not a positive integer.")
    if bond_dim > BOND_MAX:
        raise MatrixProductStateDimensionError("Bond dimension is in excess of maximal allowed value.")        
    if (not isinstance(site_dims, list)) or len(site_dims) < 3:
        raise MatrixProductStateDimensionError("Site dimensions not given by list with at least 3 entries.")
    for i in range(1, len(site_dims)):
        if (not isinstance(site_dims[i], int)) or site_dims[i] < 1:
            raise MatrixProductStateDimensionError("Site dimension is not a positive integer.")
    
    length = len(site_dims) - 1   # the dummy subtracted
    site_dims[0] = 1
    row_dims = [1]   # site 0
    col_dims = []   
    row_dim = 1
    for i in range(0, length):
        row_dim *= site_dims[i]
        if row_dim < bond_dim:
            row_dims.append(row_dim)
            col_dims.append(row_dim)   # one site to the left!
        else:
            row_dims = row_dims + (length - i) * [bond_dim]
            col_dims = col_dims + (length - i + 1) * [bond_dim]
            break
    row_dim = 1
    for i in range(length, 0, -1):
        row_dim *= site_dims[i]
        if row_dim < row_dims[i]:
            row_dims[i] = row_dim
            col_dims[i-1] = row_dim
        else:
            break
    col_dims[length] = 1
    return row_dims, col_dims


def mps_from_site_dimensions(site_dims, bond_dim, init_random=False):
    """
    generates random- or zero-initialized MPS from site dimensions and maximum bond length
    
    arguments:
    site_dims: list of local physical dimensions (int), implies length, with dummy entry on 0
    bond_dim: maximal allowed bond dimension
    init_random: False (default), True: initialize with zeros (default) or random complex numbers
    
    output:
    returns MPS class object
    """
            
    # actual construction of MPS
        
    row_dims, col_dims = mps_matrix_sizes(site_dims, bond_dim)
    psi = MPS()
    for i in range(1,len(site_dims)):
        psi.m.append(MPSTensor(site_dims[i], row_dims[i], col_dims[i], init_random, 'U'))
    return psi


def mps_from_tensor_list(tensor_list, prefactor=1.0):
    """
    generates MPS from list of tensors and a prefactor
    
    arguments:
    tensor_list: list of MPS tensors, with arbitrary first tensor (will be ignored, site 0!)
    prefactor: prefactor to store normalization
    
    output:
    returns MPS class object
    """
    
    # error handling
    
    if not isinstance(tensor_list, list):
        raise MPSTensorError("Tensor list must be list.")
    for i in range(1, len(tensor_list)):
        if not isinstance(tensor_list[i], MPSTensor):
            raise MPSTensorError("Tensor list must be made from MPSTensor objects.")
        
# consistency checks of matrices in list

    if not isobccompatible(tensor_list[1], 'L'):
        raise MatrixProductStateDimensionError("OBC violated at first site")
    if not isobccompatible(tensor_list[-1], 'R'):
        raise MatrixProductStateDimensionError("OBC violated at last site")
    for i in range(1, -2):
        if not mps_tensor_ismatched(tensor_list[i], tensor_list[i+1]):
            raise MatrixProductStateDimensionError("dimensional inconsistency in tensor list")
        
    psi = MPS()
    psi.prefactor = prefactor
    for i in range(1, len(site_dims)):
        psi.m.append(tensor_list[i])
    return psi


# generation of specific, simple MPS that are frequently useful

def mps_aklt(length, deg1, deg2):
    """
    generates OBC AKLT state (MPS state of bond dimension 2) 
    
    arguments:
    length: length of the AKLT state
    deg1: 0, 1: left boundary state 
    deg2: 0, 1: right boundary state 
        the four resulting states are orthonormal in the TD limit
        
    output:
    returns AKLT state in MPS form
    """

    # error handling
    
    if (not isinstance(length, int)) or length < 2:
        raise MatrixProductStateDimensionError("AKLT state length must be 2 or more.")
    if not deg1 in [0, 1]:
        raise MatrixProductStateDimensionError("Left boundary state must be fixed by values 0, 1.")
    if not deg2 in [0, 1]:
        raise MatrixProductStateDimensionError("Right boundary state must be fixed by values 0, 1.")

    site_dims = [1] + length * [3]   # dummy on site 0
    psi = mps_from_site_dimensions(site_dims, 2)

# the singlet bond and the (S=1/2)x(S=1/2)->(S=1) projection, Eq. (83), (85)
# M^{Sz=0} modified for the first and last site from (85) to have bulk conserved 
# while being properly normalized

    sigma = array([[0.0, +1.0/sqrt(2.0)], [-1.0/sqrt(2.0), 0.0]])
    m=array([[[1.0, 0.0], [0.0, 0.0]],\
             [[0.0, 1.0], [1.0, 0.0]],\
             [[0.0, 0.0], [0.0, 1.0]]])
    for k in range(3):
        for j in range(2):
            for l in range(2):
                psi.m[1].tensor[k, 0, j] += m[k, deg1, l] * sigma[l, j]
                psi.m[length].tensor[k, j, 0] += sigma[j, l] * m[k, l, deg2]
    for i in range(2, length):
        psi.m[i].tensor[0, 0, 1] = +sqrt(2.0) / sqrt(3.0)
        psi.m[i].tensor[1, 0, 0] = -1.0 / sqrt(3.0)
        psi.m[i].tensor[1, 1, 1] = +1.0 / sqrt(3.0)            
        psi.m[i].tensor[2, 1, 0] = -sqrt(2.0) / sqrt(3.0)        
    return psi


def mps_ghz(length):
    """
    generates GHZ state of size 'length' as MPS of dimension 2
    
    argument:
    length: length of the GHZ state
    
    output:
    returns GHZ state as MPS
    """
    
    # error handling

    if (not isinstance(length, int)) or length < 2:
        raise MatrixProductStateDimensionError("GHZ needs length 2 or more.")

    site_dims = [1] + length * [2]
    psi = mps_from_site_dimensions(site_dims, 2)

    b = array([[[1.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [0.0, 1.0]]])
    psi.m[1].tensor = array([[[1.0, 0.0]], [[0.0,1.0]]])
    psi.m[length].tensor = array([[[1.0], [0.0]], [[0.0], [1.0]]])
    for i in range(2, length):
        psi.m[i].tensor = b
    return psi


def mps_neel(length, startspin=0):
    """
    generates Neel state as MPS of dimension 1
    
    arguments:
    length: length of Neel state
    startspin: 0 (default), 1: first spin up (default), down
    
    output:
    returns Neel state as MPS
    """

    # error handling
    
    if (not isinstance(length,int)) or length < 2:
        raise MatrixProductStateDimensionError("Neel state needs length 2 or more.")
    if not startspin in [0, 1]:
        raise MatrixProductStateDimensionError("First spin must be set by 0 (up) or 1 (down).")
    
    site_dims = [1] + length * [2]
    psi = mps_from_site_dimensions(site_dims, 1)
    for pos in range(1, length+1):
        if pos % 2 == 1:
            psi.m[pos].tensor = array([[[1.0]], [[0.0]]])
        else:
            psi.m[pos].tensor = array([[[0.0]], [[1.0]]])
    return psi


def mps_domain(length, right=2, startspin=0):
    """
    generates domain wall state as MPS of bond dimension 1
    
    arguments:
    length: length of domain wall MPS
    right (default 2): location of first spin in right domain
    startspin 0 (default), 1: left part up (0), down (1)
    
    output:
    returns domain wall state as MPS
    """

    # error handling
    
    if (not isinstance(length, int)) or length < 2:
        raise MatrixProductStateDimensionError("Domain wall state needs length 2 or more.")
    if not startspin in [0, 1]:
        raise MatrixProductStateDimensionError("Left domain orientation must be 0 (up) or 1 (down).")
    if (not isinstance(right,int)):
        raise MatrixProductStateDimensionError("Starting point of right domain must be integer.")

    site_dims = [1] + length * [2]
    psi=mps_from_site_dimensions(site_dims, 1)
    for pos in range(1, length+1):
        if pos < right:
            psi.m[pos].tensor[startspin, 0, 0] = 1.0
        else:
            psi.m[pos].tensor[1-startspin, 0, 0] = 1.0
    return psi


def mps_singlet_ladder(length):
    """
    generates spin ladder from spin-1/2 as MPS of minimal bond dimensions
    singlets sit on sites (1,2), (3,4), (5,6), ...
    hence bond dimensions 2 on bond (1,2), ... but 1 on (2,3), ...
    i.e. row_dims 1 on odd and 2 on even sites, col_dims vice versa
    prepares finite-T evolution
    
    argument:
    length: number of spins, i.e. twice the number of rungs, must therefore be even
    
    output:
    spin ladder as an MPS
    """

    # error handling
    
    if (not isinstance(length,int)) or length < 2 or length % 2 != 0:
        raise MatrixProductStateDimensionError("Singlet ladder needs to be of even length.")

    site_dims = [1] + length * [2]
    row_dims = [1] + (length // 2) * [1, 2]
    col_dims = [1] + (length // 2) * [2, 1]
    psi = mps_from_dimensions(site_dims, row_dims, col_dims)

    sqrt2 = 1.0 / sqrt(sqrt(2.0))
    for i in range(1, length + 1):
        if i % 2 == 1:
            psi.m[i].tensor = array([[[sqrt2, 0.0]],[[0.0, -sqrt2]]])
        else:
            psi.m[i].tensor = array([[[0.0], [sqrt2]], [[sqrt2], [0.0]]])
    return psi

