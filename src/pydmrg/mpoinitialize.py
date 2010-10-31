from numpy import *
from numpy import tensordot
from constants import *
from myexceptions import *
from mps import *
from mpo import *
import physop

# This file contains general routines to initialize matrix product operators.
# Three cases are covered:
#
#     1) mpo_from_dimensions(site_dims, row_dims, col_dims, init_random=False)
#     ------------------------------------------------------------------------ 
#     Dimensions of the bra and ket local physical state spaces and the bond state spaces
#     are available, corresponding to the row and column dimensions of the local matrices.
#     Initialization is to 0.0, unless init_random=True.
#
#     2) mpo_from_tensor_list(tensor_list, prefactor=1.0)
#     -----------------------------------------------
#     A list of MPO tensors is combined to form an MPO. A prefactor to the MPO can be given.
#
#     IMPORTANT: All lists (site_dims, row_dims, col_dims, tensor_list) must contain a 
#     dummy element list[0], to get in line with the standard physical notation of labelling
#     the first site of a chain by 1, and the output by other routines. This avoids counterintuitive
#     notation elsewhere.
#
#     3) mpo_from_dictionary()
#     ------------------------
#
#     Moreover, a number of further generation routines for specific useful quantum states are
#     provided.


def mpo_from_dimensions(site_bra_dims, site_ket_dims, row_dims, col_dims, init_random=False):
    """
    generates random- or zero-initialized MPS from dimensional lists
    
    arguments:
    site_bra_dims: list of local bra physical dimensions (int), implies length, with dummy entry on 0
    site_ket_dims: list of local ket physical dimensions (int), implies length, with dummy entry on 0
    row_dims:      list of local matrix row dimensions (int)
    col_dims:      list of local matrix column dimensions (int)
    init_random:   False (default), True: initialize with zeros (default) or random complex numbers
    
    output:
    returns MPO class object
    """
    
    # consistency and error checks
    
    if not isinstance(site_bra_dims, list):
        raise MatrixProductOperatorDimensionError("Site bra dimensions do not form a list.")
    if not isinstance(site_ket_dims, list):
        raise MatrixProductOperatorDimensionError("Site ket dimensions do not form a list.")
    if not isinstance(row_dims, list):
        raise MatrixProductOperatorDimensionError("Row dimensions do not form a list.")
    if not isinstance(col_dims, list):
        raise MatrixProductOperatorDimensionError("Column dimensions do not form a list.")
    if len(site_bra_dims) != len(site_ket_dims):
        raise MatrixProductOperatorDimensionError("Lengths of bra and ket site dimension lists do not match.")
    if len(site_bra_dims) != len(row_dims):
        raise MatrixProductOperatorDimensionError("Lengths of site and row dimension lists do not match.")
    if len(site_bra_dims) != len(col_dims):
        raise MatrixProductOperatorDimensionError("Lengths of site and column dimension lists do not match.")
    for i in range(1, len(site_bra_dims)):
        if not isinstance(site_bra_dims[i], int):
            raise MatrixProductOperatorDimensionError("Site bra dimension is not an integer.")
        if not isinstance(site_ket_dims[i], int):
            raise MatrixProductOperatorDimensionError("Site ket dimension is not an integer.")
        if not isinstance(row_dims[i], int):
            raise MatrixProductOperatorDimensionError("Row dimension is not an integer.")
        if not isinstance(col_dims[i], int):
            raise MatrixProductOperatorDimensionError("Column dimension is not an integer.")
    for i in range(1, len(site_dims) - 1):
        if col_dims[i] != row_dims[i+1]:
            raise MatrixProductOperatorDimensionError("Dimensional mismatch of row and column dimensions.")
    if row_dim[1] != 1:
            raise MatrixProductOperatorDimensionError("Left boundary not consistent with OBC.")
    if col_dim[len(site_dims) - 1] != 1:
            raise MatrixProductOperatorDimensionError("Right boundary not consistent with OBC.")
        
    # actual construction of MPS
        
    psi = MPO()
    for i in range(1, len(site_dims)):
        psi.m.append(MPOTensor(site_bra_dims[i], site_ket_dims[i], row_dims[i], col_dims[i], 'U', init_random))
    return psi


# from dictionaries of operators {1 : physop.sps, 4: physop.sms} - equaling S^+_1 S^-_4 for generic spin -
# generate MPO; dictionaries can be unordered

# shortcuts to often-used n-point correlators

# expectation values are coded by dictionaries:
# key is site, the result is the operator matrix


def sz(i):
# Dictionary for Sz(i) with S=1/2
    result={}
    result[i]=physop.sz
    return result


def szsz(i,j):
# Dictionary for Sz(i)Sz(j) with S=1/2
    result={}
    result[i]=physop.sz
    result[j]=physop.sz
    return result


def spsm(i,j):
# Dictionary for S+(i)S-(j) with S=1/2
    result={}
    result[i]=physop.sp
    result[j]=physop.sm
    return result


def sz1sz1(i,j):
# Dictionary for Sz(i)Sz(j) with S=1
    result={}
    result[i]=physop.sz1
    result[j]=physop.sz1
    return result


def sp1sm1(i,j):
# Dictionary for S+(i)S-(j) with S=1
    result={}
    result[i]=physop.sp1
    result[j]=physop.sm1
    return result


def strcor(i,j):
# Dictionary for Sz(i) exp(j*pi*(Sz(i+1)+...+Sz(j-1)))Sz(j) with S=1
    result={}
    result[i]=physop.sz1
    result[j]=physop.sz1
    for k in range(i+1,j):
        result[k]=physop.fl1
    return result


# generation of MPO, need the dictionary of operators, but also bra, ket site_dims, otherwise identities
# undefined; error if on non-trivial sites dimensions do not match

def mpo_from_dictionary(opdict, site_bra_dims, site_ket_dims):
    """
    generates an MPO from a dictionary of operators
    -----------------------------------------------
    
    arguments:
    opdict: dictionary of operators of type {1 : physop.sps, 4: physop.sms} 
        equaling S^+_1 S^-_4 for generic spin
    site_bra_dims: list of local bra physical dimensions 
    site_ket_dims: list of local ket physical dimensions 
    
    output:
    MPO containing operators and matching bra and ket physical dimensions
    """
    
    # error handling
    
    if not isinstance(site_bra_dims, list):
        raise MatrixProductOperatorDimensionError("Site bra dimensions do not form a list.")
    if not isinstance(site_ket_dims, list):
        raise MatrixProductOperatorDimensionError("Site ket dimensions do not form a list.")
    if len(site_bra_dims) != len(site_ket_dims):
        raise MatrixProductOperatorDimensionError("Lists for bra and ket state local dimensions of inconsistent length.")
    for i in range(1, len(site_bra_dims)):
        if not isinstance(site_bra_dims[i], int):
            raise MatrixProductOperatorDimensionError("Site bra dimension is not an integer.")
        if not isinstance(site_ket_dims[i], int):
            raise MatrixProductOperatorDimensionError("Site ket dimension is not an integer.")
    
    col_dims = row_dims = len(site_bra_dims) * [1]     # MPO factorizes
    
    ops = mpo_from_dimensions(site_bra_dims, site_ket_dims, row_dims, col_dims)
    
    for i in range(1, len(site_bra_dims)):
        if i in opdict.keys():
            if size(opdict[i], 0) != site_bra_dims[i] or size(opdict[i],1) != site_ket_dims[i]:
                raise MatrixProductOperatorDimensionError("bra and ket local dimensions do not match operator.")
            ops.w[i].tensor = expand_dims(expand_dims(opdict[i], 2), 3) # make op(s,t) -> w(s,t,0,0)
        else:
            if site_bra_dims[i] != size(opdict[i], 1):
                raise MatrixProductOperatorDimensionError("bra and ket local dimensions must be equal on identity sites.")
            ops.w[i].tensor = expand_dims(expand_dims(identity(site_bra_dims[i],complex), 2), 3)
    
    return ops


def tensor_list_hafm(length, J, Jz, h=0.0):

    opsarray1 = zeros((1, 5, 2, 2), complex)
    opsarray = zeros((5, 5, 2, 2), complex)
    opsarrayL = zeros((1, 5, 2, 2), complex)
    
    opsarray1[0, 0] = h * physop.sz.copy()
    opsarray1[0, 1] = physop.sm.copy()
    opsarray1[0, 2] = physop.sp.copy()
    opsarray1[0, 3] = physop.sz.copy()
    opsarray1[0, 4] = physop.one2.copy()
    
    opsarrayL[0, 0] = physop.one2.copy()
    opsarrayL[1, 0] = J * 0.5 * physop.sp.copy()
    opsarrayL[2, 0] = J * 0.5 * physop.sm.copy()
    opsarrayL[3, 0] = Jz * physop.sz.copy()
    opsarrayL[4, 0] = h * physop.sz.copy()

    opsarray[4, 0] = h * physop.sz.copy()
    opsarray[4, 1] = physop.sm.copy()
    opsarray[4, 2] = physop.sp.copy()
    opsarray[4, 3] = physop.sz.copy()
    opsarray[4, 4] = physop.one2.copy()
    opsarray[0, 0] = physop.one2.copy()
    opsarray[1, 0] = J * 0.5 * physop.sp.copy()
    opsarray[2, 0] = J * 0.5 * physop.sm.copy()
    opsarray[3, 0] = Jz * physop.sz.copy()
    
    opsarray1 = opsarray1.transpose(2, 3, 0, 1)
    opsarray = opsarray.transpose(2, 3, 0, 1)
    opsarrayL = opsarrayL.transpose(2, 3, 0, 1)
    
    tensor_list = []
    tensor_list.append(mpo_tensor_dummy)
    tensor_list.append(MPOTensor(2, 2, 1, 5))
    tensor_list.w[1].tensor = opsarray1.copy()
    for i in range(2, length):
        tensor_list_append(MPOTensor(2, 2, 5, 5))
        tensor_list.w[i].tensor = opsarray.copy()
    tensor_list.append(MPOTensor(2, 2, 5, 1))
    tensor_list.w[length].tensor = opsarrayL.copy()
        
    return tensor_list


def mpo_from_tensor_list(tensor_list):

    
    if tensor_list[1].tensor.row_dim() != 1:
        raise MatrixProductOperatorDimensionError("MPO dimension on site 1 not compatible with OBC.")
    if tensor_list[len(tensor_list)].tensor.col_dim() != 1:
        raise MatrixProductOperatorDimensionError("MPO dimension on site L not compatible with OBC.")
    for i in range(1, len(tensor_list) - 1):
        if tensor_list[i].tensor.col_dim() != tensor_list[i + 1].tensor.row_dim():
            raise MatrixProductOperatorDimensionError("MPO dimensions on bonds do not match.")
                   
    mpo = MPO()
    for i in range(1, len(tensor_list)):
        mpo.w.append(tensor_list[i])
    return mpo


# bond Hamiltonians coded as explicit complex matrices


def bond_hafm(J, Jz, h):

    # Heisenberg (an)isotropic AFM: sum_i J/2(S+(i)S-(i+1)+Hc)+ Jz Sz(i)Sz(i+1) + h Sz(i) for S=1/2
    # field is split on 'inner' sites into two
    
    hamiltonian_list=[]
    ham = zeros((4, 4), complex)
    ham1 = zeros((4, 4), complex)
    hamL = zeros((4, 4), complex)
    ham += Jz * outer(physop.sz, physop.sz)
    ham += J * 0.5 * outer(physop.sp, physop.sm)
    ham += J * 0.5 * outer(physop.sm, physop.sp)
    ham1 = ham.copy()
    hamL = ham.copy()
    ham += h * 0.5 * (outer(physop.sz, physop.ones2) + outer(physop.ones2, physop.sz))
    ham1 += h * (outer(physop.sz, physop.ones2) + 0.5 * outer(physop.ones2, physop.sz))
    hamL += h * (0.5 * outer(physop.sz, physop.ones2) + outer(physop.ones2, physop.sz))
    ham = ham.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    ham1 = ham1.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    hamL = hamL.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    hamiltonian_list.append(ham1)
    hamiltonian_list.append(ham)
    hamiltonian_list.append(hamL)
    return hamiltonian_list


def bond_hafm_spin1(J, Jz, h):
# Heisenberg (an)isotropic AFM: sum_i J/2(S+(i)S-(i+1)+Hc)+Sz(i)Sz(i+1) + h Sz(i) for S=1

    hamiltonian_list = []
    ham = zeros((9,9),complex)
    ham1 = zeros((9,9),complex)
    hamL = zeros((9,9),complex)
    ham += Jz * outer(physop.sz1, physop_sz1)
    ham += J * 0.5 * outer(physop_sp1, physop_sm1)
    ham += J * 0.5 * outer(physop_sm1, physop_sp1)
    ham1 = ham.copy()
    hamL = ham.copy()
    ham += h * 0.5 * (outer(physop.sz1, physop.ones3) + outer(physop.ones3, physop.sz1))
    ham1 += h * (outer(physop.sz1, physop.ones3) + 0.5 * outer(physop.ones3, physop.sz1))
    hamL += h * (0.5 * outer(physop.sz1, physop.ones3) + outer(physop.ones3, physop.sz1))
    ham = ham.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
    ham1 = ham1.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
    hamL = hamL.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
    hamiltonian_list.append(ham1)
    hamiltonian_list.append(ham)
    hamiltonian_list.append(hamL)
    return hamiltonian_list



def bond_bilinbiquad1(lin, quad):
# bilinear-biquadratic spin chain for S=1: sum_i lin*S(i)S(i+1)+ quad*(S(i)S(i+1))**2

    hamiltonian_list = []
    ham = zeros((9, 9),complex)
    ham1 = zeros((9, 9),complex)
    hamL = zeros((9, 9),complex)
    ham += Jz * outer(physop_sz1, physop_sz1)
    ham += J * 0.5 * outer(physop_sp1, physop_sm1)
    ham += J * 0.5 * outer(physop_sm1, physop_sp1)
    ham = ham.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
    ham2 = dot(ham, ham)
    ham = ham * lin + ham2 * quad
    ham1 = ham.copy()
    hamL = ham.copy()
    hamiltonian_list.append(ham1)
    hamiltonian_list.append(ham)
    hamiltonian_list.append(hamL)
    return hamiltonian_list


def unitary_from_bond(hamiltonian_list, timestep):
    
    unitary_list = []
    for ham in hamiltonian_list:
        size2 = ham.shape[0]
        size = int(sqrt(size2))
        eigvals, eigvecs = eigh(ham)
        expmat = diag(exp(-1.0j * timestep * eigvals))
        evolution_operator = dot(dot(eigvecs, expmat), eigvecs.T.conj())
        u, s, vh = svd(evolution_operator.reshape(size, size, size, size).\
                       transpose(0,2,1,3).reshape(size2, size2),0)
        compress_list = [item > SV_TOL for item in list(s)]
        if len(compress_list) < size2:
            s = compress(compress_list, s, 0)
            u = compress(compress_list, u, 1)
            vh = compress(compress_list, vh, 0)
        s = sqrt(s)
        ops1 = MPOTensor(size, size, 1, len(compress_list))
        ops2 = MPOTensor(size, size, len(compress_list), 1)
        ops1.tensor = tensordot(u, s, (1, 0)).reshape(size, size, 1, len(compress_list)).copy()
        ops2.tensor = tensordot(s, vh, (0, 0)).reshape(size, size, len(compress_list), 1).copy()
        unitary_list.append(ops1, ops2)
        
    return unitary_list

