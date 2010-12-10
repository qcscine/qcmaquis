from numpy import *
from mps import *
from mpsinitialize import *
from mpo import *
from numpy import tensordot
import physop

# 2 states (maybe identical) and MPO sandwiched in between
# overlap <bra|MPO|ket> by zipping from left or right ->
# single steps of carrying zip one step forward from left or from right
# test for identity as local MPO
# routines: overlap_MPO_left_step, overlap_MPO_right_step

def overlap_mpo(bra, ket, mpo, withprefactors=True, direction='L'):        
    """
    calculates overlap <bra|mpo|ket>
    --------------------------------
    arguments:
    self: bra MPS
    ket: ket MPS
    mpo: MPO
    withprefactors=True: take *.prefactors into account
    direction='L': zip from 'L'eft or 'R'ight
    
    output:
    <bra|mpo|ket> 
    """
    if bra.length() != ket.length():
        raise MatrixProductStateDimensionError("Lengths of the operators do not match.")

    c=ones((1, 1),complex)
    if direction=='L':
        for i in range(1, ket.length()+1):
            c = overlap_mpo_step_left(bra.w[i], ket.w[i], c).copy() # is this copy() necessary? check!
    else:
        for i in range(ket.length(), 0, -1):
            c = overlap_mpo_step_right(bra.w[i], ket.w[i], c).copy() # is this copy() necessary? check!    
        c[0, 0] = conj(c[0, 0]) # complex conjugate because of construction
    if withprefactors == True:
        c[0, 0] *= conj(bra.prefactor) * ket.prefactor
    return c[0, 0]


# expectation values with proper normalization

def expectationvalue(psi, mpo):
    """
    expectation value
    -----------------

    arguments:
    psi: MPS
    mpo: MPO

    output: <psi|mpo|psi>/<psi|psi>
    """

    denominator = psi.norm2()
    if abs(denominator) == 0.0:
        raise NormalizationError("state not normalizable")
    numerator = overlap_mpo(psi, psi, mpo)
    result = numerator / denominator
    return result


# addition of two matrix product states


def mpo_times_mps_tensor(mpo_tensor,mps_tensor):
    """
    local MPO tensor applied from left to (ket) MPS tensor
    ------------------------------------------------------

    arguments:
    mpo_tensor: local MPO tensor
    mps_tensor: local MPS tensor

    output:
    product of them
    """

    site_dim = mpo_tensor.site_bra_dim()
    row_dim = mpo_tensor.row_dim() * mps_tensor.row_dim()
    col_dim = mpo_tensor.col_dim() * mps_tensor.col_dim()
    mpomps = MPSTensor(site_dim, row_dim, col_dim, normalization='U')
    mpomps.tensor = tensordot(mpo_tensor.tensor, mps_tensor.tensor, (1, 0)). \
          transpose(0, 1, 3, 2, 4).reshape(site_dim, row_dim, col_dim).copy()
    return mpomps


def mpo_times_mps(mpo, mps):
    """
    MPO applied from left to (ket) MPS
    ----------------------------------

    arguments:
    mpo: MPO
    mps: MPS

    output:
    resulting mps
    """

    # error handling

    if mpo.length() != mps.length():
        raise MatrixProductStateDimensionError("Length of MPO and MPS do not match.")
    if mpo.site_ket_dims() != mps.site_dims():
        raise MatrixProductStateDimensionError("Local (ket) dimensions of MPO and MPS do not match.")

    site_dims = mpo.site_bra_dims()
    row_dims = list(array(mpo.row_dims()) * array(mps.row_dims()))
    col_dims = list(array(mpo.col_dims()) * array(mps.col_dims()))
    psi=mps_from_dimensions(site_dims, row_dims, col_dims, prefactor = mpo.prefactor * mps.prefactor)
    for i in range(1, mps.length()+1):
        psi.m[i] = mpo_times_mps_local(mpo.w[i], mps.m[i]).copy()
    return psi


def mpo_trotter_from_unitaries(length, unitary_list, oddeven):

    opsunit = MPOTensor(unitary_list(1).site_bra_dim(), unitary_list(1).site_ket_dim(), 1, 1)
    for k in range(opsunit.site_bra_dim())::
        opsunit.tensor[k, k, 0, 0] = 1.0
    
    opslist = []
    opslist.append(MPOTensor(1, 1, 1, 1))  # dummy for site 0
    if oddeven == 1:
        for i in range(1, length, 2):
            opslist.append(ops1)
            opslist.append(ops2)
        if length % 2 == 1:
            opslist.append(opsunit) # at position length if odd length
    else:
        opslist.append(opsunit) # at position 1
        for i in range(2,length,2):
            opslist.append(ops1)
            opslist.append(ops2)
        if length % 2 == 0:
            opslist.append(opsunit)

    trottermpo = mpo_from_tensor_list(opslist)
    return trottermpo


def trotter_1_step(psiold, mpo_odd, mpo_even, sv_min=0.0, bond_dim=BOND_MAX):
    """
    first order Trotter time step
    ------------------------------
    
    arguments:
    psiold: MPS at beginning of time step
    mpo_odd:  MPO for odd bond time evolution for time step tau
    mpo_even: MPO for odd bond time evolution for time step tau
    sv_min: cutoff of small singular values in decomposition
    bond_dim: maximally allowed bond dimension
    """

    psi = mpo_times_mps(mpo_odd, psiold)
    psi = mpo_times_mps(mpo_even, psi)
    psi.canonize_left()  # now L-canonical, if before AND real time evolution, omittable!
    psi.compress_svd_left(sv_min, bond_dim, test=False)
    return psi


def time_evolution_trotter_2order(hamiltonian_list, psi_initial, steps, eval_every, timestep, sv_min, bond_dim):

    unitary_list = unitary_from_bond(hamiltonian_list, timestep)
    unitary_list_half = unitary_from_bond(hamiltonian_list, 0.5 * timestep)
    mpo_odd = mpo_trotter_from_unitaries(psi_initial.length(),unitary_list, 1)
    mpo_even = mpo_trotter_from_unitaries(psi_initial.length(),unitary_list, 0)
    mpo_odd_half = mpo_trotter_from_unitaries(psi_initial.length(),unitary_list, 1)
    timecount = 0
    psi = psi_initial
    while timecount < steps:
        psi = mpo_times_mps (mpo_odd_half, psi)
        psi.canonize_left()  # now L-canonical, if before AND real time evolution, omittable!
        psi.compress_svd_left(sv_min, bond_dim, test=False)
        for i in range(eval_every-1):
            psi = trotter_1_step(psi, mpo_odd, mpo_even, sv_min,bonddim)
        psi = mpo_times_mps (mpo_even, psi)
        psi = mpo_times_mps (mpo_odd_half, psi)
        psi.canonize_left()  # now L-canonical, if before AND real time evolution, omittable!
        psi.compress_svd_left(sv_min, bond_dim, test=False)
        timecount += eval_every
        print "Trotter step, time: ", timecount, timecount * timestep
        energy = 0.0
        for i in range(1, psi.length):
            bondenergy = expectationvalue(psi, szsz(i, i+1)) + expectationvalue(psi, spsm(i, i+1))
            energy += bondenergy
            print round(real(bondenergy),5),
        print " "
        print "energy: ",energy
        print "current maximal bond dimension: ",max(psi.row_dims)
    return psi


def hamiltonian_to_state(lc, local_state, rc, local_operator):
    
    m_new = MPSTensor(local_operator.site_bra_dim(), size(lc, 0), size(rc, 0))
    for s in range(local_operator.site_ket_dim()):
        for bl in range(local_operator.row_dim()):
            for br in range(local_operator.col_dim()):
                calc_done = False
                for t in range(local_operator.site_bra_dim()):
                    if local_operator.tensor[t, s, bl, br] != 0.0:
                        if calc_done == False:
                            if br != 0:
                                mr = tensordot(local_state.tensor[s], rc.conj().transpose(1, 0, 2)[br], (2, 2))
                            else:
                                mr = local_state.tensor[s]
                            if bl != local_operator.row_dim()-1:
                                lmr = tensordot(lc.transpose(1, 0, 2)[bl], mr, (2, 0))
                            else:
                                lmr = mr
                            calc_done = True
                        m_new.tensor[t] += local_operator.tensor[t, s, bl, br] * lmr
    return m_new
                        
                            

def extremize(lc, local_state, rc, local_operator):
    
    alpha = []
    beta = []
    lanczos_list = []
    norm_factor = scalar_norm(local_state)
    multiply_by_scalar(local_state, 1.0 / norm_factor)
    lanczos_list.append(local_state)  # Krylov vector v_0 - MPSTensor
    lanczos_converged = False
    eval_min_old = 1.0e10
    while lanczos_converged == False:
        lanczos_list.append(hamiltonian_to_state(lc, lanczos_list[-1], rc, local_operator))
        alpha.append(scalar_overlap(lanczos_list[-2], lanczos_list[-1]))
        # test for convergence
        if alpha != []:     # at least 1x1 matrix
            tridiag = diag(alpha) + diag(beta, -1) + diag(beta, +1)
            evals, evecs = eigh(tridiag)
            improve = abs (eval_min_old - min(evals))
            eval_min_old = min(evals)
            if improve <= 1.0e-14:
                lanczos_converged == True            
        lanczos_list[-1] -= multiply_by_scalar(lanczos_list[-2], alpha[-1])
        if beta != []:
            lanczos_list[-1]-= multiply_by_scalar(lanczos_list[-3], beta[-1])
        beta.append(scalar_norm(lanczos_list[-1]))
        multiply_by_scalar(lanczos_list[-1], 1.0 / beta[-1])
        
    result = MPSTensor(local_state.site_dim(), local_state.row_dim(), local_state.col_dim())
    pos = argmin(evals)
    for i in range(size(evals, 0)):
        result += multiply_by_scalar(lanczos_list[i], evecs[i, pos]) 
        
    return local_state, eval_min_old


def ground_state(psi_guess, mpo_ham):


    if not (psi_guess.isrightcanonical(test) or psi_guess.ismixedcanonical(1,test)):
        psi_guess.canonize_right()
            
    list_lc = []
    list_rc = []
        
    rc = array((1, 1, 1),complex)
    rc[0, 0, 0] = 1.0
    lc = array((1, 1, 1),complex)
    lc[0, 0, 0] = 1.0
    list_rc.append(rc)
    list_lc.append(lc)
        
    convergence = False
        
# create list of right-contractions        

    for i in range(psi_guess.length(), 1, -1):
        rc = overlap_mpo_left_step(psi_guess.m[i], psi_guess.m[i], mpo_ham, rc)
        list_rc.append(rc)
        
# now dynamical sweeping forth and back

    while convergence == False:
        
# from left to right ...

        for i in range(1, self.length()):
            rc = list_rc.pop()
            psi_guess.m[i].tensor, eigenvalue = extremize(list_lc[-1], psi_guess.m[i].tensor, rc.conj(), mpo_ham.w[i].tensor).copy()
            psi_guess.m[i].normalization = 'U'
            psi_guess.canonize_left_step(i, False)
            list_lc.append(overlap_mpo_left_step(psi_guess.m[i], psi_guess.m[i], mpo_ham.w[i].tensor, list_lc[-1]))
                            
# from right to left ...

        for i in range(self.length(), 1, -1):
            lc = list_lc.pop()
            psi_guess.m[i].tensor, eigenvalue = extremize(lc, psi_guess.m[i].tensor, list_rc[-1], mpo_ham.w[i].tensor).copy()
            psi_guess.m[i].normalization = 'U'
            psi_guess.canonize_right_step(i, False)
            list_lc.append(overlap_mpo_right_step(psi_guess.m[i], psi_guess.m[i], mpo_ham.w[i].tensor, list_rc[-1]))

    return eigenstate, eigenvalue