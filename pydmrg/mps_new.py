import math
import cmath
from numpy import *
from numpy.dual import svd, eigh

from numpy.linalg import LinAlgError

# global constants 

# maximal number of bosons per site, local dimension is this plus 1
max_bosons_per_site = 4
# length of spin: 0.5, 1, 1.5, 2, ...
spin_length=0.5

# deal with the SVD issue

def SchmidtGram(mat,startind):

    aux=mat.T.copy()
    for i in range(startind,aux.shape[0]):
        aux[i]=random.rand(aux.shape[1])+1.0j*\
               random.rand(aux.shape[1])
    for i in range(startind,aux.shape[0]):
        for k in range(i):
            aux[i]-=vdot(aux[k],aux[i])*aux[k]
        factor=sqrt(vdot(aux[i],aux[i]))
        aux[i]/=factor
    mat=aux.T.copy()
    del aux
    return
    
def sort_eigenpair(evals,evecs):
    evecs_hold=evecs.copy()
    evals_hold=evals.copy()
    sortlist=[]
    for i in range(evals.shape[0]):
        sortlist.append((evals_hold[i],i))
    sortlist.sort()
    sortlist.reverse()
    for i in range(evals.shape[0]):
        evals[i]=evals_hold[sortlist[i][1]]
        evecs[:,i]=evecs_hold[:,sortlist[i][1]]
    del evecs_hold
    return 

def dumbsvd(a):
    tol=10e-8
    noOK=0

    # A(m x n): m>=n or not?
    if a.shape[0]>=a.shape[1]: awork=a.copy()
    else: awork=a.T.conj().copy()

    a2=dot(awork.T.conj(),awork)
    evals,evecs=eigh(a2)
    sort_eigenpair(evals,evecs) # eigenvalues come unordered
    for i in range(evals.shape[0]):
        if evals[i]<0.0: evals[i]=0.0 
        if evals[i]>=tol*tol: noOK=noOK+1
    s=sqrt(evals)
    vh=evecs.T.conj()
    umod=dot(awork,evecs)
    u=umod.copy()
    for k in range(noOK): u[:,k]=umod[:,k]/s[k]
    if noOK !=evals.shape[0]: SchmidtGram(u,noOK+1)
    del evals,evecs,awork,umod
    
    if a.shape[0]>=a.shape[1]: return u,s,vh
    else: return vh.T.conj(),s,u.T.conj()
    
def mysvd(a):
    try:
        u,s,vh=svd(a,full_matrices=0,compute_uv=1)
    except LinAlgError:
        print "info from mysvd: numpy svd failed, dumb svd called"
        u,s,vh=dumbsvd(a)
    return u,s,vh



class PlainQuantumState:
    """General quantum states with exponentially many coefficients"""
    
# initialization

    def __init__(self,length,localdimensions):

# error catches:
# length must be an integer >=2
# localdimensions must be an array of integers of length length+1 
# (because localdimensions[0] is not used, but forced to be 1 for consistency)
# on sites, localdimensions[i] must be at least 1 
    
        if (not isinstance(length,int)) or length<2:
            print "error in PlainQuantumState.__init__: \
            length must be integer and at least 2"
            return
        
        if (not isinstance(localdimensions,list)):
            print "error in PlainQuantumState.__init__: \
            localdimensions must be list"
            return

        if (not len(localdimensions) == length+1):
            print "error in PlainQuantumState.__init__: \
            localdimensions must be list of length length+1"
            return

        for i in localdimensions:
            if (not isinstance(i,int)) or i<1:
                print "error in PlainQuantumState.__init__: \
                local dimension must be integer and at least 1"
                return
        
        self.length=length
        self.localdimensions=localdimensions.copy()
        self.localdimensions[0]=1
        self.size=1
        for i in localdimensions:
            self.size*=i
        self.c=zeros(self.size,complex)
        
    def copy(self):
        
        _psicopy_=PlainQuantumState(self.length,self.localdimensions)
        _psicopy_.c=self.c.copy()
        return _psicopy_
                 
def randomPlainQuantumState(length,localdimensions):
    
    rPQS=PlainQuantumState(length,localdimensions)
    rPQS.c=random.rand(rPQS.size)+1.0j*random.rand(rPQS.size)
### ??? normalize, is there more compact notation?
    norm_square=0
    for i in range(rPQS.size):
        norm_square+=rPQS.c[i]*rPQS.c[i].conj
    rPQS.c/=sqrt(norm_square)    
    del norm_square
    return rPQS

# MatrixProductState is the class for general matrix product states
# Python arrays start at index 0, sites are labeled from 1:
# convention to set A^0 to be a scalar (1x1 matrix)
# we allow a self.prefactor, which can be used to give canonical form to non-normalized MPS
# row and column dimensions are kept at top level and at the operator level for convenience
# this is informational overkill, but no strain on memory and time


# matrix product states in the Gamma/Lambda notation introduced by Vidal (2003)
# for inheritance, I call "gamma" "tensor"; "lambda" becomes "lbda" 
# I take it as a matrix, although often only the diagonal will be filled in view
# of McCulloch's iDMRG

class MatrixProductStateVidal(MatrixProductState):
    """Matrix product state in the Vidal notation"""

# initialization

    def __init__(self,length,localdimensions,rowdimensions,\
                 boundarycondition='OBC',prefactor=1.0+0.0j):
        MatrixProductState.__init__(self,length,localdimensions,rowdimensions,\
                 boundarycondition='OBC',prefactor=1.0+0.0j)
        self.bonddimensions=self.columndimensions.copy()
        self.lbda=[]
        for i in range(length+1):
            self.lbda.append(zeros((self.bonddimensions[i],self.bonddimensions[i]),complex))
        self.lbda[0][0,0]=1.0

# copy method

    def copy(self):
        _psicopy_=MatrixProductStateVidal(self.length,\
                                     self.localdimension,self.rowdimensions,\
                                     self.boundarycondition,self.prefactor)
        _psicopy_.rowdimensions=self.rowdimensions.copy()
        for i in range(self.length+1):
            _psicopy_.ldba[i]=self.lbda[i].copy()
            _psicopy_.gamma[i].tensor=self.a[i].tensor.copy()
        return _psicopy_

# conversion of plain quantum state into left normalized MPS

def Plain2MPS_L(plainpsi):

    MPSpsi=Plain2MPS(plainpsi,plainpsi.length+1)
    return MPSpsi

# conversion of plain quantum state into right normalized MPS

def Plain2MPS_R(plainpsi):

    MPSpsi=Plain2MPS(plainpsi,0)
    return MPSpsi

# conversion of plain quantum state into mixed normalized MPS
# other conversions default to this one

def Plain2MPS(plainpsi,mixedsite):

    tol=10e-8

# create dummy MPS with scalar matrices (to be changed dynamically)

    rowdims=ones(plainpsi.length+1,int)
    coldims=ones(plainpsi.length+1,int)

    MPSpsi=MatrixProductState(plainpsi.length,
                              plainpsi.localdimension,
                              rowdims,
                              coldims,
                              'OBC')
    coefficients=plainpsi.c.copy()
    local_rowdim=1
    svd_rowdim=plainpsi.localdimension
    svd_coldim=plainpsi.localdimension**(plainpsi.length-1)
    for i in range(1:mixedsite):
        u,s,vh=svd(coefficients.reshape(svd_rowdim,svd_coldim))
        svd_rowdim=0
        for j in range(min(svd_rowdim,svd_coldim)):
            if s[j,j]<tol:
                break
            svd_rowdim+=1
        
# generation routines for states 


def groundstate_1site(psi,hamiltonian,precision=1.e-14):

    # check whether psiguess is right-canonical or left-canonical to start with

    right_canonized_flag=True
    left_canonized_flag=True
    for i in range(1,psi.length+1):
        if psi.a[i].canonized !='R':
            right_canonized_flag=False
        if psi.a[i].canonized !='L':
            left_canonized_flag=False
    if right_canonized_flag==False and left_canonized_flag==False:
        print "error in groundstate_1site: only a fully canonized guess state is allowed"
        return
    
    # this being OK we initialise
    
    converged_flag=0

    left_contraction=[]
    right_contraction=[]
    left_contraction.append(LocalMatrixSet(1,1,1,'U'))
    right_contraction.append(LocalMatrixSet(1,1,1,'U'))
    left_stackpos=0
    right_stackpos=0
    
    if right_canonized_flag==True:
        for i in range(psi.length,1,-1):
            right_contraction.append(LocalMatrixSet(psi.a[i].localdimension,
                                                      psi.a[i].rowdimension,psi.a[i].rowdimension,'U'))
            right_stackpos+=1
            right_contraction[right_stackpos].tensor[0]=
            identity(right_contraction[right_stackpos].rowdimension,complex) 
            for k1 in range(1,right_contraction[right_stackpos].localdimension):
                for k2 in range(0,right_contraction[right_stackpos-1].localdimension):
                    for s1 in range(psi.a[i].localdimension):
                        for s2 in range(psi.a[i].localdimension):
                            if abs(hamiltonian.a[i].tensor[s1,s2,k1,k2]) !=0.0:
                                hamfactor=hamiltonian.a[i].tensor[s1,s2,k1,k2]
                                right_contraction[right_stackpos].tensor[k1]+=
                                hamfactor*dot(psi.a[i].tensor[s1].conj(),
                    dot(right_contraction[right_stackpos-1].tensor[k2],psi.a[i].tensor[s2].T))
    else:
        for i in range(1,psi.length):
            left_contraction.append(LocalMatrixSet(psi.a[i].localdimension,
                                                      psi.a[i].columndimension,psi.a[i].columndimension,'U'))
            left_stackpos+=1
            left_contraction[left_stackpos].tensor[left_contraction[left_stackpos].localdimension]=diag(1.0) # look up proper syntax, I want diagonal
            # notation to be amended
        
            for k1 in range(0,left_contraction.a[left_stackpos-1].localdimension):
                for k2 in range(0,left_contraction.a[left_stackpos].localdimension-1):
                    for s1 in range(psi.a[i].localdimension):
                        for s2 in range(psi.a[i].localdimension):
                            if abs(hamiltonian.a[i].tensor[s1,s2,k1,k2]) !=0.0:
                                hamfactor=hamiltonian.a[i].tensor[s1,s2,k1,k2]
                                left_contraction[left_stackpos].tensor[k2]+=
                                hamfactor*dot(psi.a[i].tensor[s1].T.conj(),
                    dot(left_contraction[left_stackpos-1].tensor[k1],psi.a[i].tensor[s2]))
 
#####################                                
                                
    while converged_flag==0:

        # this order only good for right-normalized; simple copy? or two routines? after debugging ...
        
        # sweep from left to right ...
        
        for i in range(1,psiguess.length+1):
            current_right_contraction=right_contraction.pop()
            right_stackpos-=1
            for sigma in range(psiguess.a[i].localdimension):
                psiguess.a[i].tensor[sigma]=dot(left_contraction[left_stackpos], 
                                            dot(psi.a[i].tensor[sigma],current_right_contraction))
            MatrixProductState.canonize_left_step(psiguess,i)
            if i != psiguess.length:
                left_contraction.append(zeros(psiguess.a[i].columndimension,psi.a[i].columndimension),complex)
                left_stackpos+=1
                for sigma in range(psiguess.a[i].localdimension):
                    left_contraction[left_stackpos]+=dot(psiguess.a[i].tensor[sigma].conj(),
                                dot(left_contraction[left_stackpos-1],psi.a[i].tensor[sigma].T))
            
        right_contraction.append(ones(1,1),complex)
        right_stackpos+=1
        
        # ... and from right to left ...
        
        for i in range(psiguess.length,0,-1):
            current_left_contraction=left_contraction.pop()
            left_stackpos-=1
            for sigma in range(psiguess.a[i].localdimension):
                psiguess.a[i].tensor[sigma]=dot(current_left_contraction,
                                                dot(psi.a[i].tensor[sigma],right_contraction[right_stackpos]))
            MatrixProductState.canonize_right_step(psiguess,i)
            if i != psiguess.length:
                right_contraction.append(zeros(psiguess.a[i].rowdimension,psi.a[i].rowdimension),complex)
                right_stackpos+=1
                for sigma in range(psiguess.a[i].localdimension):
                    right_contraction[right_stackpos]+=dot(psiguess.a[i].tensor[sigma].conj(),
                                dot(right_contraction[right_stackpos-1],psi.a[i].tensor[sigma].T))
                    
        left_contraction.append(ones(1,1),complex)
        left_stackpos+=1
                
        if abs(overlap(psiold,psiguess))>1.0-conv:
            converged_flag=1
        psiold=psiguess
##########
    
    return energy,psi


# ##################### # ################### # ######################


#################################

def NNunitary(ham,timestep):
    size2=ham.shape[0]
    size=int(sqrt(size2))
    eigvals,eigvecs=eigh(ham)
    expmat=diag(exp(-1.0j*timestep*eigvals))
    evolutionoperator=dot(dot(eigvecs,expmat),eigvecs.T.conj())
    NNevolop=zeros((size,size2,size,size2),complex)
    for i in range(size):
        NNevolop[i,:,i,:]=evolutionoperator[:,:]
    
    u2,s2,vh2=mysvd(NNevolop.reshape(size,size,size,size,size,size).\
                 transpose(1,0,4,3,2,5).reshape(size2*size2,size2))
    singvals2=0
    tol=10e-8
    for i in range(size2):
        if real(s2[i])>tol:
            singvals2+=1
    ops3_set=LocalMatrixOperatorSet(size,size,singvals2,1)
    for k in range(singvals2):
        ops3_set.tensor[:,:,k,0]=sqrt(real(s2[k]))*vh2.reshape(size2,size,size)[k,:,:]
    u1,s1,vh1=mysvd(u2.reshape(size,size,size,size,size2).\
                    transpose(0,2,1,3,4).\
                    reshape(size2,size2*size2))
    singvals1=0
    tol=10e-8
    for i in range(size2):
        if real(s1[i])>tol:
            singvals1+=1
    ops1_set=LocalMatrixOperatorSet(size,size,1,singvals1)
    ops2_set=LocalMatrixOperatorSet(size,size,singvals1,singvals2)
    for k in range(singvals1):
        ops1_set.tensor[:,:,0,k]=sqrt(real(s1[k]))*u1.reshape(size,size,size2)[:,:,k]
    for k1 in range(singvals1):
        for k2 in range(singvals2):
            ops2_set.tensor[:,:,k1,k2]=\
            sqrt(real(s1[k1]))*sqrt(real(s1[k2]))*\
            vh1.reshape(size2,size,size,size2).transpose(1,2,0,3)[:,:,k1,k2]
    del u1,u2,s1,s2,vh1,vh2,eigvals,eigvecs,expmat
    return singvals1,singvals2,ops1_set,ops2_set,ops3_set



def NNmakeTrotterMPO(length,ops1,ops2,ops3,oddeven):
    opslist=[]
    rowdimensions=zeros((length+1),int)
    opsunit=LocalMatrixOperatorSet(ops1.localdimension1,ops1.localdimension2,1,1)
    for k in range(opsunit.localdimension1):
        opsunit.tensor[k,k,0,0]=1.0
    opslist.append(opsunit) #at position 0
    rowdimensions[0]=1
    if oddeven==1:
        for i in range(1,length+1,4):
            if i+3<=length:
                opslist.append(ops1)
                opslist.append(ops2)
                opslist.append(ops3)
                opslist.append(opsunit)
                rowdimensions[i]=1
                rowdimensions[i+1]=ops2.rowdimension
                rowdimensions[i+2]=ops3.rowdimension
                rowdimensions[i+3]=1
            else:
                opslist.append(opsunit) # at position L if unpaired
                opslist.append(opsunit) # at position L+1 if unpaired
                rowdimensions[i]=1
                rowdimensions[i+1]=1
    else:
        opslist.append(opsunit) # at position 1
        opslist.append(opsunit) # at position 2
        rowdimensions[1]=1
        rowdimensions[2]=1
        for i in range(3,length+1,4):
            if i+3<=length:
                opslist.append(ops1)
                opslist.append(ops2)
                opslist.append(ops3)
                opslist.append(opsunit)
                rowdimensions[i]=1
                rowdimensions[i+1]=ops2.rowdimension
                rowdimensions[i+2]=ops3.rowdimension
                rowdimensions[i+3]=1
            else:
                opslist.append(opsunit)
                opslist.append(opsunit)
                rowdimensions[i]=1     
                rowdimensions[i+1]=1     
    trottermpo=makeMPOfromUnitaries(length,ops1.localdimension1,ops1.localdimension2,rowdimensions,opslist)
    del opslist,rowdimensions,opsunit
    return trottermpo



def finiteT_timeevolution(ham,psi,steps,timestep,bonddim):
    singvals1,singvals2,ops1,ops2,ops3=NNunitary(ham,timestep)
    singvals1,singvals2,ops1_half,ops2_half,ops3_half=NNunitary(ham,timestep*0.5)
    mpo_odd=NNmakeTrotterMPO(psi.length,ops1,ops2,ops3,1)
    mpo_even=NNmakeTrotterMPO(psi.length,ops1,ops2,ops3,0)
    mpo_odd_half=NNmakeTrotterMPO(psi.length,ops1_half,ops2_half,ops3_half,1)
    mpo_even_half=NNmakeTrotterMPO(psi.length,ops1_half,ops2_half,ops3_half,0)
    timecount=0
    while timecount<steps:
        print "timestep: ",timecount
        psi=Trotter2step(psi,mpo_odd,mpo_even,mpo_odd_half,mpo_even_half,bonddim)
        energy=0.0
        for i in range(1,psi.length):
            bondenergy=expectationvalue(psi,szsz(i,i+1))+expectationvalue(psi,spsm(i,i+1))
            energy+=bondenergy
            print round(real(bondenergy),5),
        print " "
        print "energy: ",energy
        print "current maximal bond dimension: ",max(psi.rowdimensions)
        timecount+=1    
    return psi


# attempt a timeevolution:

psi=neel(32,30,1)
psi=timeevolution(HAFM(1.0,1.0),psi,80,-0.05j,30)
psi=timeevolution(HAFM(1.0,1.0),psi,80,-0.01j,40)
psi=timeevolution(HAFM(1.0,1.0),psi,80,-0.005j,50)
psi=timeevolution(HAFM(1.0,1.0),psi,80,-0.001j,60)
psi=timeevolution(HAFM(1.0,1.0),psi,80,-0.0005j,70)
psi=timeevolution(HAFM(1.0,1.0),psi,80,-0.0001j,70)

    
    
    
