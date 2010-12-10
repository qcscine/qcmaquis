from constants import *
from numpy import *

# definitions of physical operators 
# in explicit arrays of suitable dimensions

one2 = ones((2,2), complex)
one3 = ones((3,3), complex)
one4 = ones((4,4), complex)


# spin-1/2

# 0: Sz=+1/2, 1: Sz=-1/2
# sz=S^{z}, sp=S^{+}, sm=S^{-} for spin-1/2

sz=zeros((2,2),float)
sp=zeros((2,2),float)
sm=zeros((2,2),float)
sz[0,0]=+0.5
sz[1,1]=-0.5
sp[0,1]=+1.0
sm[1,0]=+1.0

# spin-1
# 0: Sz=+1, 1: Sz=0, 2: Sz=-1
# sz1=S^{z}, sp1=S^{+}, sm1=S^{-}, fl1=e^{iS^z\pi} for spin-1

sz1=zeros((3,3),float)
sp1=zeros((3,3),float)
sm1=zeros((3,3),float)
fl1=zeros((3,3),float)
sz1[0,0]=+1.0
sz1[1,1]=0.0
sz1[2,2]=-1.0
sp1[0,1]=sqrt(2.0)
sp1[1,2]=sqrt(2.0)
sm1[1,0]=sqrt(2.0)
sm1[2,1]=sqrt(2.0)
fl1[0,0]=-1.0
fl1[1,1]=+1.0
fl1[2,2]=-1.0

# for spin-S
# 0: Sz=S, 1: Sz=S-1, 2: Sz=S-2, and so on up to Sz=-S
# szs=S^{z}, sps=S^{+}, sms=S^{-} for spin-SPIN_LENGTH

szs=zeros(2*SPIN_LENGTH+1,2*SPIN_LENGTH+1),float)
sps=zeros(2*SPIN_LENGTH+1,2*SPIN_LENGTH+1),float)
sms=zeros(2*SPIN_LENGTH+1,2*SPIN_LENGTH+1),float)
for i in range(2*SPIN_LENGTH+1):
    szs[i,i]=SPIN_LENGTH-i
for i in range(2*SPIN_LENGTH):
    sps[i,i+1]=sms[i+1,i]=\
    sqrt(SPIN_LENGTH*(SPIN_LENGTH+1)-(SPIN_LENGTH-i)*(SPIN_LENGTH-i+1))

# spinless bosons
# 0: zero bosons, 1: 1 boson, ..., up to MAX_BOSONS_PER_SITE
# bcr=b^{\dagger}, ban=b, bn=n=b^{\dagger}b, bnn=n(n+1)

bcr=zeros((MAX_BOSONS_PER_SITE+1,MAX_BOSONS_PER_SITE+1),float)
ban=zeros((MAX_BOSONS_PER_SITE+1,MAX_BOSONS_PER_SITE+1),float)
bn =zeros((MAX_BOSONS_PER_SITE+1,MAX_BOSONS_PER_SITE+1),float)
bnn1=zeros((MAX_BOSONS_PER_SITE+1,MAX_BOSONS_PER_SITE+1),float)
for i in range(MAX_BOSONS_PER_SITE):
    bcr[i+1,i]=ban[i,i+1]=sqrt(i+1)
for i in range(MAX_BOSONS_PER_SITE+1):
    bn[i,i]=i
    bnn1[i,i]=i*(i-1)

# spinless fermions
# 0: empty, 1: occupied
# fcr=c^{\dagger}, fan=c, fn=n=c^{\dagger}c, fp=e^{i n \pi}

fcr=zeros((2,2),float)
fan=zeros((2,2),float)
fn=zeros((2,2),float)
fp=zeros((2,2),float)
fcr[1,0]=fan[0,1]=+1.0
fn[1,1]=+1.0
fp[0,0]=+1.0
fp[1,1]=-1.0

# spin-1/2 fermions (electrons)
# 0: empty, 1: spin-up, 2: spin-down, 3: doubly occupied
# fermionic sign: ordering up-1 down-1 up-2 down-2 ... 
# sign in down-i depending on up-i empty or occupied
# fermionic signs due to previous sites will be treated explicitly
# cucr=c_{\uparrow}^{\dagger}, cuan=c_{\uparrow}
# cdcr=c_{\downarrow}^{\dagger}, cdan=c_{\downarrow}
# cn=n_{\uparrow}+n_{\downarrow}
# csz=0.5*(n_{\uparrow}-n_{\downarrow})
# cun=n_{\uparrow}, cdn=n_{\downarrow}
# cnn=n_{\uparrow}n_{\downarrow}
# cp=e^{i\pi n}

cucr=zeros((4,4),float)
cuan=zeros((4,4),float)
cdcr=zeros((4,4),float)
cdan=zeros((4,4),float)
cn=zeros((4,4),float)
csz=zeros((4,4),float)
cun=zeros((4,4),float)
cdn=zeros((4,4),float)
cnn=zeros((4,4),float)
cp=zeros((4,4),float)

# creators & annihilators

cucr[1,0]=cuan[0,1]=1.0
cucr[3,2]=cuan[2,3]=1.0
cdcr[2,0]=cdan[0,2]=1.0
cdcr[3,1]=cdan[1,3]=-1.0

# densities and magnetizations

cn[1,1]=cn[2,2]=1.0
cn[3,3]=2.0
csz[1,1]=+0.5
csz[2,2]=-0.5
cun[1,1]=cun[3,3]=+1.0
cdn[2,2]=cdn[3,3]=+1.0
cnn[3,3]=+1.0

# fermionic sign

cp[0,0]=cp[3,3]=+1.0
cp[1,1]=cp[2,2]=-1.0



