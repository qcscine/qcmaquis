# dummy executable

from numpy import *
from mps import *
from mpo import *
from mpsinitialize import *

a=MPSTensor(3,4,5)
print MPSTensor.site_dim(a)
print a.site_dim()
print MPSTensor.isnormalized(a)
print a.isnormalized()

a=MPSTensor(3,4,10,init_random=True)

print a.tensor
print a.tensor[2,0,0]
print a.isleftnormalized(test=True)
print a.isrightnormalized(test=True)
print a.isnormalized(test=True)

r=a.normalize_right('qr', truncation=0.0)

print a.isleftnormalized(test=True)
print a.isrightnormalized(test=True)
print a.isnormalized(test=True)

site_dims=[1]+ 20*[2]
row_dims, col_dims=mps_matrix_sizes(site_dims, 40)
psi=mps_from_dimensions(site_dims, row_dims, col_dims, init_random=True)
phi=mps_from_dimensions(site_dims, row_dims, col_dims, init_random=True)
print psi.overlap(phi)
print phi.overlap(psi)
print psi.overlap(phi,direction='R')
print phi.overlap(psi,direction='R')

print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
psi.canonize(11)
print "after:"
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
print psi.ismixedcanonical(10,test=True)
print psi.ismixedcanonical(11,test=True)
print psi.ismixedcanonical(12,test=True)
psi.canonize_right()
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
psi.normalize(7)
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
print psi.isnormalized(7,test=True)
print psi.isnormalized(8,test=True)

print psi.site_dims()
print psi.row_dims()

psi=mps_aklt(30,0,0)
print psi.length()
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
psi.canonize(11)
print "after:"
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
print psi.ismixedcanonical(10,test=True)
print psi.ismixedcanonical(11,test=True)
print psi.ismixedcanonical(12,test=True)
psi.canonize_left()
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
psi.normalize(7)
print psi.isleftcanonical(test=True)
print psi.isrightcanonical(test=True)
print psi.isnormalized(7,test=True)
print psi.isnormalized(8,test=True)

print psi.site_dims()
print psi.row_dims()

psi=mps_ghz(30)
psi=mps_neel(30,0)
psi=mps_domain(30,3,0)
psi=mps_singlet_ladder(16)



