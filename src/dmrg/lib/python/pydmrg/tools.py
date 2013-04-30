#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import numpy as np


def labels_argsort(x):
    # sorting observables
    x = np.array(x)
    if len(x.shape) > 1:
        x = x.reshape( x.shape[0], np.prod(x.shape[1:]) )
        # keys = transpose(x)[::-1]
        keys = np.rot90(x, 2)
        ind = np.lexsort(keys)
    else:
        ind = np.argsort(x)
    
    return ind

