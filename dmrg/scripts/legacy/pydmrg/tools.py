#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import numpy as np
import pyalps

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


def mergeXY(sets, foreach=[]):
    foreach_sets = {}
    for iset in pyalps.flatten(sets):
        fe_par_set = tuple((iset.props[m] for m in foreach))
        
        if fe_par_set in foreach_sets:
            foreach_sets[fe_par_set].append(iset)
        else:
            foreach_sets[fe_par_set] = [iset]
    
    for k,v in foreach_sets.items():
        common_props = pyalps.dict_intersect([q.props for q in v])
        res = pyalps.DataSet()
        res.props = common_props
        for im in range(0,len(foreach)):
            m = foreach[im]
            res.props[m] = k[im]
        for data in v:
            if len(res.x) > 0 and len(res.y) > 0:
                res.x = np.concatenate((res.x, data.x ))
                res.y = np.concatenate((res.y, data.y))
            else:
                res.x = data.x
                res.y = data.y
        
        order = np.argsort(res.x, kind = 'mergesort')
        res.x = res.x[order]
        res.y = res.y[order]
        res.props['label'] = ''
        for im in range(0,len(foreach)):
            res.props['label'] += '%s = %s ' % (foreach[im], k[im])
        
        foreach_sets[k] = res
    return foreach_sets.values()
