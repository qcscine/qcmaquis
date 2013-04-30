#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import pyalps
from pyalps.load import *
from pyalps.ngs import archive

from copy import deepcopy


def ReadDMRGSweeps(ar, measurements, props, path='/simulation', selector=None):
    ret = []
    for sweep_name in ar.list_children(path):
        if 'sweep' in sweep_name:
            sweep = eval(sweep_name.strip('sweep'))
            ret_sweep = []
            p_sweep = path+'/'+sweep_name
            props_sweep = deepcopy(props)
            props_sweep['sweep'] = sweep
            if 'parameters' in ar.list_children(p_sweep):
                props_sweep.update( dict(ar[p_sweep+'/parameters']) )
            if selector==None or selector(props_sweep):
                for meas in measurements:
                    if meas in ar.list_children(p_sweep+'/results'):
                        try:
                            d = pyalps.DataSet()
                            d.props = deepcopy(props_sweep)
                            d.props['observable'] = meas
                            d.props['hdf5_path'] = p_sweep+'/results/'+meas
                            d.y = ar[ d.props['hdf5_path']+'/mean/value' ]
                            if 'labels' in ar.list_children(d.props['hdf5_path']):
                                d.x = parse_labels(ar[ d.props['hdf5_path']+'/labels' ])
                            else:
                                d.x = range(len(d.y))
                            ret_sweep.append(d)
                        except Exception, e:
                            print "Could not load " + meas
                            print e
                            pass
            if len(ret_sweep) > 0:
                ret.append(ret_sweep)
    return ret

def LoadDMRGSweeps(files, what, selector=None):
    """ loads measurements in each DMRG iteration from HDF5 result files
        
        The parameters are:
        
        files:    a list of DMRG result files in HDF5 format.
        what:     an optional argument that is either a string or list of
                  strings, specifying the names of the observables which
                  should be loaded
        selector: (if specified) function taking a dict of properties and
                  returning a boolean of the dataset has to be loaded
        
        The function returns a list of list of lists of DataSet objects. 
        The elements of the outer list each correspond to the file names
        specified as input.
        The elements of the next level are different sweeps
        The elements of the inner-most list are each for a different observable.
        The y-values of the DataSet objects is an array of the measurements
        in all eigenstates calculated in this sector, and the x-values
        optionally the labels (indices) of array-valued measurements
    """
    
    if isinstance(what,str):
      what = [what]
    
    ret = []
    for fname in files:
        ar = archive(fname)
        if 'simulation' in ar.list_children('/'):
            props = dict(ar['/parameters'])
            props['filename'] = fname
            tmp = ReadDMRGSweeps(ar, what, props, path='/simulation', selector=selector)
            if len(tmp) > 0:
                ret.append(tmp)
    return ret

def LoadLastTruncation(sets, pname='TruncatedWeight'):
    last_sweep_selector = lambda props: props['sweep'] == props['nsweeps']-1
    for d in pyalps.flatten(sets):
        truncs = LoadDMRGSweeps([d.props['filename']], what=pname, selector=last_sweep_selector)
        if truncs:
            d.props[pname] = max(pyalps.flatten(truncs)[0].y)

