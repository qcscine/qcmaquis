#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import pyalps
from pyalps.load import *
try:
    from pyalps.ngs import archive
except:
    from pyalps.hdf5 import archive

from copy import deepcopy

def ReadMeasurements(ar, measurements, path, props):
    ret = []
    if measurements is None:
        measurements = ar.list_children(path)

    for meas in measurements:
        if meas in ar.list_children(path):
            try:
                d = pyalps.DataSet()
                d.props = deepcopy(props)
                d.props['observable'] = meas
                d.props['hdf5_path'] = path+'/'+meas
                d.y = ar[ d.props['hdf5_path']+'/mean/value' ]
                if 'labels' in ar.list_children(d.props['hdf5_path']):
                    d.x = parse_labels(ar[ d.props['hdf5_path']+'/labels' ])
                else:
                    d.x = range(len(d.y))
                ret.append(d)
            except Exception, e:
                print "Could not load " + meas
                print e
                pass
    return ret

def ReadDMRGSweeps(ar, measurements, props, path='/simulation', selector=None):
    ret = []
    for sweep_name in ar.list_children(path):
        if sweep_name == 'iteration':  ## new iteration format
            # add sweeps in order
            sweeps = sorted(map(int, ar.list_children(path+'/iteration')))
            for sweep in sweeps:
                sweep_num = str(sweep)
                ret_sweep = []
                p_sweep = path+'/iteration/'+sweep_num
                props_sweep = deepcopy(props)
                props_sweep['sweep'] = sweep
                if 'parameters' in ar.list_children(p_sweep):
                    props_sweep.update( dict(ar[p_sweep+'/parameters']) )
                if selector==None or selector(props_sweep):
                    ret_sweep += ReadMeasurements(ar, measurements, p_sweep+'/results', props_sweep)
                if len(ret_sweep) > 0:
                    ret.append(ret_sweep)
        elif 'sweep' in sweep_name:  ## old iterarion format
            sweep = eval(sweep_name.strip('sweep'))
            ret_sweep = []
            p_sweep = path+'/'+sweep_name
            props_sweep = deepcopy(props)
            props_sweep['sweep'] = sweep
            if 'parameters' in ar.list_children(p_sweep):
                props_sweep.update( dict(ar[p_sweep+'/parameters']) )
            if selector==None or selector(props_sweep):
                ret_sweep += ReadMeasurements(ar, measurements, p_sweep+'/results', props_sweep)
            if len(ret_sweep) > 0:
                ret.append(ret_sweep)
        elif sweep_name == 'results' and 'simulation' in path:
            props_sweep = deepcopy(props)
            props_sweep['sweep'] = -1
            if 'parameters' in ar.list_children(path):
                props_sweep.update( dict(ar[path+'/parameters']) )
            if selector==None or selector(props_sweep):
                tmp = ReadMeasurements(ar, measurements, path+'/results', props_sweep)
                if len(tmp) > 0:
                    ret.append(tmp)
    return ret

def LoadDMRGSweeps(files, what=None, selector=None):
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
        if 'spectrum' in ar.list_children('/') and 'iteration' in  ar.list_children('/spectrum'):
            props = dict(ar['/parameters'])
            props['filename'] = fname
            tmp = ReadDMRGSweeps(ar, what, props, path='/spectrum', selector=selector)
            if len(tmp) > 0:
                ret.append(tmp)
    return ret

def LoadLastTruncation(sets, pname='TruncatedWeight'):
    last_sweep_selector = lambda props: props['sweep'] == props['nsweeps']-1
    for d in pyalps.flatten(sets):
        truncs = LoadDMRGSweeps([d.props['filename']], what=pname, selector=last_sweep_selector)
        if truncs:
            d.props[pname] = max(pyalps.flatten(truncs)[0].y)

