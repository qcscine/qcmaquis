
# -*- coding: utf-8 -*-
# ****************************************************************************
# 
# MAQUIS DMRG Project 
# 
# Copyright (C) 2013 by Sebastian Keller
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************


# common utility functions

import pyalps
import pyalps.hdf5 as h5
import string, re
import os.path
from copy import deepcopy

def ReadCpuTime(fname, ALPS=False):
    f = open(fname, 'r')
    content = f.read()
    val = content[content.rfind('Task'):content.rfind('seconds')].translate(None, string.letters+' ')
    return float(val)

def readParameters(h5ar, fname=None, proppath='/parameters'):
    if fname != None:
        d = {'filename' : fname}
    else:
        d = {}
    LOP=h5ar.list_children(proppath)
    for m in LOP:
            try:
                d[m] = h5ar.xml(proppath+'/'+m)
                try:
                    d[m] = float(d[m])
                except:
                    d[m] = map(float,d[m])
            except ValueError:
                pass
    return d 

def loadProperties(fname):
    ar = h5.archive(fname)
    return dict(ar['/parameters'])

#partly adapted
def ReadDMRGSweeps(ar, measurements, props, path='/spectrum/iteration'):
    ret = []
    print sorted(ar.list_children(path), key=lambda x: int(x.translate(None, string.letters)) )
    for sweep in sorted(ar.list_children(path), key=lambda x: int(x.translate(None, string.letters)) ):
        sweep_nr = int(sweep)
        ret_sweep = []
        p_sweep = path+'/'+sweep
        props_sweep = deepcopy(props)
        props_sweep['sweep'] = sweep_nr
        for meas in measurements:
            if meas in ar.list_children(p_sweep+'/results'):
                try:
                    d = pyalps.DataSet()
                    d.props = deepcopy(props_sweep)
                    d.props['observable'] = meas
                    d.props['hdf5_path'] = p_sweep+'/results/'+meas
                    d.y = ar[d.props['hdf5_path']+'/mean/value']
                    if 'labels' in ar.list_children(d.props['hdf5_path']):
                        d.x = pyalps.load.parse_labels(ar.read(d.props['hdf5_path']+'/labels'))
                    else:
                        d.x = range(len(d.y))
                    ret_sweep.append(d)
                except Exception, e:
                    print "Could not load " + meas
                    print e
                    pass
        ret.append(ret_sweep)
    return ret

#partly adapted
def LoadDMRGSweeps(files, measurements):
    ret = []
    for fname in files:
        ar = h5.archive(fname)
        if 'spectrum' in ar.list_children('/'):
            props = readParameters(ar, fname = fname)
            ret.append( ReadDMRGSweeps(ar, measurements, props) )
    return ret
