#!/usr/bin/env python2
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

import sys
import pyalps.hdf5 as h5

def loadProperties(fname):
    ar = h5.archive(fname)
    return ar['/parameters']


if __name__ == "__main__":
    f = sys.argv[1]
    props = loadProperties(f)
    for k,v in props.iteritems():
        if isinstance(v, str):
            print k, " = ", "\"" + v + "\""
        else:
            print k, " = ", v
