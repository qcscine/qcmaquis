#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import sys
import inspect
import argparse

from testbase import DMRGTestBase

def find_subclasses(module_name='__main__', clazz=DMRGTestBase):
    module = __import__(module_name)
    return [ cls for name, cls in inspect.getmembers(module) if inspect.isclass(cls)
                                                                and issubclass(cls, clazz) ]


def main():
    parser = argparse.ArgumentParser(description='Testing tool for MAQUIS DMRG')
    parser.add_argument('--create', '-c', dest='action', action='store_const',
                      const='create', default='run', help='create test reference results')
    parser.add_argument('dmrg_app', help='path to `dmrg` application')
    parser.add_argument('meas_app', help='path to `dmrg_meas` application')
    
    args = parser.parse_args()
    
    for test in find_subclasses():
        t = test()
        if args.action == 'create':
            t.create(dmrg_app=args.dmrg_app, meas_app=args.meas_app)
        else:
            t.run(dmrg_app=args.dmrg_app, meas_app=args.meas_app)
            t.check()
