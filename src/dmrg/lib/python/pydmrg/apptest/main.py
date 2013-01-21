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
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-c', '--create', dest='action', action='store_const',
                      const='create', default='run', help='create test reference results')
    group.add_argument('-p', '--params', dest='action', action='store_const',
                       const='params', default='run', help='create parameter files')
    parser.add_argument('dmrg_app', help='path to `dmrg` application')
    parser.add_argument('meas_app', help='path to `dmrg_meas` application')
    
    args = parser.parse_args()
    
    for test in find_subclasses():
        t = test()
        if args.action == 'create':
            t.create(dmrg_app=args.dmrg_app, meas_app=args.meas_app)
        elif args.action == 'params':
            t.write_parameters()
        else:
            t.run(dmrg_app=args.dmrg_app, meas_app=args.meas_app)
            t.check()
