#############################################################################
#
# MAQUIS DMRG Project
# Vistrails package
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

identifier = 'org.comp-phys.maquis'
version = '2.0.0'
name = 'MAQUIS - DMRG'


##############################################################################

def package_requirements():
    import core.requirements
    if not core.requirements.python_module_exists('pydmrg'):
        raise core.requirements.MissingRequirement('pydmrg')

def package_dependencies():
    return ['org.comp-phys.alps']
