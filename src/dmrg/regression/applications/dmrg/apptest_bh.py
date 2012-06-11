#!/usr/bin/env python

from pyalps import apptest
# Explicitly specify "compMethod=..." and "outputs=..." if needed
apptest.runTest( './bh.testin.xml', outputs='auto', compMethod='auto' )
