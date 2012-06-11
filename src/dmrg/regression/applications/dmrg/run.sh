#!/bin/bash

tests="bh ff_middle fermiladder"

for tt in $tests; do
  echo "Running $tt... "
  python apptest_$tt.py 
done


