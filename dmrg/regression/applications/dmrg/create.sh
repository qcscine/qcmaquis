#!/bin/bash

tests="bh ff_middle fermiladder"

for tt in $tests; do
  echo "Creating $tt... "
  python $tt.py --create
done


