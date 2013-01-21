#!/bin/bash

DMRGAPP=$1
MEASAPP=$2

for tt in `ls test_*py`; do
  echo "Running $tt... "
  ./$tt $DMRGAPP $MEASAPP
done

