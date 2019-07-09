#!/bin/bash

DMRGAPP=$1
MEASAPP=$2

for tt in `ls test_*py`; do
  echo "Creating $tt... "
  ./$tt --create $DMRGAPP $MEASAPP
done

