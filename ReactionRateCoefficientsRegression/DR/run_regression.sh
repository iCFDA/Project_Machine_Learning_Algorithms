#!/bin/bash
cd src
# vibrational levels of N2
for i in {0..46}; do \
  python3 regression_dis.py $i
  #python3 regression_rec.py $i
done
cd ..
