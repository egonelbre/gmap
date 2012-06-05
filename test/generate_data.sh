#!/bin/bash

SIZE=10000x1000

../bin/gmapgenerate -s $SIZE -g data/test.gt data/test
../bin/gmapconvert -i data/test data/test
../bin/gmaprandpheno -i data/test data/pheno.txt
rm data/test.gmap data/test.midx data/test.gidx data/test.pidx
