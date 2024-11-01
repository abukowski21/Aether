#!/bin/sh


rm -rf run.2d
 
cp -r ../../share/run ./run.2d
cd run.2d
cp ../aether.json.2d ./aether.json

./aether

cd UA/output

~/bin/postAether.py -rm
~/bin/aether_plot.py -var=Temperature_neutral -alt=2 3DALL_20110320_010000.nc
