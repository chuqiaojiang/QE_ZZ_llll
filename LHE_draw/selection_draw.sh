#!/bin/bash

mkdir plots
rm -r plots/*

paraname="NoAQGC"
rwstart=-100
rwstep=10
rwnstep=20

root -l -q -b 'selection_draw.C("'${paraname}'",'${rwstart}','${rwstep}','${rwnstep}')' > selection_draw.log
