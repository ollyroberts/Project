#!/bin/bash

i=0; 
for f in *.pdb; 
do 
    d=dir_$(printf %03d $((i/100+1))); 
    mkdir -p $d; 
    mv "$f" $d; 
    let i++; 
done
