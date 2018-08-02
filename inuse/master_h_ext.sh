#!/bin/sh
pdbsecstr 1ct5.pdb 1ct5.pdb.sec
helix_extractor.py 1ct5.pdb.sec > resfile
pdbgetresidues resfile 1ct5.pdb
