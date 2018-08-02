#!/bin/python
import re
test =
"""ATOM     47  N   GLU A   9      25.070  11.305  -6.635  1.00 23.65           N  
ATOM     48  CA  GLU A   9      23.943  11.996  -7.247  1.00 27.47           C  
ATOM     49  C   GLU A   9      22.623  11.558  -6.624  1.00 21.91           C  
ATOM     50  O   GLU A   9      21.712  12.407  -6.493  1.00 21.71           O  
ATOM     51  CB  GLU A   9      23.927  11.759  -8.758  1.00 33.99           C  
ATOM     52  CG  GLU A   9      25.173  12.252  -9.476  1.00 42.66           C  
ATOM     53  CD  GLU A   9      25.128  11.990 -10.969  1.00 46.64           C  
ATOM     54  OE1 GLU A   9      24.118  11.430 -11.443  1.00 51.29           O  
ATOM     55  OE2 GLU A   9      26.103  12.345 -11.664  1.00 48.91           O  
ATOM     56  N   ASP A  10      22.384  10.275  -6.272  1.00 20.79           N  
ATOM     57  CA  ASP A  10      21.123   9.959  -5.590  1.00 19.91           C  
ATOM     58  C   ASP A  10      21.044  10.774  -4.296  1.00 21.04           C  
ATOM     59  O   ASP A  10      19.923  11.195  -3.966  1.00 20.30           O  
ATOM     60  CB  ASP A  10      20.942   8.479  -5.285  1.00 19.58           C  
ATOM     61  CG  ASP A  10      20.508   7.707  -6.545  1.00 21.82           C  
ATOM     62  OD1 ASP A  10      20.145   8.291  -7.582  1.00 22.06           O  
ATOM     63  OD2 ASP A  10      20.575   6.481  -6.454  1.00 24.91           O  
ATOM     64  N   ARG A  11      22.162  10.914  -3.578  1.00 18.00           N  
ATOM     65  CA  ARG A  11      22.125  11.786  -2.389  1.00 18.52           C  
ATOM     66  C   ARG A  11      21.798  13.230  -2.777  1.00 17.66           C  
ATOM     67  O   ARG A  11      21.051  13.864  -2.036  1.00 17.03           O  
ATOM     68  CB  ARG A  11      23.384  11.791  -1.549  1.00 18.18           C  
ATOM     69  CG  ARG A  11      23.442  12.729  -0.352  1.00 19.32           C  
ATOM     70  CD  ARG A  11      22.307  12.504   0.668  1.00 19.28           C  
ATOM     71  NE  ARG A  11      22.507  13.279   1.907  1.00 16.77           N  
ATOM     72  CZ  ARG A  11      21.574  13.432   2.827  1.00 19.25           C  
ATOM     73  NH1 ARG A  11      20.375  12.833   2.649  1.00 17.83           N  
ATOM     74  NH2 ARG A  11      21.808  14.108   3.963  1.00 18.74           N  
ATOM     75  N   LYS A  12      22.377  13.713  -3.845  1.00 18.11           N  
ATOM     76  CA  LYS A  12      22.083  15.058  -4.327  1.00 18.05           C  
ATOM     77  C   LYS A  12      20.567  15.239  -4.542  1.00 19.80           C  
ATOM     78  O   LYS A  12      19.969  16.155  -4.008  1.00 18.10           O  
ATOM     79  CB  LYS A  12      22.794  15.310  -5.643  1.00 20.78           C  
ATOM     80  CG  LYS A  12      22.477  16.705  -6.156  1.00 24.94           C  
ATOM     81  CD  LYS A  12      22.180  16.775  -7.621  1.00 29.14           C  
ATOM     82  CE  LYS A  12      23.295  17.360  -8.442  1.00 32.80           C"""

residue_list = []

#complete to extract the residue,chain, res_num and xyz cords
# e.g LYS 12 23.295 17.360 -8.442
# group1. = AAname  group2 =Chain group3 = residue.no group4 = x_cord
# group5 = y_cord group6 = z cord
p = re.compile(r'CA\s+(\w+?)\s(\w)\s+?(\d+?)\s+?(\S+?)\s+?(\S+?)\s+?(\S+?)\s')
for line in test:
    if p.search(test):
        
