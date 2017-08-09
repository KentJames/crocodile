#!/bin/env python3

import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

import argparse
import h5py
import numpy
import math
import scipy.optimize

from arl.test_support import import_visibility_from_oskar, export_visibility_to_fits
from crocodile.synthesis import w_kernel

# Parse arguments
parser = argparse.ArgumentParser(description='Generate w kernel file for gridding a dataset')
parser.add_argument('output', type=str, help='output file')
parser.add_argument('--theta', dest='theta', type=float, required=True,
                    help='Field of view size')
parser.add_argument('--ff', dest='ff', type=int, required=True,
                    help='Far field size')
parser.add_argument('--wstep', dest='wstep', type=float, required=True,
                    help='Step length for w coordinates')
parser.add_argument('--wcount', dest='wcount', type=int, required=True,
                    help='Number of w-planes to generate kernels for')
parser.add_argument('--error', dest='error', type=float, default=0.01,
                    help='Estimated Error of kernels to generate')
parser.add_argument('--oversample', dest='oversample', type=int, default=8,
                    help='Amount of oversampling to use for kernel (default 8)')
parser.add_argument('--overwrite', dest='overwrite', const=True, default=False, action='store_const',
                    help='Overwrite existing kernels in output file?')
args = parser.parse_args()

# Open output file
output = h5py.File(args.output, "a")

# Calculate recommended w-kernel size
wmax = args.wstep * args.wcount
def recommended_size(wval,eps=0.01):
    usupp = numpy.sqrt((wval*args.theta/2)**2 + (wval**1.5 * args.theta / 2 / numpy.pi / eps))
    return math.ceil(2 * args.theta * usupp) // 2 * 2 +1 #numpy.ceil doesn't work???
def recommended_diff(eps):
    return abs(recommended_size(eps) - args.size)

print("w range:        %.1f - %.1f lambda" % (-wmax, wmax))

# Generate kernels
for iw in range(-args.wcount, args.wcount+1):
    w = iw * args.wstep

    npixkern = recommended_size(w,eps=args.error)

    if npixkern < 3:
        npixkern = 3
    
    print("Npixkern: %f " %(npixkern),end="", flush=True)
    
    print("W-Plane: %f " % w, flush=True)

    # Check whether it already exists
    kern_name = 'wkern/%s/%s/kern' % (args.theta, w)
    if kern_name in output:
        if args.overwrite:
            del output[kern_name]
        else:
            continue

    # Make kernel
    kern = w_kernel(args.theta, w, NpixFF=args.ff, NpixKern=npixkern, Qpx=args.oversample)
    output[kern_name] = kern
    output.flush()

print("done")
output.close()
