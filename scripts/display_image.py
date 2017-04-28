import sys
import os
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

import argparse
import h5py
import numpy
import itertools
import pylru
import ctypes
from multiprocessing import Process, Array, Queue

import arl.test_support
from crocodile.synthesis import *
import util.visualize

parser = argparse.ArgumentParser(description='Display Image')
parser.add_argument('image', metavar='image', type=argparse.FileType('r'),
        help='input visibilities')
parser.add_argument('--lambda', dest='lam', type=float, required=True, default=300000, help='Grid size')
parser.add_argument('--theta', dest='theta', type=float, required=True, default=0.08, help='Field of view size')
parser.add_argument('--grid', dest='grid', type=bool, required=False, default=False, help='Show grid?')
args = parser.parse_args()

if args.grid is True:
    img = numpy.fromfile(args.image, dtype=complex)
    print(img.shape)
    img.reshape(24000,24000)
    print(img.shape)
    img = img[6000:18000,6000:18000]
    print(img.shape)
    util.visualize.show_image(img,"grid",args.theta)
else:

    img = numpy.fromfile(args.image, dtype=float)
    print(img.shape)
    img = img.reshape(24000,24000)
    img = img[6000:18000,6000:18000]
    print(img.shape)
    print(img[12000:12010,12000:12010])
    util.visualize.show_image(img,"result",args.theta)

