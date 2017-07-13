#!/usr/bin/python
######################################################################
# Name:         interferometry.py
# Author:		A. Marocchino, F. Filippi
# Date:			2017-07-13
# Purpose:      synthetic inteferometry from binary ALaDyn files
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, sys
import numpy as np
import matplotlib
#matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes','Code_ALaDyn','tools-ALaDyn','pythons'))
from read_ALaDyn_bin import *
### --- ###


# - #
#bin_path = '/'
#n,x,y,z = read_ALaDyn_bin(bin_path,'rho_out3.bin','grid')
