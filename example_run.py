# example script to run PSG to get telluric spectrum
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.time import Time
from scipy import interpolate

import glob,os,sys

plt.ion()

sys.path.append('psg/')
from run_psg import run_psg
from tools import *


if __name__=='__main__':
	# define everything
	output_path = './outputs/'
	obs_time = '2019/11/18 10:07:16.000'
	lon, lat, elev = 33.1504, 242.8173, 1.871 # palomar 
	site = [lon,lat,elev] # pull from header info?
	l0, l1, res    = 700, 720, 100000
	line_list = 'HIT'

	# run psg, save telluric spectra to file
	outfile = run(
				l0,
				l1,
				res,
				obs_time,
				site,
				line_list,
				output_path='./outputs/', # path to save final spectrum
				config_path='./configs/', # path to save intermediate config files
				plot_path = './outputs/', # path so save plot
				mode='fits')              # save as fits



