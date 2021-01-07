# example script to run PSG to get telluric spectrum
#http://bianca.palomar.caltech.edu:8000/maintenance/weather/user_gen_200_file.tcl
import sys

sys.path.append('psg/')
from run_psg import run_psg
from tools import *

if __name__=='__main__':
	# define everything
	output_path    = './outputs/'               # where to dump final spectrum
	config_path    = './configs/'               # where to dump PSG intermediate config files
	plot_path      = './outputs/'               # where to dump PSG final plot
	obs_time       = '2020/09/27 10:46:01.000'  # time of observation, must be in this format 
	#lat, lon, pres = 19.820664, -155.468066, 4.084  #Mauna Kea
	#lat,lon, pres  = 33.1504, 242.8173, 0.85   # palomar  lat lon (deg), surface pressure in bar
	lon, lat, pres  = 9.9158, 51.5413, 0.99     # gottingen
	l0, l1, res    = 500, 520, 100000           # wavelength range and resolving power
	line_list      = 'HIT'                      # 'HIT' or 'GEISA' - must have lines package, HIT more up to date

	# run psg, save telluric spectra to file
	site = [lon,lat,pres] # has to be this order
	outfile = run(
				l0,
				l1,
				res,
				obs_time,
				site,
				line_list,
				output_path=output_path,   # path to save final spectrum
				config_path=config_path,   # path to save intermediate config files
				plot_path = plot_path,     # path so save plot
				extension='fits')          # save as fits (currently only option)
