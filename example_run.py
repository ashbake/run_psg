# example script to run PSG to get telluric spectrum
#http://bianca.palomar.caltech.edu:8000/maintenance/weather/user_gen_200_file.tcl
path = '/Users/ashbake/Documents/Research/Projects/PSG/psg_main/'
import sys

sys.path.append(path + 'psg/')
from run_psg import run_psg
from tools import *



if __name__=='__main__':
	# define everything
	output_path    = './outputs/'               # where to dump final spectrum
	config_path    = './configs/'               # where to dump PSG intermediate config files
	plot_path      = './outputs/'               # where to dump PSG final plot
	obs_time       = '2020/12/18 10:07:16.000'  # time of observation, must be in this format 
	lon, lat, elev = 33.1504, 242.8173, 1.871   # palomar lon lat (deg), elevation in km
	l0, l1, res    = 3000, 3030, 1000000        # wavelength range and resolving power
	line_list      = 'HIT'                      # line list to use. 'HIT' or 'GEISA'
	
	site = [lon,lat,elev]
	# run psg, save telluric spectra to file
	outfile = run(
				l0,
				l1,
				res,
				obs_time,
				site,
				line_list,
				output_path=path + output_path,   # path to save final spectrum
				config_path=path + config_path,   # path to save intermediate config files
				plot_path = path + plot_path,     # path so save plot
				extension='fits',                 # save as fits (currently only option)
				cleanup=False)

