# example script to run PSG to get telluric spectrum
# higher resolutions now fail - contact Geronimo about it
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
	obs_time       = '2015/06/12 16:17:49.661'  # time of observation, must be in this format 
	#lon, lat, pres = 204.53, 19.82, 0.5826     # mauna kea lon lat (deg), pres in km
	lon, lat, pres = 9.9158, 51.5413, 0.98468   # goettingen
	l0, l1, res    = 510, 520, 0.004          # wavelength range and resolving power
	line_list      = 'HIT'                      # line list to use. 'HIT' or 'GEISA'
	
	site = [lon,lat,pres]
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

