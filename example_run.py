# example script to run PSG to get telluric spectrum
import sys

sys.path.append('psg/')
from run_psg import run_psg
from tools import *

if __name__=='__main__':
	# define everything
	output_path    = './outputs/'              # where to dump final spectrum
	config_path    = './configs/'              # where to dump PSG intermediate config files
	obs_time       = '2019/11/18 10:07:16.000' # time of observation, must be in this format 
	lon, lat, elev = 33.1504, 242.8173, 1.871  # palomar lon lat (deg), elevation in km
	l0, l1, res    = 700, 720, 100000          # wavelength range and resolving power
	line_list      = 'HIT'                     # 'HIT' or 'GEISA' - must have lines package, HIT more up to date

	# run psg, save telluric spectra to file
	outfile = run(
				l0,
				l1,
				res,
				obs_time,
				[lon,lat,elev],
				line_list,
				output_path=output_path,  # path to save final spectrum
				config_path='./configs/', # path to save intermediate config files
				plot_path = output_path,  # path so save plot
				mode='fits')              # save as fits (currently only option)



