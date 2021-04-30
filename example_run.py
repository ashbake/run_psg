# example script to run PSG to get telluric spectrum
# higher resolutions now fail - contact Geronimo about it
#http://bianca.palomar.caltech.edu:8000/maintenance/weather/user_gen_200_file.tcl
path = '/Users/ashbake/Documents/Research/Projects/PSG/psg_main/'
import sys

sys.path.append(path + 'psg/')
from tools import pick_site, run



if __name__=='__main__':
	# define everything
	output_path    = './outputs/'  		        # where to dump final spectrum
	config_path    = './configs/'               # where to dump PSG intermediate config files
	plot_path      = output_path                # where to dump PSG final plot
	obs_time       = '2021/04/20 10:50:49.661'  # time of observation, must be in this format 
	site           = pick_site(sitename='palomar')
	l0, l1, res    = 1100, 1900, 0.001          # wavelength range and resolving power


	# run psg, save telluric spectra to file
	outfile = run(
				l0,
				l1,
				res,
				obs_time,
				site,
				'HIT',
				output_path=path + output_path,   # path to save final spectrum
				config_path=path + config_path,   # path to save intermediate config files
				plot_path = path + plot_path,     # path so save plot
				extension='fits',                 # save as fits (currently only option)
				cleanup=True,					  # deletes intermediate files
				run_atm=True) 					  # do or don't regenerate atm
