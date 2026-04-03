# example script to run PSG to get telluric spectrum
# higher resolutions now fail - contact Geronimo about it
#http://bianca.palomar.caltech.edu:8000/maintenance/weather/user_gen_200_file.tcl
import sys

# add psg to path
sys.path.append('./psg/')

from tools import pick_site, run


if __name__=='__main__':
	# define everything
	output_path    = './outputs/'  		        # where to dump final spectrum
	config_path    = './configs/'               # where to dump PSG intermediate config files
	plot_path      = output_path                # where to dump PSG final plot
	obs_time       = '2016/02/27 11:46:07.8'    # time of observation, must be in this format 
	site           = pick_site(sitename='gottingen')
	l0, l1, res    = 1000, 1200, 0.001           # wavelength range and resolving power
	config_name    = 'psg_cfg_20160227_obsnum0025_IAG_solar.txt'#'psg_cfg_20150617_obsnum0041_IAG_solar.txt'# name of config to load. if none, will load default. Must be in config_path folder

	# run psg, save telluric spectra to file
	# TODO make config file an input and make run part of a class
	outfile = run(
				l0,
				l1,
				res,
				obs_time,
				site,
				'HIT',
				output_path= output_path,  # path to save final spectrum
				config_path=config_path,   # path to save intermediate config files
				plot_path = plot_path,     # path so save plot
				extension='fits',          # save as fits (currently only option)
				cleanup=True,			   # deletes intermediate files
				run_atm=False,             # do or don't regenerate atm
				config_name = config_name) 	   # name of config to load. if none, will load default		   

