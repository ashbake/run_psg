# tools for running a psg sequence
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.time import Time
from scipy import interpolate
import inspect

import glob,os,sys

plt.ion()

sys.path.append('../scripts/')
from run_psg import run_psg


def psg_sequence(
	l0,
	l1,
	res,
	obs_time,
	site,
	line_list,
	config_path='./outputs/',
	output_path='./outputs/',
	plot_path='./outputs/',
	cleanup=False):
	"""
	"""
	psgpath = os.path.abspath(inspect.getfile(run_psg)) 
	config_basic = psgpath.strip('run_psg.py') + 'config_gen.txt' # base configuration file downloaded from PSG - code modifies this to user's specs

	if l1-l0 < 10:
		psg = run_psg('generate',config_basic,
					data=[l0,l1,res],
					obs_time=obs_time,
					line_list=line_list,
					site=site,
					config_path=config_path,
					output_path=output_path,
					plot_path=plot_path,
					ploton=False,
					run_atm=True)
		
		if i==0: pwv  = calc_pwv(site[2], psg.config_to_run)
		filelist = [psg.output_name]

		if cleanup:
			os.system('rm %s' %psg.temp_config)
			os.system('rm %s' %psg.new_atm_config)
		
	else:
		# step through wavelengths in 10nm chunks
		wls = np.arange(l0,l1+10,10)

		# save filenames to combine them later
		filelist = []

		for i in range(len(wls)-1):
			run_atm = True if i==0 else False
			config_input = config_basic if i==0 else psg.config_to_run
			psg = run_psg('generate',config_input,
						data=[wls[i],wls[i+1],res],
						obs_time=obs_time,
						line_list=line_list,
						site=site,
						config_path=config_path,
						output_path=output_path,
						plot_path=plot_path,
						ploton=False,
						run_atm=run_atm)
			
			if i==0: pwv  = calc_pwv(site[2], psg.config_to_run)
			if cleanup:
				if i>0:
					os.system('rm %s' %config_input) # remove old configtorun file
				elif i==0:
					os.system('rm %s' %psg.temp_config) # these are only made when i==0 when run_atm=true
					os.system('rm %s' %psg.new_atm_config)

			filelist.append(psg.output_name)
			print(i)

	if cleanup:
		os.system('rm %s' %psg.config_to_run) # rm last config to run file for either case

	return filelist, pwv

def read_output(file):
	"""
	read output from generating PSG telluric file
	"""
	f = open(file,'r')
	lines = f.readlines()	
	f.close()

	i = 0
	line = lines[i]
	while line.startswith('#'):
		line = lines[i] # will end after last # row, which is varnames row
		i+=1

	varnames = lines[i-2].strip('\n').strip('#').split() # could replace with row #...too lazy to count

	# make data_out dic for lambda + all output columns
	data_out = {}

	f = np.loadtxt(file)
	for i in range(len(f[0])):
		data_out[varnames[i]] = f[:,i]

	return data_out

def compile_segments(filelist):
	"""
	load all generated telluric spectra from run_sequence and save to one file
	"""

	# initiate array
	for i, file in enumerate(filelist):
		if i == 0:
			dat = read_output(file)
		else:
			temp_dic = read_output(file)
			for key in temp_dic.keys():
				if key != 'CIA': # CIA not in all columns so skip it
					dat[key] = np.concatenate((dat[key], temp_dic[key]))

	return dat

def extract_atm(config_file):
	"""
	read atm info from config file
	"""
	# extract water vapor info from config file
	f = open(config_file,'r')
	lines = f.readlines()	
	f.close()

	i = 0
	while lines[i].startswith('<ATMOSPHERE-LAYERS>') == False:
		i+=1

	weight = float(lines[i-4].strip('\n') .split('>')[1])
	nlayers = int(lines[i].strip('\n').split('>')[-1])
	molecules = lines[i-1][29:].strip('\n').split(',')

	atm = {}
	atm['weight'] = weight
	atm['pres'] = []
	atm['temp'] = []
	for molecule in molecules: atm[molecule] = []

	for j in range(i,i+nlayers):
		layer = np.array(lines[j+1].strip('\n').split('>')[1].split(','),dtype='float')
		atm['pres'].append(layer[0])
		atm['temp'].append(layer[1])
		for k,molecule in enumerate(molecules): atm[molecule].append(layer[k + 2])

	# pull out ZA
	iza = 0
	while lines[iza].startswith('<GEOMETRY-OBS-ANGLE>') == False:
		iza+=1

	ZA =  float(lines[iza].strip('\n').split('>')[1])
	atm['ZA'] = ZA

	return atm

def pres_at_height(elevation, m=28.97):
	"""
	get height in earth's atm given pressure 

	assuming constant temperature so pressure will be overestimated

	m: g/mol
	elevation: km
	"""
	atm_to_pascal = 101325 # N/m2/atm = kg /s2/m/atm
	g = 9.81 #m/s2
	pres_0 = 1
	T = 260 # K
	R = 8.31 #joule per kelvin (K)

	pres = pres_0 * np.exp(-1 * (m/1000) * g * elevation*1000 / R / T )	

	return pres

def calc_pwv(pres, config_file):
	"""
	calculate the PWV from the config file

	https://glossary.ametsoc.org/wiki/Precipitable_water#:~:text=The%20total%20atmospheric%20water%20vapor,the%20same%20unit%20cross%20section.

	check why this doesnt match pwv_frac_tot * column_mas (kg/m2)
	"""
	atm     = extract_atm(config_file)

	#p_surf = pres_at_height(elevation) # gives all pressures in file, but only want to integrate from altitube of site, so much calculate surface pressure
	p_surf=pres
	g = 9.81 #m/s2

	# locate index of surface pressure at elevation
	i_pres_surf = np.where(np.array(atm['pres']) < p_surf)[0][0]

	# integrate h2o column density starting at i_pres_surf
	int_h2o = np.abs(np.trapz(atm['H2O'][i_pres_surf:], x=atm['pres'][i_pres_surf:])) # integrated h2o column (atm)

	atm_to_pascal = 101325 # N/m2/atm = kg /s2/m/atm
	rho_g   = 9.81 * 1000 # m/s2 * kg/m3

	# weather equation approach from link in fxn header
	pwv_m = int_h2o * atm_to_pascal / rho_g
	pwv   = pwv_m * 1000 # convert m to mm

	airmass = airmass_from_ZA(atm['ZA'])

	# column mass approach -tbd go through derivation
	#col_mass = p_surf * atm_to_pascal/g
	#pwv = col_mass * int_h2o

	return round(pwv*airmass,2)

def airmass_from_ZA(ZA):
	"""
	calc airmass from ZA given in deg
	"""
	return 1/np.cos(ZA * np.pi/180)

def plot_tellurics(l0,l1,res,outfile, plot_path):
	"""
	quick plot of results

	CIA is weird for some reason
	"""

	f = fits.open(outfile)
	fdat = f[1].data
	f.close()

	wl_key = f[1].columns[0].name
	mol_keys = f[1].columns[1:]

	plt.figure()
	for mol_key in mol_keys:
		if mol_key.name == 'CIA': continue
		plt.plot(fdat[wl_key], fdat[mol_key.name], label=mol_key.name)

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')

	plt.savefig(plot_path + 'plot_%s_%s_R%s.png'%(l0,l1,res))

def gen_final_spec_name(date, l0, l1, lon, lat, pres):
	"""
	spectrum save name fxn
	date =  obs_time[0:10] (just yyyy/mm/dd)

	"""
	filename = 'psg_out_%s_l0_%snm_l1_%snm_lon_%s_lat_%s_pres_%s.fits' %(date.replace("/", '.'),l0,l1,lon,lat,pres)
	return filename

def save_dat(
	dat, 
	pwv, 
	l0,
	l1,
	res,
	obs_time,
	site,
	line_list,
	output_path='outputs/',
	extension='fits'):
	"""
	save compiled data dictionary as either 'fits' or 'text' or 'tapas'
	'tapas' will save them as close to TAPAS format as possible
	"""
	date =  obs_time[0:10]
	lon, lat, pres = site
	final_spec_name = gen_final_spec_name(date, l0, l1, lon, lat, pres)
	filename = output_path + final_spec_name

	if extension=='fits':
		cols = []
		for key in dat.keys():
			cols.append(fits.Column(name=key, format='E', array=np.array(dat[key])))

		tbhdu = fits.BinTableHDU.from_columns(cols)

		# Header info
		hdr = fits.Header()
		hdr['CV']   = pwv
		hdr['PWV']  = pwv
		hdr['PWV_unit']  = 'mm'
		hdr['OBSTIME']   = obs_time
		hdr['LINELIST']  = line_list
		hdr['L0']   = l0
		hdr['LF']   = l1
		hdr['R']    = res
		hdr['LAMUNIT']    = 'nm'
		hdr['LON']  = lon
		hdr['LAT']  = lat
		hdr['PRES'] = pres

		primary_hdu = fits.PrimaryHDU(header=hdr)

		hdu = fits.HDUList([primary_hdu, tbhdu])
		if os.path.exists(filename):
			os.system('rm %s' % filename)
			hdu.writeto(filename)
		else:
			hdu.writeto(filename)
	else:
		pass
		# tbd text file option

	return filename

def run(
	l0,
	l1,
	res,
	obs_time,
	site,
	line_list,
	config_path='./outputs/',
	output_path='./outputs/',
	plot_path='./outputs/',
	extension='fits',
	cleanup=False):
	"""
	run full fit
	"""
	# run psg sequence - save to files
	filelist, pwv = psg_sequence(
					l0,
					l1,
					res,
					obs_time,
					site,
					line_list,
					config_path=config_path,
					output_path=output_path,
					plot_path=plot_path,
					cleanup=cleanup)

	# open files, compile all results to dat dictionary
	dat = compile_segments(filelist)

	# save, pick save mode
	outfile = save_dat(
			dat, 
			pwv, 
			l0,
			l1,
			res,
			obs_time,
			site,
			line_list,
			output_path=output_path,
			extension=extension)

	if cleanup:
		for file in filelist:
			os.system('rm %s' %file)

	plot_tellurics(l0,l1,res,outfile,plot_path)

	return outfile

