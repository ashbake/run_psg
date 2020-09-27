##############################################################
# RUN PSG
#to. do
# figure out retrieval option, add more sites to store long/lat/alt for each observatory?

# Usage Notes
# The run_psg class currently takes a basic PSG config file (made with tellurics in mind) and can
# (1) return PSG config files with updated atmospheric profile sample
# (2) generate telluric transmission profile for given observing time
# (3) and will be able to run fits to data (not completed with this version tho)
###############################################################

import matplotlib
import os
import numpy as np
import matplotlib.pylab as plt
from datetime import datetime

plt.ion()

font = {'size'   : 14}
matplotlib.rc('font', **font)

def load_data(lam1=770,lam2=780,instrument='IAG',default=True,target=None,obs_num=None):
	"""
	load earth transmission data - currently picking one file from IAG telluric atlas

	lam1- lower wavelength limit (nm except uses microns for Parvi and Kpic- to do: add units)
	lam2 - upper wavelength limit
	mode - whether to load KittPeak, IAG, KPIC, or PARVI data
	default - if true, will load some Ashley-picked file from all the options, if False, must set target and obs_num
	target - target name and folder name for some instruments - check data
	obs_num - observation number or keyword to select the desired observation - see code and data

	returns [v,s,e] , tt - wavelength, flux, errors, time of observation


	define site info here


	"""
	if instrument =='IAG':
		if default:
			target = 1
			obs_num = '0163'

		filenames = glob.glob('%s/%s/*.fits' %(data_dir,instrument))
		f = fits.open('%s/%s/telluric_spectra_%s.fits' %(data_dir,instrument,target)) # target 1
		hdr = f[0].header
		jd = hdr['JD_%s'%obs_num] # obs_num: 0163

		# convert jd to time to put in config file
		tt = Time(jd,format='jd')

		v = 1e7/f[1].data['v'][::-1]
		s = f[1].data['telluric_%s'%obs_num][::-1]
		e = 10*np.abs(f[1].data['res_med'][::-1]-1)
		e[np.where(e > 0.1)] = 0.1

		isub = np.where((v > lam1) & (v < lam2))[0]
		
		return np.vstack([v[isub],s[isub],e[isub]]).T, tt.ymdhms
	
	if instrument =='IAG_NIR':
		pass

	if instrument=='KittPeak':
		f = np.loadtxt('%s/%s/transdata_0.5_1_mic'%(data_dir,instrument))

		v = 1e7/f[:,0][::-1]
		s = f[:,1][::-1]
		e = np.ones_like(s) * 0.01
		e[np.where(e > 0.1)] = 0.1

		isub = np.where((v > lam1) & (v < lam2))[0]

		tt =  Time(2451545,format='jd') # dummy time

		return np.vstack([v[isub],s[isub],e[isub]]).T, tt.ymdhms

	if instrument=='KPIC':
		if default:
			target = 'HR8799'
			obs_num = '200702_0100'

		filename = '%s/%s/%s/nspec%s_fluxes.fits' %(data_dir,instrument,target,obs_num)
		hdulist = fits.open(filename)
		header  = hdulist[0].header

		fluxes  = hdulist[0].data
		errors  = hdulist[1].data 
		slits   = hdulist[2].data
		darks   = hdulist[3].data

		# determine which fiber was on the star
		ibad = np.where(np.isnan(fluxes))
		flux_mod = fluxes * 1.0
		flux_mod[ibad] = 0
		flux_levels = np.median(np.median(flux_mod,axis=1),axis=1)
		istar = np.where(flux_levels == np.max(flux_levels))[0][0]

		wl = 1000*fits.getdata('%s/%s/%s/calib/20200702_HIP_81497_wvs.fits' %(data_dir,instrument,target))

		isub = np.where((wl[istar] > lam1) & (wl[istar] < lam2))

		w, s, e, _, _ = wl[istar][isub],fluxes[istar][isub],errors[istar][isub],slits[istar][isub],darks[istar][isub]
		if len(w) > 0:
			order = np.unique(isub[0])

		tt = Time(header['MJD'],format='mjd')

		bad = np.where(np.isnan(s))
		wgd   = np.delete(w,bad)
		sgd   = np.delete(s,bad)
		egd   = np.delete(e,bad)

		return np.vstack([wgd, sgd/np.max(sgd), egd/np.max(sgd)]).T, tt.ymdhms

	if instrument=='PARVI':
		if default:
			target = 'GJ229'
			obs_num = '20191119102854'

		filename = '%s/%s/%s/%s_R01_%s_deg0_sp.fits' %(data_dir,instrument, target, target, obs_num)

		hdulist = fits.open(filename)
		header  = hdulist[0].header

		science = hdulist[1].data
		sky     = hdulist[2].data
		lfc     = hdulist[3].data
		test    = hdulist[4].data

		wl, raw, rawerr, fluxflattened, errflattened, flat, flaterr, flux, err = science

		timeisot = obs_num[0:4] + '-' + obs_num[4:6] + '-' + obs_num[6:8] + 'T' + obs_num[8:10] + ':' + obs_num[10:12] + ':' + obs_num[12:]
		tt = Time(timeisot,format='isot')

		isub = np.where((wl > lam1) & (wl < lam2))
		w, s, e = wl[isub], raw[isub], rawerr[isub]
		if len(w) > 0:
			order = np.unique(isub[0])

		bad = np.where(np.isnan(s))
		wgd   = np.delete(w,bad)
		sgd   = np.delete(s,bad)
		egd   = np.delete(e,bad)

		return np.vstack([wgd, sgd/np.max(sgd), egd/np.max(sgd)]).T, tt.ymdhms



class run_psg():
	def __init__(self,
				 mode,
				 config,
				 data=None,
				 obs_time=None,
				 site=None,
				 line_list='HIT',
				 data_dir='../data',
				 config_path='./psg_files/',
				 output_path='./psg_files/',
				 plot_path='./psg_files/',
				 savedat=False,
				 ploton=False):
		"""
		data loading and config editing tools

		Inputs:
		-------
		mode: 'generate' or 'retrieve'

		Attributes:
		-----------

		"""
		self.mode = mode 
		self.config_path = config_path
		self.date = datetime.now().isoformat().replace(':','_') 
		self.output_name = output_path + 'run_%s.txt' %self.date

		self.gen_config(config, mode, obs_time, data, site) # creates self.config_to_run

		if mode=='generate':
			self.generate(self.config_to_run, self.output_name)

		elif mode=='retrieve':
			self.retrieve(self.config_to_run, self.output_name)
			if savedat:
				np.savetxt(output_path + 'datachunk_%s.txt' %self.date, 
								data,
								header='wavelength (nm) F_norm error')
		else:
			print('Mode not an option')

		if ploton:
			self.plot_from_file(self.output_name,mode=mode)

	def define_args(self,mode='generate',data=None, date=None, site=None):
		"""
		define the args dictionary containing data/key words to change in PSG config file
		before running - must hand edit for other variables

		Inputs:
		-------
		data - default 'same', or can be None (no data in config file), or 3 by N array of x, y, err as outputted by load_data fxn. the latter will replace data section in config
				or if mode=='generate' then data should be [lam1, lam2, unit, Resolving Power] ...[float, float, str, float]
		
		site - name of observing site location, currently can choose from pre-sets in PSG e.g. "Maunakea (Keck)" or "Mt Graham (LBT)" or "USA California"
				spelling and spaces must match exactly and 
				Future: Could upgrade in future to enable "User-Defined", but would need to pass long/lat/elevation
				Future: Could allow positioning to object by inputing zenith angle and azimuth, but haven't checked how much that matters (i.e. if earth PT profile actually changes or if just a simple airmass scaling)

		date - date of observation, used to pull satellite-informed PT profile. 
			   code will automatically run Merra2 (satellite query) to update atmosphere parameters if date is provided

		Output
		------
		args - dictionary to feed edit_psg_cong

		Note:
		-----
		this could use work,, currently most things are defaulted and if want to change things from base config referencing (like more molecules) will take extra work...tbd for easier solution
		ideas: have defaults here with args entry not repeating the key name (do formatting in second step)
		have each psg variable group (i.e. retrieval, geometry) be an input that the user can define 
		more simply as a dictionary that feeds into the args dictionary here that overrides defaults 
		so this fxn would mostly just do that formatting
		"""
		args = {}

		if date != None:
			args['OBJECT-DATE'] = ['<OBJECT-DATE>%s\n'%date]

		if site != None:
			# add site info
			#args['GEOMETRY-REF'] = ['<GEOMETRY-REF>%s\n'%site] # this doesn't work in API mode
			lon, lat, alt = site
			args['OBJECT-OBS-LONGITUDE']   = ['<OBJECT-OBS-LONGITUDE>%s\n'%lon]
			args['OBJECT-OBS-LATITUDE']    = ['<OBJECT-OBS-LATITUDE>%s\n'%lat]
			args['GEOMETRY-OBS-ALTITUDE']  = ['<GEOMETRY-OBS-ALTITUDE>=%s\n'%alt]
			args['GEOMETRY-ALTITUDE-UNIT'] = ['<GEOMETRY-ALTITUDE-UNIT>km\n']

		# generate mode is pretty decent as is
		if mode == 'generate':
			args['GENERATOR-INSTRUMENT'] = ['<GENERATOR-INSTRUMENT>user\n']
			args['GENERATOR-RANGE1']     = ['<GENERATOR-RANGE1>%s\n'%data[0]]
			args['GENERATOR-RANGE2']     = ['<GENERATOR-RANGE2>%s\n'%data[1]]
			args['GENERATOR-RANGEUNIT']  = ['<GENERATOR-RANGEUNIT>nm\n']
			args['GENERATOR-RESOLUTION']   = ['<GENERATOR-RESOLUTION>%s\n'%data[2]]
			args['GENERATOR-RESOLUTIONUNIT'] = ['<GENERATOR-RESOLUTIONUNIT>RP\n']

		# hand edit retrieval variables for now, but fix to user inputs in future?
		if mode=='retrieve':
			if np.all(data !=None):
				# update lambda range
				lam1, lam2 = data[0][0], data[-1][0] 
				args['GENERATOR-RANGE1'] = ['<GENERATOR-RANGE1>%s\n'%lam1] 
				args['GENERATOR-RANGE2'] = ['<GENERATOR-RANGE2>%s\n'%lam2] 

				# format data for config
				args['DATA'] = ['<DATA>\n']
				for i in range(len(data)):
					args['DATA'].append('%s %s %s\n'%(data[i,0],data[i,1],data[i,2]))
			#args['DATA'].append('</DATA>')
			args['RETRIEVAL-RESOLUTION']    = ['<RETRIEVAL-RESOLUTION>400000\n']
			args['RETRIEVAL-FITTELLURIC']   = ['<RETRIEVAL-FITTELLURIC>N\n']
			args['RETRIEVAL-FITSTELLAR']    = ['<RETRIEVAL-FITSTELLAR>N\n']
			args['RETRIEVAL-FITFREQ']       = ['<RETRIEVAL-FITFREQ>N\n']
			args['RETRIEVAL-FITGAIN']       = ['<RETRIEVAL-FITGAIN>N\n']
			args['RETRIEVAL-FITRESOLUTION'] = ['<RETRIEVAL-FITRESOLUTION>N\n']
			args['RETRIEVAL-NVARS']         = ['<RETRIEVAL-NVARS>1\n']
			args['RETRIEVAL-FLUXLABELS']    = ['<RETRIEVAL-FLUXLABELS>Data,Noise\n']
			args['RETRIEVAL-VARIABLES']     = ['<RETRIEVAL-VARIABLES>ATMOSPHERE-H2O\n']
			args['RETRIEVAL-VALUES']        = ['<RETRIEVAL-VALUES>0.4\n']
			args['RETRIEVAL-MIN']			= ['<RETRIEVAL-MIN>0.8\n']
			args['RETRIEVAL-MAX']			= ['<RETRIEVAL-MIN>6\n']
			args['RETRIEVAL-STATUS']        = ['<RETRIEVAL-STATUS>OK\n']

			#args['ATMOSPHERE-NGAS']   = ['<ATMOSPHERE-NGAS>2\n']
			#args['ATMOSPHERE-GAS']    = ['<ATMOSPHERE-GAS>H2O,O2\n']
			#args['ATMOSPHERE-TYPE']   = ['<ATMOSPHERE-TYPE>HIT[1],HIT[7]\n']
			#args['ATMOSPHERE-ABUN']   = ['<ATMOSPHERE-ABUN>0.4,1.8267\n']
			#args['ATMOSPHERE-UNIT']   = ['<ATMOSPHERE-UNIT>scl,scl\n']
			#args['ATMOSPHERE-TAU']    = ['<ATMOSPHERE-TAU>1,1e5\n']


		return args

	def edit_config(self,filename,outname=None,args=None,return_data=False,line_list='HIT'):
		"""
		Read the PSG config file to edit

		inputs:
		-------
		filename - name of config file to load and edit
		outname = if defined, the name of config file to save to
		args - if defined, dictionary with entries corresponding to config parameter data to replace
		return_data - default True, will read the data in old config file and return

		outputs:
		--------
		none, saves new config file to psg.new_config which is set to config_CURRENTTIME.txt
		"""
		f = open(filename,'r')
		lines = f.readlines()	
		f.close()

		if args != None:
			keys = args.keys()
			newlines = []
			#step through config file lines, replace if variable is defined in args
			for i,line in enumerate(lines):
				varname = line.split('>')[0].strip('<')
				if varname in keys:
					print('Editing ~~ %s ~~' %varname)
					for new_entry in args[varname]:
						newlines.append(new_entry) # new entries must match depending on the variable
				else:
					if line.startswith('<'):
						newlines.append(line) # use old line if not changing it
				if line.startswith('<ATMOSPHERE-TYPE>'): line_atm_type = i

		# swap line_list catalog
		if line_list=='HIT':
			newlines[line_atm_type] = newlines[line_atm_type].replace('GEISA',line_list)
		elif line_list=='GEISA':
			newlines[line_atm_type] = newlines[line_atm_type].replace('HIT',line_list)

		# save new config
		if outname != None:
			f = open(outname,'w')
			for line in newlines:
				f.write(line)
			f.close()
		
	def gen_config(self, config, mode, obs_time, data, site):
		"""
		generate config file in steps because watm=y will change everything to the defaults
		so must:
		1 - create config from base config with new date and location (temp config)
		2 - run watm=y to generate config with atmospheric profile with  (new atm config)
		3 - create config from new_atm_config with line_list choice update and same date,site options (config_)
		
		inputs: 
		config - name of base config file to edit 
		mode - [str] 'generate' or 'retrieval'
		obs_time - user input observation date/time -(see run_psg class requirements)
		data - user input data (see run_psg class requirements)
		site - user site info (long, lat, alt) in degrees and km (see run_psg class requirements)

		outputs
		None, saves config files, config_DATE.txt will be used to run psg on
		"""
		# define default new config names (need several temps bc of how psg runs)
		self.temp_config = self.config_path  + 'temp_config_%s.txt'%self.date
		self.new_atm_config = self.config_path  + 'new_atm_config_%s.txt'%self.date
		self.config_to_run = self.config_path  + 'config_%s.txt'%self.date

		# update basic general config file with new date, site, data ranges
		args = self.define_args(mode=mode,date=obs_time,data=data, site=site)
		self.edit_config(config, outname=self.temp_config, args=args)
		self.config(self.temp_config, self.new_atm_config)

		# run psg config generator to get new atm for that date, site, then update config with params
		args = self.define_args(mode=mode,date=obs_time, data=data, site=site)
		self.edit_config(self.new_atm_config, outname=self.config_to_run, args=args, line_list=line_list)

	def config(self,config_file,output_name):
		"""
		run merra2 and save results to config file

		move all config things here
		"""
		success = os.system('curl -d type=cfg -d wephm=y -d watm=y -d wgeo=y --data-urlencode file@%s http://localhost:3000/api.php > %s' %(config_file,output_name))

	def generate(self,config_file,output_name):
		"""
		generate spectrum from psg with config file

		must fix- aug 27
		"""
		# run job
		# can amke wephm, wgeo y to calc ephem or sky position- make this work plz
		success = os.system('curl -d type=trn -d watm=n -d weph=y --data-urlencode file@%s http://localhost:3000/api.php > %s' %(config_file,output_name))
		#success = os.system('curl -d type=cfg -d wephm=y -d watm=y -d wgeo=n --data-urlencode file@%s http://localhost:3000/api.php > %s' %(config_file,output_name))

		return success

	def retrieve(self,config_file,output_name):
		"""
		fit spectrum with config file

		**add watm=y? or assume watm was run before?
		"""
		# run job
		# can amke wephm, wgeo y to calc ephem or sky position- make this work plz
		success = os.system('curl -d type=ret --data-urlencode file@%s http://localhost:3000/api.php > %s' %(config_file,output_name))

		return success

	def read_retrieval_output(self,outfile):
		"""
		read output from PSG retrieval routine
		"""
		f = open(outfile,'r')
		lines = f.readlines()	
		f.close()

		i = 0
		line = lines[i]
		while line.startswith('results_dat.txt') is False:
			if line.startswith('<RETRIEVAL-FLUXLABELS>'):
				labels = line.strip('\n').split('>')[1].split(',')
				labels = ['lam'] + labels
			i+=1
			line = lines[i]
		# iterate past results_dat.txt line
		i+=1
		line = lines[i]

		# make data_out dic for lambda + all output columns
		data_out = {}
		for label in labels:
			data_out[label] = []

		while line.startswith('results_log.txt') is False:
			s = line.strip('\n').strip(' ').split('  ')
			for j, label in enumerate(labels):
				data_out[label].append(float(s[j].strip(' ')))
			i+=1
			line = lines[i]

		return data_out

	def plot_from_file(self,filename, mode='config'):
		"""
		plot data either from config file or from online data file

		input
		-----
		source: str
				options are 'config' if must pull data from config file, 
				'retrieval' if pulling from a retrieval output file, or
				'data' if pulling from text file with just data in it

		"""
		if mode=='retrieve':
			data = self.read_retrieval_output(filename)
			x, y, e, model = np.array(data['lam']), np.array(test['Data']), np.array(data['Noise']), np.array(data['Model'])

			f, (ax1, ax2) = plt.subplots(2, 1, sharex=True,figsize=(9,6))
			ax1.errorbar(x,y,e,label='data')
			ax1.plot(x,model,'k--',lw=1,zorder=100,label='model')
			
			ax2.plot(x,y-model,'r')
			ax2.plot(x,np.zeros_like(x),'k--')

			ax1.set_ylim(np.min(y)-0.05,1.1)
			ax2.set_ylim(-.5,0.5)

			ax2.set_xlabel('Wavelength (nm)')
			ax1.set_ylabel('Transmission')
			ax2.set_ylabel('O-M')

			f.subplots_adjust(bottom=0.15,left=0.15,hspace=0,right=0.9,top=0.9)
		
			self.data = data

		if mode=='generate':
			#filename = 'data_chunks/psg_retrieval_output_parvi_data_lam1_1474_lam2_1480_2019-11-18T10.16.01.000.txt'
			data = np.loadtxt(filename)
			n_entries = np.shape(data)[1]

			x, y = data.T[0:2] # make this flexible - load all entries (n_entries can change, use header info)
			
			self.x, self.y = x,y

			f, ax1 = plt.subplots(1, 1,figsize=(8,4))
			ax1.plot(x,y)

			ax1.set_ylim(np.min(y)-0.05,1.1)

			ax1.set_xlabel('Wavelength (nm)')
			ax1.set_ylabel('Transmission')

			f.subplots_adjust(bottom=0.15,left=0.15,hspace=0,right=0.9,top=0.9)


if __name__=='__main__':
	#lon, lat, alt = 33.1504, 242.8173, 1.871      #Palomar
	lon, lat, alt = -155.468066, 19.820664, 4.084 #Mauna Kea
	#### EXAMPLE USAGE FOR GENERATE MODE
	config_file = './config_gen.txt' # base configuration file downloaded from PSG - code modifies this to user's specs
	data     = [900,920,800000]      # lam1 (nm), lam2 (nm), resolving power of spectrum to retrieve (if R~10^6, must keep lambda range to ~10nm, but depends on computing power probs)
	obs_time = '2020/05/02 10:50'    # observing date and time in format 'YYYY/MM/DD HH:MM' - UT I think
	site     = [lon, lat, alt]       # Longitude (deg), latitude(deg), and altitude (km) of observer
	line_list= 'GEISA'                 # line list to use - either 'HIT' for hitran 2019 or 'GEISA' for geisa database (less updated)
	psg2 = run_psg('generate',config_file,
				data=data,
				obs_time=obs_time,
				line_list=line_list,
				site=site,
				ploton=True)

	
	#### EXAMPLE USAGE FOR RETRIEVE MODE
	config_file = './config_ret.txt'
	instrument = 'KittPeak' #options: KittPeak, IAG
	lam1, lam2 = 700,702
	data_dir = '../data'
	data, obs_time = load_data(lam1=lam1,lam2=lam2,instrument=instrument,default=True,target=None,obs_num=None)
	
	psg2 = run_psg('retrieve',config_file,
				data=data,
				obs_time=obs_time,
				line_list=line_list,
				site=site,
				ploton=True)



