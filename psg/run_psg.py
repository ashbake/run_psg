##############################################################
# RUN PSG
#to. do
# figure out retrieval option
#add more sites to store long/lat/alt for each observatory? 
# make README helpful

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
				 run_atm=True,
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
		self.run_atm=run_atm

		if mode not in ['retrieve', 'generate']: raise TypeError

		self.config_path = config_path
		self.date = datetime.now().isoformat().replace(':','_') 
		self.output_name = output_path + 'run_%s.txt' %self.date
		self.config_to_run = self.config_path  + 'config_%s.txt'%self.date #name of config to run on

		self.gen_config(config, mode, obs_time, data, site, line_list, run_atm=run_atm) # creates self.config_to_run

		if mode=='generate':
			self.generate(self.config_to_run, self.output_name)

		elif mode=='retrieve':
			self.retrieve(self.config_to_run, self.output_name)
			if savedat:
				np.savetxt(output_path + 'datachunk_%s.txt' %self.date, 
								data,
								header='wavelength (nm) F_norm error')
		else:
			raise TypeError

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
		
		site - presets dont work in docker framework - input lon, lat, pres as tuple (deg, deg, bar)

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
			lon, lat, pres = site
			args['OBJECT-OBS-LONGITUDE']   = ['<OBJECT-OBS-LONGITUDE>%s\n'%lon]
			args['OBJECT-OBS-LATITUDE']    = ['<OBJECT-OBS-LATITUDE>%s\n'%lat]
			if self.run_atm:
				args['GEOMETRY-OBS-ALTITUDE']  = ['<GEOMETRY-OBS-ALTITUDE>=%s\n'%0.0] # set alt to 0 (ground level)
				args['GEOMETRY-ALTITUDE-UNIT'] = ['<GEOMETRY-ALTITUDE-UNIT>km\n']
				args['ATMOSPHERE-PRESSURE']    = ['<ATMOSPHERE-PRESSURE>%s\n'%pres] #update w/ user pressure
				args['ATMOSPHERE-PUNIT']       = ['<ATMOSPHERE-PUNIT>bar\n']

		# generate mode is pretty decent as is
		if mode == 'generate':
			args['GENERATOR-INSTRUMENT'] = ['<GENERATOR-INSTRUMENT>user\n']
			args['GENERATOR-RANGE1']     = ['<GENERATOR-RANGE1>%s\n'%data[0]]
			args['GENERATOR-RANGE2']     = ['<GENERATOR-RANGE2>%s\n'%data[1]]
			args['GENERATOR-RANGEUNIT']  = ['<GENERATOR-RANGEUNIT>nm\n']
			args['GENERATOR-RESOLUTION']   = ['<GENERATOR-RESOLUTION>%s\n'%data[2]]
			args['GENERATOR-RESOLUTIONUNIT'] = ['<GENERATOR-RESOLUTIONUNIT>nm\n']

		# hand edit retrieval variables for now, but fix to user inputs in future?
		if mode=='retrieve':
			if np.all(data !=None):
				# update lambda range
				x, y, e, res = data[0], data[1], data[2], data[3]
				lam1, lam2 = np.min(x), np.max(x)
				args['GENERATOR-RANGE1'] = ['<GENERATOR-RANGE1>%s\n'%lam1] 
				args['GENERATOR-RANGE2'] = ['<GENERATOR-RANGE2>%s\n'%lam2] 

				# format data for config
				args['DATA'] = ['<DATA>\n']
				for i in range(len(y)):
					args['DATA'].append('%s %s %s\n'%(x[i],y[i],e[i]))
			#args['DATA'].append('</DATA>')
			args['RETRIEVAL-RESOLUTION']    = ['<RETRIEVAL-RESOLUTION>%s\n'%res]
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

	def edit_config(self,filename,outname,args=None,return_data=False,line_list='HIT'):
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

		if len(lines)==0:
			raise ValueError

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
		
		else:
			np.savetxt(outname, lines)

	def gen_config(self, config, mode, obs_time, data, site, line_list, run_atm=True):
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
		# update basic general config file with new date, site, data ranges for atm generation
		if run_atm:
			self.temp_config = self.config_path  + 'temp_config_%s.txt'%self.date
			self.new_atm_config = self.config_path  + 'new_atm_config_%s.txt'%self.date

			args = self.define_args(mode=mode,date=obs_time,data=data, site=site)
			self.edit_config(config, self.temp_config, args=args)
			# run psg config generator to get new atm for that date, site, then update config with params
			self.config(self.temp_config, self.new_atm_config)

			# re-edit config to user-settings that were erased in psg run (e.g. like line list pref)
			args = self.define_args(mode=mode,date=obs_time, data=data, site=site)
			self.edit_config(self.new_atm_config, self.config_to_run, args=args, line_list=line_list)
		
		else:
			# re-edit config to user-settings that were erased in psg run (e.g. like line list pref)
			args = self.define_args(mode=mode, date=obs_time, data=data, site=site)
			self.edit_config(config, self.config_to_run, args=args, line_list=line_list)


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
		#success = os.system('curl -d type=trn -d watm=n -d weph=y --data-urlencode file@%s https://psg.gsfc.nasa.gov/api.php > %s' %(config_file,output_name))
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
			x, y, e, model = np.array(data['lam']), np.array(data['Data']), np.array(data['Noise']), np.array(data['Model'])

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





