Instructions for running PSG using DOCKER
-----------------------------------------
PSG API webpage: https://psg.gsfc.nasa.gov/helpapi.php#calling 

Basic instructions:
1. Setup PSG docker container (See PSG API webpage for instructions)
2. install any extra packages that you want (see note on packages below)
3. Edit example_run.py with parameters that you want
4. At command line, run from run_psg folder:
docker start psg
python example_run.py
docker stop psg   (if want to stop docker)
5. Use your new telluric spectrum


Note on Packages Needed:
-----------------------
To install packages, open terminal and type:

"curl --data-urlencode file@config.txt http://localhost:3000/api.php" with 'config.txt' as some local psg config file, then can go to http://localhost:3000 in browser to install new packages

OR to do it just with command line can run:

"curl http://localhost:3000/index.php?install=package"
Where 'package' is the name of the package to install

(See PSG API webpage for full instructions)

***Will need atmospheres package if want to access MERRA2 - e.g. if requesting telluric spectra at different observing time/location than in the default config file and want to update atmospheric profile succesfully***

***Will need lines package if want to access GEISA database - e.g. if comparing HITRAN line lists to GEISA***

Notes on Running:
-----------------------
- Running O3 seems to be too computationally expensive (at high res?) at low lambda<550nm. I take this molecule out
- All wavelengths and the resolution element (defined as d_lambda) are in nanometers. This can be changed if you wish by editing run_psg.py where I edit the config file to define the units to be nm
- If the resolution is too high, it fails, although in theory (by how it was described to me) it shouldn't. This is something I will follow up on. Generally a minimum d_lambda of 0.004 works across most of the spectrum. There is now a Resolution Kernel key word that seems like it should be turned on. Will double check this
- Rayleigh scattering does some weird stuff around 2 microns. I therefore generally don't include Rayleigh in the NIR since it's not crucial anyhow

