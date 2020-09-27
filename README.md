# run_psg

Instructions for running PSG using DOCKER
-----------------------------------------
See https://psg.gsfc.nasa.gov/helpapi.php#calling for more info

/Users/ashbake/Documents/Research/Projects/PSG/README.txt
See PSG API webpage to setup docker

Once docker setup:
------------------
docker start psg
(docker stop psg)

curl -d -type=ret --data-urlencode file@config.txt http://localhost:3000/api.php

Can go to http://localhost:3000 to install new packages

OR

curl http://localhost:3000/index.php?install=package

***Will need atmospheres package if want to access MERRA2 - e.g. if requesting telluric spectra at different observing times and want to update atmospheric profile***

***Will need lines package if want to access GEISA database - e.g. if comparing HITRAN line lists to GEISA***

***have not implemented user changes to target sky position key words and retrieval option currently not working***

Running run_psg
---------------
-Make sure psg_files/ folder exists or feed/edit paths to run_psg input

-See example at bottom of run_psg.py

