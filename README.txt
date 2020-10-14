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

(Can go to http://localhost:3000 to install new packages)
OR
curl http://localhost:3000/index.php?install=package

***Will need atmospheres package if want to access MERRA2 - e.g. if requesting telluric spectra at different observing times and want to update atmospheric profile***

***Will need lines package if want to access GEISA database - e.g. if comparing HITRAN line lists to GEISA***


Calculating Ephem 
-----------------
This keyword enables the calculation of ephemeris parameters based on the object name and the provided date:
wephm=y computes ephemeris for the date provided.
wephm=N computes ephemeris for the current date/time.
wephm=n turns it off

Calculating Geometry 
-----------------
The geometry module computes the observational angles (GEOMETRY-SOLAR-ANGLE, GEOMETRY-OBS-ANGLE) and the beam/planet ratio (GEOMETRY-PLANET-FRACTION) employing the geometry information provided by the user in the XML configuration file. For exoplanets, it also computes the planet transit factor (GEOMETRY-STAR-FRACTION), and the planet-star distance (GEOMETRY-STAR-DISTANCE). This information is then saved into the configuration file to be used by the other modules. The default is wgeo=y, yet one can disable this computation with wgeo=n.

How to run these:
----------------
a) Create a text file config.txt with this content:
<OBJECT-DATE>2017/01/15 14:30
<OBJECT>Mars
<GEOMETRY-REF>MRO
 
b) Call the API with this command:
curl -d type=cfg -d wephm=y -d watm=y --data-urlencode file@config.txt https://psg.gsfc.nasa.gov/api.php
 