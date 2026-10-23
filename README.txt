Instructions for running PSG using DOCKER
-----------------------------------------
PSG API webpage: https://psg.gsfc.nasa.gov/helpapi.php#calling


Overview:
-----------------------
This code generates telluric (Earth atmosphere) transmission spectra with the
Planetary Spectrum Generator (PSG) running locally in a Docker container. Given a
wavelength range, wavelength sampling, observation time, and observing site, it:
  - edits a base PSG config file with your settings
  - optionally queries MERRA2 (through PSG) for the atmospheric profile at that time/site
  - calls the local PSG API (http://localhost:3000/api.php) with curl to compute the
    transmittance, stepping through the wavelength range in 5 nm chunks
  - stitches the chunks together, computes PWV and airmass from the atmospheric
    profile, and saves a single FITS file plus a quick-look plot


Directory Layout:
-----------------------
example_run.py          - script to edit and run; sets all user parameters and calls tools.run()
psg/tools.py            - high-level functions: pick_site(), run(), chunking, stitching,
                          PWV/airmass calculation, saving to FITS, plotting
psg/run_psg.py          - run_psg class: edits PSG config files and sends them to the
                          PSG API via curl ('generate' mode works; 'retrieve' mode is unfinished)
psg/example_psg_customatm.txt - default base config used if no config_name is given
psg/example_psg_cfg.txt - another example PSG config file
configs/                - user config files to load, and where intermediate configs are written
outputs/                - final FITS spectra, plots, and intermediate PSG output files
make_linelist/          - tools for building custom HITRAN line lists for PSG
                          (see make_linelist/README.txt)


Requirements:
-----------------------
- Docker with the PSG container installed and named "psg", serving on localhost:3000
- curl (used to talk to the PSG API)
- Python with numpy, scipy, matplotlib, astropy


Basic instructions:
1. Setup PSG docker container (See PSG API webpage for instructions)
2. install any extra packages that you want (see note on packages below)
3. Edit example_run.py with parameters that you want
4. At command line, run from run_psg folder:
docker start psg
python example_run.py
docker stop psg   (if want to stop docker)
5. Use your new telluric spectrum


Parameters in example_run.py:
-----------------------
output_path  - folder for the final spectrum (default ./outputs/)
config_path  - folder for intermediate PSG config files, and where config_name is loaded from
plot_path    - folder for the quick-look plot
obs_time     - observation time, format 'YYYY/MM/DD HH:MM:SS.S' (UT)
site         - [lon, lat, pressure] in degrees, degrees, bar. Use pick_site(sitename=...)
               with one of 'maunakea', 'palomar', 'gottingen', 'kittpeak' (pres can be
               overridden), or pick_site(sitename=None, lon=..., lat=..., pres=...) for a custom site
l0, l1       - start and end wavelength in nm
res          - resolution element d_lambda in nm (not resolving power)
config_name  - name of a config file in config_path to use as the base config.
               If None, psg/example_psg_customatm.txt is used

Options passed to run():
line_list    - 'HIT' (HITRAN) or 'GEISA'; swaps the line list in the config's ATMOSPHERE-TYPE
extension    - 'fits' (currently the only option)
cleanup      - True deletes intermediate config files and per-chunk PSG outputs
run_atm      - True has PSG query MERRA2 for the atmospheric profile at obs_time/site
               (requires the atmospheres package, see below). False uses the atmosphere
               already in the base config file as-is


How it works:
-----------------------
tools.run() calls:
1. psg_sequence() - splits l0..l1 into 5 nm chunks and creates a run_psg object for each.
   For each chunk run_psg:
     a. edits the config (date, site, wavelength range, resolution in nm, line list)
     b. if run_atm=True (first chunk only), sends the config to PSG with watm=y to pull
        the MERRA2 atmosphere, then re-applies the user settings that PSG overwrote
     c. sends the config to PSG with type=trn and saves the transmittance output
        to outputs/run_<timestamp>.txt
   Later chunks reuse the config from the previous chunk, so the atmosphere is only
   generated once. PWV and airmass are computed from the first chunk's config.
2. compile_segments() - reads and concatenates the per-chunk outputs
3. save_dat() - writes the combined spectrum to FITS
4. plot_tellurics() - plots the transmittance of each component and saves a PNG


Outputs:
-----------------------
outputs/psg_out_<date>_l0_<l0>nm_l1_<l1>nm_res_<res>nm_lon_<lon>_lat_<lat>_pres_<pres>.fits
  - extension 1: table with a wavelength column (nm) and transmittance columns
    for the total and each molecule/component returned by PSG
  - primary header: PWV (mm), AIRMASS, OBSTIME, LINELIST, L0, LF, R, LAMUNIT, LON, LAT, PRES
outputs/plot_<l0>_<l1>_R<res>.png - quick-look plot of the transmittance


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
- To change which molecules are included (e.g. removing O3 or Rayleigh), edit the base config file you are loading (config_name or psg/example_psg_customatm.txt)
- With run_atm=False, the observation time and site only set the geometry; the atmospheric profile comes from the base config file


Custom Line Lists:
-----------------------
See make_linelist/README.txt for how to convert a HITRAN .par file with make_linelist/convert.py
and copy the resulting line list into the PSG docker container.
