Instructions for Creating New Line List
---------------------------------------


1. Download HITRAN .par file according to the display of selections shown in convert.png - I have this format saved as PSG on my hitran account

2. If desired, edit line parameters in .par file (by hand or write a script - hapi.py useful if going the script route)

3. Run convert.py on the .par download file (edit filename at top of convert.py)

4. Move the output file to /var/www/html/data/linelists/hitran (tbd on this - maybe make own folder) by running: 
""docker cp foo.dat psg:/var/www/html/data/linelists/hitran/data/foo.dat""

foo.dat will be in format "hitran##-##.dat" with [molecule#]_[isotope#] e.g. hitran01_01.dat for 1 isotopologue of water. The script will generate this filename automatically regardless of input name

5. Run: ""docker exec -t -i psg /bin/bash"" to access the Docker container - check date of "/var/www/html/data/linelists/hitran/data/foo.dat" to confirm that file updated as expected

6. Rerun with new line lists! :) 




Test (Oct 19, 2020):
--------------------
Edit this line:
Change from 6.89e-27 to 6.89e-22
1 1 11790.68986 6.89e-27 326.6255 0.3 # # # 0.0962 0.67 -0.012986 # # # # # # # # # # # # 9 11 ElecStateLabel=X;v1=2;v2=3;v3=0;J=4;Ka=2;Kc=2 ElecStateLabel=X;v1=0;v2=0;v3=0;J=5;Ka=1;Kc=5

docker cp hitran01-01.dat psg:/var/www/html/data/linelists/hitran/data/hitran01-01.dat
docker cp hitran01.dat psg:/var/www/html/data/linelists/hitran/data/hitran01.dat


Did it change anything? -Compare to PSG online spectrum. Wavelength around 848.1268nm

