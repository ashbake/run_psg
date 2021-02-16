import numpy as np

fr = open('hitran/hitran01.par')
while True:
    line = fr.readline()
    if not line : break
    if 1:
        # Extraction of information from HITRAN file with special PSG format
        line = line.replace('#','0')
        sstr = line.split()
        mol=int(sstr[0]); iso=int(sstr[1]); wn=float(sstr[2]); lstr=float(sstr[3]);  ener   = float(sstr[4])
        self_g= float(sstr[5]);  self_n = float(sstr[6]);  self_d = float(sstr[7]);  self_t = float(sstr[8])
        air_g = float(sstr[9]);  air_n  = float(sstr[10]); air_d  = float(sstr[11]); air_t  = float(sstr[12])
        h2_g  = float(sstr[13]); h2_n   = float(sstr[14]); h2_d   = float(sstr[15]); h2_t   = float(sstr[16])
        he_g  = float(sstr[17]); he_n   = float(sstr[18]); he_d   = float(sstr[19])
        co2_g = float(sstr[20]); co2_n  = float(sstr[21]); co2_d  = float(sstr[22])
        sw1=float(sstr[23]); sw2=float(sstr[24]);
    else:
        # Extraction of information from old/classical HITRAN file
        delim = [2,1,12,10,10,5,5,10,4,8,15,15,15,15,19,7,7]; ipos=0; i=0
        sstr = np.empty(30,dtype='U30')
        for dl in delim:
            sstr[i]=line[ipos:ipos+dl]
            ipos=ipos+dl; i=i+1
        #End extracting values
        mol=int(sstr[0]); iso=int(sstr[1]); wn=float(sstr[2]); lstr=float(sstr[3]); ener=float(sstr[7])
        if lstr<1e-25: continue
        self_g= float(sstr[6]); self_n=0; self_d=0; self_t=0
        air_g = float(sstr[5]); air_n = float(sstr[8]); air_d = float(sstr[9]); air_t = 0
        h2_g = 0; h2_n = 0; h2_d = 0; h2_t = 0
        he_g = 0; he_n = 0; he_d = 0
        co2_g = 0; co2_n = 0; co2_d = 0
        sw1=float(sstr[15]); sw2=float(sstr[16]);
    # End readtype

    # Assign default
    if iso==0: iso=10
    if co2_g==0.0: co2_g = air_g
    if co2_n==0.0: co2_n = air_n
    if co2_d==0.0:
        co2_d = air_d
        co2_t = air_t
    if self_g==0.0: self_g = air_g
    if self_n==0.0: self_n = air_n
    if self_d==0.0:
        self_d = air_d
        self_t = air_t
    if h2_g==0.0: h2_g = air_g
    if h2_n==0.0: h2_n = air_n
    if h2_d==0.0:
        h2_d = air_d
        h2_t = air_t
    if he_g==0.0: he_g = air_g
    if he_n==0.0: he_n = air_n
    if he_d==0.0:
        he_d = air_d
        he_t = air_t

    # Save the results
    for k in [0,1]:
        if k==0: fout = 'data/hitran%02d' % mol
        else: fout = 'data/hitran%02d-%02d' % (mol,iso)
        fw=open('%s.txt' % fout,'a+'); fw.write('%12.6f %11.4e %10.4f\n' % (wn,lstr,ener)); fw.close()
        fw=open('%s.dat' % fout,'ab+')
        fw.write(np.ubyte(mol)); fw.write(np.ubyte(iso))
        fw.write(np.float64(wn)); fw.write(np.float64(lstr)); fw.write(np.float64(ener))
        fw.write(np.float32(self_g)); fw.write(np.float32(self_n)); fw.write(np.float32(self_d)); fw.write(np.float32(self_t))
        fw.write(np.float32(air_g));  fw.write(np.float32(air_n));  fw.write(np.float32(air_d));  fw.write(np.float32(air_t))
        fw.write(np.float32(h2_g));   fw.write(np.float32(h2_n));   fw.write(np.float32(h2_d));   fw.write(np.float32(h2_t))
        fw.write(np.float32(he_g));   fw.write(np.float32(he_n));   fw.write(np.float32(he_d));   fw.write(np.float32(he_t))
        fw.write(np.float32(co2_g));  fw.write(np.float32(co2_n));  fw.write(np.float32(co2_d));  fw.write(np.float32(co2_t))
        fw.write(np.float32(sw1));    fw.write(np.float32(sw2));
        fw.write(np.int16(0));     fw.write(np.int16(0));
        fw.close()
    #End save files
#End hitran line
fr.close()
