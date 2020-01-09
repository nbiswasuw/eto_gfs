# -*- coding: utf-8 -*-
"""

@author: Nishan Kumar Biswas
PhD Student and Graduate Research Assistant
Dept. of CEE, University of Washington

"""
import datetime, urllib, os, math
import numpy as np

class eto_from_gfs():
    def __init__(self):
        self.demfile = 'DEM_Pakistan.txt'
        self.invars = ['APCP', 'TMAX', 'TMIN', 'UGRD', 'VGRD', 'RH']
        self.row = 134
        self.col = 169
        self.xll = 60.9    
        self.yll = 23.7
        self.cell = 0.1
        self.currentdir = os.getcwd()
        self.todaydate = datetime.date.today()
        self.noforedays = 7
        self.forestartdate = datetime.datetime.combine(self.todaydate, datetime.time(0,0)) + datetime.timedelta(days=-1)
        self.Gsc = 0.0820		# Solar constant (MJ m^-2 min^-1)
        self.alpha = 0.23		# Albedo of reference crop (see FAO56)
        self.sigma = 4.903e-9		# Stefan-Boltzmann constant (MJ K^-4 m^-2 day^-1)
        self.Krs = 0.16		# Adjustment coeffecient for Rs (near coast ~0.19 far from coast ~0.16)

    def gfs_download(self):
        strdate = self.forestartdate.strftime("%Y%m%d")
        todayyr = self.forestartdate.strftime("%Y")
        todaymonth = self.forestartdate.strftime("%m")
        xmin = self.xll
        xmax = self.xll + self.cell*self.col
        ymin = self.yll
        ymax = self.yll + self.cell*self.row
        
        time = '00'
        print todayyr, todaymonth, strdate
        for forecasthr in range(0, (self.noforedays+1)*24+1, 3):
            fdate = strdate + '.f' + str(forecasthr).zfill(3)
            for vartype in self.invars:
                ## Data from NCEP NOaa
                path = r'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t' + time + 'z.pgrb2.0p25.f' + str(forecasthr).zfill(3) + '&var_' + vartype + '=on&subregion=&leftlon=' + str(xmin) + '&rightlon=' + str(xmax) + '&toplat=' + str(ymax) + '&bottomlat=' + str(ymin) + '&dir=%2Fgfs.' + strdate+ '%2F' + time
                print strdate, forecasthr, fdate
                fpath = self.currentdir + r'/Data/' + vartype + '/' +  fdate + '.' + vartype.lower() +'.gfs.raw.pakistan.tif'
                urllib.urlretrieve(path, fpath)
                urllib.urlcleanup()
    
    def gfs_gtiff_resample_ascii(self):
        strdate = self.forestartdate.strftime("%Y%m%d")
        xmin = self.xll
        xmax = self.xll + self.cell*self.col
        ymin = self.yll
        ymax = self.yll + self.cell*self.row
        
        for forecasthr in range(0, (self.noforedays+1)*24+1, 3):
            fdate = strdate + '.f' + str(forecasthr).zfill(3)
            for vartype in self.invars:
                fpath = self.currentdir + r'/Data/' + vartype + '/' +  fdate + '.' + vartype.lower() +'.gfs.raw.pakistan.tif'
                resamplepath =  self.currentdir + r'/Data/' + vartype + '/' +  fdate + '.' + vartype.lower() +'.gfs.pakistan.tif'
                asciipath =  self.currentdir + r'/Data/' + vartype + '/' +  fdate + '.' + vartype.lower() +'.gfs.pakistan.txt'
                print('C:\OSGeo4W64\OSGeo4W.bat gdalwarp -of GTiff -tr ' + str(self.cell) + ' ' + str(self.cell) + ' -te ' + str(xmin) + ' ' + str(ymin) + ' ' + str(xmax) + ' ' + str(ymax) + ' ' + fpath + ' ' + resamplepath)
                os.system('C:\OSGeo4W64\OSGeo4W.bat gdalwarp -of GTiff -tr ' + str(self.cell) + ' ' + str(self.cell) + ' -te ' + str(xmin) + ' ' + str(ymin) + ' ' + str(xmax) + ' ' + str(ymax) + ' ' + fpath + ' ' + resamplepath)
                os.system('C:\OSGeo4W64\OSGeo4W.bat gdal_translate -of AAIGrid -a_nodata -9999 ' + resamplepath + ' -b 1 ' + asciipath)
    
    def rh_max_min(self):
        strdate = self.forestartdate.strftime("%Y%m%d")
        minarray = np.zeros([self.row,self.col])
        maxarray = np.zeros([self.row,self.col])
        for forecasthr in range(3, (self.noforedays+1)*24+1, 24):
            print forecasthr
            for i in xrange(8):
                fpath1 = self.read_ascii(self.currentdir + r'/Data/RH/' +  strdate + '.f' + str(forecasthr + i*3).zfill(3) + '.rh.gfs.pakistan.txt')
                maxarray = np.maximum(np.asarray(fpath1), maxarray)
                minarray = np.minimum(np.asarray(fpath1), minarray)
            outfile1 = self.currentdir + r'/Data/RHMAX/' +  strdate + '.L' + str(int((forecasthr-3)/24.0)) + '.rhmax.gfs.pakistan.txt'
            outfile2 = self.currentdir + r'/Data/RHMIN/' +  strdate + '.L' + str(int((forecasthr-3)/24.0)) + '.rhmin.gfs.pakistan.txt'
            self.write_ascii(outfile1,maxarray)
            self.write_ascii(outfile2,minarray)
    
    def daily_mean(self):
        strdate = self.forestartdate.strftime("%Y%m%d")
        for var in self.invars:
            if var !='RH':
                meanarray = np.zeros([self.row,self.col])
                for forecasthr in range(3, (self.noforedays+1)*24+1, 24):
                    print forecasthr
                    meanarray = np.zeros([self.row,self.col])
                    for i in xrange(8):
                        fpath = self.read_ascii(self.currentdir + r'/Data/' + var + '/' +  strdate + '.f' + str(forecasthr + i*3).zfill(3) + '.' + var.lower() + '.gfs.pakistan.txt')
                        print var, forecasthr, len(fpath)
                        meanarray = np.asarray(fpath) + meanarray
                    outfile = self.currentdir + r'/Data/' + var + '/' +  strdate + '.L' + str(int((forecasthr-3)/24.0)) + '.' + var.lower() + '.gfs.pakistan.txt'
                    self.write_ascii(outfile,meanarray/8)
            

    def wind_calc(self):
        strdate = self.forestartdate.strftime("%Y%m%d")
        for i in range(0, self.noforedays + 1):
            print i
            uwind = self.read_ascii(r'Data/UGRD/'+ strdate + '.L' + str(i) + '.ugrd.gfs.pakistan.txt')
            vwind = self.read_ascii(r'Data/VGRD/'+ strdate + '.L' + str(i) + '.vgrd.gfs.pakistan.txt')
            print len(uwind), len(uwind[0])
            wind = np.asarray(uwind)*np.asarray(uwind) + np.asarray(vwind)*np.asarray(vwind)
#            print len(wind), len(wind[0])
            outfile = self.currentdir + r'/Data/WIND/' +  strdate + '.L' + str(i) + '.wind.gfs.pakistan.txt'
            self.write_ascii(outfile,np.sqrt(wind))

    def read_ascii(self, filepath):
        var = []
        fn = open(filepath, 'r')
        lines = fn.readlines()
        for inum in range(0, self.row):
            var1 = []
            line = lines[inum+6].split()
            for jnum in range(0, self.col):
                var1.append(float(line[jnum]))
            var.append(var1)
        return var
        
    def write_ascii(self, filepath, data):
        with open(filepath,'wb') as txt:
            txt.write("ncols    \t"+str(self.col)+"\n") 
            txt.write("nrows    \t"+str(self.row)+"\n")
            txt.write("xllcorner\t"+str(self.xll)+"\n")
            txt.write("yllcorner\t"+str(self.yll)+"\n")
            txt.write("cellsize \t"+str(self.cell)+"\n")
            txt.write("NODATA_value  -9999\n")
            np.savetxt(txt,data,fmt='%f')

    def julian_day(self, dt):
        tt = dt.timetuple()
        return tt.tm_yday

    def calc_grid_eto(self, Lat,Z,TH,TL,RHH,RHL,WS,J):
        M_PI = np.pi
        P = 101.3 * (((293.0 - 0.0068 * Z) / 293.0)**5.26)    # Pressure in kPa
        T = (TH + TL) / 2.0																				# Mean air temperature in C
        gamma = 0.665e-3 * P																				# Psychrometric constant (kPa C^-1)
        delta = 4098*(0.6108*np.exp(17.27*T/(T+237.3)))/((T+237.3)**2)										# Slope of sat. vapor pressure (kPa C^-1)
        eth = 0.6108 * np.exp(17.27 * TH / (TH + 237.3))														# Sturated vapor pressure in highest temp. (kPa)
        etl = 0.6108 * np.exp(17.27 * TL / (TL + 237.3))														# Sturated vapor pressure in lowest temp. (kPa)
        es = (eth + etl) / 2.0																				# Mean sturated vapor pressure in kPa
        ea = (eth * RHL + etl * RHH) / (2.0 * 100)															# Actual vapor pressure in kPa
        	
        Phi = M_PI * Lat / 180																				# Latitue in radian
        dr = 1 + 0.033 * np.cos(2 * M_PI * J / 365)															# Inverse relative Earth-Sun distance
        delta2 = 0.409 * np.sin(2 * M_PI * J/365 - 1.39)														# Solar declination
        ws = np.arccos(-np.tan(Phi) * np.tan(delta2))																	# Sunset hour angle
        Ra = 24*60*self.Gsc*dr*(ws*np.sin(Phi)*np.sin(delta2)+np.cos(Phi)*np.cos(delta2)*np.sin(ws))/M_PI						# Extraterrestrial radiation (MJ m^-2 day^-1)
        Rs = self.Krs * np.sqrt(TH - TL) * Ra																		# Solar or shortwave radiation (MJ m^-2 day^-1)
        Rso = (0.75 + 2e-5 * Z) * Ra																		# Clear sky solar radiation (MJ m^-2 day^-1)
        if (Rs > Rso):
            Rs = Rso
        Rns = (1 - self.alpha) * Rs																				# Net solar radiation (MJ m^-2 day^-1)
        Rnl = self.sigma*(pow((TH+273.16),4)+pow((TL+273.16),4))*(0.34-0.14*np.sqrt(ea))*(1.35*Rs/Rso-0.35)/2		# Net outgoing longwave radiation (MJ m^-2 day^-1)
        Rn = Rns - Rnl																						# Net radiation (MJ m^-2 day^-1)
        ETo = (0.408*delta*Rn+gamma*900*WS*(es-ea)/(T+273))/(delta+gamma*(1+0.34*WS))						# Evapotranspiration for reference crop in standard condition (mm day^-1)
        return(ETo)

    def eto_iteration(self):
        strdate = self.forestartdate.strftime("%Y%m%d")
        VAR =["tmax", "tmin", "wind", "rhmax", "rhmin"]
        dayn = self.last_week(self.forestartdate)
        Elev = self.read_ascii(self.demfile)
        Evap = np.zeros([self.row,self.col])
        EvapW = np.zeros([self.row,self.col])
        ## Main call for Daily ETo
        for j in range(self.noforedays+1):
            
            date = self.forestartdate + datetime.timedelta(days=j)
            print date
            julday = self.julian_day(date)
        
            for m in range(len(VAR)):   #for each variable
                filename = self.currentdir + "/Data/" + VAR[m].upper() + '/' + strdate + ".L" + str(j) + '.'+ VAR[m] + ".gfs.pakistan.txt"
                
                if(m==0):
                   Tmax = self.read_ascii(filename)
                elif(m==1):
                    Tmin = self.read_ascii(filename)
                elif(m==2):
                    WSP = self.read_ascii(filename)
                elif(m==3):
                    RHmax = self.read_ascii(filename)
                elif(m==4):
                    RHmin = self.read_ascii(filename)
                        
            # Call for ETo
            for m in range(self.row):
                lat = self.yll + (self.row - m)*self.cell - self.cell/2;
                for n in range(self.col):
                    
                    if Elev[m][n] == -9999:
                        Evap[m][n] = -9999
                        EvapW[m][n] = -9999
                    else:
                        Evap[m][n] = self.calc_grid_eto(lat,Elev[m][n],Tmax[m][n],Tmin[m][n],RHmax[m][n],RHmin[m][n],WSP[m][n],julday)
                        EvapW[m][n] = EvapW[m][n] + Evap[m][n]
            # Write daily ETo
            outfile = self.currentdir + "/Data/ET/" + strdate + '.L' + str(j) + ".pakistan.eto.forecast.txt"
            self.write_ascii(outfile,Evap)
            
        # Write weekly ETo
        date = datetime.datetime.strftime(dayn[0],"%Y%m%d")   # Name weekly file with the last day of forecasts (Lead 7)
        outfileweekly = self.currentdir + "/Data/ET/" + strdate + '.L' + str(j) + ".pakistan.eto.forecast.weekly.txt"
        self.write_ascii(outfileweekly,EvapW)

if __name__ == '__main__':
    forecast = eto_from_gfs()
#    forecast.gfs_download()
    forecast.gfs_gtiff_resample_ascii()
    forecast.rh_max_min()
    forecast.daily_mean()
    forecast.wind_calc()
    forecast.eto_iteration()
