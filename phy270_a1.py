import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii
from scipy.stats import linregress
from astropy import constants
def q2():
# Read in asigment data	and plot RA, Dec
	data=ascii.read('phy270_ass1_parallax.txt')
	fig,ax=plt.subplots()
	ax.scatter(data['RA'],data['Dec'])
	plt.axis([-430,400,-550,280])
	plt.xlabel('RA (mas)')
	plt.ylabel('Dec (mas)')

# Choose a few points that are separated by a year. 
	tsample=np.arange(3)*365.25+56
# Fit a line to the RA, Dec as a function of time to calculate the proper motion
	yRA=np.interp(tsample,data['Days'],data['RA'])
	yDec=np.interp(tsample,data['Days'],data['Dec'])
	plt.plot(yRA,yDec)
	rafit=linregress(tsample,yRA)
	decfit=linregress(tsample,yDec)
	pm=np.sqrt(rafit[0]*rafit[0]+decfit[0]*decfit[0])*365.25
	print 'Proper motion is ',str(pm), 'mas per year'
# Fit a line to RA as a function of Dec for the annually-separated points.  Then calculaate 
# the maximum perpendicular distance from this line to get the parallax	
	posfit=linregress(yRA,yDec)
	xmin=(posfit[0]*(data['Dec']-posfit[1])+data['RA'])/(1+posfit[0]*posfit[0])
	ymin=xmin*posfit[0]+posfit[1]
	r=np.sqrt((xmin-data['RA'])*(xmin-data['RA'])+(ymin-data['Dec'])*(ymin-data['Dec']))
	print 'Parallax is ',np.max(r)
	plt.show()

def q3(Tarr=np.array([8000,9000,9500,10000,11000])):
	#8190
	data=ascii.read('phy270_ass1_spectrum.txt')
	fig,ax=plt.subplots()
	plt.loglog(data['Wavelength'],data['Flux'])
#	ax.scatter()	
	dlambda=data['Wavelength'][1]-data['Wavelength'][0]
	irradiance=np.sum(data['Flux'])*dlambda
	print irradiance
	#T=8190
	h=constants.h.value
	k=constants.k_B.value
	c=constants.c.value
	lam=np.arange(1000,20000,5)*1.e-10
	for T in Tarr:
		B=1.e-10*2*h*c*c/(lam**5)/(np.exp((h*c)/(lam*k*T))-1) # in W/m^2/A
		renorm=np.average(B[np.where((lam*1.e10>10000)&(lam*1.e10<20000))])/np.average(data['Flux'][np.where((data['Wavelength']>10000)&(data['Wavelength']<20000))])
		plt.loglog(lam*1.e10,B/renorm,label='T='+str(T))
	distance=np.sqrt(20.*constants.L_sun.value/(4.*np.pi*irradiance))/3.086e16
	print 'Distance is ',distance,' pc'
	plt.axis([1000,20000,1e-15,1e-11])
	plt.xlabel('Wavelength $(\AA)$')
	plt.ylabel('Flux ($W/m^2/\AA$)')
	plt.legend(loc='best')
	plt.show()
    
print q3()