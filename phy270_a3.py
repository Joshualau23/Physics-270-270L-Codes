import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii,fits
from scipy.stats import linregress
from astropy import constants
import matplotlib.ticker as mticker
h=constants.h.value
k=constants.k_B.value
c=constants.c.value
def q3():
	wavelength=np.arange(300,900)
	Diameters=np.array([3.5,0.2,30,6.5])
	r0array=np.array([15,10,0.,0])
	labels=np.array(['A','B','C','D'])
	fig,ax=plt.subplots()
	for D,r0,txt in zip(Diameters,r0array,labels):
		print D,r0
		resolution=getres(D,r0,wavelength,0.)
		plt.semilogy(wavelength,resolution,label=txt)
		plt.legend(loc='best')
	ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
	ax.yaxis.get_major_formatter().set_scientific(False)
	ax.yaxis.get_major_formatter().set_useOffset(False)
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Resolution (arcsec)')
	plt.show()
def getres(D,r0,l,zeta):
	diffraction_limit=2*1.22*l*1.e-9/D*206265
	if r0>0:
		seeing_limit=1.0/(r0/10.)*(l/500.)**(-0.2)#*(np.cos(zeta))**(-3./5)
	else:
		seeing_limit=l*0.
	resolution=np.maximum(diffraction_limit,seeing_limit)
	return resolution
def q4():
	Tarr=np.array([293.,243.])
	dl=5
	ll=np.arange(10000,100000,dl)
	fig,ax=plt.subplots()
	for T in Tarr:
		AB=getBB(ll,T)
		plt.plot(ll,AB,label='T='+str(T)+'K')
		plt.xlabel('Wavelength (A)')
		plt.ylabel('AB magnitude')
		plt.legend(loc='best')
	plt.show()
	print getBB(22000,50),getBB(100000,50)
def getBB(ll,T):
	llm=ll*1.e-10
	B=1.e-10*2*h*c*c/(llm**5)/(np.exp((h*c)/(llm*k*T))-1) # in W/m^2/A
	nu=c/llm
	log10Bnu=np.log10(2*h)+3.*np.log10(nu)-2.*np.log10(c)-np.log10(np.exp(h*nu/(k*T))-1) # in W/m^2/Hz
	log10Bnu=log10Bnu-10.628 # in W/m^2/Hz/arcsec2
	AB=-2.5*log10Bnu-56.1
	return AB

def calcJK5(l,f):
	Jlim=np.array([10000,12000])
	Klim=np.array([20000,25000])
	Fivelim=np.array([50000,51000])
	Jrange=np.where((l>Jlim[0])&(l<Jlim[1]))
	Krange=np.where((l>Klim[0])&(l<Klim[1]))
	Fiverange=np.where((l>Fivelim[0])&(l<Fivelim[1]))
	dlambda=l[1]-l[0]
	Jmag=-2.5*np.log10(np.sum(l[Jrange]*f[Jrange])*dlambda)
	Kmag=-2.5*np.log10(np.sum(l[Krange]*f[Krange])*dlambda)
	Fivemag=-2.5*np.log10(np.sum(l[Fiverange]*f[Fiverange])*dlambda)
	return (Jmag,Kmag,Fivemag)

print q3()