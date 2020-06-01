import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii,fits
from scipy.stats import linregress
from astropy import constants
def q1():
	h=constants.h.value
	k=constants.k_B.value
	c=constants.c.value
# 2D array that has (lambda_min,lambda_max) for each filter under consideration
	band_definitions=np.array([[4000,5400],[5500,6850],[6700,8250]])
	data=ascii.read('phy270_ass1_spectrum.txt')
	lam=data['Wavelength'] 
	flux=data['Flux']
	fluxes=calc_photon_flux(lam,flux,band_definitions)
	print 'Photon flux: g=',fluxes[0],'r=',fluxes[1],'i=',fluxes[2]

	hdulist=fits.open("alpha_lyr_mod_002.fits")
	Vega=hdulist[1].data
	Vega_lambda=Vega['WAVELENGTH']
	Vega_flambda=Vega['FLUX']*0.001 # convert to W/m^2/s
	Vega_fluxes=calc_photon_flux(Vega_lambda,Vega_flambda,band_definitions)
	print Vega_fluxes
	gmag=-2.5*np.log10(fluxes[0]/Vega_fluxes[0])
	rmag=-2.5*np.log10(fluxes[1]/Vega_fluxes[1])
	imag=-2.5*np.log10(fluxes[2]/Vega_fluxes[2])
	print gmag,rmag,imag
	Tarr=np.arange(5000,20000,500)
	dl=5
	ll=np.arange(1000,20000,dl)
	llm=ll*1.e-10
	BBgmag=np.array([])
	BBrmag=np.array([])
	BBimag=np.array([])
	for T in Tarr:
		B=1.e-10*2*h*c*c/(llm**5)/(np.exp((h*c)/(llm*k*T))-1) # in W/m^2/A
		BB_fluxes=calc_photon_flux(ll,B,band_definitions)
		BBgmag=np.append(BBgmag,-2.5*np.log10(BB_fluxes[0]/Vega_fluxes[0]))
		BBrmag=np.append(BBrmag,-2.5*np.log10(BB_fluxes[1]/Vega_fluxes[1]))
		BBimag=np.append(BBimag,-2.5*np.log10(BB_fluxes[2]/Vega_fluxes[2]))
	fig,ax=plt.subplots()
	ax.scatter(BBgmag-BBrmag,BBrmag-BBimag,facecolor='None')
	ax.scatter(gmag-rmag,rmag-imag,c='r')
	plt.xlabel('(g-r)')
	plt.ylabel('(r-i)')
	plt.show()

def calc_photon_flux(l,f,band_definitions):
#Input: l = array of wavelengths, in Angstroms
#       f = corresponding array of fluxes, in W/m2/A
#	    band_definitions - lower and upper wavelength limits for N filters.
	h=constants.h.value
	c=constants.c.value
	output=np.array([])
# Determine number of fluxes to compute based on shape of band_definitions.	
	nfilter=np.shape(band_definitions)[0]
# Calculate wavelength interval.  
	for i in np.arange(0,nfilter):
		lamrange=np.where((l>band_definitions[i][0]) & (l<band_definitions[i][1]))
		ll=l[lamrange]
		dlambda=ll[1::]-ll[0:-1]
		dlambda=np.append(dlambda,dlambda[-1])
		# Calcualte photon_flux in photons/m2/s
		photon_flux=1./(h*c)*np.sum(l[lamrange]*f[lamrange]*dlambda)*1.e-20
		output=np.append(output,photon_flux)
	return output
