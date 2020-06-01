import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.io import ascii
from astropy import constants
h=constants.h.value
k=constants.k_B.value
c=constants.c.value
import random
import scipy
import matplotlib.ticker as mticker

def q1a(name='fgd71',nspec=5):
	for i in np.arange(nspec):
		data=ascii.read(name+'_'+str(i+1)+'.dat')
		try:
			fluxes=np.dstack((fluxes,data['Flux']))
			errors=np.dstack((errors,data['Uncertainty']))
		except:
			fluxes=data['Flux']
			errors=data['Uncertainty']
			wavelength=data['Wavelength']
	fluxes=fluxes[0]
	errors=errors[0]
	average=np.average(fluxes,axis=1)
	weight=1./(errors*errors)
	waverage=np.average(fluxes,weights=weight,axis=1)
	filtered_data=sigma_clip(fluxes,sigma=2.0,iters=None,axis=1)
	clipaverage=np.ma.average(filtered_data,axis=1)
	clipwaverage=np.ma.average(filtered_data,weights=weight,axis=1)
	uncertainty_on_stack=np.sqrt(1./(np.sum(1./(errors*errors),axis=1)))
	fig,ax=plt.subplots()
	plt.loglog(wavelength,wavelength*average,label='Average')
	plt.loglog(wavelength,wavelength*waverage,label='Weighted Average')
	plt.loglog(wavelength,wavelength*clipaverage,label='Clipped Average')
	plt.loglog(wavelength,wavelength*clipwaverage,label='Clipped, weighted Average')
	plt.legend(loc='best')
	plt.xlabel('Wavelength (A)')
	plt.ylabel('$\lambda F (10^{13} W/m^2$)')
	fig.show()
	fig2,ax2=plt.subplots()
	plt.plot(wavelength,uncertainty_on_stack/clipwaverage)
	plt.xlabel('Wavelength (A)')
	plt.ylabel('Relative Uncertainty')
	fig2.show()

def q1b(name='fgd71',nspec=10):
	for i in np.arange(nspec):
		data=ascii.read(name+'_'+str(i+1)+'.dat')
		try:
			fluxes=np.dstack((fluxes,data['Flux']))
			errors=np.dstack((errors,data['Uncertainty']))
		except:
			fluxes=data['Flux']
			errors=data['Uncertainty']
			wavelength=data['Wavelength']
	wavelengthSI=wavelength*1.e-10
	fluxes=fluxes[0]
	errors=errors[0]
	weight=1./(errors*errors)
	filtered_data=sigma_clip(fluxes,sigma=2.0,iters=None,axis=1)
	clipwaverage=np.ma.average(filtered_data,weights=weight,axis=1)
	uncertainty_on_stack=np.sqrt(1./(np.sum(1./(errors*errors),axis=1)))
	Temperatures=np.arange(5000,100000,500)
	dof=np.size(wavelength)-2
	print 'Degrees of freedom=',dof
	chisqarray=np.array([])
	for T in Temperatures:
		B=1.e-10*2*h*c*c/(wavelengthSI**5)/(np.exp((h*c)/(wavelengthSI*k*T))-1) # in W/m^2/A
		select=np.where((wavelength>5000)&(wavelength<6000))
		renorm=np.average(B[select])/np.average(clipwaverage[select])
		chisq=np.sum((B/renorm-clipwaverage)**2./uncertainty_on_stack**2.)
		redchisq=chisq/dof
		chisqarray=np.append(chisqarray,redchisq)
	minchisq=np.min(chisqarray)
	print 'Minimum reduced $\chi^2$ and resulting probability: ', minchisq,scipy.stats.chi2.pdf(minchisq*dof,dof)
	Tbest=Temperatures[np.where(chisqarray==minchisq)][0]
	fig,ax=plt.subplots()
	plt.plot(Temperatures,chisqarray)
	plt.axis([30000,60000,100,400])
	plt.xlabel('Temperature (K)')
	plt.ylabel(r'$\chi^2/\nu$')
	fig.show()
	fig2,ax2=plt.subplots()
	plt.loglog(wavelength,wavelength*clipwaverage,label='Observed data')
	B=1.e-10*2*h*c*c/(wavelengthSI**5)/(np.exp((h*c)/(wavelengthSI*k*Tbest))-1)
	renorm=np.average(B[select])/np.average(clipwaverage[select])
	plt.loglog(wavelength,wavelength*B/renorm,label='Blackbody T='+str(Tbest))
	plt.xlabel('Wavelength (A)')
	plt.ylabel('$\lambda F (10^{13} W/m^2$)')
	plt.legend()
	fig2.show()


def q2():
	Gstars=ascii.read('phy270_ass5_stars1.txt')
	Fstars=ascii.read('phy270_ass5_stars2.txt')
	fig,ax=plt.subplots(3,2)
	xtest=Fstars['B']
	ytest=Gstars['B']
	testdistribution(xtest,ytest,ax[0][0],xlabel='B',ylabel='N')
	xtest=Fstars['V']
	ytest=Gstars['V']
	testdistribution(xtest,ytest,ax[1][0],xlabel='V',ylabel='N')
	xtest=Fstars['B']-Fstars['V']
	ytest=Gstars['B']-Gstars['V']
	testdistribution(xtest,ytest,ax[2][0],xlabel='(B-V)',ylabel='N')
	xtest=np.concatenate((Fstars['B'],Gstars['B']))
	ytest=np.concatenate((Fstars['V'],Gstars['V']))
	testcorrelation(xtest,ytest,ax[0,1],xlabel='B',ylabel='V')
	testcorrelation(xtest,xtest-ytest,ax[1,1],xlabel='B',ylabel='(B-V)') 
	ax[2,1].axis('off')
	thismanager=plt.get_current_fig_manager()
	thismanager.resize(1200, 800)
	fig.show()

def testcorrelation(x,y,ax,xlabel,ylabel):
	ax.scatter(x,y,c='k')
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	rho,sp=scipy.stats.spearmanr(x,y)
	print 'Spearman\'s rank and probability: ',rho,sp


def testdistribution(x,y,ax,xlabel,ylabel):
	min=np.min((np.min(x),np.min(y)))
	max=np.max((np.max(x),np.max(y)))
	bsize=(max-min)/10
	bins=np.arange(min,max,bsize)
	Nx,xbins=np.histogram(x,bins)
	Ny,ybins=np.histogram(y,bins)
	bmid=0.5*(xbins[1:]+xbins[:-1])
	ax.scatter(bmid,Nx,c='k',s=100)
	ax.scatter(bmid,Ny,edgecolors='r',facecolors='none',s=100)
	ax.errorbar(bmid,Nx,yerr=np.sqrt(Nx),fmt=',',c='k')
	ax.errorbar(bmid,Ny,yerr=np.sqrt(Ny),fmt=',',c='r')
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	t,tp=scipy.stats.ttest_ind(x,y)
	D,kp=scipy.stats.ks_2samp(x,y)
	print 'Student\'s T-test, probability: ',t,tp
	print 'KS test, probability: ',D,kp

def q3a(sigma=0.3):
	imf=ascii.read('phy270_ass5_imf.txt')['Mass']
	imf_obs=imf+np.random.normal(scale=sigma,size=np.size(imf))
	bins=np.arange(0.1,20,.1)
	(ytrue,x)=np.histogram(imf,bins)
	(yobs,x)=np.histogram(imf_obs,bins)
	xmid=0.5*(x[1:]+x[:-1])
	fig,ax=plt.subplots()
	plt.loglog(xmid,ytrue,label='True')
	plt.loglog(xmid,yobs,label='Observed')        
	ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax.xaxis.get_major_formatter().set_scientific(False)
	ax.xaxis.get_major_formatter().set_useOffset(False)
	plt.legend()
	ax.axis([0.1,5,1,20000])
	plt.xlabel('$M/M_\odot$')
	plt.ylabel('Number of stars')
	fig.show()

def q3c(nbins=10):
	imf=ascii.read('phy270_ass5_imf.txt')['Mass']
	alpha_obs=np.array([])
	sigma_array=np.arange(0.1,0.5,.01)
	for sigma in sigma_array:
		imf_obs=imf+np.random.normal(scale=sigma,size=np.size(imf))
		slope,intercept=fitimf(imf_obs,1.,3.,Nbins=nbins)
		alpha_obs=np.append(alpha_obs,slope)
	fig,ax=plt.subplots()
	plt.plot(sigma_array,alpha_obs)
	plt.xlabel(r'Uncertainty $\sigma$')
	plt.ylabel(r'Slope $\alpha$')
	fig.show()

def fitimf(imf,m1,m2,Nbins=10.):
	bins=np.arange(m1,m2,(m2-m1)/Nbins)
	y,x=np.histogram(imf,bins)
	xmid=0.5*(x[1:]+x[:-1])
	slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(np.log10(xmid),np.log10(y))
	return (slope,intercept)



