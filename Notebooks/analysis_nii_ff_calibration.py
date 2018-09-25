# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import numpy.ma as ma #For Masking the wrong or null entry and avoiding them in array calculations later on
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table,join #http://docs.astropy.org/en/stable/table/

import matplotlib.patches as mpatches
from astropy.cosmology import WMAP9
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import matplotlib.patheffects as PathEffects
from astropy import constants
from matplotlib.patches import Rectangle
from matplotlib import colors

#Define a cosmology
cosmoWMAP=WMAP9
cosmoFLAT=FlatLambdaCDM(H0=70,Om0=0.3) #http://docs.astropy.org/en/v0.3/api/astropy.cosmology.core.FlatLambdaCDM.html


oiii_nii_ff_metal_location='/Users/qe5106hi/Dropbox/ASTROPHYSICS/oiii_nii_ff_metallicity.csv'
#L_FIR is lensing corrected for lensed galaxies
#Setting up the appropriate Unites for the appropriate data columns
data_table = ascii.read(oiii_nii_ff_metal_location,data_start=2) #Return Table Probably starting a column2?
for column in data_table.columns[2,3,5,6,8,9,11,12,14,15,17,18,20,21]:
    data_table['%s' %column].unit = '1e-18 W m^-2'
for column in data_table.columns[24,25]:
    data_table['%s' %column].unit='12+log(O/H)'
for column in data_table.columns[27,28,32,33]:
    data_table['%s' %column].unit='mJy'
for column in data_table.columns[30,34]:
    data_table['%s' %column].unit='GHz'
data_table.columns[30,34].unit = 'GHz'
data_table.columns[37].unit = 'Mpc'
#
#data_table.columns[3].mask[np.where(ferkinhoff2015.columns[3]==-1)[0]]=True #mask galaxies with unknown magnification

#-----------------------------------------------------------------------------#
#---------------------Define some useful function-----------------------------#
#-----------------------------------------------------------------------------#

def make_masked_quantities(data_table,convert_unit=False,desired_unit=None):
    #print unit
    if convert_unit == False:
        values = (data_table.data.data.astype(np.float)*data_table.unit)
    else:
        values = (data_table.data.data.astype(np.float)*data_table.unit).to(desired_unit)
    mask = data_table.mask
    return ma.array(values,mask=mask)

def linear_to_log_errors(data_table,linear_errors): #Function to return the error
    upper_bound = data_table+linear_errors
    lower_bound = data_table-linear_errors
    log_data = np.log10(data_table)
    log_upper_bound = np.log10(upper_bound)
    log_lower_bound = np.log10(lower_bound)
    positive_error = log_upper_bound-log_data
    negative_error = log_data-log_lower_bound
    return positive_error,negative_error

def linear_to_log_errors_log_data(log_data_table,linear_errors):
    upper_bound = 10**log_data_table+linear_errors
    lower_bound = 10**log_data_table-linear_errors
    log_data = log_data_table
    log_upper_bound = np.log10(upper_bound)
    log_lower_bound = np.log10(lower_bound)
    positive_error = log_upper_bound-log_data
    negative_error = log_data-log_lower_bound
    return positive_error,negative_error
#-----------------------------------------------------------------------------#
#---------------------Define plotting variables ------------------------------#
#-----------------------------------------------------------------------------#

#cii = make_masked_quantities(data_table['CII'],True,'W m^-2')
cii_data = make_masked_quantities(data_table['CII'],True,'W m^-2')
ERRcii = make_masked_quantities(data_table['errCII'],True,'W m^-2')
ULcii = np.zeros(len(cii_data))
ULcii[np.where(data_table['ulCII']==1)[0]]=True

#nii = make_masked_quantities(data_table['N122'],True,'W m^-2')
nii_data = make_masked_quantities(data_table['N122'],True,'W m^-2')
ERRnii = make_masked_quantities(data_table['errN122'],True,'W m^-2')
ULnii = np.zeros(len(nii_data))
ULnii[np.where(data_table['ulN122']==1)[0]]=True

#oiii= make_masked_quantities(data_table['OIII88'],True,'W m^-2')
oiii_data= make_masked_quantities(data_table['OIII88'],True,'W m^-2')
ERRoiii = make_masked_quantities(data_table['errOIII88'],True,'W m^-2')
ULoiii = np.zeros(len(oiii_data))
ULoiii[np.where(data_table['ulOIII88']==1)[0]]=True #?????????

ff  = make_masked_quantities(data_table['FreeFree'])
ERRff = make_masked_quantities(data_table['errFreeFree'])
#For Free Free values without an uncertainty I assume an error 30% (i.e they are ~3 sigma detections)
ERRff[np.where(data_table['errFreeFree'].data.data==0.)[0]]=0.3*ff[np.where(data_table['errFreeFree'].data.data==0.)[0]].data
ULff = np.zeros(len(ff))
ULff[np.where(data_table['ulFreeFree']==1)[0]]=True

metallicity = ma.array(np.copy(data_table['Metallicity'].data.data),mask=data_table['Metallicity'].data.mask)
ERRmetallicity = data_table['errMetallicity'].data.data
#For the metallicity values without a given uncertainty I give a 0.1 dex error
ERRmetallicity[np.where(data_table['errMetallicity'].data.data<0.1)[0]]=0.1

cloudy_o_solar_abundance = 8.69
cloudy_n_solar_abundance = 7.93
#abundance = metallicity - cloudy_o_solar_abundance
abundance_data = metallicity - cloudy_o_solar_abundance


#-----------------------------------------------------------------------------#
#---------------------Plot Data  -- Line, FF, Metalicity----------------------#
#-----------------------------------------------------------------------------#

fig = plt.figure(1,figsize=(4.5,7.5),dpi=150)
fig.subplots_adjust(hspace=0,wspace=0)

#----------------NII/FF versus Metalicity ------------------------------------#
ax1  = fig.add_subplot(3,1,1)
plt.setp(ax1.get_xticklabels(),visible=False)
#plt.xlim(5e7,1e12)
plt.ylim(10,13.5)
plt.yticks([10.0,10.5,11,11.5,12.0,12.5,13.0])

ydata = ma.array(np.log10((nii/ff).data.value*10.0**29),mask=(nii/ff).mask) #Nii_data?
xdata = metallicity
zdata = ma.array(np.log10((oiii/nii).data.value),mask=(oiii/nii).mask) #???????

#cmap = plt.cm.jet
#norm = colors.Normalize(zdata.min(),zdata.max())
#ec_colors = plt.cm.jet(norm(zdata))


yerr_linear = np.sqrt((ERRnii/nii).data.value**2+(ERRff/ff).data.value**2)*(nii/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULnii

xerr = ERRmetallicity

plt.scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
plt.errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)


plt.ylabel(r'$[NII]/S_{\nu}^{ff}$')
plt.grid(True,alpha=0.3)

#----------------OIII/FF versus Metalicity ------------------------------------#

ax2 = fig.add_subplot(3,1,2)
plt.setp(ax2.get_xticklabels(),visible=False)
plt.ylim(10,13.5)
plt.yticks([10.0,10.5,11,11.5,12.0,12.5,13.0])


ydata = ma.array(np.log10((oiii/ff).data.value*10.0**29),mask=(nii/ff).mask)
xdata = metallicity
zdata = ma.array(np.log10((oiii/nii).data.value),mask=(oiii/nii).mask)

yerr_linear = np.sqrt((ERRoiii/oiii).data.value**2+(ERRff/ff).data.value**2)*(oiii/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULoiii

xerr = ERRmetallicity

plt.scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
plt.errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)
             
plt.ylabel(r'$[OIII]/S_{\nu}^{ff}$')
plt.grid(True,alpha=0.3)

#----------------CII/FF versus Metalicity ------------------------------------#

ax3 = fig.add_subplot(3,1,3)
plt.ylim(10,13.5)
plt.yticks([10.0,10.5,11,11.5,12.0,12.5,13.0])

ydata = ma.array(np.log10((cii/ff).data.value*10.0**29),mask=(nii/ff).mask)
xdata = metallicity
zdata = ma.array(np.log10((oiii/nii).data.value),mask=(oiii/nii).mask)

yerr_linear = np.sqrt((ERRcii/cii).data.value**2+(ERRff/ff).data.value**2)*(cii/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULcii

xerr = ERRmetallicity

plt.scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
plt.errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)
             
plt.ylabel(r'$[CII]/S_{\nu}^{ff}$')
plt.xlabel(r'$12+log[O/H]$')
plt.grid(True,alpha=0.3)

#cb = plt.colorbar(ax=[ax1,ax2,ax3],orientation="horizontal",fraction=0.04,pad=0.07,
#                  anchor)
CBposition=fig.add_axes([0.35,0.16,0.50,0.02])
cb = plt.colorbar(cax=CBposition,ticks=[0.25,.75,1.25,1.75],orientation="horizontal")
cbax = cb.ax
cbax.text(-0.02,0.4,r'$[OIII]/[NII]$',horizontalalignment='right')
#cb.set_label(r'$[OIII]/[NII]$')

plt.savefig('Proposals/ff_fir_metalicity/nii_oii_cii_ff_OH.png',dpi=300,bbox_inches='tight',pad_inches=0.5)
plt.savefig('Proposals/ff_fir_metalicity/nii_oii_cii_ff_OH.pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)


#-----------------------------------------------------------------------------#
#---------------------Plot Data - Square - Line, FF, Metalicity_--------------#
#-----------------------------------------------------------------------------#

fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(6,8.75),dpi=150)
fig.subplots_adjust(hspace=.175,wspace=.375)
#plt.setp([a.get_xticklabels() for a in axes[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axes[:, 1]], visible=False)

#----------------NII/FF versus Metalicity ------------------------------------#
axes[0,0].set_ylabel(r'$log([NII]/S_{\nu}^{radio})$')
axes[0,0].set_xlabel(r'$12+log[O/H]$')

ydata = ma.array(np.log10((nii_data/ff).data.value*10.0**29),mask=(nii_data/ff).mask)
xdata = metallicity
zdata = ma.array(np.log10((oiii_data/nii_data).data.value),mask=(oiii_data/nii_data).mask)

#cmap = plt.cm.jet
#norm = colors.Normalize(zdata.min(),zdata.max())
#ec_colors = plt.cm.jet(norm(zdata))


yerr_linear = np.sqrt((ERRnii/nii_data).data.value**2+(ERRff/ff).data.value**2)*(nii_data/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yuplims = ULnii
#yerr_log_pos.mask = yuplims
#yerr_log_neg.mask = yuplims
where_UL=np.where(yuplims==1.0)[0]
yerr_log_pos[where_UL]=0

yerr=np.vstack([yerr_log_neg,yerr_log_pos])
#yerr.mask=yuplims
xerr = ERRmetallicity

axes[0,0].scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
axes[0,0].errorbar(xdata,ydata,xerr=xerr,yerr=yerr,lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)


#----------------OIII/FF versus Metalicity ------------------------------------#

axes[0,1].set_ylabel(r'$log([OIII]/S_{\nu}^{radio})$')
axes[0,1].set_xlabel(r'$12+log[O/H]$')

ydata = ma.array(np.log10((oiii_data/ff).data.value*10.0**29),mask=(nii_data/ff).mask)
xdata = metallicity
zdata = ma.array(np.log10((oiii_data/nii_data).data.value),mask=(oiii_data/nii_data).mask)

yerr_linear = np.sqrt((ERRoiii/oiii_data).data.value**2+(ERRff/ff).data.value**2)*(oiii_data/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULoiii

xerr = ERRmetallicity

axes[0,1].scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
axes[0,1].errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)
             


#----------------CII/FF versus Metalicity ------------------------------------#
axes[1,0].set_ylabel(r'$log([CII]/S_{\nu}^{radio})$')
axes[1,0].set_xlabel(r'$12+log[O/H]$')

ydata = ma.array(np.log10((cii_data/ff).data.value*10.0**29),mask=(nii_data/ff).mask)
xdata = metallicity
zdata = ma.array(np.log10((oiii_data/nii_data).data.value),mask=(oiii_data/nii_data).mask)

yerr_linear = np.sqrt((ERRcii/cii_data).data.value**2+(ERRff/ff).data.value**2)*(cii_data/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULcii

xerr = ERRmetallicity

cii_plot=axes[1,0].scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
cii_error=axes[1,0].errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)
             


#cb = plt.colorbar(ax=[ax1,ax2,ax3],orientation="horizontal",fraction=0.04,pad=0.07,
#                  anchor)
CBposition=fig.add_axes([0.15,0.16,0.28,0.0125])
cb = fig.colorbar(cii_plot,cax=CBposition, ticks=[0.25,.75,1.25,1.75,2.25],orientation="horizontal")
cbax = cb.ax
cbax.text(0.0,1.5,r'$[OIII]/[NII]$',horizontalalignment='left')
#cb.set_label(r'$[OIII]/[NII]$')

for ax in axes.flat:
    ax.set_ylim([10.5,13.5])
    ax.set_xlim([7.5,9.5])
    ax.set_yticks([10.5,11,11.5,12.0,12.5,13.0])
    ax.set_xticks([8.0,8.5,9.0,9.5])
    ax.grid(True,alpha=0.3)


plt.savefig('/home/carl/Dropbox/Proposals/ff_fir_metalicity/nii_oii_cii_ff_OH_SQAURE.png',dpi=300,bbox_inches='tight',pad_inches=0.5)
plt.savefig('/home/carl/Dropbox/Proposals/ff_fir_metalicity/nii_oii_cii_ff_OH_SQAURE.pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)



#-----------------------------------------------------------------------------#
#---------------------Plot Data - Square - Line, FF, log(SOLAR)--------------#
#-----------------------------------------------------------------------------#

fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(6,8.75),dpi=150)
fig.subplots_adjust(hspace=.175,wspace=.375)
#plt.setp([a.get_xticklabels() for a in axes[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axes[:, 1]], visible=False)

#----------------NII/FF versus Metalicity ------------------------------------#
axes[0,0].set_ylabel(r'$log([NII]/S_{\nu}^{radio})$')
axes[0,0].set_xlabel(r'$log(Z/Z_{\odot})$')

ydata = ma.array(np.log10((nii_data/ff).data.value*10.0**29),mask=(nii_data/ff).mask)
xdata = abundance_data
zdata = ma.array(np.log10((oiii_data/nii_data).data.value),mask=(oiii_data/nii_data).mask)

#cmap = plt.cm.jet
#norm = colors.Normalize(zdata.min(),zdata.max())
#ec_colors = plt.cm.jet(norm(zdata))


yerr_linear = np.sqrt((ERRnii/nii_data).data.value**2+(ERRff/ff).data.value**2)*(nii_data/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
#yerr_log_pos.mask = ydata.mask
#yerr_log_neg.mask = ydata.mask
yuplims = ULnii
#yerr_log_pos.mask = yuplims
#yerr_log_neg.mask = yuplims
where_UL=np.where(yuplims==1.0)[0]
yerr_log_pos[where_UL]=0
xerr = ERRmetallicity

axes[0,0].scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
axes[0,0].errorbar(xdata,ydata,xerr=xerr,yerr=yerr,lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)


#----------------OIII/FF versus Metalicity ------------------------------------#

axes[0,1].set_ylabel(r'$log([OIII]/S_{\nu}^{radio})$')
axes[0,1].set_xlabel(r'$log(Z/Z_{\odot})$')

ydata = ma.array(np.log10((oiii_data/ff).data.value*10.0**29),mask=(nii_data/ff).mask)
xdata = abundance_data
zdata = ma.array(np.log10((oiii_data/nii_data).data.value),mask=(oiii_data/nii_data).mask)

yerr_linear = np.sqrt((ERRoiii/oiii_data).data.value**2+(ERRff/ff).data.value**2)*(oiii_data/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULoiii

xerr = ERRmetallicity

axes[0,1].scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
axes[0,1].errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)
             


#----------------CII/FF versus Metalicity ------------------------------------#
axes[1,0].set_ylabel(r'$log([CII]/S_{\nu}^{radio})$')
axes[1,0].set_xlabel(r'$log(Z/Z_{\odot})$')

ydata = ma.array(np.log10((cii_data/ff).data.value*10.0**29),mask=(nii_data/ff).mask)
xdata = abundance_data
zdata = ma.array(np.log10((oiii_data/nii_data).data.value),mask=(oiii_data/nii_data).mask)

yerr_linear = np.sqrt((ERRcii/cii_data).data.value**2+(ERRff/ff).data.value**2)*(cii_data/ff).data.value*10.0**29
yerr_log_pos,yerr_log_neg=linear_to_log_errors_log_data(ydata,yerr_linear)
yerr_log_pos.mask = ydata.mask
yerr_log_neg.mask = ydata.mask
yuplims = ULcii

xerr = ERRmetallicity

cii_plot=axes[1,0].scatter(xdata,ydata,s=40,c=zdata,cmap='jet',alpha=0.7,zorder=1.0,linewidth=0.5)
cii_error=axes[1,0].errorbar(xdata,ydata,xerr=xerr,yerr=[yerr_log_neg,yerr_log_pos],lolims=yuplims,
             ecolor='k',elinewidth=1,capthick=1,fmt=None,zorder=-1.0,alpha=0.5)
             


#cb = plt.colorbar(ax=[ax1,ax2,ax3],orientation="horizontal",fraction=0.04,pad=0.07,
#                  anchor)
CBposition=fig.add_axes([0.15,0.16,0.28,0.0125])
cb = fig.colorbar(cii_plot,cax=CBposition, ticks=[0.25,.75,1.25,1.75,2.25],orientation="horizontal")
cbax = cb.ax
cbax.text(0.0,1.5,r'$[OIII]/[NII]$',horizontalalignment='left')
#cb.set_label(r'$[OIII]/[NII]$')

for ax in axes.flat:
    ax.set_ylim([10.5,13.5])
    ax.set_xlim([-1.25,0.75])
    ax.set_yticks([10.5,11,11.5,12.0,12.5,13.0])
    ax.set_xticks([-1.0,-0.5,0.0,0.5])
    ax.grid(True,alpha=0.3)


plt.savefig('/home/carl/Dropbox/Proposals/ff_fir_metalicity/nii_oii_cii_ff_abundance_SQAURE.png',dpi=300,bbox_inches='tight',pad_inches=0.5)
plt.savefig('/home/carl/Dropbox/Proposals/ff_fir_metalicity/nii_oii_cii_ff_abundance_SQAURE.pdf',dpi=300,bbox_inches='tight',pad_inches=0.5)

























def gaussian(value,mean,sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1/2*((value-mean)/sigma)**2)
plt.figure(2)
mean=3
sigma = 0.5
x=np.arange(-1,10,0.1)
ax1=plt.subplot(311)
ax1.plot(x,gaussian(x,mean,sigma))
ax1.set_xlim(0.0,10.0)
ax2=plt.subplot(312)
ax2.semilogx(x,gaussian(x,mean,sigma))
#plt.plot(x,gaussian(np.log10(x),np.log10(mean),sigma))
ax2.set_xlim(0.0,10.0)
ax3=plt.subplot(313)
ax3.plot((np.log10(x)),gaussian(x,(mean),(sigma)))
ax3.set_xlim(0.0,1.0)
ax3=plt.subplot(313)
log_x = np.range(0.0,1.0,)
ax3.plot((np.log10(x)),gaussian(x,(mean),(sigma)))
ax3.set_xlim(0.0,1.0)
