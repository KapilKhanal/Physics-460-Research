# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 20:16:50 2015
@author: carl

collections of scripts that I use often for astronomy relate things
"""
import numpy as np
from astropy import units as u
from astropy.cosmology import WMAP9
from astropy.cosmology import FlatLambdaCDM
from astropy import constants
from astropy.cosmology import Planck13


#Define a cosmology
cosmoWMAP=WMAP9
cosmoFLAT=FlatLambdaCDM(H0=70,Om0=0.3)
cosmoPLANCK=Planck13



def luminosity(flux,z):
    'enter the flux in W m^-2 and redshift of the sorce to return the luminosity in L_solar'
    lum_dist = cosmoWMAP.luminosity_distance(z)
    area     = (4*np.pi*(lum_dist)**2.).to(u.m**2)
    luminosity = (flux*area).to(u.solLum)
    return luminosity


   
def luminoisty_prime(flux,frequency,z):
    'enter the flux in jy km/s and the rest wavelegnth of the line, and will return the'
    'line linumonsity in K pc^2 km/s'
    lum_dist = cosmoWMAP.luminosity_distance(z)
    sky_frequency = frequency/(z+1)
    luminosity_prime = 3.25e7*flux*(lum_dist.value)**2/((1+z)**3*sky_frequency**2)
    return luminosity_prime*u.Unit("K pc^2 km/(s)")


def Lnu_to_Fnu(Lnu,z):
    '''Given the speicifc luminosity in erg s^-1 Hz^-1, returns the flux density
    in erg cm^-2 s^-1 Hz^-1 and mJy. Do not enter values with astropy units'''
    Dl = cosmoWMAP.luminosity_distance(z).to('cm')
    Fnu = Lnu/(4*np.pi*Dl.value**2) #erg s^-1 cm^-2 Hz^-1
    Fnu_mJy = Fnu/10**-23*1000.
    return Fnu, Fnu_mJy

def Lsol_to_cgs(Lsol):
    '''Given a luminosity in solor luminosity units it returns the luminosity
    in cgs, erg s^-1. Do not enter values with astropy units'''
    L_cgs = (Lsol*u.Lsun).to('erg/(s)').value 
    return L_cgs

def Lum(flux,z):
    'enter the flux in W m^-2 and redshift of the sorce to return the luminosity in L_solar'
    lum_dist = cosmoWMAP.luminosity_distance(z)
    area     = (4*np.pi*(lum_dist)**2.).to(u.m**2)
    luminosity = (flux*area).to(u.solLum)
    return luminosity

def Lum_from_LumPrime(l_prime,line_freq):
    'Enter the line luminosity in K pc^2 km/s and rest frequency of line in GHz will return the line luminosity in L_solar'
    Lum = 3.0e-11*(1/(u.K*u.km/u.s*u.pc**2*u.GHz**3))*line_freq**3*l_prime*u.Lsun
    return Lum

def LumPrime(line_flux,z,freq_obs):
    'given the line flux in Jy km/s, the redshift and the frequency of observation in'
    'GHz, will return the luminosity in K km/s pc^2'
    lum_dist = cosmoWMAP.luminosity_distance(z)
    LumPrime = 3.25e3*u.K/u.Jy*(100*u.GHz)**2*(u.pc/u.Mpc)**2*line_flux*lum_dist**2/((1+z)**3*(freq_obs)**2)
    return LumPrime.to(u.K*u.km/u.s*u.pc**2)

def jy_kms_watts_m(flux,frequency,z):
    'enter the line flux in Jy km/s, the line rest frequency in GHz, and redshift and will return the line flux in W m^-2'
    watts_m = (flux*frequency/(z+1)/constants.c).to(u.W/u.m**2)
    return watts_m


def watts_m_jy_kms(flux,frequency,z):
    'enter the line flux in W m^-2, the rest frequency of the line in GHz, and redshift and will return the line flux in Jy km/s'
    jy_kms = (flux/frequency*(z+1)*constants.c).to(u.Jy*u.km/u.s)
    return jy_kms

def cont_luminosity_density(flux_density,z):
    'enter the flux density in mJy and redshift of the sorce to return the luminosity density in L_solar/Hz'
    lum_dist = cosmoWMAP.luminosity_distance(z)
    area     = (4*np.pi*(lum_dist)**2.).to(u.m**2)
    luminosity = (flux_density*area).to(u.solLum/u.Hz)
    return luminosity

def cont_luminosity(flux_density,z,v_obs):
    'enter the flux density in mJy and redshift and observed frequency of the sorce to return'
    'the luminosity in L_solar at the rest frequency of freq*(1+z)'
    lum_dist = cosmoWMAP.luminosity_distance(z)
    area     = (4*np.pi*(lum_dist)**2.).to(u.m**2)
    luminosity = (flux_density*area*v_obs).to(u.solLum)
    return luminosity

    #Flux ratios in Jansky Units
#	Bournaud et al. 2014 Model	
#	High-z disk	SB Merger
#CO(2-1)/CO(1-0)	2.5	3.5
#CO(3-2)/CO(1-0)	3	6.25
#CO(4-3)/CO(1-0)	3.75	8.75
#CO(5-4)/CO(1-0)	3.5	9.5
#CO(6-5)/CO(1-0)	2.75	8
#CO(7-6)/CO(1-0)	1.75	5.25
    


ALMA_Bands =np.array([('Band 3', 84,    116,    'GHz',  3.6,    2.6 ,  'mm'),
                     ('Band 4', 125,    163,    'GHz',  2.4,    1.8 ,  'mm'),
                     ('Band 6', 211,    275,    'GHz',  1.4,    1.1 ,  'mm'),
                     ('Band 7', 275,    373,    'GHz',  1.1,    0.8 ,  'mm'),
                     ('Band 8', 385,    500,    'GHz',  0.78,   0.6 ,  'mm'),
                     ('Band 9', 602,    720,    'GHz',  0.5,    0.42,  'mm'),
                     ('Band10', 787,    950,    'GHz',  0.38,   0.32 ,  'mm')],
                    dtype={'names'  :['Band','loF',  'upF','freqUnit', 'loW', 'upW','waveUnit'],
                           'formats':['|S10','float64','float64','|S10','float64','float64','|S10']})

FSL_lines  = np.array([('[CI] 610', 492.16065,  'GHz',  609.13537,  'um'),
                       ('[CI] 370', 809.34197,  'GHz',  370.41506,  'um'),
                       ('[NII] 205',1461.13141, 'GHz',  205.17830,  'um'),
                       ('[CII] 158',1900.53690, 'GHz',  157.74093,  'um'),
                       ('[OI] 145', 2060.06886, 'GHz',  145.52545,  'um'),
                       ('[NII] 122',2459.38010, 'GHz',  121.89757,  'um'),
                       ('[OIII] 88',3393.00624, 'GHz',  88.35600,   'um'),
                       ('[OI] 63'  ,4744.77,    'GHz',  63.18,      'um'),
                       ('[NIII] 57',5230.43,    'GHz',  57.32,      'um'),
                       ('[OIII] 52',5786.89659, 'GHz',  51.81450,   'um')],
                    dtype={'names'  :['Line','freq', 'freqUnit', 'wave','waveUnit'],
                           'formats':['|S10','float64','|S10','float64','|S10']})  

def CO_line_frequency(co_name):
    CODict = np.array([(115.2712,230.5380,345.79599,461.04077,576.26793,691.47308,806.65180)],
                    dtype={'names'  :['CO(1-0)','CO(2-1)','CO(3-2)','CO(4-3)','CO(5-4)','CO(6-5)','CO(7-6)'],
                               'formats':['float64','float64','float64','float64' ,'float64','float64','float64' ]})
    return CODict[co_name]

def FSL_line_frequency(line_name):
    '''Given a line name returns the frequency of the line in GHz
    Possible line names include [CI] 610, [CI] 370,[NII] 205, [CII] 158, [OI] 145
    [NII] 122, [OIII] 88, [OI] 63, [OIII] 52'''
    
    FSL_lines  = np.array([('[CI] 610', 492.16065,  'GHz',  609.13537,  'um'),
                           ('[CI] 370', 809.34197,  'GHz',  370.41506,  'um'),
                           ('[NII] 205',1461.13141, 'GHz',  205.17830,  'um'),
                           ('[CII] 158',1900.53690, 'GHz',  157.74093,  'um'),
                           ('[OI] 145', 2060.06886, 'GHz',  145.52545,  'um'),
                           ('[NII] 122',2459.38010, 'GHz',  121.89757,  'um'),
                           ('[OIII] 88',3393.00624, 'GHz',  88.35600,   'um'),
                           ('[OI] 63'  ,4744.77,    'GHz',  63.18,      'um'),
                           ('[NIII] 57',5230.43,    'GHz',  57.32,      'um'),
                           ('[OIII] 52',5786.89659, 'GHz',  51.81450,   'um')],
                        dtype={'names'  :['Line','freq', 'freqUnit', 'wave','waveUnit'],
                               'formats':['|S10','float64','|S10','float64','|S10']})   
    line = np.where(FSL_lines['Line']==line_name.encode())[0][0]
    return FSL_lines['freq'][line]


def redshifted_wavelength(z,line_name):
    line = np.where(FSL_lines['Line']==line_name.encode())[0][0]
    return (z+1)*FSL_lines['wave'][line]*u.Unit(FSL_lines['waveUnit'][line])

def redshifted_frequency(z,line_name):
    line = np.where(FSL_lines['Line']==line_name.encode())[0][0]
    return FSL_lines['freq'][line]*u.Unit(FSL_lines['freqUnit'][line])/(z+1)   

def linear_to_log_errors(data_table,linear_errors):
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