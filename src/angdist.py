import numpy as np
from cosmolopy import fidcosmo, magnitudes, luminosityfunction,cc, cd
import numpy, math
import sys
from astropy import units as u
from astropy import constants as c

cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}

def angular_diameter_distance(z1,z2):
  ''' Return the angular diameter distance between z1 and z2 in Mpc '''
  dds = cd.angular_diameter_distance(z2, z0=z1, **cosmo)
  return dds

def inf_to_Dds_ds1(z_lens):
  ''' Return the factor to scale a map at z=inf to Dds_Ds=1'''
  dm2= cd.comoving_distance_transverse(np.inf, **cosmo)
  dm1= cd.comoving_distance_transverse(z_lens, **cosmo)
  dh = cd.hubble_distance_z(1.0, **cosmo)*cd.e_z(1.0, **cosmo)
  ddsminf = (dm2*np.sqrt(1+cosmo['omega_k_0']*np.square(dm1/dh)) - dm1*np.sqrt(1+cosmo['omega_k_0']*np.square(dm2/dh)))
  return dm2/ddsminf

def redshift_weight(z_lens,z_source):
  ''' get the weight you need to apply to a z=inf map to get the map at z=z_source for a lens at z=z_lens '''
  dd = cd.angular_diameter_distance(z_lens, **cosmo)*u.parsec*1E6 # angular diameter distance between observer and lens
  ds = cd.angular_diameter_distance(z_source, **cosmo)*u.parsec*1E6 # angular diameter distance between observer and source
  dds = cd.angular_diameter_distance(z_source, z0=z_lens, **cosmo)[0]*u.parsec*1E6
  dm2= cd.comoving_distance_transverse(np.inf, **cosmo)
  dm1= cd.comoving_distance_transverse(z_lens, **cosmo)
  dh = cd.hubble_distance_z(1.0, **cosmo)*cd.e_z(1.0, **cosmo)
  ddsminf = (dm2*np.sqrt(1+cosmo['omega_k_0']*np.square(dm1/dh)) - dm1*np.sqrt(1+cosmo['omega_k_0']*np.square(dm2/dh)))
  Z = (dm2/ddsminf) / (ds/dds)
  return Z.value

def Dds_ds(z_lens,z_source):
  ''' Get the factor that you multiply a map scaled to Dds_Ds=1 by to get a map at z=z_source'''
  ds = cd.angular_diameter_distance(z_source, **cosmo)*u.parsec*1E6 # angular diameter distance between observer and source
  dds = cd.angular_diameter_distance(z_source, z0=z_lens, **cosmo)[0]*u.parsec*1E6
  return dds/ds
