#=======================================================================
#+
# NAME:
#   calc_magnification_glass.py
#
# PURPOSE:
#   Main script to calculate the magnification of an object from the GLASS lensing maps
#
#	
# Known BUGS:
# 
#
# Needed Fixes:
# 		
# 		
# REVISION HISTORY:
#   2016-05-30  started by Austin Hoag (UC Davis)
#	2018-02-20  fixed bug where single line did not work
#-
#=======================================================================

# default python libraries
import numpy as np

# Additionally-installed python libraries
from astropy.wcs import WCS
from astropy.io import fits 

# Homemade python libraries
from angdist import redshift_weight, inf_to_Dds_ds1, Dds_ds

parent_dir = '/data/external/GLASS/GLASS21/lensmodels' # main directory where lens maps live
cluster_z_dict = {'A2744':0.308,'MACS0416':0.396,'MACS0717':0.545,'MACS1149':0.543,'RXJ2248':0.348,'A370':0.375,'MACS2129':0.59,'MACS1423':0.545,'RXJ1347':0.451,'MACS0744':0.686}
HFF_name_dict = {'A2744':'abell2744','MACS0416':'macs0416','MACS0717':'macs0717','MACS1149':'macs1149','RXJ2248':'abells1063','A370':'abell370','MACS2129':'macs2129','MACS1423':'macs1423','RXJ1347':'rxj1347','MACS0744':'macs0744'}
version_dict = {'A2744':'v2','MACS0416':'v3','MACS0717':'v1','MACS1149':'v1','RXJ2248':'v1','A370':'v1','MACS2129':'v1','MACS1423':'v2','RXJ1347':'v2','MACS0744':'v1'}

def get_magnification_list_singlefits_vectorized(cluster,ras,decs,zs,errors=False):
	"""
	-----PURPOSE-----
	Look up the magnification from a list of objects in a single cluster
	using the single fits file for errors and the vectorized approach 
	-----INPUT-------
	cluster 		shortname for cluster, e.g. 'A2744'
	ras 			list of right ascensions in degrees (e.g. [3.5713877])
	decs 			list of declinations in degrees (e.g. [-30.3770812])
	zs 				List of source redshifts (e.g. [1.05])
	errors 			If True, use the additional reconstructions in range/ to determine confidence interval
	-----OUTPUT------
	magstrs 		list of strings of the format: 'mag_best/mag_median [mag 68% conf. limit]'' 
					or '-99' if object not in lensing field of view
	"""
	HFF_name = HFF_name_dict[cluster]
	scale_dir = parent_dir + '/' + HFF_name
	ffconv_pref = 'hlsp_glass_model_%s_bradac' % HFF_name
	version = version_dict[cluster]
	source_pref = scale_dir + '/' + ffconv_pref + '_' + version + '_'
	kappa_Dds_Ds1 = source_pref + 'kappa.fits'
	gamma1_Dds_Ds1 = source_pref + 'gamma1.fits'
	gamma2_Dds_Ds1 = source_pref + 'gamma2.fits'
	W_best = WCS(kappa_Dds_Ds1) 
	''' First, see if ra and dec are in the lensing field-of-view '''
	kappa_Dds_Ds1_array = fits.getdata(kappa_Dds_Ds1)
	xmax,ymax = np.shape(kappa_Dds_Ds1_array)
	gamma1_Dds_Ds1_array = fits.getdata(gamma1_Dds_Ds1)
	gamma2_Dds_Ds1_array = fits.getdata(gamma2_Dds_Ds1)
	if errors: 
		bootdir =  scale_dir + '/range'
		bootfile = '%s/all_range_%s_%s.fits' % (bootdir,HFF_name,version)
		bootdata = fits.getdata(bootfile)

	magstrs = []
	for ii in range(len(ras)):
		ra = float(ras[ii])
		dec = float(decs[ii])
		z = float(zs[ii])
		x,y = W_best.wcs_world2pix(ra,dec,1)
		# print x,y
		if x < 0 or y < 0 or x> xmax or y> ymax:
			print "Point is outside lensing field of view"
			magstr = '-99'
			magstrs.append(magstr)
			continue
		else:
			pass
		
		kappa_Dds_Ds1 = bilinear_interpolate(kappa_Dds_Ds1_array,x,y)
		gamma1_Dds_Ds1 = bilinear_interpolate(gamma1_Dds_Ds1_array,x,y)
		gamma2_Dds_Ds1 = bilinear_interpolate(gamma2_Dds_Ds1_array,x,y)
		# print kappa_Dds_Ds1, gamma1_Dds_Ds1, gamma2_Dds_Ds1
		z_lens = cluster_z_dict[cluster]
		# zweight = lm.Dds_Ds(zl=z_lens,zs=z) # what you mulitply Dds_Ds=1 map by to get z=z_source map
		zweight = Dds_ds(z_lens=z_lens,z_source=z)
		kappa_source = kappa_Dds_Ds1*zweight
		gamma1_source = gamma1_Dds_Ds1*zweight
		gamma2_source = gamma2_Dds_Ds1*zweight
		det_source = (1-kappa_source)**2 - gamma1_source**2 - gamma2_source**2
		mag_source = np.abs(1/det_source)

		if errors: 
			kappa_Dds_Ds1_boot_array = np.array([bilinear_interpolate(bootdata[jj],x,y) for jj in range(100)])
			gamma_Dds_Ds1_boot_array = np.array([bilinear_interpolate(bootdata[jj+100],x,y) for jj in range(100)])
			kappavals_boot = kappa_Dds_Ds1_boot_array*zweight 
			gammavals_boot = gamma_Dds_Ds1_boot_array*zweight
			magvals_boot = 1/np.abs((1-kappavals_boot)**2-(gammavals_boot**2))
			med,low,high = confidence_interval(magvals_boot)
			magstr = '%.2f/%.2f [%.2f-%.2f]' % (mag_source,med,low,high)
		else:
			magstr = '%.2f' % mag_source
		magstrs.append(magstr)
	return magstrs

def get_magnification_list_smart(infile):
	"""
	-----PURPOSE-----
	Return the magnifications for the objects in an input catalog
	by looking up magnifications of objects in the same cluster simultaneously
	each time. This is the method I plan to use in the actual webtool.
	-----INPUT-------
	infile 		The name of the comma-separated input file formatted cluster,id,ra,dec,z
	-----OUTPUT------
	cluster_mag_dict 		A dictionary where the keys are cluster names 
							and the values are dictionaries where the keys are object IDs
							and the values are magnification strings 
	"""
	all_clusters,ids,ras,decs,zs = np.genfromtxt(infile,unpack=True,usecols=(0,1,2,3,4),dtype='S20',delimiter=',')
	if type(ras) == np.string_: # if there is only one line, then genfromtxt sets type to string, not ndarray
		all_clusters = np.array([all_clusters])
		ids = np.array([ids])
		ras = np.array([float(ras)])
		decs = np.array([float(decs)])
		zs = np.array([float(zs)])
	else:
		ras = map(float,ras)
		decs = map(float,decs)
		zs = map(float,zs)

	''' If one cluster appears more than once, only open the fits file once '''
	unique_cluster_list = list(set(all_clusters))
	cluster_id_dict = dict((m,[]) for m in unique_cluster_list)
	cluster_mag_dict = dict((m,{}) for m in unique_cluster_list)
	for ii in range(len(all_clusters)):
		cluster = all_clusters[ii]
		ID = ids[ii]
		ra = ras[ii]
		dec = decs[ii]
		z = zs[ii]
		cluster_id_dict[cluster].append([ID,ra,dec,z])	
	for cluster in cluster_id_dict.keys():	
		idlist = np.array(cluster_id_dict[cluster])[:,0]
		ralist = np.array(cluster_id_dict[cluster])[:,1]
		declist = np.array(cluster_id_dict[cluster])[:,2]
		zlist = np.array(cluster_id_dict[cluster])[:,3]

		if cluster == 'MACS0744': # don't have lens models for these clusters yet
			magstrs = ['-99' for x in range(len(idlist))]
		else:
			magstrs = get_magnification_list_singlefits_vectorized(cluster,ralist,declist,zlist,errors=True)
		for jj in range(len(idlist)):
			ID = idlist[jj]
			magstr = magstrs[jj]
			cluster_mag_dict[cluster][ID] = magstr
	''' Return the magnifications of objects in the same order as appeared in the input catalog '''
	outstr = '# Cluster,ID,RA,DEC,z,mu (best/med [lower68%,upper68%])<br/> '
	for kk in range(len(all_clusters)):
		cluster = all_clusters[kk]
		ID = ids[kk]
		ra = ras[kk]
		dec = decs[kk]
		z = zs[kk]
		mag = cluster_mag_dict[cluster][ID]
		outstr += '%s,%s,%f,%f,%f,%s<br/>' % (cluster,ID,ra,dec,z,mag) 
	return outstr

def confidence_interval(data,interval=68.3):
	"""
	-----PURPOSE-----
	Calculate the confidence interval of a 1D data array 
	-----INPUT-------
	im  	A 2D numpy array
	x 		x coordinate at which to calculate the value of im
	y 		y coordinate at which to calculate the value of im
	"""
	halfint = interval/2.
	percent_low = 50-halfint
	percent_high = 50+halfint
	lower=np.percentile(data,percent_low)
	upper=np.percentile(data,percent_high)
	med=np.percentile(data,50)
	return med,lower,upper

def bilinear_interpolate(im, x, y):
	"""
	-----PURPOSE-----
	Calculate the value of im at an exact x,y 
	not necessarily on the grid nodes of im by
	bilinear interpolation 
	-----INPUT-------
	im  	A 2D numpy array
	x 		x coordinate at which to calculate the value of im
	y 		y coordinate at which to calculate the value of im
	-----OUTPUT-------
	The value of im evaluated at the input x,y
	"""
	x = np.asarray(x)
	y = np.asarray(y)

	x0 = np.floor(x).astype(int)
	x1 = x0 + 1
	y0 = np.floor(y).astype(int)
	y1 = y0 + 1

	x0 = np.clip(x0, 0, im.shape[1]-1);
	x1 = np.clip(x1, 0, im.shape[1]-1);
	y0 = np.clip(y0, 0, im.shape[0]-1);
	y1 = np.clip(y1, 0, im.shape[0]-1);

	Ia = im[ y0, x0 ]
	Ib = im[ y1, x0 ]
	Ic = im[ y0, x1 ]
	Id = im[ y1, x1 ]

	wa = (x1-x) * (y1-y)
	wb = (x1-x) * (y-y0)
	wc = (x-x0) * (y1-y)
	wd = (x-x0) * (y-y0)

	return wa*Ia + wb*Ib + wc*Ic + wd*Id
