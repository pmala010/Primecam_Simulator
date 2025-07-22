# ----------------------------------------------------------------------------------------
# utility functions for doing pseudo-Cl with actpol maps
# v2: nominal
# v3: fixing "one pixel offset" for projecting maps to full sky for clenshaw_curtis
# v4: minor updates for doing CEA and CAR more conveniently
# ----------------------------------------------------------------------------------------

from __future__ import print_function
import sys
import numpy as np
import json
from astropy.io import fits as pf
from scipy import ndimage
import copy
import os
# import astropy.wcs as wcs
from pixell import enmap
from pixell import fft as enfft

def write_map(m,header,outname):
    d_hdu = pf.PrimaryHDU(m)
    d_hdu.header = header
    d_hdu.writeto(outname,clobber=True)

def zip_folder(folder_path, output_path):
    import zipfile
    """Zip the contents of an entire folder (with that folder included
    in the archive). Empty subfolders will be included in the archive
    as well.
    """
    parent_folder = os.path.dirname(folder_path)
    # Retrieve the paths of the folder contents.
    contents = os.walk(folder_path)
    try:
        zip_file = zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED)
        for root, folders, files in contents:
            # Include all subfolders, including empty ones.
            for folder_name in folders:
                absolute_path = os.path.join(root, folder_name)
                relative_path = absolute_path.replace(parent_folder,'')
                # print "Adding '%s' to archive." % absolute_path
                zip_file.write(absolute_path, relative_path)
            for file_name in files:
                absolute_path = os.path.join(root, file_name)
                relative_path = absolute_path.replace(parent_folder,'')
                # print "Adding '%s' to archive." % absolute_path
                zip_file.write(absolute_path, relative_path)
        print("'%s' created successfully." % output_path)
    except IOError:#, message:
        print(message)
        sys.exit(1)
    except OSError:#, message:
        print(message)
        sys.exit(1)
    except zipfile.BadZipfile:#, message:
        print(message)
        sys.exit(1)
    finally:
        zip_file.close()

def string_to_binary(input):
   import hashlib
   a = hashlib.md5(input)
   b = a.hexdigest()
   as_int = int(b, 16)
   return bin(as_int)[2:]

def fill_Cl_to_zero(cl):
   return np.append([0,0],cl)

def convolve_Cl(cl,bl):
   ret = copy.copy(cl)
   sh = ret.shape
   assert np.max(sh) == len(bl)
   if len(sh) == 2:
      for i in range(np.min(sh)):
         if sh[0]>sh[1]:
            ret[:,i] *= bl**2
         else:
            ret[i] *= bl**2
      return ret

   return cl*bl**2

def beam_select(path,season,field,ar):
   flist = os.listdir(path)
   for f in flist:
      if season in f and ar in f:
         return os.path.join(path,f)

def beam_select_new(path,season,field,ar):
   if season == 's12':
      return beam_select(path,season,field,ar)
   else:
      flist = os.listdir(path)
      for f in flist:
         if season in f and ar in f and field in f:
            return os.path.join(path,f)

# misc functions for run_ps.py

def get_season_from_path(fpath):
   flist = fpath.split('/')
   for f in flist:
      if 's1' in f:
         if len(f) == 3:
            return f


def mkdir(dirpath):
   try:
      os.makedirs(dirpath)
   except:
      pass

# misc functions for pseudo-Cl

def get_tfunc(kx,ky,lmax,lbands=None):
   cut = (ky+kx)*4
   ell = np.arange(lmax+1)
   tfunc = np.zeros(lmax+1)
   tfunc[1:] = 1 - cut / (2*np.pi*ell[1:])
   if lbands is not None:
      flat_bin = get_flat_binning(lbands)
      return np.dot(flat_bin,tfunc)
   return tfunc

def get_window_from_fname(fname):
   window = None
   for n in fname.split('_'):
      if 'window' in n:
         try:
            window = int(n[6:])
         except:
            try:
               if 'dat' in n[6:] or 'txt' in n[6:] or 'mmcm' in n[6:] :
                  window = int(n[6:].split('.')[0])
            except:
               print("no window")
               return None
         return window
      else:
         continue
   return window

def get_fsky_from_fname(fname):
   for n in fname.split('_'):
      if 'fsky' in n:
         try:
            fsky = float(n[4:])
         except ValueError:
            if n.count('.')>=2:
               fsky = float('.'.join(n[4:].split('.')[:2]))
   return fsky

def get_lmax_from_fname(fname):
   for n in fname.split('_'):
      if 'lmax' in n:
         try:
            lmax = int(n[4:])
            return lmax
         except ValueError:
            for i in range(len(n.split('.'))):
               if 'lmax' in n.split('.')[i]:
                  try:
                     lmax = int(n.split('.')[i][4:])
                     return lmax
                  except ValueError:
                     if i == len(n.split('.'))-1:
                        print("huh?!?!")


def get_lmaxcut_from_fname(fname):
   for n in fname.split('_'):
      if 'lmaxcut' in n:
         try:
            lmax = int(n[7:])
            return lmax
         except ValueError:
            for i in range(len(n.split('.'))):
               if 'lmaxcut' in n.split('.')[i]:
                  try:
                     lmax = int(n.split('.')[i][7:])
                     return lmax
                  except ValueError:
                     if i == len(n.split('.'))-1:
                        print("huh?!?!")


def get_lbands_all(lbands):
   nbins = len(lbands)
   lmax = lbands[-1, -1]
   ell = np.arange(lmax + 1)
   dl = ell*(ell+1)/(2*np.pi)

   band = np.zeros(nbins)
   ell_bin = np.zeros(nbins)

   for i in range(nbins):
      band[i] = (lbands[i,1]-lbands[i,0])/2.
      ell_bin[i] = (lbands[i,0]+lbands[i,1])/2.

   return (ell_bin, band, nbins, lmax, ell, dl)

def get_flat_binning(lbands):
   def bin_op(b, l, lbands):
      if 2 <= lbands[b,0] and lbands[b,0] <= l and l <= lbands[b,1]:
         return 1.0/(lbands[b,1] - lbands[b,0] + 1.0)
      else:
         return 0
   (ell_bin, band, nbins, lmax, ell, dl) = get_lbands_all(lbands)
   flat_bin = np.zeros((nbins, lmax + 1))
   for l in range(lmax + 1):
      for b in range(nbins):
         flat_bin[b, l] = bin_op(b, l, lbands)

   return flat_bin

def get_flat_binning_inv(lbands):
   def bin_op(l, b, lbands):
      if 2 <= lbands[b,0] and lbands[b,0] <= l and l <= lbands[b,1]:
         return 1.0
      else:
         return 0
   (ell_bin, band, nbins, lmax, ell, dl) = get_lbands_all(lbands)
   flat_bin = np.zeros((lmax + 1,nbins))
   for l in range(lmax + 1):
      for b in range(nbins):
         flat_bin[l, b] = bin_op(l,b,lbands)
   return flat_bin


def get_lbands(nbins, interval, lcut):
   ret = np.zeros([nbins-1, 2])
   ret[0,0] = 2
   ret[0,1] = 1 - ret[0,0] + interval + lcut
   for i in range(1, nbins - 1):
      ret[i,0] = ret[i-1,1] + 1
      ret[i,1] = ret[i,0] + interval - 1
   ret[0,0] = lcut
   ret = np.append([[2, lcut - 1]], ret, axis = 0)
   return np.array(ret, dtype = 'int')

# def get_hivon_fsky(mask):
#    """ Hivon formula for fsky which works for equal area projections """
#    npix = len(mask.reshape(-1))
#    non0 = np.where(mask!=0)[0]
#    fs = len(non0)/float(npix)
#    w2 = np.sum(mask[non0]**2)/(npix*fs)
#    w4 = np.sum(mask[non0]**4)/(npix*fs)
#    return fs*w2**2/w4

def get_hivon_fsky(mask,area_in_sr=None):
   """ Hivon formula for fsky which works edited to take arbitrary pixelization
   with area_in_sr
   e.g.
   equal area pixelizations (EAP), such as Healpix and CEA, do not need area_in_sr 
   but non-EAP, such as CAR, can take area per pixel (in sr) map to compute w2 and w4 

   edit 170530: Healpix with no area_in_sr and CAR with area_in_sr give consistent results
   """
   sh = mask.shape
   mask = mask.reshape(-1)
   npix = len(mask)
   non0 = np.where(mask!=0)
   # fs = len(non0[0])/float(npix)
   if area_in_sr is None:
      # this is for equal area pixels
      w2 = np.sum(mask[non0]**2)/(npix)
      w4 = np.sum(mask[non0]**4)/(npix)
      mask = mask.reshape(sh)
      return w2**2/w4
   else:
      ar_sh = area_in_sr.shape
      area_in_sr = area_in_sr.reshape(-1)
      # print area_in_sr.shape

      # w2 = np.sum(area_in_sr[non0]**2*mask[non0]**2)/(npix*fs*(4*np.pi)**2)
      # w4 = np.sum(area_in_sr[non0]**4*mask[non0]**4)/(npix*fs*(4*np.pi)**4)
      
      # w2 = np.sum(area_in_sr[non0]*mask[non0]**2)/(np.sum(area_in_sr))
      # w4 = np.sum(area_in_sr[non0]*mask[non0]**4)/(np.sum(area_in_sr))
      w2 = np.sum(area_in_sr[non0]*mask[non0]**2)/(4*np.pi)
      w4 = np.sum(area_in_sr[non0]*mask[non0]**4)/(4*np.pi)
      mask = mask.reshape(sh)
      area_in_sr = area_in_sr.reshape(ar_sh)
      # print area_in_sr.shape
      return w2**2/w4

def get_hivon_fsky_enmap(mask):
   """ Hivon formula for fsky which works edited to take arbitrary pixelization
   with area_in_sr
   e.g.
   equal area pixelizations (EAP), such as Healpix and CEA, do not need area_in_sr 
   but non-EAP, such as CAR, can take area per pixel (in sr) map to compute w2 and w4 

   edit 170530: Healpix with no area_in_sr and CAR with area_in_sr give consistent results

   from blake
   """
   imap = mask/np.max(mask)
   w2w4 = (np.mean(imap**2.)**2./np.mean(imap**4.))
   skyarea = imap.area()/(4.*np.pi)*w2w4

   return skyarea

def get_pix_area(d_th,pixel,return_1d=True):
   """returns pixel area map for CEA and CAR with ordering that goes from
   north to south (Pi to -Pi; or 0 to Pi)

   edit 170530: CAR has been tested to give consistent result
   """
   from pixell import sharp
   if pixel == 'CEA':
      N_th = int(round(360.0/(np.pi*d_th)))
      N_phi = int(round(360.0/d_th))
      # minfo = sharp.map_info_gauss_legendre(N_th, N_phi)
      ret = (d_th*np.pi/180.0)**2*np.ones(N_th)
      return np.repeat(ret.reshape(len(ret),1),N_phi,axis=1)

   if pixel == 'CAR':
      N_th = int(round(180.0/d_th))+1
      N_phi = int(round(360.0/d_th))
      minfo = sharp.map_info_clenshaw_curtis(N_th,N_phi)

      theta = minfo.theta
      theta = np.pi/2.0-theta
      grad = -np.gradient(theta)
      th1 = theta-grad/2.0
      th2 = theta+grad/2.0
      # print np.where(np.sin(th2)-np.sin(th1)<0)
      ret = (np.sin(th2)-np.sin(th1))*d_th*np.pi/180.0
      return np.repeat(ret.reshape(len(ret),1),N_phi,axis=1)



def read_binfile(path,lcut=100000):
   lbands = np.array(np.genfromtxt('%s'%path)[:,:2],dtype='int')
   if len(np.where(lbands[:,1]>lcut)[0]) != 0:
      lbands = lbands[:len(lbands)-len(np.where(lbands[:,1]>lcut)[0])]
   if lbands[0,0] == 0:
      lbands[0,0] = 2
   return lbands


def get_analytic_var_auto(fsky, ell_bin, lband_width):
   """ compute total variance for auto spectra:
       get sig(C_l) from C_l (Knox 1995)
       watch the lband_width input which is half the actual width
   """
   return np.sqrt(2.0 / ((2.0 * ell_bin + 1.0) * fsky * (lband_width * 2 + 1)))


def get_analytic_var_cross(fsky, ell_bin, lband_width):
   """ compute total variance for cross spectra:
       get sig(C_l) from C_l (Knox 1995)
       watch the lband_width input which is half the actual width
   """
   return np.sqrt(2.0 / ((2.0 * ell_bin + 1.0) * fsky * (lband_width * 2 + 1)))/np.sqrt(2)

# make sin apodized fields for simulations
# the convention goes [0, N_RA] for [0, 2pi] in RA and [0, N_dec] for [0, pi] in dec

# stuff below is wrong:

# def dec_ang2pix_not_use(d_th, dec):
#    """ provide dec angle in degrees get 0th pixel index in 2d array"""
#    n_dec = int(180.0/d_th)
#    i_dec = int((1.0-np.cos(np.radians(dec))) / 2.0 * n_dec)
#    return i_dec


# def dec_ang2pix(d_th, dec):
#    """ provide dec angle in degrees get 0th pixel index in 2d array"""
#    n_dec = int(180.0/d_th)
#    i_dec = int(n_dec/2.0 - np.cos(np.radians(dec)) * n_dec/2.0)
#    return i_dec


# def ra_ang2pix(d_th, ra):
#    """ provide RA angle in degrees get 1st pixel index in 2d array"""
#    n_ra = int(360.0/d_th)
#    i_ra = int(ra/360.0*n_ra)
#    return i_ra


# def ang2pix(d_th, dec, ra):
#    """ provide RA, dec angles in degrees get 0th, and 1st pixel indices in 2d array"""
#    return dec_ang2pix(d_th, dec), ra_ang2pix(d_th, ra)


# def pix2ang(d_th, i_dec, i_ra):
#    """ provide indices in 2d array and get corresponding RA, dec in degrees"""
#    n_dec, n_ra = int(180.0/d_th), int(360.0/d_th)
#    ra = i_ra/float(n_ra)*360.0
#    dec = np.degrees(np.arccos(-(i_dec/float(n_dec) * 2.0 - 1.0)))
#    return dec, ra


# def ang_dist(d_th, i_dec1, i_ra1, i_dec2, i_ra2):
#    """ compute angular distance between two pixels on 2d array, assuming flat geometry"""
#    dec1, ra1 = pix2ang(d_th, i_dec1, i_ra1)
#    dec2, ra2 = pix2ang(d_th, i_dec2, i_ra2)
#    return np.sqrt((dec1-dec2)**2 + (ra1-ra2)**2)

# stuff below is probably correct:

class pix_ang_trans:
   def __init__(self,mpath):
      hdulist = pf.open(mpath)
      # m = hdulist[0].data
      self.w = wcs.WCS(hdulist[0].header, hdulist)
      hdulist.close()

   def ang2pix(self, ra, dec):
      """ provide RA, dec angles in degrees and get pixel index for the 2d map input
      RA goes from 0 to 360, dec goes from 90 to -90
      """
      ra, dec = self.w.all_world2pix(ra, dec, 1)
      return [np.array([ra],'int').reshape(-1), np.array([dec],'int').reshape(-1)]

   def pix2ang(self, ra, dec):
      """ provide RA, dec pixel index (or indices) for the 2d map input
      RA goes from 0 to 360, dec goes from 90 to -90
      """
      return self.w.all_pix2world(ra, dec, 1)



def get_sin_apomask(d_th, radius, dec, ra, apo_length):
   """ get sin apodized mask by providing, pixel length, radius, center dec and RA coordinates
       and length over which apodization happens
   """
   n_dec, n_ra = int(180.0/d_th), int(360.0/d_th)
   i_dec, i_ra = ang2pix(d_th, dec, ra)
   l_ra = radius/360.0 * n_ra
   l_dec_n = dec_ang2pix(d_th, dec-radius)
   l_dec_s = dec_ang2pix(d_th, dec+radius)

   ret = np.zeros([n_dec, n_ra])
   xtry = np.array(np.arange(i_ra-l_ra, i_ra+l_ra),dtype='int')

   for y in range(l_dec_n, l_dec_s):
      ret[y, xtry[np.where(ang_dist(d_th, i_dec, i_ra, y, xtry) < radius-apo_length)[0]]] = 1.0
      keep = np.where((ang_dist(d_th, i_dec, i_ra, y, xtry) <= radius) & (ang_dist(d_th, i_dec, i_ra, y, xtry) >= radius-apo_length))[0]
      for i in xtry[keep]:
         dist = ang_dist(d_th, i_dec, i_ra, y, i)
         ret[y, i] = -1.0 / (2.0 * np.pi) * np.sin(2.0 * np.pi * (radius - dist) / apo_length) + (radius - dist) / apo_length

   return ret

def get_sin_apomask_enmap(m_apo,apo_dist=120,cut_dist=10):
   """ get sin apodized mask for ndmaps.
   assumes 0.5 arcmin CAR maps.
   because it uses pixel counts, the output is not quite in degrees
   apo_dist=120 e.g., is 1 meant to be 1 deg

   """
   import scipy as sp
   import scipy.ndimage

   # apo_dist = 120
   # cut_dist = 10

   outdeg = apo_dist*0.5/60 # deg

   smooth_noise = 0
   remove_holes = 0

   boundary = copy.copy(m_apo)
   boundary[boundary!=0] = 1.
   boundary[boundary!=1] = 0.
   # get distance from zero
   dist = scipy.ndimage.distance_transform_edt(boundary)
   # throw out cut_dist (15 arcmin) 
   boundary[dist<cut_dist] = 0.
   if remove_holes:
       boundary[scipy.ndimage.distance_transform_edt(1-boundary)<=remove_holes] = 1
       boundary[scipy.ndimage.distance_transform_edt(boundary)<=remove_holes] = 0

   dist = scipy.ndimage.distance_transform_edt(boundary)
   dist[dist>apo_dist] = np.min(dist[dist>apo_dist])
   # make linear apodization into sin
   dist_max = np.max(dist)
   keep = np.where((dist>0.)&(dist<dist_max))
   edge = dist[keep]
   edge = dist_max - ((dist_max-edge)/dist_max-1./(2*np.pi)*np.sin(2*np.pi*(dist_max-edge)/dist_max))
   dist[keep] = edge
   dist[dist>0] -= np.min(edge)
   dist /= np.max(dist)
   # dist *= mask
   # dist[dist==0] = np.inf
   # im = plt.imshow(dist, origin='lower',interpolation='none')
   # plt.colorbar(im, orientation='horizontal')
   # plt.show()
   if smooth_noise:
       smooth_noise_sigma = 5
       print("smoothing noise maps with sigma = %i"%smooth_noise_sigma)
       m_apo = sp.ndimage.filters.gaussian_filter(m_apo, smooth_noise_sigma, mode='mirror')

   window = dist*m_apo
   apo_mask = enmap.samewcs(dist,m_apo)

   return apo_mask


def get_enmap_bounding_box(m,offset=5,cross=True):
    """get biggest bounding box + offset arcmin
    offset in arcmin
    cross = True crosses 180 deg in RA
    tested to separate DR6 region into two
    """
    off = np.radians(offset/60.)
    zeros = np.where(m!=0)
    if cross:
        top = m.pix2sky([np.min(zeros[0]),np.min(zeros[1])])-np.array([off,-off])
        bot = m.pix2sky([np.max(zeros[0]),np.max(zeros[1])])+np.array([off,-off])
        if top[1]<bot[1] and top[1]<0:
            top[1] = 2*np.pi+top[1]
        box = np.array([top,bot])
    else:
        top = m.pix2sky([np.min(zeros[0]),np.min(zeros[1])])-np.array([off,-off])
        bot = m.pix2sky([np.max(zeros[0]),np.max(zeros[1])])+np.array([off,-off])
        box = np.array([top,bot])        
    return box
# functions for act spectra


class split_null:
   def __init__(self,cross_index,cov,cl):
      self.cross_index = cross_index
      self.cov = cov
      self.cl = cl

   def get_bigger(self,i,j):
      if i>j:
         return '%i%i'%(i,j)
      else:
         return '%i%i'%(j,i)

   def get_null_inds(self,ind,i,j,k,l):
      null_inds = []
      null_inds.append(np.where(ind==self.get_bigger(i,k))[0][0])
      null_inds.append(np.where(ind==self.get_bigger(j,l))[0][0])
      null_inds.append(np.where(ind==self.get_bigger(i,l))[0][0])
      null_inds.append(np.where(ind==self.get_bigger(j,k))[0][0])
      return null_inds

   def get_null(self,outname,spec,fpath,ell_bin,i,j,k,l,bs,be,tfunc=None,plot_on=False):
      from scipy.stats import chisqprob
      outdir,oname = os.path.split(outname)
      null_inds = self.get_null_inds(self.cross_index,i,j,k,l)
      P = np.array([1,1,-1,-1])
      nbins = self.cl.shape[-1]
      null_list = self.cl[null_inds]
      null = np.zeros(nbins)
      null_err = np.zeros(nbins)
      for n in range(nbins):
         null[n] = np.dot(P,null_list[:,n])
         null_err[n] = np.sqrt(np.dot(P,np.dot(self.cov[null_inds][:,null_inds][:,:,n],P)))

      if tfunc is not None:
         null_err /= np.sqrt(tfunc)

      bin_start = bs
      bin_end = be

      chi2 = np.sum((null/null_err)[bs:be]**2)
      dof = len(null[bs:be])
      pte = chisqprob(chi2,dof)
      print('chi^2 = %.3f PTE =%.3f'%(chi2,pte))

      if plot_on:
         if spec == 'TT' and 'f150' in fpath:
            plt.errorbar(ell_bin[bs:be],(ell_bin*null)[bs:be]/1e6,(ell_bin*null_err)[bs:be]/1e6,
               fmt='o',label='$\chi^2=%.3f/%i$\n$\mathrm{PTE}=%.3f$'%(chi2,dof,pte),
               color='RoyalBlue',alpha=0.9,markersize=4,elinewidth=1.5)
            plt.ylabel('$\ell\mathcal{D}_\ell^{%s}/1e6 \, [\mu\mathrm{K}^2]$'%spec,fontsize=16)
            plt.ylim(-1,1)
         elif spec == 'EE':
            plt.errorbar(ell_bin[bs:be],(null/np.sqrt(ell_bin))[bs:be],(null_err/np.sqrt(ell_bin))[bs:be],
               fmt='o',label='$\chi^2=%.3f/%i$\n$\mathrm{PTE}=%.3f$'%(chi2,dof,pte),
               color='RoyalBlue',alpha=0.9,markersize=4,elinewidth=1.5)
            plt.ylabel('$\mathcal{D}_\ell^{%s}/\sqrt{\ell} \, [\mu\mathrm{K}^2]$'%spec,fontsize=16)
         elif spec == 'BB':
            plt.errorbar(ell_bin[bs:be],(null/np.sqrt(ell_bin))[bs:be],(null_err/np.sqrt(ell_bin))[bs:be],
               fmt='o',label='$\chi^2=%.3f/%i$\n$\mathrm{PTE}=%.3f$'%(chi2,dof,pte),
               color='RoyalBlue',alpha=0.9,markersize=4,elinewidth=1.5)
            plt.ylabel('$\mathcal{D}_\ell^{%s}/\sqrt{\ell} \, [\mu\mathrm{K}^2]$'%spec,fontsize=16)
         else:
            plt.errorbar(ell_bin[bs:be],null[bs:be],null_err[bs:be],
               fmt='o',label='$\chi^2=%.3f/%i$\n$\mathrm{PTE}=%.3f$'%(chi2,dof,pte),
               color='RoyalBlue',alpha=0.9,markersize=4,elinewidth=1.5)
            plt.ylabel('$\mathcal{D}_\ell^{%s} \, [\mu\mathrm{K}^2]$'%spec,fontsize=16)

         plt.xlabel('$\ell$',fontsize=16)
         plt.legend(numpoints=1,loc=0)

         plt.axhline(y=0,color='k')
         plt.grid()
         plt.xscale('log')
         plt.title('%s\nsplit null (%i-%i)x(%i-%i)'%(oname,i,j,k,l))
         plt.savefig('%s_%s_split_null_%i-%ix%i-%i.png'%(outname,spec,i,j,k,l))
         plt.clf()
         # plt.show()

         # plt.plot(ell_bin, null_list[0])
         # plt.plot(ell_bin, null_list[1])
         # plt.plot(ell_bin, null_list[2])
         # plt.plot(ell_bin, null_list[3])
         # plt.show()
         # plt.plot(ell_bin, null_list[0]-null_list[1])
         # plt.plot(ell_bin, null_list[2]-null_list[3])
         # plt.plot(ell_bin, null_list[0]-null_list[3])
         # plt.plot(ell_bin, null_list[1]-null_list[2])
         # plt.plot(ell_bin, null/null_err)
         # plt.show()


         # from scipy.stats import chi2
         # def chi2_vs_df(x, df):
         #    return chi2(df).pdf(x)

         # N = len((null/null_err)[bs:be])
         # plt.hist((null/null_err)[bs:be]**2,normed=True,alpha=0.3,bins=np.arange(30)-0.5)
         # plt.plot(np.arange(30),chi2_vs_df(np.arange(30),1))
         # plt.show()

         # exit()


# functions for actpol maps


def return_paths(pathfile):
   ret = []
   for l in open(pathfile):
      if len(l.split()) > 0:
         if l.split()[1] == '=':
            ret.append(l.split()[2])

   return ret


def load_dict(dictpath):
  with open(dictpath, 'r') as file:
    return json.load(file)

def load_js_map(dpath):
   return pf.open(dpath)[0].data

def load_skn_map(dpath,comp):
   return pf.open(dpath)[0].data[comp]

def load_skn_map_enmap(dpath,comp):
   return enmap.read_map(dpath)[comp]

# older version
def get_js_map_pos_not_use(dpath):
   d_th = abs(pf.open(dpath)[0].header['CDELT1'])
   ra1 = int(180/d_th)-pf.open(dpath)[0].header['CRPIX1']
   dec1 = int(90/d_th)-pf.open(dpath)[0].header['CRPIX2']
   D_ra = pf.open(dpath)[0].header['NAXIS1']
   D_dec = pf.open(dpath)[0].header['NAXIS2']
   return np.array([d_th, dec1, dec1+D_dec, ra1, ra1+D_ra])

def extend_full_sky_not_use(m,x,return_1d=True):
   (d_th,dec1,dec2,ra1,ra2) = x
   ret = np.zeros([int(180.0/d_th), int(360.0/d_th)])
   ret[dec1:dec2,ra1:ra2] = m
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret


def get_sht(wpar, lmax, alm_order='triangular', zmask=None):
   """ 
   alm_order can be 'triangular' or 'rectangular' when doing mcut
   pass mpos for wpar if CEA or CAR; pass nside for healpix
   """
   from pixell import sharp

   if isinstance(wpar,int):
      minfo = sharp.map_info_healpix(nside=wpar)
      # print "[psc_note] maps being analyzed are HEALPIX"
   else:
      # print "[psc_note] mpos in pcl_actpool_utility hard coded"
      d_th = wpar[0]
      N_th = wpar[-2]
      N_phi = wpar[-1]

      if N_th == N_phi*0.5:
         minfo = sharp.map_info_gauss_legendre(N_th, N_phi)
         # print "[psc_note] maps being analyzed are CEA (but assuming Gauss Legendre)"

      if N_th == N_phi*0.5+1:
         minfo = sharp.map_info_clenshaw_curtis(N_th, N_phi)
         # print "[psc_note] maps being analyzed are CAR"

   if zmask is not None:
      # print "[psc_note] using zmask to redefine minfo"
      # print minfo.theta[zmask]
      minfo = sharp.map_info(theta=minfo.theta[zmask], nphi=minfo.nphi[zmask],
                        phi0=minfo.phi0[zmask], offsets=minfo.offsets[zmask]-minfo.offsets[zmask][0],
                        stride=minfo.stride[zmask],  weight=minfo.weight[zmask])

   ainfo = sharp.alm_info(lmax=lmax,layout=alm_order)
   return sharp.sht(minfo, ainfo)


def get_sht_enmap(m, lmax, alm_order='triangular', zmask=None):
   """ 
   alm_order can be 'triangular' or 'rectangular' when doing mcut
   pass mpos for wpar if CEA or CAR; pass nside for healpix
   """
   from pixell import sharp

   if isinstance(m,int):
      minfo = sharp.map_info_healpix(nside=m)
      # print "[psc_note] maps being analyzed are HEALPIX"
   else:
      print("[psc_note] mpos in pcl_actpool_utility hard coded")
      d_th = abs(m.wcs.to_header()['CDELT1'])
      sh, wcs = enmap.fullsky_geometry(res=d_th*np.pi/180.,proj="car")
      N_th = sh[0]
      N_phi = sh[1]

      if N_th == N_phi*0.5:
         minfo = sharp.map_info_gauss_legendre(N_th, N_phi)
         # print "[psc_note] maps being analyzed are CEA (but assuming Gauss Legendre)"

      if N_th == N_phi*0.5+1:
         minfo = sharp.map_info_clenshaw_curtis(N_th, N_phi)
         # print "[psc_note] maps being analyzed are CAR"

   if zmask is not None:
      # print "[psc_note] using zmask to redefine minfo"
      # print minfo.theta[zmask]
      minfo = sharp.map_info(theta=minfo.theta[zmask], nphi=minfo.nphi[zmask],
                        phi0=minfo.phi0[zmask],
                        offsets=minfo.offsets[zmask]-minfo.offsets[zmask][0],
                        stride=minfo.stride[zmask],weight=minfo.weight[zmask])

   ainfo = sharp.alm_info(lmax=lmax,layout=alm_order)
   return sharp.sht(minfo, ainfo)


# after realizing maps are saved with decreasing RA increasing dec
# e.g. 0 to N_phi-1 would be 50 to -19 for d56 in RA
#  and 0 to N_th-1 would be -9 to 7 for d56 in dec

# works for d56 and boss north but probably not for d8 and d9
def get_ACT_map_pos(dpath,field=0):
   """when ACT maps are loaded the y axis (dec) is in increasing order 
   (south [180 or -Pi] to north [0 or Pi]) and RA is in decreasing 
   order(2Pi to Pi or Pi to -Pi)
   """
   pfobj = pf.open(dpath)
   assert np.sign(pfobj[field].header['CDELT1'])==-1 # making sure RA is in decreasing order 2pi->0
   assert np.sign(pfobj[field].header['CDELT2'])==1 # making sure dec is in increasing order 0->pi
   N_ra = pfobj[field].header['NAXIS1'] # map size in x
   N_dec = pfobj[field].header['NAXIS2'] # map size in y
   d_th = abs(pfobj[field].header['CDELT1']) # pixel size
   ra_ref = pfobj[field].header['CRVAL1'] # reference pixel RA position in deg
   dec_ref = pfobj[field].header['CRVAL2'] # reference pixel dec position in deg
   pixtype = pfobj[field].header['CTYPE1'] # to get pixelization type
#  
   if 'CAR' in pixtype:
      # print "pixel type for mapinfo: CAR!"
      N_th = int(round(180.0/d_th))+1
      N_phi = int(round(360.0/d_th))
   elif 'CEA' in pixtype:
      # print "pixel type for mapinfo: CEA!"
      N_th = int(round(180.0/d_th))
      N_phi = int(round(360.0/d_th))
   else:
      print("what pixelization?!?!")
      sys.exit(2)
#
   ra_ref_pix = ra_ref/360.0 * N_phi
   dec_ref_pix = dec_ref/180.0 * N_th
#
   ra1 = int(round(ra_ref_pix + (pfobj[field].header['CRPIX1']-1))) # this is the bigger RA corner unless it crosses 0; it's the positive RA if the other side is negative
   ra2 = int(round(ra_ref_pix - (N_ra - (pfobj[field].header['CRPIX1']-1)))) # this goes negative when RA crosses 0

   if ra1 > N_phi:
      if ra2 > ra1 - N_phi:
         ra1 -= N_phi
         ra2 -= N_phi
#
   dec1 = int(round(N_th/2.0 - (pfobj[field].header['CRPIX2']-1))) # this is the dec corner closer to 0
   dec2 = int(round(N_th/2.0 + (N_dec - (pfobj[field].header['CRPIX2']-1)))) # this is the dec corner closer to pi
#
   return np.array([d_th, dec1, dec2, ra1, ra2, N_th, N_phi])

def get_hers_map_pos(dpath,field=0):
   """when ACT maps are loaded the y axis (dec) is in increasing order 
   (south [180 or -Pi] to north [0 or Pi]) and RA is in decreasing 
   order(2Pi to Pi or Pi to -Pi)
   """
   pfobj = pf.open(dpath)
   assert np.sign(pfobj[field].header['CD1_1'])==-1 # making sure RA is in decreasing order 2pi->0
   assert np.sign(pfobj[field].header['CD2_2'])==1 # making sure dec is in increasing order 0->pi
   N_ra = pfobj[field].header['NAXIS1'] # map size in x
   N_dec = pfobj[field].header['NAXIS2'] # map size in y
   d_th = abs(pfobj[field].header['CD1_1']) # pixel size
   ra_ref = pfobj[field].header['CRVAL1'] # reference pixel RA position in deg
   dec_ref = pfobj[field].header['CRVAL2'] # reference pixel dec position in deg
   pixtype = pfobj[field].header['CTYPE1'] # to get pixelization type
#  
   if 'CAR' in pixtype:
      # print "pixel type for mapinfo: CAR!"
      N_th = int(round(180.0/d_th))+1
      N_phi = int(round(360.0/d_th))
   elif 'CEA' in pixtype:
      # print "pixel type for mapinfo: CEA!"
      N_th = int(round(180.0/d_th))
      N_phi = int(round(360.0/d_th))
   else:
      print("what pixelization?!?!")
      sys.exit(2)
#
   ra_ref_pix = ra_ref/360.0 * N_phi
   dec_ref_pix = dec_ref/180.0 * N_th
#
   ra1 = int(round(ra_ref_pix + (pfobj[field].header['CRPIX1']-1))) # this is the bigger RA corner unless it crosses 0; it's the positive RA if the other side is negative
   ra2 = int(round(ra_ref_pix - (N_ra - (pfobj[field].header['CRPIX1']-1)))) # this goes negative when RA crosses 0

   if ra1 > N_phi:
      if ra2 > ra1 - N_phi:
         ra1 -= N_phi
         ra2 -= N_phi
#
   dec1 = int(round(N_th/2.0 - (pfobj[field].header['CRPIX2']-1))) # this is the dec corner closer to 0
   dec2 = int(round(N_th/2.0 + (N_dec - (pfobj[field].header['CRPIX2']-1)))) # this is the dec corner closer to pi
#
   return np.array([d_th, dec1, dec2, ra1, ra2, N_th, N_phi])



def get_pix_coord(mpos,RA,Dec):
   """
   wrong !!! use astropy.wcs instead!!!!
   """
   (d_th,dec1,dec2,ra1,ra2,N_th,N_phi) = mpos
   ra = copy.copy(RA)
   dec = copy.copy(Dec)
   d1 = -90.0+dec1/N_th*180.0
   d2 = -90.0+dec2/N_th*180.0
   r1 = ra1/N_phi*360.0
   r2 = ra2/N_phi*360.0
   # print r1, r2, d1, d2, N_th, N_phi
   if r2 < 0: # RA crossing 0
      print("RA crossing 0")
      ra[ra>r1] = ra[ra>r1] - 360.0
   RApix = np.array(np.round((r1-ra)*N_phi/360.0),'int')
   decpix = np.array(np.round((dec-d1)*N_th/180.0),'int')
   return RApix, decpix


def extend_full_sky(m,x,return_1d=True):
   """
   takes an act map m and the mpos variables then puts it on the full sky
   in the ordering ready for sharp map2alm 
   (2d array with y axis (dec) going 0 (Pi) to 180 (-Pi) 
   and x (RA) axis 0 to 360 (increasing))
   """
   (d_th,dec1,dec2,ra1,ra2,N_th,N_phi) = x
   ret = np.zeros([int(N_th), int(N_phi)])
   if ra2 > 0:
      ret[int(dec1):int(dec2),int(ra2):int(ra1)] = np.fliplr(m) # fliplr is applied to since ACT maps saved in RA decreasing order
   else:
      if ra1>0:
         ret[int(dec1):int(dec2),:int(ra1)] = np.fliplr(m[:,:int(ra1)])
         ret[int(dec1):int(dec2),int(N_phi)+int(ra2):] = np.fliplr(m[:,int(ra1):])
      else:
         ret[int(dec1):int(dec2),int(ra2):int(ra1)] = np.fliplr(m) # fliplr is applied to since ACT maps saved in RA decreasing order

   ret = np.flipud(ret) # flipud is applied to be compatible with sharp harmonic transform i.e. dec going from north to south
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret

def get_patch_from_fs(m_fs,x,return_1d=True):
   """
   takes an full sky act map m_fs and the mpos varaiables then just grabs the 
   patch given by mpos variables and reorders to save in the ACT map ordering
   (2d array with y axis going 180 to 0 and x axis 360 0)
   """
   (d_th,dec1,dec2,ra1,ra2,N_th,N_phi) = x
   Nx = int(ra1-ra2)
   Ny = int(dec2-dec1)
   # print dec2, dec1
   ret = np.zeros([Ny,Nx])
   Ny1 = int(N_th - dec2)
   Ny2 = int(N_th - dec1)
   if ra2 > 0:
      ret = np.fliplr(m_fs[Ny1:Ny2,int(ra2):int(ra1)]) # fliplr is applied since ACT maps saved in RA decreasing order
   else:
      if ra1>0:
         ret[:,:int(ra1)] = np.fliplr(m_fs[Ny1:Ny2,:int(ra1)])
         ret[:,int(ra1):] = np.fliplr(m_fs[Ny1:Ny2,int(N_phi)+int(ra2):])
      else:
         ret = np.fliplr(m_fs[Ny1:Ny2,int(ra2):int(ra1)]) # fliplr is applied since ACT maps saved in RA decreasing order

   ret = np.flipud(ret) # flipud is applied since ACT maps saved in dec going from 180 to 0
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret

def enmap_to_sharp(d):
   ret = copy.copy(d)
   N_phi = ret.shape[-1]
   ret[:,:int(N_phi/2)] = np.fliplr(ret[:,:int(N_phi/2)])
   ret[:,int(N_phi/2):] = np.fliplr(ret[:,int(N_phi/2):])
   ret = np.flipud(ret)
   return ret

def extend_full_sky_enmap(m,return_1d=True):
   """need to test!!"""
   d_th = abs(m.wcs.to_header()['CDELT1'])
   sh, wcs = enmap.fullsky_geometry(res=d_th*np.pi/180.,proj="car")
   ret = enmap.extract(m, sh, wcs)
   ret = enmap_to_sharp(ret)
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret

def get_patch_from_fs_enmap(m_fs,m_patch,return_1d=True):
   """need to test!!"""
   ret = enmap_to_sharp(m_fs)
   if return_1d:
      return enmap.extract(ret,m_patch.shape,m_patch.wcs).reshape(-1)
   else:
      return enmap.extract(ret,m_patch.shape,m_patch.wcs)

def get_zcut(mask,mpos=None):
   """ returns the nonzero region of the mask 
   """
   if mpos is None:
      keep = np.where(mask>0)
      ylim = (np.min(keep[0]),np.max(keep[0]))
      return range(ylim[0],ylim[1]+1)
   else:
      (d_th,dec1,dec2,ra1,ra2,N_th,N_phi) = mpos
      return range(int(N_th)-int(mpos[2]),int(N_th)-int(mpos[1]))

def get_zcut_enmap(mask):
   """ input mask here is cut sky
   """
   d_th = abs(mask.wcs.to_header()['CDELT1'])
   sh, wcs = enmap.fullsky_geometry(res=d_th*np.pi/180.,proj="car")
   N_th = sh[0]
   coord0 = enmap.sky2pix(sh,wcs,mask.box()[0])
   coord1 = enmap.sky2pix(sh,wcs,mask.box()[1])
   return range(N_th-int(coord1[0]+0.5),N_th-int(coord0[0]-0.5))


def reproject(m,mpos,outpos):
   """ input map m, output to have the footprint of another """
   tmp = extend_full_sky(m,mpos,return_1d=False)
   return get_patch_from_fs(tmp,outpos,return_1d=False)


def extend_full_sky_zcut(m,x,zmask,return_1d=True):
   """ zmask works on the 2d array to keep only the desired theta
   """
   (d_th,dec1,dec2,ra1,ra2,N_th,N_phi) = x
   ret = np.zeros([int(N_th), int(N_phi)])
   if ra2 > 0:
      ret[int(dec1):int(dec2),int(ra2):int(ra1)] = np.fliplr(m)
   else:
      if ra1>0:
         ret[int(dec1):int(dec2),:int(ra1)] = np.fliplr(m[:,:int(ra1)])
         ret[int(dec1):int(dec2),int(N_phi)+int(ra2):] = np.fliplr(m[:,int(ra1):])
      else:
         ret[int(dec1):int(dec2),int(ra2):int(ra1)] = np.fliplr(m)

   ret = np.flipud(ret)[zmask]
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret

def extend_full_sky_zcut_enmap(m,zmask,return_1d=True):
   """need to test!!"""
   ret = extend_full_sky_enmap(m,return_1d=False)[zmask]
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret

def get_patch_from_fs_zcut(m_fs,x,return_1d=True):
   """
   takes an zcut sky map m_fs and the mpos varaiables then just grabs the 
   patch given by mpos variables and reorders to save in the ACT map ordering
   (2d array with y axis going 180 to 0 and x axis 360 0)
   """
   (d_th,dec1,dec2,ra1,ra2,N_th,N_phi) = x
   Nx = int(ra1-ra2)
   Ny = int(dec2-dec1)
   # print dec2, dec1
   ret = np.zeros([Ny,Nx])
   Ny1 = int(N_th - dec2)
   Ny2 = int(N_th - dec1)
   if ra2 > 0:
      ret = np.fliplr(m_fs[:,int(ra2):int(ra1)]) # fliplr is applied since ACT maps saved in RA decreasing order
   else:
      ret[:,:int(ra1)] = np.fliplr(m_fs[:,:int(ra1)])
      ret[:,int(ra1):] = np.fliplr(m_fs[:,int(N_phi)+int(ra2):])

   ret = np.flipud(ret) # flipud is applied since ACT maps saved in dec going from 180 to 0
   if return_1d:
      return ret.reshape(-1)
   else:
      return ret

def fractional_shift(map, off=[0.5,0], keepwcs=False, nofft=False):
  """Shift map cyclically by a non-integer amount off [{y_off,x_off}]
  to deal with CAR maps have 0.5 pixel shift with DR6v4 and on"""
  omap = enmap.samewcs(enfft.shift(map, off, nofft=nofft), map)
  if not keepwcs: omap.wcs.wcs.crpix += off[::-1]
  return omap

def combine_two_maps(m_a,mpos_a,m_b,mpos_b,favor=None,sum=False,save=False,mpath=None,outname=None,clobber=True):
   """ this sums or multiplies two maps that have different footprints
   the resulting footprints include the boundaries of both maps
   favor selects footprint of a or b or combined for options {None,a,b}
   """

   (d_th,dec1a,dec2a,ra1a,ra2a,N_th,N_phi) = mpos_a
   (d_th,dec1b,dec2b,ra1b,ra2b,N_th,N_phi) = mpos_b
   dec1 = min(dec1a,dec1b)
   dec2 = max(dec2a,dec2b)
   ra1 = max(ra1a,ra1b)
   ra2 = min(ra2a,ra2b)
   ret_a = np.zeros([dec2-dec1, ra1-ra2])
   ret_b = np.zeros([dec2-dec1, ra1-ra2])
   ret_a[dec1a-dec1:dec2a-dec1,ra1-ra1a:ra1-ra2a] = m_a # dec1a-dec1+(dec2a-dec1a) ra1-ra1a+(ra1a-ra2a)
   ret_b[dec1b-dec1:dec2b-dec1,ra1-ra1b:ra1-ra2b] = m_b # dec1a-dec1+(dec2a-dec1a) ra1-ra1a+(ra1a-ra2a)
   if sum:
      ret = ret_a+ret_b
   else:
      ret = ret_a*ret_b
   if favor == "a":
      ret = ret[dec1a-dec1:dec2a-dec1,ra1-ra1a:ra1-ra2a]
      ra1 = ra1a
      ra2 = ra2a
      dec1 = dec1a
      dec2 = dec2a
   if favor == "b":
      ret = ret[dec1b-dec1:dec2b-dec1,ra1-ra1b:ra1-ra2b]
      ra1 = ra1b
      ra2 = ra2b
      dec1 = dec1b
      dec2 = dec2b
   if save:
      N_ra = len(ret[0])
      N_dec = len(ret)
      d_hdu = pf.PrimaryHDU(ret)
      d_hdu.header = pf.open(mpath)[0].header
      d_hdu.header['NAXIS1'] = N_ra
      d_hdu.header['NAXIS2'] = N_dec

      d_th = abs(d_hdu.header['CDELT1']) # pixel size
      ra_ref = d_hdu.header['CRVAL1'] # reference pixel RA position in deg
      dec_ref = d_hdu.header['CRVAL2'] # reference pixel dec position in deg
      pixtype = d_hdu.header['CTYPE1'] # to get pixelization type
     
      if 'CAR' in pixtype:
         N_th = int(round(180.0/d_th))+1
         N_phi = int(round(360.0/d_th))
      elif 'CEA' in pixtype:
         N_th = int(round(180.0/d_th))
         N_phi = int(round(360.0/d_th))
      else:
         print("what pixelization?!?!")
         sys.exit(2)

      ra_ref_pix = int(round(ra_ref/360.0 * N_phi))
      dec_ref_pix = int(round(dec_ref/180.0 * N_th))

      d_hdu.header['CRPIX1'] = ra1-ra_ref_pix+1
      assert d_hdu.header['CRPIX1'] == ra2-ra_ref_pix+N_ra+1
      d_hdu.header['CRPIX2'] = int(round(N_th/2.0))-dec1+1
      assert d_hdu.header['CRPIX2'] == int(round(N_th/2.0))-dec2+N_dec+1
      d_hdu.writeto(outname,clobber=clobber)
      return ret
   else:
      return ret


def block_mean(ar, fact):
   assert isinstance(fact, int), type(fact)
   sx, sy = ar.shape
   X, Y = np.ogrid[0:sx, 0:sy]
   regions = sy/fact * (X/fact) + Y/fact
   res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
   res.shape = (sx/fact, sy/fact)
   return res


def downgrade(emap, factor):
   # from sigurd
   """Returns enmap "emap" downgraded by the given integer factor
   (may be a list for each direction, or just a number) by averaging
   inside pixels."""
   fact = np.full(2, 1, dtype=int)
   fact[:] = factor
   tshape = emap.shape[-2:]/fact*fact
   res = np.mean(np.reshape(emap[...,:tshape[0],:tshape[1]],emap.shape[:-2]+(tshape[0]/fact[0],fact[0],tshape[1]/fact[1],fact[1])),(-3,-1))
   return res

# stuff below not tested

def gen_white_noise_map(uK_arcmin,area_map):
   # area_map in sr
   noise = np.random.normal(0,uK_arcmin,size=area_map.shape)
   area_arcmin = (1/60.0)**2*(np.pi/180.0)**2
   return noise*np.sqrt(area_arcmin/area_map)


def gen_noise_map(Nl,d_th,area_map):
   import matplotlib.pyplot as plt

   lmax = len(Nl)-1
   ell = np.arange(lmax+1)
   dl = ell*(ell+1)/(2*np.pi)
   Nl[2:] /= dl[2:]

   Ny,Nx = area_map.shape

   ell_scale_factor = 2*np.pi/(d_th*np.pi/180.)

   inds  = [(np.arange(Nx)+.5 - Nx/2.) /(Nx-1.)]
   X = np.repeat(inds,Ny,axis=0) * ell_scale_factor

   inds  = (np.arange(Ny)+.5 - Ny/2.) /(Ny-1.)
   Y = np.repeat(inds.reshape(Ny,1),Nx,axis=1) * ell_scale_factor

   ell2d = np.sqrt(X**2. + Y**2.)
   Nl2d = Nl[ell2d.astype(int)]

   random_array = np.random.normal(0,1.0,(Ny,Nx))
   random_array2 = np.random.normal(0,1.0,(Ny,Nx))
   print(np.var(random_array))

   normalization = None

   FT_random_array = np.fft.fft2(random_array,norm=normalization)
   print(np.var(FT_random_array.real))
   print(np.var(FT_random_array.imag))

   # FT_random_array /= np.std(FT_random_array.real)

   # plt.imshow(random_array.real)
   # plt.show()
   # plt.imshow(FT_random_array.real)
   # plt.show()
   # plt.imshow(FT_random_array.imag)
   # plt.show()

   # exit()
   # random_array = np.random.normal(0,1,(Ny,Nx))

   FT_2d = np.sqrt(Nl2d) * FT_random_array
   # FT_2d = np.sqrt(Nl2d) * (random_array+random_array*1j)

   # plt.imshow(ell2d)
   # plt.show()

   # plt.imshow(Nl2d)
   # plt.show()

   # plt.imshow(FT_2d.real)
   # plt.show()

   noise_map = np.fft.ifft2(np.fft.fftshift(FT_2d),norm=normalization)/(d_th*np.pi/180.)
   # noise_map = np.real(noise_map)

   plt.imshow(noise_map.real)
   plt.show()

   plt.imshow(noise_map.imag)
   plt.show()

   # print (d_th*np.pi/180.0)/np.sqrt(area_map)
   noise_map *= (d_th*np.pi/180.0)/np.sqrt(area_map)

   Nl[2:] *= dl[2:]

   return np.real(noise_map)-np.imag(noise_map)


# N = 2**10.  # this is the number of pixels in a linear dimension
#             ## since we are using lots of FFTs this should be a factor of 2^N
# pix_size  = 0.5 # size of a pixel in arcminutes

# ## variables to set up the map plots
# c_min = -400  # minimum for color bar
# c_max = 400   # maximum for color bar
# X_width = N*pix_size/60.  # horizontal map width in degrees
# Y_width = N*pix_size/60.  # vertical map width in degrees


# ell_scale_factor = 2.*np.pi/(d_th*np.pi/180.)

# Ny,Nx = np.real(m).shape

# inds  = [(np.arange(Nx)+.5 - Nx/2.) /(Nx-1.)]
# X = np.repeat(inds,Ny,axis=0) * ell_scale_factor

# inds  = (np.arange(Ny)+.5 - Ny/2.) /(Ny-1.)
# Y = np.repeat(inds.reshape(Ny,1),Nx,axis=1) * ell_scale_factor
# # plt.imshow(X)
# # plt.show()
# # plt.imshow(Y)
# # plt.show()
# ell2d = np.sqrt(X**2. + Y**2.)

# def make_CMB_T_map(N,pix_size,ell,DlTT):
#     "makes a realization of a simulated CMB sky map"

#     # convert Dl to Cl
#     ClTT = DlTT * 2 * np.pi / (ell*(ell+1.))
#     ClTT[0] = 0.
#     ClTT[1] = 0.

#     # make a 2d coordinate system
#     ones = np.ones(N)
#     inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
#     X = np.outer(ones,inds)
#     Y = np.transpose(X)
#     R = np.sqrt(X**2. + Y**2.)
    
#     # now make a 2d CMB power spectrum
#     ell_scale_factor = 2. * np.pi / (pix_size/60. * np.pi/180.)
#     ell2d = R * ell_scale_factor
#     ClTT_expanded = np.zeros(ell2d.max()+1)
#     ClTT_expanded[0:(ClTT.size)] = ClTT
#     CLTT2d = ClTT_expanded[ell2d.astype(int)]
#     ## make a plot of the 2d cmb power spectrum, note the x and y axis labels need to be fixed
#     #Plot_CMB_Map(CLTT2d**2. *ell2d * (ell2d+1)/2/np.pi,0,np.max(CLTT2d**2. *ell2d * (ell2d+1)/2/np.pi)/10.,ell2d.max(),ell2d.max())  ###
 
#     # now make a realization of the CMB with the given power spectrum in fourier space
#     ramdomn_array_for_T = np.fft.fft2(np.random.normal(0,1,(N,N)))    
#     FT_2d = np.sqrt(CLTT2d) * ramdomn_array_for_T
#     ## make a plot of the 2d cmb simulated map in fourier space, note the x and y axis labels need to be fixed
#     #Plot_CMB_Map(np.real(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)/2/np.pi),0,np.max(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)/2/np.pi),ell2d.max(),ell2d.max())  ###
#     CMB_T = np.fft.ifft2(np.fft.fftshift(FT_2d)) /(pix_size /60.* np.pi/180.)
#     CMB_T = np.real(CMB_T)

#     ## return the map
#     return(CMB_T)




# ClTT = DlTT * 2 * np.pi / (ell*(ell+1.))
# ClTT[0] = 0.
# ClTT[1] = 0.

# # make a 2d coordinate system
# ones = np.ones(N)
# inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
# X = np.outer(ones,inds)
# Y = np.transpose(X)
# R = np.sqrt(X**2. + Y**2.)

# # now make a 2d CMB power spectrum
# ell_scale_factor = 2. * np.pi / (pix_size/60. * np.pi/180.)
# ell2d = R * ell_scale_factor
# ClTT_expanded = np.zeros(ell2d.max()+1)
# ClTT_expanded[0:(ClTT.size)] = ClTT
# CLTT2d = ClTT_expanded[ell2d.astype(int)]
# ## make a plot of the 2d cmb power spectrum, note the x and y axis labels need to be fixed
# #Plot_CMB_Map(CLTT2d**2. *ell2d * (ell2d+1)/2/np.pi,0,np.max(CLTT2d**2. *ell2d * (ell2d+1)/2/np.pi)/10.,ell2d.max(),ell2d.max())  ###

# # now make a realization of the CMB with the given power spectrum in fourier space
# ramdomn_array_for_T = np.fft.fft2(np.random.normal(0,1,(N,N)))    
# FT_2d = np.sqrt(CLTT2d) * ramdomn_array_for_T
# ## make a plot of the 2d cmb simulated map in fourier space, note the x and y axis labels need to be fixed
# #Plot_CMB_Map(np.real(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)/2/np.pi),0,np.max(np.conj(FT_2d)*FT_2d*ell2d * (ell2d+1)/2/np.pi),ell2d.max(),ell2d.max())  ###
# CMB_T = np.fft.ifft2(np.fft.fftshift(FT_2d)) /(pix_size /60.* np.pi/180.)
# CMB_T = np.real(CMB_T)
