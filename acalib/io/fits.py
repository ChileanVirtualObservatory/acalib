from astropy.io import fits
from astropy import log
import numpy as np
import astropy.units as u
from astropy.wcs import wcs
import astropy.nddata as ndd

def HDU_to_NDData(hdu):
   data=hdu.data
   meta=hdu.header
   mask=np.isnan(data)
   # Hack to correct wrong uppercased units generated by CASA
   try:
     bscale=meta['BSCALE']
   except KeyError:
     bscale=1.0
   try:
     bzero=meta['BZERO']
   except KeyError:
     bzero=0.0
   try:
     bsu=meta['BUNIT']
     bsu=bsu.lower()
     bsu=bsu.replace("jy","Jy")
     bunit=u.Unit(bsu,format="fits")
   except KeyError:
     bunit=u.Unit("u.Jy/u.beam")
   rem_list=[]
   for i in meta.items():
      if i[0].startswith('PC00'):
         rem_list.append(i[0])
   for e in rem_list:
      meta.remove(e)

   mywcs=wcs.WCS(meta)
   # Create astropy units
   if len(data.shape) == 4:
       # Put data in physically-meaninful values, and remove stokes
       # TODO: Stokes is removed by summing (is this correct? maybe is averaging?)
       log.info("4D data detected: assuming RA-DEC-FREQ-STOKES (like CASA-generated ones), and dropping STOKES")
       data=data.sum(axis=0)*bscale+bzero
       mywcs=mywcs.dropaxis(3)
   elif len(data.shape) == 3:
       log.info("3D data detected: assuming RA-DEC-FREQ")
       data=data*bscale+bzero
   else:
       log.error("Only 3D data allowed (or 4D in case of polarization)")
       raise TypeError
   return ndd.NDData(data, uncertainty=None, mask=mask,wcs=wcs, meta=meta, unit=unit)

def HDU_to_Table(hdu):
   log.warning("FITS Table ---> AstroPy Table not implemented Yet")
   #return atable.ATable(data=hdu.data,meta=hdu.header)


def save_fits_from_cont(filepath,acont):
   if acont.primary == None:
      phdu=fits.PrimaryHDU()
   else:
      phdu=acont.primary.get_hdu(True)
   nlist=[phdu]
   count=0
   for elm in acont.adata:
       count+=1
       hdu=elm.get_hdu()
       hdu.header['EXTNAME'] = 'SCI'
       hdu.header['EXTVER'] = count
       nlist.append(hdu)
   count=0
   for elm in acont.atable:
       count+=1
       hdu=elm.get_hdu()
       hdu.header['EXTNAME'] = 'TAB'
       hdu.header['EXTVER'] = count
       nlist.append(hdu)
   hdulist = fits.HDUList(nlist)
   hdulist.writeto(filepath,clobber=True)

def load_fits_to_cont(filePath,name,acont):
   hdulist = fits.open(filePath)
   for counter,hdu in enumerate(hdulist):
           if isinstance(hdu,fits.PrimaryHDU) or isinstance(hdu,fits.ImageHDU):
                   log.info("Processing HDU "+str(counter)+" (Image)")
                   try:
                           ndd=HDU_to_adata(hdu)
                           if isinstance(hdu,fits.PrimaryHDU):
                                   acont.primary = ndd
                           acont.adata.append(ndd)
                   except TypeError:
                           log.info(str(counter)+" (Image) wasn't an Image")

           if isinstance(hdu, fits.BinTableHDU):
                   table = HDU_to_atable(hdu)
                   acont.atable.append(table)
   if acont.primary is None:
           if len(acont.adata)==0:
                   acont.primary = acont.atable[0]
           else:
                   acont.primary = acont.adata[0]
