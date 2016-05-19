from  astropy.coordinates.name_resolve import get_icrs_coordinates
import urllib2
import os

for i in xrange(1,111):

	not_downloaded = []
        fnd = open('NOTDOWNLOADED', 'a')
        dss2b = os.path.exists("./M"+str(i)+'-dss2b.fits')
        dss2r = os.path.exists("./M"+str(i)+'-dss2r.fits')
        dss2ir = os.path.exists("./M"+str(i)+'-dss2ir.fits')
        if not (dss2b and dss2r and dss2ir):
		try:
			radec = get_icrs_coordinates("M"+str(i))
			ra = radec.ra.value
			dec = radec.dec.value
			print "Downloading M"+str(i)
			for survey in ['dss2b','dss2r','dss2ir']:
				url = 'http://skyview.gsfc.nasa.gov/cgi-bin/images?position='+str(ra)+'%2C'+str(dec)+'&survey='+survey+'&pixels=30%2C30&sampler=Clip&size=1.0%2C1.0&projection=Tan&coordinates=J2000.0&requestID=skv1434491073066&return=FITS'
				response = urllib2.urlopen(url)
				content = response.read()
				f = open('M'+str(i)+'-'+survey+'.fits' , 'wb' )
				f.write( content )
				f.close()
			print "M"+str(i)+" is ready"
		except:
			print "M"+str(i)+" not Downloaded"
			not_downloaded.append(i)

print "Not downloaded: " + " ".join(not_downloaded) 
