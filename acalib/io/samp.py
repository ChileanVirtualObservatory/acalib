from astropy.vo.samp import SAMPIntegratedClient
from urlparse import urlparse, urljoin


#TODO Make a SAMP protocol interface to core components
def send_aladin(fit,name):

    """
    Sends data to aladin
    fits only supported
    """
    #fitPath = os.path.splitext(fit)[0]
    fitExtension = os.path.splitext(fit)[1]

    if (fitExtension != ".fits"): # for now, next with exceptions
        print "Please use a file with .fits extension"
        return

    client = SAMPIntegratedClient()
    client.connect()

    params = {}
    params["url"] = 'file://' + fit
    params["name"] = name

    message = {}
    message["samp.mtype"] = "image.load.fits"
    message["samp.params"] = params

    client.notify_all(message) # the easiest way is notify all

    client.disconnect()

