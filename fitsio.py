# test script for reading in fits random group (i.e., vis) files
try:
    from astropy.io import fits
    import pyfits
except: 
    print("No astropy on this system. Check python path / installation.")
    exit(-1)


import pprint 
import numpy as np
import sys
import os
try:
    import psutil     
    HAVE_PSUTIL = True
except: 
    print("No psutil on this system. Memory check will not work.")

CC   = 299792458.0
HIGH = 5 
debug = HIGH

def get_alma_info(file):
    hdul = fits.open(file)

#    u = hdul[0].data['UU']
#    v = hdul[0].data['VV']

    if(debug == HIGH):
        print("Printing FITS info for debug: ")
        pprint.pprint(hdul)
    data   = hdul[0].data['DATA']
    pheader= hdul[0].header
    sheader= hdul[1].data

    ncorr  = np.shape(data)[5]
    nspw   = np.shape(data)[3]
    nchan  = np.shape(data)[4]
    nint   = np.shape(data)[0]
    
    info={}

    info['nints']  = nint
    info['ncorr']  = ncorr
    info['nchan']  = nchan
    info['nspw']   = nspw

    # grab frequency info
    if(pheader['ctype4'].rstrip() == 'FREQ'):
        info['reffreq'] = pheader['crval4']
        info['refchan'] = pheader['crpix4']
        info['ifoffsets']=[]
        info['chwidth']=[]
        for i in range(nspw):
            info['ifoffsets'].append(sheader['IF FREQ'][0][i])
            info['chwidth'].append(sheader['CH WIDTH'][0][i])
    if(debug == HIGH):
        pprint.pprint(info)    
    hdul.close()
    return info

def load_data(file):
    info = get_alma_info(file)
    # now grab the u / v and freq info 
    hdul = fits.open(file)
    u = hdul[0].data['UU']
    v = hdul[0].data['VV']
    data = hdul[0].data['DATA']
    if(debug == HIGH and HAVE_PSUTIL):
        BUFFER = 2 # something low
        mem = psutil.virtual_memory()
        print 
#        assert(mem.available > BUFFER*(info['nspw']*info['nchan']*np.size(u))*sys.getsizeof(u[0])),"Insufficient memory. Is the channel range you selected appropriate?"

#    u= np.repeat(u,info['nspw']*info['nchan'])
#    v= np.repeat(v,info['nspw']*info['nchan'])
    chan_freq = np.array([])
    for i in range(info['nspw']):
        chan_freq = np.append(chan_freq,np.repeat(info['ifoffsets'][i],info['nchan']) + (np.arange(info['nchan'])-info['refchan'])*info['chwidth'][i]+info['reffreq'])
        
    u_dimless  = np.array([])
    v_dimless  = np.array([])
    vis_freq   = np.array([])
    vis_re     = np.empty(shape=(0,info['ncorr']))
    vis_im     = np.empty(shape=(0,info['ncorr']))
    vis_wgt    = np.empty(shape=(0,info['ncorr']))

# maybe do all the appending in list form and then convert to numpy later to save valuable time?
# eg http://stackoverflow.com/questions/22392497/how-to-add-a-new-row-to-an-empty-numpy-array

    for i in xrange(len(u)):
        print (i/(len(u)+0.0))
        u_dimless = np.append(u_dimless,chan_freq*u[i])
        v_dimless = np.append(v_dimless,chan_freq*v[i])
        vis_freq= np.append(vis_freq,chan_freq)
        # so this is for the data that is sorted in chan / spw order.
        # blaggggg can't use np.reshape() since there is no guarantee of the memory layout 
        # of the returned array. mother fuck.
        for j in xrange(info['nspw']):
            for k in xrange(info['nchan']):
                vis_re  = np.append(vis_re ,np.reshape(data[i,0,0,j,k,0:info['ncorr'],0],(1,info['ncorr'])),axis=0)
                vis_im  = np.append(vis_im ,np.reshape(data[i,0,0,j,k,0:info['ncorr'],1],(1,info['ncorr'])),axis=0)
                vis_wgt = np.append(vis_wgt,np.reshape(data[i,0,0,j,k,0:info['ncorr'],2],(1,info['ncorr'])),axis=0)

    
                # clean up
    del data
    hdul.close()
    data = {}
    data['u'] = u_dimless
    data['v'] = v_dimless
    data['re'] = vis_re
    data['im'] = vis_im
    data['wgt'] = vis_wgt
    data['freq'] = vis_freq
    #data['iffreq'] = [for i in range(info['nspw']): np.replicate(info['ifoffset']+info['reffreq'],info['n
    print(np.shape(data['u']))
    print(np.shape(data['v']))
    print(np.shape(data['re']))
    print(np.shape(data['im']))
    print(np.shape(data['wgt']))
    print(np.shape(data['freq']))
    return data
    
    
def make_fits_model(filename,array,delta,rpix):
    os.system('rm ' + filename)
    hdu = pyfits.PrimaryHDU(array)
    hdulist = pyfits.HDUList([hdu])

#   setup the output fits file coordinate system. 
    hdulist[0].header.append(('cdelt1', delta/3600.))
    hdulist[0].header.append(('cdelt2', delta/3600.))
    hdulist[0].header.append(('ctype1', 'RA---TAN'))
    hdulist[0].header.append(('ctype2', 'DEC--TAN'))
    hdulist[0].header.append(('crpix1', rpix))
    hdulist[0].header.append(('crpix2', rpix))
    hdulist[0].header.append(('crval1', 0))
    hdulist[0].header.append(('crval2', 90.0))
    hdulist.writeto(filename)
    return 

def main():
    file = './data/IRAS16293_Band9.fixed.rebin.rest.uvfits'
    blah = load_data(file)

    return blah

