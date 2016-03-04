#cython: cdivision=True 
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

import numpy as np
cimport numpy as cnp


ctypedef unsigned int uint
ctypedef cnp.float64_t float64_t
ctypedef cnp.uint8_t uint8_t
ctypedef cnp.ndarray ndarray


"""
*  Description:
*     This function contracts (erodes) or expands (dilates) the pixels
*     which are marked as "on" in the supplied input array using a cellular
*     automata, and returns the result in an output array of the same shape
*     and size as the input array.
*
*     Each output pixel value is created in turn as follows: If the
*     corresponding input value has value equal to or greater than "magic",
*     then the output pixel is set equal to "magic". Otherwise, a
*     3x3x3 cube (or a 3x3 square for 2D arrays) is defined within the input
*     array which is centred on the position of the current output pixel. The
*     fraction of pixels within this input cube which are flagged as "on" is
*     then counted. If this fraction is larger than "thresh" then the output
*     pixel value is set to "on", otherwise it is set to "off". Optionally,
*     an additional requirement for an output pixel to be set on is that
*     the corresponding input pixel must be on.
"""

cdef ndarray _remove_isolate(float64_t[:,:,::1] inp, float frac, int on, int off, int centre):
   cdef:
       Py_ssize_t ox, oy, oz
       Py_ssize_t ix, iy, iz
       Py_ssize_t nx = inp.shape[0]
       Py_ssize_t ny = inp.shape[1]
       Py_ssize_t nz = inp.shape[2]
       float s, tot
       ndarray[float64_t, ndim=3] out = np.empty_like(inp, dtype='np.float64')

   for ox in range(nx):
       for oy in range(ny):
           for oz in range(nz):
              
               if centre and inp[ox,oy,oz] == off:
                   """
                   If the corresponding input pixel is off, then the output must be also be
                   off if "centre" is true.
                   """
                   out[ox,oy,oz] = off

               else:
                   """
                   Otherwise, loop round all input pixels in the neighbourhood of the current
                   output pixel, this is a cube of 3x3x3 input pixels, centred on the current
                   output pixel. Count how many of these input pixels are set to "on". If
                   the current output pixel is close to an edge of the array, there will be
                   fewer than 3x3x3 pixels in the cube. Count the total number of pixels
                   in the cube.
                   """
                   s = 0.
                   tot = 0.
                   for ix in range(ox-1, ox+2):
                       if ix < 0 or ix >= nx : continue
                       for iy in range(oy-1, oy+2):
                           if iy < 0 or iy >= ny: continue
                           for iz in range(oz-1, oz+2):
                               if iz <0 or iz >= nz: continue
                               if ix==ox and iy==oy and iz==oz: continue 
                               tot += 1
                               if inp[ix,iy,iz] == on: s += 1
                   if s/tot > frac:
                       out[ox,oy,oz] = on
                   else: 
                       out[ox,oy,oz] = off
   return out


# cdef ndarray _smooth_boundary(float64_t[:,:,::1] inp):
#    out = np.empty_like(inp)
#    shape = inp.shape

#    for ox in range(shape[0]):
#        for oy in range(shape[1]):
#            for oz in range(shape[2]):

#                party = np.zeros((26,2))
#                maxvotes = 0
#                reached = False

#                for ix in range(ox-1, ox+2):
#                    if ix < 0 or ix >= shape[0] : continue
#                    for iy in range(oy-1, oy+2):
#                        if iy < 0 or iy >= shape[1]: continue
#                        for iz in range(oz-1, oz+2):
#                            if iz < 0 or iz >= shape[2]: continue

#                            inp_value = inp[ix, iy, iz]
#                            #update party
#                            if inp_value not in party:
#                                party[inp_value] = 1
#                            else:
#                                party[inp_value] += 1

#                            #update winner up to now
#                            if party[inp_value] > maxvotes:
#                                winner = inp_value
#                                maxvotes = party[inp_value]

#                            #verify if winner has been reached
#                            if party[inp_value] > 13:
#                                winner = inp_value
#                                reached = True
#                                continue
#                        if reached: continue
#                    if reached: continue

#                #assign the winner
#                out[ox, oy, oz] = winner
#    return out



"""
Wraper Functions
"""

def remove_isolate(inp, frac, on, off, centre):
   return _remove_isolate(inp, frac, on, off, centre)
