# distutils: language = c
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# -*- coding: utf-8 -*-


import numpy as np
cimport numpy as np
#from libc.math cimport sqrt,pi,log,fmax,floor,ceil
#from libc.stdlib cimport malloc, free
# import scipy.special.wofz as wofz
#cimport scipy.special.cython_special as scs
#from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
import cython
cimport cython
from cython.parallel import prange
from cython.parallel cimport prange

#cdef extern from "complex.h" nogil:
#    double complex CMPLX(double,double)
#    double creal(double complex)
    
cdef extern from "Faddeeva.h" nogil:
    double complex Faddeeva_w(double complex z, double relerr) nogil

ctypedef np.float_t DTYPE_t
ctypedef np.complex128_t DTYPE_t_c
DTYPE = np.float
DTYPE_c = np.complex128


@cython.nonecheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
cdef double complex wofz_single(double complex zz) nogil:
    
    cdef double complex ww
    cdef double rel = 0.0
    
    ww = Faddeeva_w(zz,rel)
        
    return ww

#@cython.nonecheck(False)
#@cython.cdivision(True)
#@cython.profile(True)
#@cython.linetrace(True)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
def wofz(double complex [:] z):
    
    cdef int N = z.shape[0]
    cdef double complex [:] w = np.zeros(N,dtype=DTYPE_c)
    cdef int n
    
    with nogil:
        for n in range(N):
            w[n] = wofz_single(z[n])
        
    return w


#@cython.wraparound(False)
#@cython.nonecheck(False)
#@cython.cdivision(True) 
#@cython.boundscheck(False)
#cdef double VOIGT_DIST(double lnctr, double sgmd, double gm0, double nu):
#    
#    #cdef <double complex*> *z=malloc(N*sizeof(double complex))
#    #cdef int i
#    cdef double y
#    cdef double complex y_c,z
#    cdef double sqrt2 = sqrt(2)
#    cdef double sqrt_pi = sqrt(pi)
#    
#    z = CMPLX(nu-lnctr,gm0)
#    z /= CMPLX(sgmd*sqrt2,0)
#    #y_c = Faddeeva_w(z,0)
#    y_c = scs.wofz(z)
#    y = creal(y_c) 
#    y /= (sgmd*sqrt2*sqrt_pi)
#    
#    #for i in range(N):
#    #    z[i] = CMPLX(nu_grid[i]-lnctr,gm0)
#        
#    return y
#
#@cython.wraparound(False)
#@cython.nonecheck(False)
#@cython.boundscheck(False)
#def VOIGT_PROFILE(double lnctr, double sgmd, double gm0, np.ndarray[DTYPE_t,ndim=1] nu_grid):
# 
#    cdef Py_ssize_t N = nu_grid.shape[0]
#    cdef Py_ssize_t n
#    cdef double ytmp
#    cdef double nu 
#    cdef np.ndarray[DTYPE_t,ndim=1] y = np.zeros(N, dtype=DTYPE)
#    
#    for n in range(N):
#        nu = nu_grid[n]
#        #z = CMPLX(nu-lnctr,gm0)
#        #z /= CMPLX(sgmd*sqrt2,0)
#        ytmp = VOIGT_DIST(lnctr, sgmd, gm0, nu)
#        y[n] = ytmp
#
#    return y

#@cython.wraparound(False)
#@cython.boundscheck(False)
#@cython.nonecheck(False)
#cdef void VOIGT_PROFILE_malloc(double lnctr, double sgmd, double gm0, double[:] nu_grid, double* ym, int N):
# 
#    # cdef int N = nu_grid.shape[0]
#    cdef int n
#    #cdef double ytmp
#    # cdef double nu
#    # cdef double* ym = <double*>malloc(N*sizeof(double))
#    # cdef np.ndarray[DTYPE_t, ndim=1] y = np.zeros(N, dtype=DTYPE)
#    
#    for n in range(N):
#        ym[n] = VOIGT_DIST(lnctr, sgmd, gm0, nu_grid[n])
#
#    return;
#
##@cython.profile(True)
##@cython.linetrace(True)
#@cython.wraparound(False)
#@cython.boundscheck(False)
#@cython.nonecheck(False)
#@cython.cdivision(True) 
#def XSEC_VOIGT_PROFILE(np.ndarray[DTYPE_t,ndim=1] SW, np.ndarray[DTYPE_t,ndim=1] lnctrs, np.ndarray[DTYPE_t,ndim=1] gmds, np.ndarray[DTYPE_t,ndim=1] gm0s, np.ndarray[DTYPE_t,ndim=1] nu_grid, double nu_step, double wing, double wingHW):
#    
#    cdef int N = nu_grid.shape[0]
#    cdef int Nl = lnctrs.shape[0]
#    cdef int ir,il
#    cdef double ird,ild
#    cdef double nu
#    cdef int nl,nn
#    cdef double z_real,z_im
#    cdef double y
#    cdef double complex y_c,z
#    cdef double sw,gmd, sgmd,lnctr,gm0,wing_nl,wing_nl_left,wing_nl_right
#    cdef double sqrt_2log2 = sqrt(2*log(2))
#    cdef double sgmd_sqrt2
#    cdef double sqrt2 = sqrt(2)
#    cdef double sqrtpi = sqrt(pi)
#    cdef double sqrt2_sqrtpi = sqrt2*sqrtpi
#    cdef double sgmd_sqrt2_sqrt_pi
#    #cdef double* xsecm = <double*>malloc(N*sizeof(double))
#    cdef double[:] xsec
#    xsec = np.zeros(N, dtype=DTYPE)
#    
#    for nl in range(Nl):
#        gmd = gmds[nl]
#        sgmd = gmd/sqrt_2log2
#        sgmd_sqrt2 = sgmd*sqrt2
#        sgmd_sqrt2_sqrt_pi = sgmd*sqrt2_sqrtpi
#        lnctr = lnctrs[nl]
#        gm0 = gm0s[nl]
#        sw  = SW[nl]
#        
#        # first decide the indices of nu_grid: assume nu_grid is equally spaced
#        wing_nl = fmax(wing,fmax(wingHW*gm0,wingHW*gmd))
#        wing_nl_left = lnctr - wing_nl
#        wing_nl_right = lnctr + wing_nl
#        if wing_nl_left < nu_grid[0]:
#            il = 0
#        else:
#            ild = floor((wing_nl_left-nu_grid[0])/nu_step)
#            il = <int>ild
#            
#        if wing_nl_right > nu_grid[N-1]:
#            # print(wing_nl_right,nu_grid[N-1])
#            ir = N
#        else:
#            ird = ceil((wing_nl_right-nu_grid[0])/nu_step)
#            ir = <int>ird
#        
#        # print(nl)
#        # print(ir-il+1)
#        for nn in range(il,ir):
#            nu = nu_grid[nn]
#            z_real = (nu-lnctr)/sgmd_sqrt2
#            z_im = gm0/sgmd_sqrt2
#            z = CMPLX(z_real, z_im) 
#            y_c = wofz_single(z)
#            y = creal(y_c) / sgmd_sqrt2_sqrt_pi
#            xsec[nn] = sw*y
#        #xsecm[nn] += sw*VOIGT_DIST(lnctr,sgmd,gm0,nu)
#    
#        #xsec = np.asarray(xsecm)
#        #free(xsecm)
#    return xsec