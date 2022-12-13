#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 19:11:45 2022

@author: mlampert
"""
import numpy as np
import scipy

class FitGaussian:
    def __init__(self,
                 x=None,
                 y=None,
                 data=None):


        self._fwhm_to_sigma=(2*np.sqrt(2*np.log(2)))
        self.x=x
        self.y=y
        self.data=data
        self.fit_gaussian(x,y,data)


    def fit_gaussian(self,x,y,data):
        xdata=np.vstack((x.ravel(),y.ravel()))
        initial_guess=[data.max(),                                #Amplitude
                       np.sum(x*data)/np.sum(data),               #x0
                       np.sum(y*data)/np.sum(data),               #y0
                       (x.max()-x.min())/2/self._fwhm_to_sigma,   #Sigma_x
                       (y.max()-y.min())/2/self._fwhm_to_sigma,   #Sigma_y
                       0.,                                        #Angle
                       np.mean(data)                              #Offset
                       ]

        try:
            popt, pcov = scipy.optimize.curve_fit(gaussian2D_fit_function,
                                                  xdata,
                                                  data,
                                                  p0=initial_guess)
            popt[5]=np.arcsin(np.sin(popt[5]))
            self.popt=popt
        except:
            self.popt=np.zeros(7)
            self.popt[:]=np.nan
            print('Gaussian fitting failed.')


        theta=self.popt[5]
        a,b=np.abs(np.asarray([self.popt[3],
                               self.popt[4]])*self._fwhm_to_sigma)

        if a < b:
            self._axes_length=np.asarray([a,b])
            self._angle=np.arcsin(np.sin(theta))
        else:
            self._axes_length=np.asarray([b,a])
            self._angle=np.arcsin(np.sin(theta-np.pi/2))

        self._center=np.asarray([self.popt[1],
                                 self.popt[2]])

    def set_invalid(self):
        self.popt=np.zeros(7)
        self.popt[:]=np.nan

    @property
    def angle(self):
        return self._angle

    @property
    def axes_length(self):
        return self._axes_length

    @property
    def center(self):
        return self._center

    @property
    def center_of_gravity(self):
        return np.asarray([np.sum(self.x*self.data)/np.sum(self.data),
                           np.sum(self.y*self.data)/np.sum(self.data)])

    @property
    def elongation(self):
        size=self.size
        return (size[0]-size[1])/(size[0]+size[1])

    @property
    def half_level(self):
        #WARNING: NOT PRESENT IN ELLIPSE
        return (self.popt[0]-self.popt[6])/2

    @property
    def size(self):

        alfa=self.angle
        a=self.axes_length[0]
        b=self.axes_length[1]

        a0=np.cos(alfa)**2/a**2 + np.sin(alfa)**2/b**2 #*x**2
        a1=2*np.cos(alfa)*np.sin(alfa)*(1/a**2-1/b**2) #*xy
        a2=np.sin(alfa)**2/a**2 + np.cos(alfa)**2/b**2 #*y**2

        xsize=2/np.sqrt(a0)
        ysize=2/np.sqrt(a2)

        if np.imag(xsize) != 0 or np.imag(xsize) !=0:
            raise ValueError('Size is complex')

        return np.array([xsize,ysize])

def gaussian2D_fit_function(coords, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):

    x,y=coords
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = (np.sin(2*theta))/(4*sigma_x**2) - (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (+ a*(x-xo)**2
                                      + 2*b*(x-xo)*(y-yo)
                                      + c*(y-yo)**2)
                                  )
    return g.ravel()