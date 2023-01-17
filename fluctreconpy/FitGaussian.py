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
                 data=None,):


        self._fwhm_to_sigma=(2*np.sqrt(2*np.log(2)))
        self.x=x
        if y is not None:
            self.y=y
        self.data=data
        if y is not None:
            self.dimension=2
            self.fit_gaussian_2d(x,y,data)
        else:
            self.dimension=1
            self.fit_gaussian_1d()


    def fit_gaussian_1d(self):
        initial_guess=[self.data.max(),                                     #Amplitude
                       np.sum(self.x*self.data)/np.sum(self.data),          #x0
                       (self.x.max()-self.x.min())/2/self._fwhm_to_sigma,   #signma
                       np.mean(self.data),                                  #offset
                       ]
        try:
            popt, pcov = scipy.optimize.curve_fit(self.gaussian1d_fit_function,
                                                  self.x,
                                                  self.data,
                                                  p0=initial_guess)
            self.popt=popt
        except Exception as e:
            self.popt=np.zeros(7)
            self.popt[:]=np.nan
            print('Gaussian fitting failed /w exception: '+str(e))

        self._center=self.popt[1]
        self._axes_length=self.popt[2]*self._fwhm_to_sigma

    def fit_gaussian_2d(self,x,y,data):
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
            popt, pcov = scipy.optimize.curve_fit(self.gaussian2d_fit_function,
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
        if self.dimension==2:
            self.popt=np.zeros(7)
        else:
            self.popt=np.zeros(4)
        self.popt[:]=np.nan

    @property
    def angle(self):
        if self.dimension == 2:
            return self._angle
        else:
            raise ValueError('Angle is not defined for 1D Gaussian fits')

    @property
    def axes_length(self):
        return self._axes_length

    @property
    def center(self):
        return self._center

    @property
    def center_of_gravity(self):
        if self.dimension == 2:
            return np.asarray([np.sum(self.x*self.data)/np.sum(self.data),
                               np.sum(self.y*self.data)/np.sum(self.data)])
        else:
            return np.sum(self.x*self.data/np.sum(self.data))

    @property
    def elongation(self):
        if self.dimension == 2:
            size=self.size
            return (size[0]-size[1])/(size[0]+size[1])
        else:
            raise ValueError('Elongation is not defined for 1D fit.')

    @property
    def half_level(self):
        if self.dimension == 2:
            return (self.popt[0]-self.popt[6])/2
        else:
            return (self.popt[0]-self.popt[3])/2

    @property
    def size(self):
        if self.dimension == 2:
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
        else:
            return self.axes_length

    @staticmethod
    def gaussian2d_fit_function(coords, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):

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

    @staticmethod
    def gaussian1d_fit_function(x, amplitude, x0, sigma_x, offset):
        return offset+amplitude*np.exp( - (x-x0)**2/(2*sigma_x**2))