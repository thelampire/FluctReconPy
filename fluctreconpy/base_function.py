#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 01:16:20 2022

@author: mlampert
"""
import numpy as np

def base_function(r=None,
                  z=None,
                  nvec=100,
                  r_interp=None,
                  z_interp=None,
                  ):

    if len(r) != 3 or len(z) != 3:
      raise ValueError('The boundary of the pyramid needs to be defined by 3-3 r z points.')

    if r_interp is None or z_interp is None:
        raise ValueError('r_interp and z_interp need to be set.')


        r_vec=np.arange(nvec)/nvec*(r[2]-r[0])+r[0]
        z_vec=np.arange(nvec)/nvec*(z[2]-z[0])+z[0]
        base=np.zeros([nvec,nvec])
        for i_r in range(nvec):
            for j_z in range (nvec):
                if i_r < nvec/2:
                    if j_z < nvec/2:
                        base[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])*(z_vec[j_z]-z[0])/(z[1]-z[0])
                    else:
                        base[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])*(z_vec[j_z]-z[2])/(z[1]-z[2])

                else:
                    if j_z < nvec/2:
                        base[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])*(z_vec[j_z]-z[0])/(z[1]-z[0])
                    else:
                        base[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])*(z_vec[j_z]-z[2])/(z[1]-z[2])

        return base
    else:
        if r_interp < r[0] or r_interp > r[2] or z_interp < z[0] or z_interp > z[2] :
            print('r_interp and z_interp should be between the given r,z boundaries. Returning -1...')
            return -1.
        if r_interp < r[1] :
            if z_interp < z[1] :
                ret_val=(r_interp-r[0])/(r[1]-r[0])*(z_interp-z[0])/(z[1]-z[0])
            else:
                ret_val=(r_interp-r[0])/(r[1]-r[0])*(z_interp-z[2])/(z[1]-z[2])
        else:
            if z_interp < z[1]:
              ret_val=(r_interp-r[2])/(r[1]-r[2])*(z_interp-z[0])/(z[1]-z[0])
            else:
              ret_val=(r_interp-r[2])/(r[1]-r[2])*(z_interp-z[2])/(z[1]-z[2])
        return ret_val
