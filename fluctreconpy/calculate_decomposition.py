#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 01:27:56 2022

@author: mlampert
"""

import numpy as np
import pickle
import copy
import matplotlib.pyplot as plt

from fluctreconpy import undulation_matrix, reform_matrix

def calculate_decomposition(light_profile=None,
                            error_profile=None,
                            spatial_pos=None,

                            threshold = 0.001,
                            maxiter=100,
                            increment=1.25,
                            k_factor=[1e-10,1e3],

                            n_vector=None,
                            m_matrix=None,

                            floating=True,
                            iterate=True,
                            visual=True,
                            contour=True,
                            pdf=False,
                            nocalc=False,
                            save_filename='calculate_decomposition_save.pickle',

                            test=False,
                            ):

    """
    ;*******************************************************
    ;***   calculate_decomposition       by M. Lampert   ***
    ;*******************************************************
    ;
    ; This procedure calculates the density matrix from
    ; the measured BES light signal based on the arbitrary
    ; expansion method. The steps followed in the procedure
    ; are the following:
    ; 1: Arbitrary initial K value is considered.
    ; 2: Matrix N and L is calculated
    ; 3: Minimum point of X(P) is calculated
    ; 4: X is calculated
    ; 5: If X is not close enough to unity, then go back to
    ;    the second step
    ;
    ; INPUTs:
    ;   shot: shotnumber to be analyzed
    ;   time: time to be analyzed
    ;   t_win: time window around time, timerange: [time-t_win/2,time+t_win/2]
    ;   nocalibrate: No relative amplitude calibration
    ;   test: Test the code with dummy atomic physics
    ;   iterate: solve the equation with iteration
    ;   floating: the perimeter values of the undulation matrix are floating
    ;   For testing:
    ;   noise_level: the relative noise level of the signal
    ;   electronic_noise: electronic noise in mV for the calculated spectral range
    ;   blog_size: size of the simulated blob in mm
    ;   blob_pos: position of the blob in mm [R,z]
    ;   /visual: visualize the calculation
    ;   /contour: contour plot instead of simple plotting
    ; OUTPUTs:
    ;   return: calculated density fluctuation matrix
    ;   k_factor: smoothing factor for further calculation
    ;   n_vector: the original density distribution
    ;
    ;
    """

    if not nocalc:
        if iterate:
            k_factor_vector=copy.deepcopy(k_factor)
            ksi_vector=np.zeros(2)
        else:
            ksi_vector=np.zeros(maxiter)
            k_factor_vector=np.zeros(maxiter)

        size_m_mat=m_matrix.shape

        n_vert_calc=size_m_mat[0]
        n_rad_calc=size_m_mat[1]
        n_vert_meas=size_m_mat[2]
        n_rad_meas=size_m_mat[3]

        s_vector_ref=np.reshape(light_profile,
                                n_vert_meas * n_rad_meas,
                                )
        error_profile_s_ref=np.reshape(error_profile,
                                       n_vert_meas * n_rad_meas,
                                       )

        m_matrix_ref = reform_matrix(m_matrix,
                                     [n_vert_calc * n_rad_calc,n_vert_meas * n_rad_meas])

        h_matrix = undulation_matrix(spatial_pos=spatial_pos,
                                     floating=floating) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]

        h_matrix_ref = reform_matrix(h_matrix,
                                     [n_vert_calc*n_rad_calc, n_vert_meas*n_rad_meas])

        ind_iter=1
        k_factor_curr=k_factor[0]

        while True:

            m_matrix_prime_ref=copy.deepcopy(m_matrix_ref)
            for i in range(n_rad_calc*n_vert_calc-1):
                m_matrix_prime_ref[:,i]=m_matrix_ref[:,i]/error_profile_s_ref[:]

            p_vector = k_factor_curr * np.matmul(np.matmul(np.linalg.inv(h_matrix_ref.T +
                                                                         k_factor_curr * np.matmul(m_matrix_prime_ref.T,
                                                                                                   m_matrix_ref)),
                                                            m_matrix_prime_ref.T,),
                                                  s_vector_ref
                                                  )
            #Calculate ksi
            ksi=n_vert_meas**(-1)*n_rad_meas**(-1)*np.sum((s_vector_ref -
                                                          np.matmul(m_matrix_ref,p_vector))**2/error_profile_s_ref[:])

            if test:
                plt.contourf(np.reshape(p_vector, (n_vert_meas,n_rad_meas)))
                plt.pause(0.2)
            #     print(f"p_vector: {p_vector}")
            #     print(f"ksi: {ksi}")
                print(k_factor_curr)
            if ksi < 0:
                print('Ksi is lower than 0, iteration is aborted.')
                print(p_vector)
                break

            if iterate:
                if ind_iter == 1:
                    ksi_vector[0]=ksi
                    k_factor_curr=k_factor_vector[1]

                if ind_iter == 2:
                    ksi_vector[1]=ksi
                    k_factor_curr=(k_factor_vector[0]+k_factor_vector[1])/2.

                if ind_iter > 3:
                    if ksi < 1:
                        ksi_vector[1]=ksi
                        k_factor_vector[1]=k_factor_curr
                    else:
                        ksi_vector[0]=ksi
                        k_factor_vector[0]=k_factor_curr

                k_factor_curr=(k_factor_vector[0]+k_factor_vector[1])/2.

            else:
                k_factor_vector[ind_iter-1]=k_factor_curr
                k_factor_curr *= increment
                ksi_vector[ind_iter-1]=ksi

            ind_iter += 1

            if ind_iter == maxiter+1:
                break
            if (np.abs(ksi - 1) < threshold):
                break
        p_vector = np.reshape(p_vector,(n_vert_calc,n_rad_calc))

        results={'p_vector':p_vector,
                 'ksi_vector':ksi_vector,
                 'k_factor':k_factor,
                 }
        with open(save_filename,'wb') as f: pickle.dump(results,f)
    else:
        try:
            with open(save_filename,'rb') as f: results=pickle.load(f)
        except Exception as e:
            print('Raised exception: '+e)
            raise ValueError('The save file doesnt exist, please run the procedure without /nocalc')

    return results