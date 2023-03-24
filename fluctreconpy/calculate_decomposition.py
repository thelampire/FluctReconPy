#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 01:27:56 2022

@author: mlampert
"""

import numpy as np
import copy

from fluctreconpy import undulation_matrix

def calculate_decomposition(light_profile=None,
                            error_profile=None,
                            r_vec=None,
                            z_vec=None,

                            iter_method='bisect',
                            threshold = 0.001,
                            maxiter=100,
                            increment=1.25,
                            k_factor=[-1e-10,-1e3],

                            n_vector=None,
                            m_matrix=None,

                            floating=True,
                            visual=True,
                            contour=True,
                            pdf=False,

                            test=False,
                            transposed=False
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

    s_vector_ref=np.reshape(light_profile,
                            m_matrix.shape[2] * m_matrix.shape[3],
                            )
    error_profile_s_ref=np.reshape(error_profile,
                                       m_matrix.shape[2] * m_matrix.shape[3]
                                       )

    m_matrix_ref = np.reshape(m_matrix,
                              (m_matrix.shape[0] * m_matrix.shape[1],
                               m_matrix.shape[2] * m_matrix.shape[3]),
                              )
    h_matrix = undulation_matrix(r_vec=r_vec,
                                 z_vec=z_vec,
                                 floating=floating,
                                 transposed=transposed,
                                 ) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]

    h_matrix_ref = np.reshape(h_matrix,
                              (h_matrix.shape[0]*h_matrix.shape[1],
                               h_matrix.shape[2]*h_matrix.shape[3]),
                               )
    m_matrix_prime_ref=copy.deepcopy(m_matrix_ref)
    for i in range(m_matrix.shape[0]*m_matrix.shape[1]-1):
        m_matrix_prime_ref[:,i]=m_matrix_ref[:,i]/error_profile_s_ref[:]


    def function_to_iterate(k_factor_curr,
                            s_vector_ref, m_matrix_ref, m_matrix_prime_ref,
                            h_matrix_ref, error_profile_s_ref, m_matrix_shape,
                            return_p_vector):

        p_vector = k_factor_curr * np.matmul(np.matmul(np.linalg.inv(h_matrix_ref.T +
                                                                     k_factor_curr * np.matmul(m_matrix_prime_ref.T,
                                                                                               m_matrix_ref)),
                                                        m_matrix_prime_ref.T,),
                                              s_vector_ref
                                              )
        if not return_p_vector:
            return np.sum((s_vector_ref -
                           np.matmul(m_matrix_ref,p_vector))**2/error_profile_s_ref[:])/m_matrix_shape[2]/m_matrix_shape[3]-1
        else:
            return p_vector

    args=(s_vector_ref, m_matrix_ref, m_matrix_prime_ref,
          h_matrix_ref, error_profile_s_ref, m_matrix.shape,False)

    if iter_method == 'scipy_bisect' :                                             #~34s by far the fastest even though it should be the slowest.
        from scipy.optimize import bisect as root_finder
    elif iter_method == 'scipy_brenth':                                           #~48s
        from scipy.optimize import brenth as root_finder
    elif iter_method == 'scipy_brentq':                                           #~40s
        from scipy.optimize import brentq as root_finder
    elif iter_method == 'scipy_ridder':                                           #~70s
        from scipy.optimize import ridder as root_finder
    elif iter_method == 'scipy_toms748':                                          #Really slow, not recommended
        from scipy.optimize import toms748 as root_finder
    elif iter_method == 'bisect':                                                 #~40s
        root_finder=bisect_root_finder
    else:
        raise ValueError('iter_method is not in the list of methods: bisect, scipy_[bisect, brenth, brentq, ridder, toms748]')
    k_factor_opt=root_finder(function_to_iterate, k_factor[0], k_factor[1],
                             args=args, maxiter=maxiter, xtol=threshold)

    args=(s_vector_ref, m_matrix_ref, m_matrix_prime_ref,
          h_matrix_ref, error_profile_s_ref, m_matrix.shape, True)

    p_vector=function_to_iterate(k_factor_opt, *args)


    results={'p_vector':np.reshape(p_vector,(m_matrix.shape[0],m_matrix.shape[1])),
             'k_factor':k_factor,
             'k_factor_opt':k_factor_opt,
             }
    return results

def bisect_root_finder(function_to_iterate, a, b, args=None, maxiter=1000, xtol=0.001):
    ind_iter=1
    k_factor_curr=a
    while True:

        result=function_to_iterate(k_factor_curr, *args)
        if result < -1:
            print('Ksi is lower than 0, iteration is aborted.')
            break

        k_factor_vector=copy.deepcopy([a,b])
        result_vector=np.zeros(2)
        if ind_iter == 1:
            result_vector[0]=result
            k_factor_curr=k_factor_vector[1]

        if ind_iter == 2:
            result_vector[1]=result
            k_factor_curr=(k_factor_vector[0]+k_factor_vector[1])/2.

        if ind_iter > 3:
            if result < 0:
                result_vector[1]=result
                k_factor_vector[1]=k_factor_curr
            else:
                result_vector[0]=result
                k_factor_vector[0]=k_factor_curr

        k_factor_curr=(k_factor_vector[0]+k_factor_vector[1])/2.

        ind_iter += 1

        if ind_iter == maxiter+1:
            break
        if (np.abs(result) < xtol):
            break
    return k_factor_curr