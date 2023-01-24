import numpy as np
import time
import os
import pickle

from fluctreconpy import rescale_spat_pos
from fluctreconpy.test import test_fluct_decomp
from fluctreconpy import FitGaussian

import matplotlib.pyplot as plt

def test_decomp_all(nwin_x=32,
                    nwin_y=8,
                    stationary=True,
                    electronic_noise=0.001,
                    time_win=100e-6, #100 samples
                    sampling_time=0.5e-6,
                    recalc=False,
                    test=False,
                    plot=False,
                    iter_method='scipy_bisect',
                    calculate_params=True,
                    transposed=False,
                    decompose=True
                    ):

    noise_level=np.arange(10.)/((10.*2)) #0-50%
    spatial_resolution=np.asarray([10]) #Only the relative resolution counts
    blob_size=(np.arange(10.)+1.)*0.5 #times the spatial resolution
    nwin_time=int(time_win/sampling_time)

    if test:
        #blob_size=[1.5,2,2.5]
        blob_size=[2.]
        noise_level=[0.1]
        #noise_level=[0.05,0.1,0.15,0.2,0.3]
        nwin_time=5
        time_win=nwin_time*sampling_time


    pos_orig=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time,2])
    pos_calc=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time,2])
    fwhm_orig=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time,2])
    fwhm_calc=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time,2])
    density_calc=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time])
    density_samp=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time])
    density_orig=np.zeros([len(noise_level),len(blob_size),len(spatial_resolution),nwin_time])

    filename_all=f'tmp/test_fluct_decomp_all_BS_{blob_size[0]:.1f}_{blob_size[-1]:.1f}_NL_{noise_level[0]:.1f}_{noise_level[-1]:.1f}_SR_{spatial_resolution[0]:.1f}_{spatial_resolution[-1]:.1f}'
    if test:
        filename_all+='_test'
    if not decompose:
        filename_all+='_nodecomp'
    filename_all+='.pickle'

    n_noise=len(noise_level)
    n_size=len(blob_size)
    n_res=len(spatial_resolution)
    plt.figure()
    for i_noise in range(0, len(noise_level)):
        for i_size in range(0, len(blob_size)):
            for i_res in range(0, len(spatial_resolution)):
                start_time=time.time()
                spatial_pos=rescale_spat_pos(nwin_x=nwin_x,
                                             nwin_y=nwin_y,
                                             rad_res=spatial_resolution[i_res],
                                             vert_res=spatial_resolution[i_res],
                                             transposed=transposed)
                if not stationary:
                    moving={'vx':spatial_resolution[i_res]/1e3/time_win,
                            'vy':spatial_resolution[i_res]/1e3/time_win,
                            'v_size':0.,
                            'rot_freq':0.} #it only moves one spatial resolution step during the analysis
                else:
                    moving={'vx':0,'vy':0, 'v_size':0., 'rot_freq':0.}

                blob_pos=[np.mean(spatial_pos[:,:,0]),
                          np.mean(spatial_pos[:,:,1])]

                filename=f'tmp/test_fluct_decomp_N_{noise_level[i_noise]:.1f}_BS_{blob_size[i_size]:.1f}_R_{spatial_resolution[i_res]:.1f}'
                if test:
                    filename+='_test'
                if not decompose:
                    filename+='_nodecomp'
                filename+='.pickle'


                results=test_fluct_decomp(nocalib=True,
                                          #k_factor=k_factor,
                                          test=test,
                                          recalc=recalc,
                                          noise_level=noise_level[i_noise],
                                          electronic_noise=electronic_noise,
                                          blob_size=[blob_size[i_size]*spatial_resolution[i_res],
                                                     blob_size[i_size]*spatial_resolution[i_res]],
                                          blob_pos=blob_pos,
                                          moving=moving,
                                          plot=plot,
                                          time_win=time_win,
                                          spatial_pos=spatial_pos,
                                          sampling_time=sampling_time,
                                          floating=False,
                                          iter_method=iter_method,
                                          save_filename=filename,
                                          transposed=transposed,
                                          decompose=decompose
                                          )

                if calculate_params:
                    blob_density_samp=results['blob_density_samp']
                    blob_density_orig=results['blob_density_orig']
                    n_vector_calc=results['n_vector_calc']
                    #n_vector_orig=results['n_vector_orig']
                    time_vec=results['time_vec']

                    density_samp[i_noise,i_size,i_res,:]=blob_density_samp
                    density_orig[i_noise,i_size,i_res,:]=blob_density_orig

                    for k in range(nwin_time):
                        #contour, n_vector_calc[:,:,k], /fill, nlevel=21, /iso, title='N_'+strtrim(noise_level[i_noise],2)+'_BS_'+strtrim(blob_size[i_size],2)+'_R_'+strtrim(spatial_resolution[i_res],2)
                        pos_orig[i_noise,i_size,i_res,k,0]=blob_pos[0]+moving['vx'] * time_vec[k]
                        pos_calc[i_noise,i_size,i_res,k,0]=np.sum(n_vector_calc[:,:,k]*spatial_pos[:,:,0])/np.sum(n_vector_calc[:,:,k])

                        pos_orig[i_noise,i_size,i_res,k,1]=blob_pos[1]+moving['vy'] * time_vec[k]
                        pos_calc[i_noise,i_size,i_res,k,1]=np.sum(n_vector_calc[:,:,k]*spatial_pos[:,:,1])/np.sum(n_vector_calc[:,:,k])

                        fwhm_orig[i_noise,i_size,i_res,k,0]=blob_size[i_size]*spatial_resolution[i_res]
                        fwhm_orig[i_noise,i_size,i_res,k,1]=blob_size[i_size]*spatial_resolution[i_res]

                        # gauss_x=FitGaussian(x=np.mean(spatial_pos[:,:,0],axis=1),
                        #                     data=np.mean(n_vector_calc[:,:,k],axis=1))
                        # xsize=gauss_x.size
                        # fwhm_calc[i_noise,i_size,i_res,k,0]=xsize

                        # gauss_y=FitGaussian(x=np.mean(spatial_pos[:,:,1],axis=0),
                        #                     data=np.mean(n_vector_calc[:,:,k],axis=0))
                        # ysize=gauss_y.size
                        # fwhm_calc[i_noise,i_size,i_res,k,1]=ysize
                        try:
                            y_mean_pos=np.mean(spatial_pos[:,:,1],axis=0)
                            y_mean_data=np.mean(n_vector_calc[:,:,k],axis=0)
                            ind_over_half=np.where(y_mean_data > (np.max(y_mean_data)-np.min(y_mean_data))/2)

                            ysize=(np.max(y_mean_pos[ind_over_half])-
                                   np.min(y_mean_pos[ind_over_half]))
                            fwhm_calc[i_noise,i_size,i_res,k,1]=ysize

                            x_mean_pos=np.mean(spatial_pos[:,:,1],axis=0)
                            x_mean_data=np.mean(n_vector_calc[:,:,k],axis=0)
                            ind_over_half=np.where(x_mean_data > (np.max(x_mean_data)-np.min(x_mean_data))/2)
                            xsize=(np.max(x_mean_pos[ind_over_half])-
                                   np.min(x_mean_pos[ind_over_half]))
                            fwhm_calc[i_noise,i_size,i_res,k,1]=xsize
                        except Exception as e:
                            print('Raised exception: '+str(e))
                            pass


                        density_calc[i_noise,i_size,i_res,k]=np.sum(n_vector_calc[:,:,k])/(np.pi*xsize*ysize)

                finish_time=time.time()

                rem_time=(finish_time-start_time)*(n_noise * n_size * n_res - i_noise*n_size*n_res - i_size*n_res - i_res)
                print(f"One calculation took {finish_time-start_time} seconds.")
                print(f"Remaining time from the calculation: {rem_time//60} min.")

    results_all={
                'noise_level':noise_level,
                'blob_size':blob_size,
                'spatial_resolution':spatial_resolution,
                'nwin_time':nwin_time,
                'pos_orig':pos_orig,
                'pos_calc':pos_calc,
                'fwhm_orig':fwhm_orig,
                'fwhm_calc':fwhm_calc,

                'density_orig':density_orig,
                'density_samp':density_samp,
                'density_calc':density_calc,
                 }

    with open(filename_all, 'wb') as f: pickle.dump(results_all, f)
    return results_all