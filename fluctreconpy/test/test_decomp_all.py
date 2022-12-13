import numpy as np
import time
import os
import pickle

from fluctreconpy import rescale_spat_pos
from fluctreconpy.test import test_fluct_decomp
from fluctreconpy import FitGaussian

def test_decomp_all(nwin=16,
                    stationary=True,
                    electronic_noise=0.001,
                    time_win=100e-6, #100 samples
                    sampling_time=0.5e-6,
                    nocalc=False,
                    test=True,
                    plot=False,
                    ):

    n_try=10
    noise_vector=np.arange(n_try)/((n_try*2)) #0-50%
    spatial_resolution=np.asarray([10]) #Only the relative resolution counts
    blob_size=(np.arange(n_try)+1.)*0.5 #times the spatial resolution

    nwin_time=int(time_win/sampling_time)


    pos_orig=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time,2])
    pos_calc=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time,2])
    fwhm_orig=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time,2])
    fwhm_calc=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time,2])
    density_calc=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time])
    density_samp=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time])
    density_orig=np.zeros([len(noise_vector),len(blob_size),len(spatial_resolution),nwin_time])



    n_noise=len(noise_vector)
    n_size=len(blob_size)
    n_res=len(spatial_resolution)

    for i_noise in range(0, len(noise_vector)):
        for i_size in range(0, len(blob_size)):
            for i_res in range(0, len(spatial_resolution)):
                start_time=time.time()
                spatial_pos=rescale_spat_pos(nwin=nwin,
                                             rad_res=spatial_resolution[i_res],
                                             vert_res=spatial_resolution[i_res])
                if not stationary:
                    moving={'vx':spatial_resolution[i_res]/1e3/time_win,
                            'vy':spatial_resolution[i_res]/1e3/time_win,
                            'v_size':0.,
                            'rot_freq':0.} #it only moves one spatial resolution step during the analysis
                    blob_pos=[np.mean(spatial_pos[:,:,0]),
                              np.mean(spatial_pos[:,:,1])]
                else:
                    moving={'vx':150,'vy':10, 'v_size':1., 'rot_freq':1.}
                    blob_pos=[2200,0]


                filename='tmp/test_fluct_decomp_N_'+str(noise_vector[i_noise])+'_BS_'+str(blob_size[i_size])+'_R_'+str(spatial_resolution[i_res])+'.pickle'

                # if os.path.exists(filename) and (not nocalc) :
                #     with open(filename, 'rb') as f: results = pickle.load(f)
                # else:
                results=test_fluct_decomp(nocalib=True,
                                          #k_factor=k_factor,
                                          test=test,
                                          iterate=True,
                                          noise_level=noise_vector[i_noise],
                                          electronic_noise=electronic_noise,
                                          blob_size=[blob_size[i_size]*spatial_resolution[i_res],
                                                     blob_size[i_size]*spatial_resolution[i_res]],
                                          blob_pos=blob_pos,
                                          moving=moving,
                                          plot=plot,
                                          time_win=time_win,
                                          spatial_pos=spatial_pos,
                                          sampling_time=sampling_time,
                                          )

                blob_density_samp=results['blob_density_samp']
                blob_density_orig=results['blob_density_orig']
                n_vector_calc=results['n_vector_calc']
                n_vector_orig=results['n_vector_orig']
                time_vec=results['time_vec']

                density_samp[i_noise,i_size,i_res,:]=blob_density_samp
                density_orig[i_noise,i_size,i_res,:]=blob_density_orig

                for k in range(nwin):
                    #contour, n_vector_calc[:,:,k], /fill, nlevel=21, /iso, title='N_'+strtrim(noise_vector[i_noise],2)+'_BS_'+strtrim(blob_size[i_size],2)+'_R_'+strtrim(spatial_resolution[i_res],2)
                    pos_orig[i_noise,i_size,i_res,k,0]=blob_pos[0]+moving['vx'] * time_vec[k]
                    pos_calc[i_noise,i_size,i_res,k,0]=np.sum(n_vector_calc[:,:,k]*spatial_pos[:,:,0])/np.sum(n_vector_calc[:,:,k])
                    pos_orig[i_noise,i_size,i_res,k,1]=blob_pos[1]+moving['vy'] * time_vec[k]
                    pos_calc[i_noise,i_size,i_res,k,1]=np.sum(n_vector_calc[:,:,k]*spatial_pos[:,:,1])/np.sum(n_vector_calc[:,:,k])

                    fwhm_orig[i_noise,i_size,i_res,k,0]=blob_size[i_size]*spatial_resolution[i_res]
                    fwhm_orig[i_noise,i_size,i_res,k,1]=blob_size[i_size]*spatial_resolution[i_res]
                    gauss=FitGaussian(x=spatial_pos[:,:,0],
                                      y=spatial_pos[:,:,1],
                                      data=n_vector_calc[:,:,k])
                    size=gauss.size
                    fwhm_calc[i_noise,i_size,i_res,k,0]=size[0]

                    fwhm_calc[i_noise,i_size,i_res,k,1]=size[1]

                    density_calc[i_noise,i_size,i_res,k]=np.sum(n_vector_calc[:,:,k])/(np.pi*size[0]*size[1])

                finish_time=time.time()

                rem_time=(finish_time-start_time)*(n_noise * n_size * n_res - i_noise*n_size*n_res - i_size*n_res - i_res)
                print(f"Remaining time from the calculation: {rem_time//60} min.")

#TODO: maybe a script which could continue the calculation if it takes too long.
            results={'pos_orig':pos_orig,
                     'pos_calc':pos_calc,
                     'fwhm_orig':fwhm_orig,
                     'fwhm_calc':fwhm_calc,
                     'density_orig':density_orig,
                     'density_samp':density_samp,
                     'density_calc':density_calc,
                     }
            with open(filename, 'wb') as f: pickle.dump(results, f)
