import numpy as np
import matplotlib.pyplot as plt

from fluctreconpy import calculate_decomposition,reform_matrix,\
                         get_fluct_resp_matrix,get_fluct_resp_matrix,\
                         fluctuation_matrix_test, get_spatcal

def test_fluct_decomp(nocalib=True,
                      k_factor=None,
                      test=True,
                      iterate=True,
                      floating=True,
                      noise_level=0.01,
                      electronic_noise=0.003,
                      blob_size=[10,10],
                      blob_pos=[2200,40],
                      blob_density_orig=1e20,
                      visual=False,
                      contour=None,
                      spatial_pos=None,
                      pdf=True,
                      nocalc=False,
                      save_file=None,
                      maxiter=100.,
                      moving=None,
                      hole_fluct_amp=0.,
                      plot=False,
                      n_vector_calc=None,
                      n_vector_orig=None,
                      time_vec=None,
                      time_win=50e-6,
                      sampling_time=0.5e-6,
                      threshold=0.001,
                      ):

    spatial_pos = get_spatcal(shot=14110, device='KSTAR', nbi=2)['spatcal']
    n_rad_meas=  len(spatial_pos[0,:,0])
    n_vert_meas= len(spatial_pos[:,0,0])
    n_rad_calc=  len(spatial_pos[0,:,0])

    n_vert_calc=len(spatial_pos[:,0,0])

    electronic_noise=electronic_noise/np.sqrt(20) #1Mhz bandwidth --> 50kHz bandwidth for blobs, approx. power decrease
    #Fill up of the error_vector

    s_vector=np.zeros([n_vert_meas,n_rad_meas])
    n_vector=np.zeros([n_vert_calc,n_rad_calc])
    error_vector=np.zeros([n_vert_meas,n_rad_meas]) #square of the standard error
    if (moving is not None) :
        m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos, test=test) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
        m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])
        n_vector=fluctuation_matrix_test(spatial_pos=spatial_pos,
                                         rad_pos=blob_pos[0],
                                         vert_pos=blob_pos[1],
                                         rad_size=blob_size[0],
                                         vert_size=blob_size[1],
                                         blob_density=blob_density_orig)

        n_vector_ref = np.reshape(n_vector.T, n_rad_meas*n_vert_meas)
        blob_density_samp=np.sum(n_vector_ref)/(np.pi*blob_size[0]*blob_size[1])
        s_vector_ref = np.matmul(m_matrix_ref, n_vector_ref) + (np.max(np.abs(np.matmul(m_matrix_ref, n_vector_ref)))*noise_level+electronic_noise)*np.random.rand(64) #s_vector is a simulated measurement)
        error_vector_s_ref=np.abs((np.matmul(m_matrix_ref,n_vector_ref))*noise_level+electronic_noise)**2
        s_vector=(s_vector_ref, n_vert_meas, n_rad_meas)
        error_vector=np.reshape(error_vector_s_ref,(n_vert_meas, n_rad_meas))

        n_vector_calc=calculate_decomposition(light_profile=s_vector,
                                              error_profile=error_vector,
                                              m_matrix=m_matrix,
                                              spatial_pos=spatial_pos,
                                              floating=floating,
                                              iterate=iterate,
                                              k_factor=k_factor,
                                              visual=visual,
                                              contour=contour,
                                              pdf=pdf,
                                              nocalc=nocalc,
                                              save_file=save_file,
                                              n_vector_test=n_vector)

    else:
        m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos, test=test) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
        m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])
        n_vector=fluctuation_matrix_test(rad_pos=blob_pos[0],
                                         vert_pos=blob_pos[1],
                                         moving=True,
                                         vert_size=blob_size[1],
                                         rad_size=blob_size[0],
                                         vx=moving['vx'],
                                         vy=moving['vy'],
                                         v_size=moving['v_size'],
                                         rot_freq=moving['rot_freq'],
                                         time_vec=time_vec,
                                         time_win=time_win,
                                         hole_fluct_amp=hole_fluct_amp,
                                         spatial_pos=spatial_pos,
                                         sampling_time=sampling_time,
                                         blob_density=blob_density_orig)
        n_vector_calc=n_vector
        blob_density_samp=np.zeros(len(time_vec))
        for i in range[len(time_vec)]:
            value=np.max(np.abs(np.matmul(m_matrix_ref, np.reshape(((n_vector[:,:,0])).T, n_rad_meas*n_vert_meas))))
            noise=(value*noise_level+electronic_noise)*np.random.rand(n_rad_meas*n_vert_meas)
            n_vector_ref = np.reshape((n_vector[:,:,i]).T, n_rad_meas*n_vert_meas)
            blob_density_samp[i]=np.sum(n_vector_ref)/(np.pi*blob_size[0]*blob_size[1])
            s_vector_ref = np.matmul(m_matrix_ref, n_vector_ref) + noise #s_vector is a simulated measurement
            error_vector_s_ref=abs((np.matmul(m_matrix_ref, n_vector_ref))*noise_level+electronic_noise)**2
            s_vector=np.reshape(s_vector_ref,(n_vert_meas, n_rad_meas))
            error_vector=np.reshape(error_vector_s_ref,(n_vert_meas, n_rad_meas)).T

            n_vector_calc[:,:,i]=calculate_decomposition(light_profile=s_vector,
                                                         error_profile=error_vector,
                                                         m_matrix=m_matrix,
                                                         spatial_pos=spatial_pos,
                                                         floating=floating,
                                                         iterate=iterate,
                                                         k_factor=k_factor,
                                                         contour=contour,
                                                         pdf=pdf,
                                                         nocalc=nocalc,
                                                         save_file=save_file,
                                                         n_vector_test=n_vector[:,:,i])
            if plot:
                # if i == 0:
                #     zrange=zrange[[n_vector,n_vector_calc]]
                nr=len(spatial_pos[0,:,0])
                nz=len(spatial_pos[:,0,0])
                r_vec=((spatial_pos[nz/2-1,:,0]+spatial_pos[nz/2,:,0])/2)
                z_vec=((spatial_pos[:,nr/2-1,1]+spatial_pos[:,nr/2,1])/2)

                p_vector= np.reshape(n_vector_calc[:,:,i].T, n_rad_calc*n_vert_calc)
                # if i == 0:
                #     zmax_light=np.max([np.max(s_vector_ref),
                #                        np.max(np.matmul(m_matrix_ref, p_vector))])
                # if i == 0:
                #     zmin_light=np.min([np.min(s_vector_ref),
                #                        np.min(np.matmul(m_matrix_ref, p_vector))])

                plt.cla()
                figsize=(8.5/2.54,8.5/2.54*1.5)
                fig,axs=plt.subplots(5,1,figsize=figsize)

                ax=axs[0]
                plt.contourf((n_vector[:,:,i]).T,
                             r_vec,
                             z_vec,
                             nlevel=51,
                             lw=3,
                             )
                #             position=[0.05,0.8,0.45,0.95])
                ax.set_title("Original density")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")

                ax=axs[1]
                ax.contourf((n_vector_calc[:,:,i]).T,
                            r_vec,
                            z_vec, nlevel=21,
                            #position=[0.05,0.55,0.45,0.7],
                            lw=3,
                            )
                ax.set_title("Reconstructed density")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")

                ax=axs[2]
                ax.contourf(np.reshape(s_vector_ref, (n_rad_meas, n_vert_meas)),
                            r_vec,
                            z_vec,
                            nlevel=21,
                            #position=[0.5,0.8,0.95,0.95],
                            lw=3,
                            )
                ax.set_title("Original light with noise")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")
                #zrange=[zmin_light,zmax_light]

                ax=axs[2]
                ax.contourf(np.reshape(np.matmul(m_matrix_ref, p_vector),(n_rad_meas,n_vert_meas)),
                            r_vec,
                            z_vec,
                            nlevel=21,
                            #position=[0.5,0.55,0.95,0.7],
                            lw=3,
                            )
                ax.set_title("Reconstructed light")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")

                ax=axs[3]
                ax.plot(r_vec,
                        n_vector[1,:,i],
                        lw=3,
                        )
                ax.set_xlabel("R [mm]")

                ax=axs[4]
                plt.plot(r_vec,
                         n_vector_calc[1,:,i],
                         lw=3,
                         )

                ax.set_xlabel("R [mm]")
                #print, total((n_vector_calc[:,:,i]-n_vector[:,:,i])**2)/n_vert_calc/n_rad_calc
    results={'blob_density_orig':blob_density_orig,
             'blob_density_samp':blob_density_samp,
        }
    return results