import numpy as np
import matplotlib.pyplot as plt
import copy

from fluctreconpy import calculate_decomposition,reform_matrix,\
                         get_fluct_resp_matrix,get_fluct_resp_matrix,\
                         get_spatcal, rescale_spat_pos
from fluctreconpy.test import fluctuation_matrix_test

def test_fluct_decomp(nocalib=True,
                      k_factor=[1e-10,1e3],
                      test=True,
                      iterate=True,
                      floating=True,
                      noise_level=0.01, #relative
                      electronic_noise=0.003, #mV
                      fluct_amp=0.05, #percent
                      blob_size=[10.,10.], #mm
                      blob_pos=[2260.,0.], #mm#
                      blob_density_orig=1e20,
                      visual=False,
                      contour=None,
                      spatial_pos=None,
                      pdf=True,
                      nocalc=False,
                      save_filename='tmp/calculate_decomposition_save.pickle',
                      maxiter=100.,
                      moving={'vx':1000,
                              'vy':500,
                              'v_size':1e-3/0.5e-6,
                              'rot_freq':0},
                      hole_fluct_amp=0.,
                      plot=False,
                      n_vector_calc=None,
                      n_vector_orig=None,
                      time_vec=None,
                      time_win=50e-6,
                      sampling_time=0.5e-6,
                      threshold=0.001,
                      ):
    if spatial_pos is None:
        spatial_pos=rescale_spat_pos(nwin=16,
                                     rad_res=10,
                                     vert_res=10)

    n_rad_meas=  len(spatial_pos[0,:,0])
    n_vert_meas= len(spatial_pos[:,0,0])
    n_rad_calc=  len(spatial_pos[0,:,0])
    n_vert_calc=len(spatial_pos[:,0,0])

    electronic_noise=electronic_noise/np.sqrt(20) #1Mhz bandwidth --> 50kHz bandwidth for blobs, approx. power decrease
    #Fill up of the error_vector

    s_vector=np.zeros([n_vert_meas,n_rad_meas])
    n_vector=np.zeros([n_vert_calc,n_rad_calc])
    error_vector=np.zeros([n_vert_meas,n_rad_meas]) #square of the standard error

    if (moving is None) :
        m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos,
                                       test=test) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
        m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])
        results=fluctuation_matrix_test(spatial_pos=spatial_pos,

                                        rad_pos=blob_pos[0],
                                        vert_pos=blob_pos[1],
                                        rad_size=blob_size[0],
                                        vert_size=blob_size[1],

                                        fluct_amp=fluct_amp,
                                        blob_density=blob_density_orig,
                                        test=test,
                                        )
        n_vector=results['fluct_matrix']

        if test:
            print(f"n_vector: {n_vector}")
        n_vector_ref = np.reshape(n_vector, n_rad_meas*n_vert_meas)
        blob_density_samp=np.sum(n_vector_ref)/(np.pi*blob_size[0]*blob_size[1])

        s_ref_no_noise=np.matmul(m_matrix_ref, n_vector_ref.T)

        noise=(s_ref_no_noise*noise_level + electronic_noise)*(np.random.rand(spatial_pos.shape[0]*spatial_pos.shape[1])-0.5)*2 #s_vector is a simulated measurement)
        s_vector_ref = s_ref_no_noise + noise

        if test:
            print(f"s_vector_ref: {s_vector_ref}")

        error_vector_s_ref=np.abs((np.matmul(m_matrix_ref,n_vector_ref))*noise_level+electronic_noise)**2
        s_vector=np.reshape(s_vector_ref, (n_vert_meas, n_rad_meas),
                            )
        error_vector=np.reshape(error_vector_s_ref,(n_vert_meas, n_rad_meas))

        results=calculate_decomposition(light_profile=s_vector,
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
                                        save_filename=save_filename,
                                        n_vector=n_vector,
                                        test=test,
                                        )
        time_vec=None
        n_vector_calc=results['p_vector']
        if test:
            print(f"n_vector_calc: {n_vector_calc}")
    else:
        m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos,
                                       test=test) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
        m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])

        results=fluctuation_matrix_test(rad_pos=blob_pos[0],
                                        vert_pos=blob_pos[1],
                                        rad_size=blob_size[0],
                                        vert_size=blob_size[1],

                                         moving=True,
                                         vx=moving['vx'],
                                         vy=moving['vy'],
                                         v_size=moving['v_size'],
                                         rot_freq=moving['rot_freq'],

                                         time_win=time_win,
                                         hole_fluct_amp=hole_fluct_amp,
                                         spatial_pos=spatial_pos,
                                         sampling_time=sampling_time,
                                         blob_density=blob_density_orig)

        n_vector=results['fluct_matrix']
        n_vector_calc=copy.deepcopy(results['fluct_matrix'])
        time_vec=results['time_vec']
        blob_density_samp=np.zeros(len(time_vec))

        figsize=(8.5/2.54,8.5/2.54*1.5)
        fig,axs=plt.subplots(6,1,figsize=figsize)

        for i in range(len(time_vec)):
            value=np.max(np.abs(np.matmul(m_matrix_ref,
                                          np.reshape(((n_vector[:,:,0])).T,
                                                     n_rad_meas*n_vert_meas))))
            noise=(value*noise_level+electronic_noise)*np.random.rand(n_rad_meas*n_vert_meas)

            n_vector_ref = np.reshape((n_vector[:,:,i]), n_rad_meas*n_vert_meas)
            blob_density_samp[i]=np.sum(n_vector_ref)/(np.pi*blob_size[0]*blob_size[1])

            s_vector_ref = np.matmul(m_matrix_ref, n_vector_ref) + noise #s_vector is a simulated measurement
            error_vector_s_ref=abs((np.matmul(m_matrix_ref, n_vector_ref))*noise_level+electronic_noise)**2

            s_vector=np.reshape(s_vector_ref,(n_vert_meas, n_rad_meas))
            error_vector=np.reshape(error_vector_s_ref,(n_vert_meas, n_rad_meas))

            results=calculate_decomposition(light_profile=s_vector,
                                            error_profile=error_vector,
                                            m_matrix=m_matrix,
                                            spatial_pos=spatial_pos,

                                            floating=floating,
                                            iterate=iterate,
                                            k_factor=k_factor,
                                            contour=contour,
                                            pdf=pdf,
                                            nocalc=nocalc,
                                            save_filename=save_filename,
                                            )

            n_vector_calc[:,:,i]=results['p_vector']

            if plot:
                # if i == 0:
                #     zrange=zrange[[n_vector,n_vector_calc]]
                nr=len(spatial_pos[0,:,0])
                nz=len(spatial_pos[:,0,0])
                r_vec=((spatial_pos[int(nz/2-1),:,0]+spatial_pos[int(nz/2),:,0])/2)
                z_vec=((spatial_pos[:,int(nr/2-1),1]+spatial_pos[:,int(nr/2),1])/2)

                p_vector= np.reshape(n_vector_calc[:,:,i].T, n_rad_calc*n_vert_calc)
                # if i == 0:
                #     zmax_light=np.max([np.max(s_vector_ref),
                #                        np.max(np.matmul(m_matrix_ref, p_vector))])
                # if i == 0:
                #     zmin_light=np.min([np.min(s_vector_ref),
                #                        np.min(np.matmul(m_matrix_ref, p_vector))])

                ax=axs[0]
                ax.cla()
                ax.contourf(r_vec,
                            z_vec,
                            (n_vector[:,:,i]).T,
                            nlevel=51,
                            lw=3,
                            )
                #             position=[0.05,0.8,0.45,0.95])
                ax.set_title("Original density")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")

                ax=axs[1]
                ax.cla()
                ax.contourf(r_vec,
                            z_vec,
                            (n_vector_calc[:,:,i]).T,
                            nlevel=21,
                            #position=[0.05,0.55,0.45,0.7],
                            lw=3,
                            )
                ax.set_title("Reconstructed density")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")

                ax=axs[2]
                ax.cla()
                ax.contourf(r_vec,
                            z_vec,
                            np.reshape(s_vector_ref, (n_rad_meas, n_vert_meas)),
                            nlevel=21,
                            lw=3,
                            )
                ax.set_title("Original light with noise")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")
                #zrange=[zmin_light,zmax_light]

                ax=axs[3]
                ax.cla()
                ax.contourf(r_vec,
                            z_vec,
                            np.reshape(np.matmul(m_matrix_ref, p_vector),(n_rad_meas,n_vert_meas)),
                            nlevel=21,
                            lw=3,
                            )
                ax.set_title("Reconstructed light")
                ax.set_xlabel("R [mm]")
                ax.set_ylabel("z [mm]")

                ax=axs[4]
                ax.cla()
                ax.plot(r_vec,
                        n_vector[1,:,i],
                        lw=3,
                        )
                ax.set_xlabel("R [mm]")

                ax=axs[5]
                ax.cla()
                plt.plot(r_vec,
                         n_vector_calc[1,:,i],
                         lw=3,
                         )

                ax.set_xlabel("R [mm]")
                plt.pause(0.1)
                #print, total((n_vector_calc[:,:,i]-n_vector[:,:,i])**2)/n_vert_calc/n_rad_calc
    results={'blob_density_orig':blob_density_orig,
             'blob_density_samp':blob_density_samp,
             'n_vector_calc':n_vector_calc,
             'n_vector_orig':n_vector,
             'spatial_pos':spatial_pos,
             'time_vec':time_vec,
             # 's_vector':s_vector,
             # 's_vector_no_noise':np.reshape(s_ref_no_noise, (n_vert_meas, n_rad_meas))
        }
    return results