import numpy as np
import matplotlib.pyplot as plt
import copy
import pickle
import os

from fluctreconpy import calculate_decomposition,reform_matrix,\
                         get_fluct_resp_matrix,get_fluct_resp_matrix,\
                         get_spatcal, rescale_spat_pos
from fluctreconpy.test import fluctuation_matrix_test

def test_fluct_decomp(nocalib=True,
                      k_factor=[1e-10,1e3],

                      floating=True,
                      noise_level=0.01, #relative
                      electronic_noise=0.003, #mV
                      fluct_amp=0.05, #percent
                      blob_size=[10.,10.], #mm
                      blob_pos=[2200.,0.], #mm#
                      blob_density_orig=1e20,
                      visual=False,
                      contour=None,
                      spatial_pos=None,
                      pdf=True,
                      recalc=False,
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
                      test=False,
                      iter_method='scipy_bisect',
                      transposed=False,
                      decompose=True,
                      ):

    if spatial_pos is None:
        spatial_pos=rescale_spat_pos(nwin=16,
                                     rad_res=10,
                                     vert_res=20,
                                     transposed=transposed)

    electronic_noise=electronic_noise/np.sqrt(20) #1Mhz bandwidth --> 50kHz bandwidth for blobs, approx. power decrease
    #Fill up of the error_vector
    s_vector=np.zeros([spatial_pos.shape[0],
                       spatial_pos.shape[1]])
    n_vector=np.zeros([spatial_pos.shape[0],
                       spatial_pos.shape[1]])
    error_vector=np.zeros([spatial_pos.shape[0],
                           spatial_pos.shape[1]]) #square of the standard error

    m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos,
                                   test=True,
                                   transposed=transposed
                                   ) #[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]

    m_matrix_ref = np.reshape(m_matrix,
                              (spatial_pos.shape[0]*spatial_pos.shape[1],
                               spatial_pos.shape[0]*spatial_pos.shape[1]),
                                 #transposed=transposed
                                 )
    results=fluctuation_matrix_test(rad_pos=blob_pos[0],
                                    vert_pos=blob_pos[1],
                                    rad_size=blob_size[0],
                                    vert_size=blob_size[1],

                                    vx=moving['vx'],
                                    vy=moving['vy'],
                                    v_size=moving['v_size'],
                                    rot_freq=moving['rot_freq'],

                                    time_win=time_win,
                                    hole_fluct_amp=hole_fluct_amp,
                                    spatial_pos=spatial_pos,
                                    sampling_time=sampling_time,
                                    blob_density=blob_density_orig,
                                    )

    n_vector=results['fluct_matrix']

    # plt.contourf(spatial_pos[:,:,0],
    #              spatial_pos[:,:,1],
    #              np.reshape(np.matmul(m_matrix_ref, np.reshape((n_vector[:,:,0]), spatial_pos.shape[0]*spatial_pos.shape[1])),
    #                         (spatial_pos.shape[0],spatial_pos.shape[1]) ))
    # return

    n_vector_calc=copy.deepcopy(results['fluct_matrix'])
    if test:
        print(f"n_vector: {n_vector}")

    time_vec=results['time_vec']
    blob_density_samp=np.zeros(len(time_vec))

    if os.path.exists(save_filename) and not recalc:
        with open(save_filename, 'rb') as f: results = pickle.load(f)
    else:
        #figsize=(8.5/2.54,8.5/2.54*1.5)
        for i in range(len(time_vec)):

            n_vector_ref = np.reshape((n_vector[:,:,i]),
                                      spatial_pos.shape[0]*spatial_pos.shape[1])
            s_ref_no_noise=np.matmul(m_matrix_ref,n_vector_ref)
            #noise=(np.max(np.abs(s_ref_no_noise))*noise_level+electronic_noise)*np.random.rand(n_rad_meas*n_vert_meas) #oldie but not goodie
            noise=(s_ref_no_noise*noise_level + electronic_noise)*(np.random.rand(spatial_pos.shape[0]*spatial_pos.shape[1])-0.5)*2 #s_vector is a simulated measurement)

            blob_density_samp[i]=np.sum(n_vector_ref)/(np.pi*blob_size[0]*blob_size[1])

            s_vector_ref = np.matmul(m_matrix_ref, n_vector_ref) + noise #s_vector is a simulated measurement
            error_vector_s_ref=np.abs((np.matmul(m_matrix_ref, n_vector_ref))*noise_level+electronic_noise)**2

            if test:
                print(f"s_vector_ref: {s_vector_ref}")
            s_vector=np.reshape(s_vector_ref,(spatial_pos.shape[0],
                                              spatial_pos.shape[1]))
            error_vector=np.reshape(error_vector_s_ref,(spatial_pos.shape[0],
                                                        spatial_pos.shape[1]))

            if transposed:
                r_vec=np.mean(spatial_pos,axis=1)[:,0]
                z_vec=np.mean(spatial_pos,axis=0)[:,1]
            else:
                r_vec=np.mean(spatial_pos,axis=0)[:,0]
                z_vec=np.mean(spatial_pos,axis=1)[:,1]
            if not decompose:
                p_vector=np.matmul(np.linalg.inv(m_matrix_ref), s_vector_ref)
                p_vector=np.reshape(p_vector,(m_matrix.shape[0],m_matrix.shape[1]))
                n_vector_calc[:,:,i]=p_vector
            else:
                results=calculate_decomposition(light_profile=s_vector,
                                                error_profile=error_vector,

                                                m_matrix=m_matrix,
                                                r_vec=r_vec,
                                                z_vec=z_vec,

                                                floating=floating,
                                                k_factor=k_factor,
                                                contour=contour,
                                                pdf=pdf,
                                                iter_method=iter_method,
                                                transposed=transposed
                                                )

                n_vector_calc[:,:,i]=results['p_vector']
            if test:
                print(f"n_vector_calc: {n_vector_calc}")
        results={'blob_density_orig':blob_density_orig,
                 'blob_density_samp':blob_density_samp,
                 'n_vector_calc':n_vector_calc,
                 'n_vector_orig':n_vector,
                 'spatial_pos':spatial_pos,
                 'time_vec':time_vec,
                 's_vector':s_vector,
                 's_vector_no_noise':np.reshape(s_ref_no_noise, (spatial_pos.shape[0],spatial_pos.shape[1]))
                 }
        with open(save_filename, 'wb') as f: pickle.dump(results,f)
    if plot:
        #for i in range(len(results['time_vec'])):
            i=0
            plot_test_fluct_decomp(spatial_pos=results['spatial_pos'],
                                   n_vector=results['n_vector_orig'],
                                   n_vector_calc=results['n_vector_calc'],
                                   s_vector=results['s_vector'],
                                   m_matrix_ref=m_matrix_ref,
                                   ind=i,
                                   transposed=transposed)
    return results

def plot_test_fluct_decomp(spatial_pos=None,
                           n_vector=None,
                           n_vector_calc=None,
                           s_vector=None,
                           m_matrix_ref=None,
                           ind=None,
                           transposed=False,
                           order='C'):
    # if i == 0:
    #     zrange=zrange[[n_vector,n_vector_calc]]
    r_vec=spatial_pos[:,:,0]
    z_vec=spatial_pos[:,:,1]

    plt.close("all")
    fig,axs=plt.subplots(6,1,figsize=(17/2.54,17/2.54*1.5))
    ax=axs[0]
    ax.cla()
    con=ax.contourf(r_vec,
                    z_vec,
                    n_vector[:,:,ind],
                    nlevel=51,
                    lw=3,
                    )
    ax.set_title("Original density")
    ax.set_xlabel("R [mm]")
    ax.set_ylabel("z [mm]")

    fig.colorbar(con,ax=ax)

    ax=axs[1]
    ax.cla()
    con=ax.contourf(r_vec,
                    z_vec,
                    n_vector_calc[:,:,ind],
                    nlevel=21,
                    lw=3,
                    )
    ax.set_title("Reconstructed density")
    ax.set_xlabel("R [mm]")
    ax.set_ylabel("z [mm]")
    fig.colorbar(con,ax=ax)

    ax=axs[2]
    ax.cla()
    con=ax.contourf(r_vec,
                    z_vec,
                    s_vector,
                    nlevel=21,
                    lw=3,
                    )
    ax.set_title("Original light with noise")
    ax.set_xlabel("R [mm]")
    ax.set_ylabel("z [mm]")
    fig.colorbar(con,ax=ax)
    ax=axs[3]
    ax.cla()
    if transposed:
        con=ax.contourf(r_vec,
                        z_vec,
                        np.reshape(np.matmul(m_matrix_ref, np.ravel(n_vector_calc[:,:,ind])),
                                   (spatial_pos.shape[0],spatial_pos.shape[1])),
                        nlevel=21,
                        lw=3,
                        )
        fig.colorbar(con,ax=ax)
        ax.set_title("Reconstructed light")
        ax.set_xlabel("R [mm]")
        ax.set_ylabel("z [mm]")

        ax=axs[4]
        ax.cla()
        ax.plot(r_vec[:,4],
                np.mean(n_vector[:,:,ind], axis=1),
                lw=3,
                )
        ax.set_xlabel("R [mm]")

        ax=axs[5]
        ax.cla()
        plt.plot(r_vec[:,4],
                 np.mean(n_vector_calc[:,:,ind], axis=1),
                 lw=3,
                 )

        ax.set_xlabel("R [mm]")
        plt.pause(0.1)
        #print(np.sum((n_vector_calc[:,:,i]-n_vector[:,:,i])**2)/n_vert_calc/n_rad_calc)
    else:
        con=ax.contourf(r_vec,
                    z_vec,
                    np.reshape(np.matmul(m_matrix_ref, np.ravel(n_vector_calc[:,:,ind])),(spatial_pos.shape[0],spatial_pos.shape[1])),
                    nlevel=21,
                    lw=3,
                    )
        fig.colorbar(con,ax=ax)
        ax.set_title("Reconstructed light")
        ax.set_xlabel("R [mm]")
        ax.set_ylabel("z [mm]")

        ax=axs[4]
        ax.cla()
        ax.plot(r_vec[4,:],
                np.mean(n_vector[:,:,ind], axis=0),
                lw=3,
                )
        ax.set_xlabel("R [mm]")

        ax=axs[5]
        ax.cla()
        plt.plot(r_vec[4,:],
                 np.mean(n_vector_calc[:,:,ind], axis=0),
                 lw=3,
                 )

        ax.set_xlabel("R [mm]")
        plt.pause(0.1)
        #print(np.sum((n_vector_calc[:,:,i]-n_vector[:,:,i])**2)/n_vert_calc/n_rad_calc)
