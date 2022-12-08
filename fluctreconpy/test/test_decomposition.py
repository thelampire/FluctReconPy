import numpy as np
import matplotlib.pyplot as plt

from fluctreconpy import calculate_decomposition

def test_decomposition(n_noise=100.,
                       test_density=False,
                       test_blob_position=False,
                       test_blob_size=False,
                       n_error=20.,
                       blob_pos=[2200.,20.],
                       blob_size=[20.,30.],
                       ):

    noise_vec=np.arange(n_noise)/n_noise

    if (test_density):

        error_mat=np.zeros([n_noise,n_error])
        for i in range(0, n_noise):
            for j in range(0, n_error):
                print( (i*n_error+j)/(n_error*n_noise))
                results=calculate_decomposition(14110,
                                                   time=2.5,
                                                   nocalib=True,
                                                   floating=0,
                                                   blob_pos=[2200.,20.],
                                                   blob_size=[20.,30.],
                                                   electronic_noise=0.003/25.,
                                                   noise_level=noise_vec[i],
                                                   iterate=True,
                                                   )
                n_vec_sim=results['n_vec_sim']
                n_vec_calc=results['n_vec_calc']
                spatial_pos=results['spatial_pos']

                error_mat[i,j]=np.sqrt(np.mean(1./(len(n_vec_calc))*((n_vec_calc-n_vec_sim))**2))/np.mean(n_vec_sim)


            plt.plot(noise_vec,
                     np.sum(error_mat,axis=1)/n_error)
            error_vector=np.zeros([n_noise])
            for i in range(0,n_noise):
                error_vector[i]=np.sqrt(np.var(error_mat[i,:]))
            plt.errorbar(noise_vec,
                         np.sum(error_mat,axis=1)/n_error,
                         xerr=error_vector,
                         )


    if (test_blob_position is None) :
        position_mat=np.zeros([n_noise,n_error,2])
        size_mat=np.zeros([n_noise,n_error,2])
        for i in range(0, n_noise):
            for j in range(0, n_error):
                print((i*n_error+j)/(n_error*n_noise))
                #TODO:return as dict
                results=calculate_decomposition(14110,
                                                time=2.5,
                                                nocalib=True,
                                                floating=0,
                                                blob_pos=blob_pos,
                                                blob_size=blob_size,
                                                electronic_noise=0.003/25.,
                                                noise_level=noise_vec[i],
                                                iterate=True,
                                                )

                n_vec_sim=results['n_vec_sim']
                n_vec_calc=results['n_vec_calc']
                spatial_pos=results['spatial_pos']
                ind=np.argmax(n_vec_calc)
                #position_mat[i,j,*]=sqrt((spatial_pos[ind[0],ind[1],0]-blob_pos[0])**2+(spatial_pos[ind[0],ind[1],1]-blob_pos[1])**2) #maximum intensity calculation
                position_mat[i,j,:]=np.abs((np.sum(n_vec_calc*(spatial_pos[:,:,0]))/np.sum(n_vec_calc),
                                            np.sum(n_vec_calc*(spatial_pos[:,:,1]))/np.sum(n_vec_calc))-blob_pos) #mass centre calculation
                size_mat[i,j,0]=np.sqrt(np.sum(n_vec_calc*(spatial_pos[:,:,0]-position_mat[i,j,0]-blob_pos[0])**2)/np.sum(n_vec_calc))-blob_size[0]
                size_mat[i,j,1]=np.sqrt(np.sum(n_vec_calc*(spatial_pos[:,:,1]-position_mat[i,j,1]-blob_pos[1])**2)/np.sum(n_vec_calc))-blob_size[1]


            plt.plot(noise_vec,np.sum(position_mat,2)/n_error)
            error_vector=np.zeros([n_noise])
            for i in range(n_noise):
                error_vector[i]=np.sqrt(np.var(position_mat[i,:,0]+np.var(position_mat[i,:,1])))
            plt.errorbar(noise_vec,
                         np.sum(position_mat,axis=1)/n_error,
                         xerr=error_vector,
                         )