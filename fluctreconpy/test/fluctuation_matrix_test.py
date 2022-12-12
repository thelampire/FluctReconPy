import numpy as np
import matplotlib.pyplot as plt

def fluctuation_matrix_test(fluct_amp=0.05,
                            spatial_pos=None, #TODO: originally this was getcal_kstar_spat(14110)
                            plot=False,

                            moving=False,
                            rad_size=10.,
                            vert_size=10.,
                            rad_pos=2200.,
                            vert_pos=15.,
                            vx=100,
                            vy=50,
                            v_size=50.,
                            rot_freq=100.,

                            hole_size=10.,
                            hole_rad_pos=2200.,
                            hole_vert_pos=15.,
                            hole_fluct_amp=-0.02,
                            hole_vx=-50.,
                            hole_vy=100.,

                            time_win=1e-3,
                            sampling_time=0.5e-6,
                            blob_density=1e20,
                            test=False,
                            ):

    """
        #************************************************
        #*** fluctuation_matrix_test    by M. Lampert
        #************************************************
        #
        #Returns a simulated Gaussian blob with a given
        #radial & vertical size & a given radial and
        #vertical position. The amplitude of the blob
        #can be also set.
        #
        #************************************************
        #
        #INPUTs:
        #  shot: shotnumber for spatial calibration
        #  rad_size: radial size of the blob
        #  vert_size: vertical size of the blog
        #  rad_pos: radial position of the blob
        #  vert_pos: vertical position of the blob
        #  fluct_amp: amplitude of the blob
        #
        #OUTPUTs:
        #  Matrix of the blob with the same coordinates
        #  as the BES measurement.
    """



    n_rad=len(spatial_pos[0,:,0])
    n_vert=len(spatial_pos[:,0,0])
    levels=np.linspace(hole_fluct_amp, fluct_amp, 51)

    if moving:
        n_time=time_win/sampling_time
        return_matrix=np.zeros([n_vert,n_rad,n_time])
        blob_density=np.zeros(n_time) # Gaussian integral divided by the ellipse defined by the sigma_r, sigma_z of the Gaussian (surface density)
        time_vec=np.arange(n_time)/n_time*time_win
        for i in range(0, n_time):
            t=i*sampling_time
            arg=2*np.pi*rot_freq*t
            arg2=0.
            for j in range(0, n_vert):
                for k in range(0, n_rad):
                    #blob
                    x=spatial_pos[j,k,0]-(vx*1e3*t+rad_pos)
                    y=spatial_pos[j,k,1]-(vy*1e3*t+vert_pos)
                    a1=(np.cos(arg)/(rad_size+v_size*1e3*t))**2+(np.sin(arg)/(vert_size))**2
                    b1=-0.5*np.sin(2*arg)/(rad_size+v_size*1e3*t)**2+0.5*np.sin(2*arg)/(vert_size)**2
                    c1=(np.sin(arg)/(rad_size+v_size*1e3*t))**2+(np.cos(arg)/(vert_size))**2
                    return_matrix[j,k,i]=fluct_amp*np.exp(-0.5*(a1*x**2+2*b1*x*y+c1*y**2))
                    #hole
                    x=spatial_pos[j,k,0]-(hole_vx*1e3*t+hole_rad_pos)
                    y=spatial_pos[j,k,1]-(hole_vy*1e3*t+hole_vert_pos)
                    a2=(np.cos(arg2)/(hole_size))**2+(np.sin(arg2)/(hole_size))**2
                    b2=-0.5*np.sin(2*arg2)/(hole_size)**2+0.5*np.sin(2*arg2)/(vert_size)**2
                    c2=(np.sin(arg2)/(hole_size+v_size*1e3*t))**2+(np.cos(arg2)/(hole_size))**2
                    return_matrix[j,k,i]+=hole_fluct_amp*np.exp(-0.5*(a2*x**2+2*b2*x*y+c2*y**2))

                    min_axis=((a1+c1)-np.sqrt((a1+c1)**2 + 4*(b1**2/4 - a1*c1)))/2
                    maj_axis=((a1+c1)+np.sqrt((a1+c1)**2 + 4*(b1**2/4 - a1*c1)))/2
                    blob_density[i]=fluct_amp*(4*np.pi/np.sqrt(4*a1*c1-b1**2))/(np.pi*maj_axis*min_axis)

                    if plot:
                        plt.contourf(return_matrix[:,:,i],
                                     (spatial_pos[:,:,0]),
                                     (spatial_pos[:,:,1]),
                                     levels=levels)
                        plt.scatter(spatial_pos[:,:,0],
                                    spatial_pos[:,:,1],
                                    color='red',
                                    marker='x',
                                    )

    else:
        return_matrix=fluct_amp*np.exp(-0.5*(((spatial_pos[:,:,0]-rad_pos)/rad_size)**2+((spatial_pos[:,:,1]-vert_pos)/vert_size)**2))
        blob_density=fluct_amp*2/np.sqrt(rad_size*vert_size)
        levels=np.linspace(0,fluct_amp, 51)
        if (plot is None):
            plt.contourf(return_matrix,
                         spatial_pos[:,:,0],
                         spatial_pos[:,:,1],
                         levels=levels)

            plt.scatter(spatial_pos[:,:,0],
                        spatial_pos[:,:,1],
                        color='red',
                        marker='x',
                        )


    return return_matrix
