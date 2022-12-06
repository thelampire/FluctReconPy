import numpy as np
def rescale_spat_pos(shot=14110,
                     nwin=16.,
                     rad_res=10.,
                     vert_res=10,
                     ):

    spatial_pos_init=None #TODO: getcal_kstar_spat(shot)
    spat_pos_mid=[np.mean(spatial_pos_init[:,:,0]),
                  np.mean(spatial_pos_init[:,:,1])]
    spatial_pos=np.zeroes(nwin,nwin,2)


    for i in range(0, nwin):
        spatial_pos[:,i,0]=spat_pos_mid[0]-nwin/2*rad_res+i*rad_res
        spatial_pos[i,:,1]=spat_pos_mid[1]-nwin/2*vert_res+i*vert_res

    return spatial_pos