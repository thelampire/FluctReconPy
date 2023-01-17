import numpy as np

from fluctreconpy import get_spatcal

def rescale_spat_pos(shot=14110,
                     nwin=16.,
                     rad_res=10.,
                     vert_res=10,
                     ):

    spatial_pos_init=spatial_pos = get_spatcal(shot=14110, device='KSTAR', nbi=123)['spatcal']
    spat_pos_mid=[np.mean(spatial_pos_init[:,:,0]),
                  np.mean(spatial_pos_init[:,:,1])]
    spatial_pos=np.zeros([nwin,nwin,2])


    for i in range(0, nwin):
        spatial_pos[:,i,0]=spat_pos_mid[0]-nwin/2*rad_res+i*rad_res
        spatial_pos[i,:,1]=spat_pos_mid[1]-nwin/2*vert_res+i*vert_res
    spatial_pos=np.flip(spatial_pos, axis=0)
    return spatial_pos