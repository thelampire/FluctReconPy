import numpy as np

def get_fluct_resp_matrix(timerange=None, test=True, reform=None, spatial_pos=None):


    if test :
        if spatial_pos is None :
            raise ValueError('The spatial coordinates need to be defined. Returning...')

        #The matrix itself is a Lower triangular matrix
        nr=len(spatial_pos[0,:,0])
        nz=len(spatial_pos[:,0,0])
        tauvec=np.zeros([nr,nz])
        matrix=np.zeros([nz,nr,nz,nr])
        energy=90e3 #eV
        e=1.6e-19
        mh=1.6727e-27
        vel=np.sqrt(2*energy*e/mh)
        #tauvec=[10.3,10.2,10.1,10.,9.,8.,7.,6.,5.,4.,3.,2.6,2.3,2.,1.5,1.]*1e-9 #ns
        for i_rad in range(0, nr):
            for j_vert in range(0, nz):
                tauvec[i_rad,j_vert]=(np.tanh((spatial_pos[j_vert,i_rad,0]-spatial_pos[j_vert,0,0])/20.-2.)+1.)*5.*1e-9

        for k in range(0, nz):
            for j in range(0, nr):
                for i in range(j, nr):
                    matrix[k,i,k,j]=1e-9/tauvec[i,k]*np.exp(-((i-j)*0.01)/(vel*tauvec[i,k]))

        return matrix
    else:
        raise NotImplementedError('No fluctuation response matrix is available at the moment.')
        pass