import numpy as np
def undulation_matrix(spatial_pos=None,
                      r_vec=None,
                      z_vec=None,
                      floating=False,
                      onedim=True):

    """
    #*********************************************************
    #*** undulation_matrix.pro    by M. Lampert 2016-11-02 ***
    #*********************************************************
    #
    # The procedure calculates the undulation matrix for a
    # given geometry. The return matrix is a 4D matrix with
    # dimensions of [n_r,n_z,n_r,n_z]. The current calculation
    # is using bilinear base functions, which are placed on an
    # equidistant grid.
    #
    #  INPUTs:
    #      edges instead of zero edges. (default: 0)
    #    spatial_pos: [z_vec, r_vec, 2] matrix for the spatial
    #      coordinates. The center line is taken for the
    #      calculation for both radial & vertical directions.
    #    r_vec: 1D radial coordinate vector
    #    v_vec: 1D vertical coordinate vector
    #    /floating: the calculation is done with floating
    #    /onedim: one dimensional coordinate vectors
    #    errormess: error message
    #    /silent: no printing (default: 0)
    #
    #  OUTPUT:
    #    Returns the undulation matrix.
    #
    #  Caveat:
    #    The code only works for equidistant grid. Irregular
    #    grids are averaged & the middle axis is taken as
    #    radial & vertical coordinates.
    """

    if r_vec is None and z_vec is None:
        if spatial_pos is None:
            raise ValueError('spatial_pos needs to be set.')

        nr=len(spatial_pos[0,:,0])
        nz=len(spatial_pos[:,0,0])
        if onedim:
            r=(spatial_pos[int(nz/2)-1,:,0]+spatial_pos[int(nz/2),:,0])/2.
            z=(spatial_pos[:,int(nr/2)-1,1]+spatial_pos[:,int(nr/2),1])/2.
        else:
            r=spatial_pos[:,:,0]
            z=spatial_pos[:,:,1]
    elif r_vec is not None and z_vec is not None:
        nr=len(r_vec)
        nz=len(z_vec)
        r=r_vec
        z=z_vec
    else:
        raise ValueError('Either both or neither r_vec and z_vec need to be set.')


    undulation=np.zeros([nz,nr,nz,nr])

    for i in range(1, nr-2+1):
        for j in range(1, nz-2+1):


            undulation[j,i,j,i-1]  =(z[j+1]-z[j-1])/3.*(1./(r[i-1]-r[i])) + \
                                    (r[i]-r[i-1])  /6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))

            undulation[j,i,j,i]    =(z[j+1]-z[j-1])/3.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + \
                                    (r[i+1]-r[i-1])/3.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))

            undulation[j,i,j,i+1]  =(z[j+1]-z[j-1])/3.*(1./(r[i]-r[i+1])) + \
                                    (r[i+1]-r[i])  /6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))

            undulation[j,i,j-1,i-1]=(z[j]-z[j-1])  /6.*(1./(r[i-1]-r[i])) + \
                                    (r[i]-r[i-1])  /6.*(1./(z[j-1]-z[j]))

            undulation[j,i,j-1,i]  =(z[j]-z[j-1])  /6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + \
                                    (r[i+1]-r[i-1])/3.*(1./(z[j-1]-z[j]))

            undulation[j,i,j-1,i+1]=(z[j]-z[j-1])  /6.*(1./(r[i]-r[i+1])) + \
                                    (r[i+1]-r[i])  /6.*(1./(z[j-1]-z[j]))

            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])  /6.*(1./(r[i-1]-r[i])) + \
                                    (r[i]-r[i-1])  /6.*(1./(z[j]-z[j+1]))

            undulation[j,i,j+1,i]  =(z[j+1]-z[j])  /6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + \
                                    (r[i+1]-r[i-1])/3.*(1./(z[j]-z[j+1]))

            undulation[j,i,j+1,i+1]=(z[j+1]-z[j])  /6.*(1./(r[i]-r[i+1])) + \
                                    (r[i+1]-r[i])  /6.*(1./(z[j]-z[j+1]))

    if floating :
        for i in range(0, nr):
            for j in range(0, nz):

                if j == 0 :
                    if i == 0 :
                        undulation[j,i,j,i]=    (z[j+1]-z[j])/6.*(1./(r[i+1]-r[i])) + \
                                                (r[i+1]-r[i])/6.*(1./(z[j+1]-z[j]))
                        undulation[j,i,j,i+1]=  (z[j+1]-z[j])/6.*(1./(r[i+1]-r[i])) + \
                                                (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))
                        undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + \
                                                (r[i+1]-r[i])/6.*(1./(z[j+1]-z[j]))
                        undulation[j,i,j+1,i+1]=(z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + \
                                                (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))
                    else:
                        if i == nr-1 :
                            undulation[j,i,j,i-1]=  (z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j+1]-z[j]))
                            undulation[j,i,j,i]=    (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j+1]-z[j])) #corr
                            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
                            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
                        else:
                            undulation[j,i,j,i-1]=  (z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j+1]-z[j]))
                            undulation[j,i,j,i]=    (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + \
                                                    (r[i+1]-r[i-1])/3.*(1./(z[j+1]-z[j]))
                            undulation[j,i,j,i+1]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j+1]-z[j]))
                            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
                            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
                            undulation[j,i,j+1,i+1]=(z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))


                else:
                    if j == nz-1 :
                        if i == 0 :
                            undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i+1]-r[i])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))
                            undulation[j,i,j-1,i+1]=(z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))
                            undulation[j,i,j,i]=    (z[j]-z[j-1])/6.*(1./(r[i+1]-r[i])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1]))
                            undulation[j,i,j,i+1]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1]))
                        else:
                            if i == nr-1 :
                                undulation[j,i,j-1,i-1]=(z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + \
                                                        (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
                                undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])) + \
                                                        (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
                                undulation[j,i,j,i-1]=  (z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + \
                                                        (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1]))
                                undulation[j,i,j,i]=    (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])) + \
                                                        (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1])) #corr
                            else:
                                undulation[j,i,j-1,i-1]=(z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + \
                                                        (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
                                undulation[j,i,j-1,i]=  (z[j]-z[j-1])  /6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + \
                                                        (r[i+1]-r[i-1])/3.*(1./(z[j-1]-z[j]))

                                undulation[j,i,j-1,i+1]=(z[j]-z[j-1])  /6.*(1./(r[i]-r[i+1])) + \
                                                        (r[i+1]-r[i])  /6.*(1./(z[j-1]-z[j]))
                                undulation[j,i,j,i-1]=  (z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + \
                                                        (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1]))
                                undulation[j,i,j,i]=    (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + \
                                                        (r[i+1]-r[i-1])/3.*(1./(z[j]-z[j-1])) #corr
                                undulation[j,i,j,i+1]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + \
                                                        (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1]))
                    else:
                        if i == 0 :
                            undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i+1]-r[i])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))
                            undulation[j,i,j-1,i+1]=(z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))
                            undulation[j,i,j,i]=    (z[j+1]-z[j-1])/3.*(1./(r[i+1]-r[i])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            undulation[j,i,j,i+1]=  (z[j+1]-z[j-1])/3.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i+1]-r[i])) + \
                                                    (r[i+1]-r[i])/3.*(1./(z[j]-z[j+1]))
                            undulation[j,i,j+1,i+1]=(z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + \
                                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))

                        if i == nr-1 :
                            undulation[j,i,j-1,i-1]=(z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
                            undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
                            undulation[j,i,j,i-1]=  (z[j+1]-z[j-1])/3.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            undulation[j,i,j,i]=    (z[j+1]-z[j-1])/3.*(1./(r[i]-r[i-1])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
                            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + \
                                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
    return undulation