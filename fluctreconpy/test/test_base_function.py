import numpy as np
import matplotlib.pyplot as plt

def test_base_function(shot, timerange=None):

    spat_cal=None #TODO: getcal_kstar_spat(shot)
    nz=4
    nr=16
    z_inter=4
    r_inter=4

    r_vec=(spat_cal[1,:,0]+spat_cal[2,:,0])/2
    z_vec=(spat_cal[:,7,1]+spat_cal[:,8,1])/2
    si=np.zeros([(nz-1)*z_inter+1,(nr-1)*r_inter+1,2])
    si[:,:,0]=np.interpolate(spat_cal[:,:,0],
                             (np.arange((nz-1)*z_inter+1))/z_inter,
                             (np.arange((nr-1)*r_inter+1))/r_inter,
                             )
    si[:,:,1]=np.interpolate(spat_cal[:,:,1],
                             (np.arange((nz-1)*z_inter+1))/z_inter,
                             (np.arange((nr-1)*r_inter+1))/r_inter,
                             )

    r_vec_inter=np.interpolate(r_vec, (np.arange((nr-1)*r_inter+1))/r_inter)
    z_vec_inter=np.interpolate(z_vec, (np.arange((nz-1)*z_inter+1))/z_inter)

    signal=np.zeros([nz,nr])
    signal_inter=np.zeros([(nz-1)*z_inter+1,(nr-1)*r_inter+1])

    for i in range(0, nz):
        for j in range(0, nr):
            #get_rawsignal, shot, 'BES-'+strtrim(i+1,2)+'-'+strtrim(j+1,2), t, d, timerange=timerange, /nocalib
            d=None #TODO
            signal[i,j]=np.mean(d)



    for i in range(0, (nz-1)*z_inter):
        for j in range(0, (nr-1)*r_inter):
            ir=j/r_inter
            iz=i/z_inter
            p1=(r_vec_inter[j]-r_vec[ir+1])/(r_vec[ir]-r_vec[ir+1]) * (z_vec_inter[i]-z_vec[iz+1])/(z_vec[iz]-z_vec[iz+1])
            p2=(r_vec_inter[j]-r_vec[ir])  /(r_vec[ir+1]-r_vec[ir]) * (z_vec_inter[i]-z_vec[iz])/  (z_vec[iz+1]-z_vec[iz])
            p3=(r_vec_inter[j]-r_vec[ir+1])/(r_vec[ir]-r_vec[ir+1]) * (z_vec_inter[i]-z_vec[iz])/  (z_vec[iz+1]-z_vec[iz])
            p4=(r_vec_inter[j]-r_vec[ir])  /(r_vec[ir+1]-r_vec[ir]) * (z_vec_inter[i]-z_vec[iz+1])/(z_vec[iz]-z_vec[iz+1])

            signal_inter[i,j]=(p1 * signal[iz,  ir]+
                               p2 * signal[iz+1,ir+1]+
                               p3 * signal[iz+1,ir]+
                               p4 * signal[iz,  ir+1])



    for i in range(0, (nz-1)*z_inter):
        j=(nr-1)*r_inter
        ir=j/r_inter
        iz=i/z_inter
        p1=(z_vec_inter[i]-z_vec[iz+1])/(z_vec[iz]-z_vec[iz+1])
        p3=(z_vec_inter[i]-z_vec[iz])/  (z_vec[iz+1]-z_vec[iz])
        signal_inter[i,j]=p1 * signal[iz,  ir] + p3 * signal[iz+1,ir]


    for j in range(0, (nr-1)*r_inter):
        i=(nz-1)*z_inter
        ir=j/r_inter
        iz=i/z_inter
        p1=(r_vec_inter[j]-r_vec[ir+1])  /(r_vec[ir]-r_vec[ir+1])
        p4=(r_vec_inter[j]-r_vec[ir])  /(r_vec[ir+1]-r_vec[ir])
        signal_inter[i,j]=p1 * signal[iz,  ir] + p4 * signal[iz,ir+1]


    signal_inter[(nz-1)*z_inter,(nr-1)*r_inter]=signal[(nz-1),(nr-1)]


    plt.cla()
    nlevel=21
    min_range=np.min([min(signal),
                min(signal_inter)])
    max_range=np.max([max(signal),
                max(signal_inter)])
    levels=np.arange(nlevel)/nlevel*(max_range-min_range)+min_range
    plt.contourf(signal.T,
                 spat_cal[:,:,0].T,
                 spat_cal[:,:,1],
                 nlevel=nlevel,
                 levels=levels)
    #TODO: plots, transpose(reform(s[*,*,0])),transpose(reform(s[*,*,1])), psym=4, color=120
    plt.contourf(signal_inter.T,
                 si[:,:,0].T,
                 si[:,:,1],
                 nlevel=nlevel,
                 levels=levels)
    #TODO: plots, transpose(reform(s[*,*,0])),transpose(reform(s[*,*,1])), psym=4, color=120