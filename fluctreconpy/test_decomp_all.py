import numpy
pro test_decomp_all, recalc=recalc, nwin=nwin, stacionery=stacionery

default, stacionery, 1
default, nwin, 16
default, electronic_noise, 0.001

noise_vector=numpy.arange(10.)/20. #0-50%
spatial_resolution=[10] #Only the relative resolution counts
blob_size=(numpy.arange(10)+1.)*0.5 #times the spatial resolution

time_win=100e-6 #100 samples
sampling_time=0.5e-6
n=time_win/sampling_time          


pos_orig=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n,2)
pos_calc=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n,2)
fwhm_orig=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n,2)
fwhm_calc=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n,2)
density_calc=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n)
density_samp=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n)
density_orig=numpy.zeroes(len(noise_vector),len(blob_size),len(spatial_resolution),n)

time0=systime(/seconds)
for i_noise in range(0, len(noise_vector)):
    for i_size in range(0, len(blob_size)):
        for i_res in range(0, len(spatial_resolution)):
            time_left, indices=[i_noise,i_size,i_res], maxind=[len(noise_vector),len(blob_size),len(spatial_resolution)], starttime=starttime 
            spatial_pos=rescale_spat_pos(nwin=nwin,rad_res=spatial_resolution[i_res], vert_res=spatial_resolution[i_res])
            if (stacionery is None) :
                moving={vx:spatial_resolution[i_res]/1e3/time_win,vy:spatial_resolution[i_res]/1e3/time_win, v_size:0., rot_freq:0.} #it only moves one spatial resolution step during the analysis
                blob_pos=[mean(spatial_pos[*,*,0]),mean(spatial_pos[*,*,1])]
                     else begin
                        moving={vx:150,vy:10, v_size:1., rot_freq:1.}
                        blob_pos=[2120,20]
                            endelse
        
                            filename='test_fluct_decomp/test_fluct_decomp_N_'+strtrim(noise_vector[i_noise],2)+'_BS_'+strtrim(blob_size[i_size],2)+'_R_'+strtrim(spatial_resolution[i_res],2)+'.sav'
        
                            if file_test(filename) & (recalc is not None) :
                                restore, filename
                                 else begin
                                    test_fluct_decomp, /nocalib, k_factor=k_factor,$
                                        test=test, /iterate,$
                                        noise_level=noise_vector[i_noise], electronic_noise=electronic_noise, $
                                        blob_size=[blob_size[i_size]*spatial_resolution[i_res],blob_size[i_size]*spatial_resolution[i_res]], $
                                        blob_pos=blob_pos,$
                                        moving=moving,$
                                        n_vector_calc=n_vector_calc, n_vector_orig=n_vector_orig, $
                                        time_vec=time_vec, time_win=time_win,$
                                        spatial_pos=spatial_pos,sampling_time=sampling_time,$
                                        blob_density_orig=blob_density_orig,blob_density_samp=blob_density_samp

                                            save, n_vector_calc, n_vector_orig, time_vec, spatial_pos, blob_density_orig, blob_density_samp, filename=filename
                                endelse
                                density_samp[i_noise,i_size,i_res,*]=blob_density_samp
                                density_orig[i_noise,i_size,i_res,*]=blob_density_orig
                                for k in range(0, n):
                                    #contour, n_vector_calc[*,*,k], /fill, nlevel=21, /iso, title='N_'+strtrim(noise_vector[i_noise],2)+'_BS_'+strtrim(blob_size[i_size],2)+'_R_'+strtrim(spatial_resolution[i_res],2)
                                        pos_orig[i_noise,i_size,i_res,k,0]=blob_pos[0]+moving.vx*time_vec[k]
                                        pos_calc[i_noise,i_size,i_res,k,0]=total(n_vector_calc[*,*,k]*spatial_pos[*,*,0])/total(n_vector_calc[*,*,k])
                                        pos_orig[i_noise,i_size,i_res,k,1]=blob_pos[1]+moving.vy*time_vec[k]
                                        pos_calc[i_noise,i_size,i_res,k,1]=total(n_vector_calc[*,*,k]*spatial_pos[*,*,1])/total(n_vector_calc[*,*,k])
           
                                        fwhm_orig[i_noise,i_size,i_res,k,0]=blob_size[i_size]*spatial_resolution[i_res]
                                        fwhm_orig[i_noise,i_size,i_res,k,1]=blob_size[i_size]*spatial_resolution[i_res]

                                        spat_rad=reform(spatial_pos[0,*,0])
                                        n_vector_calc_rad=reform(total(n_vector_calc[*,*,k],1))
                                        a=gaussfit(spat_rad,n_vector_calc_rad, param3, nterms=3)
                                        fwhm_calc[i_noise,i_size,i_res,k,0]=param3[2]

                                        spat_vert=reform(spatial_pos[*,0,1])
                                        n_vector_calc_vert=reform(total(n_vector_calc[*,*,k],2))
                                        a=gaussfit(spat_vert,n_vector_calc_vert, param4, nterms=3)
                                        fwhm_calc[i_noise,i_size,i_res,k,1]=param4[2]
                                        density_calc[i_noise,i_size,i_res,k]=total(n_vector_calc[*,*,k])/(!pi*param3[2]*param4[2])
                                        
                                    
                        
            
            save, pos_orig, pos_calc, fwhm_orig, fwhm_calc, density_orig, density_samp, density_calc, filename="test_decomp_all_save_0-50noise_nores_dens.sav"

end
