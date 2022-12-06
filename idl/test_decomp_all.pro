pro test_decomp_all, recalc=recalc, nwin=nwin, stacionery=stacionery

default, stacionery, 1
default, nwin, 16
default, electronic_noise, 0.001

noise_vector=findgen(10.)/20. ;0-50%
spatial_resolution=[10] ;Only the relative resolution counts
blob_size=(findgen(10)+1.)*0.5 ;times the spatial resolution

time_win=100e-6 ;100 samples
sampling_time=0.5e-6
n=time_win/sampling_time          


pos_orig=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n,2)
pos_calc=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n,2)
fwhm_orig=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n,2)
fwhm_calc=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n,2)
density_calc=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n)
density_samp=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n)
density_orig=dblarr(n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution),n)

time0=systime(/seconds)
for i_noise=0,n_elements(noise_vector)-1 do begin
   for i_size=0,n_elements(blob_size)-1 do begin
      for i_res=0,n_elements(spatial_resolution)-1 do begin
        time_left, indices=[i_noise,i_size,i_res], maxind=[n_elements(noise_vector),n_elements(blob_size),n_elements(spatial_resolution)], starttime=starttime 
        spatial_pos=rescale_spat_pos(nwin=nwin,rad_res=spatial_resolution[i_res], vert_res=spatial_resolution[i_res])
        if keyword_set(stacionery) then begin
          moving={vx:spatial_resolution[i_res]/1e3/time_win,vy:spatial_resolution[i_res]/1e3/time_win, v_size:0., rot_freq:0.} ;it only moves one spatial resolution step during the analysis
          blob_pos=[mean(spatial_pos[*,*,0]),mean(spatial_pos[*,*,1])]
        endif else begin
          moving={vx:150,vy:10, v_size:1., rot_freq:1.}
          blob_pos=[2120,20]
        endelse
        
        filename='test_fluct_decomp/test_fluct_decomp_N_'+strtrim(noise_vector[i_noise],2)+'_BS_'+strtrim(blob_size[i_size],2)+'_R_'+strtrim(spatial_resolution[i_res],2)+'.sav'
        
        if file_test(filename) and not keyword_set(recalc) then begin
           restore, filename
        endif else begin
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
        for k=0,n-1 do begin
          ;contour, n_vector_calc[*,*,k], /fill, nlevel=21, /iso, title='N_'+strtrim(noise_vector[i_noise],2)+'_BS_'+strtrim(blob_size[i_size],2)+'_R_'+strtrim(spatial_resolution[i_res],2)
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
           endfor
         endfor
      endfor
   endfor
   save, pos_orig, pos_calc, fwhm_orig, fwhm_calc, density_orig, density_samp, density_calc, filename="test_decomp_all_save_0-50noise_nores_dens.sav"

end
