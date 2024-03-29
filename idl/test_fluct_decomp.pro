pro test_fluct_decomp, nocalibrate=nocalib, k_factor=k_factor, $
                       test=test, iterate=iterate, floating=floating,$
                       noise_level=noise_level, electronic_noise=electronic_noise_in, $
                       blob_size=blob_size, blob_pos=blob_pos, visual=visual, contour=contour, $
                       spatial_pos=spatial_pos, ps=ps, nocalc=nocalc, save_file=save_file,$
                       maxiter=maxiter, moving=moving, hole_fluct_amp=hole_fluct_amp,$
                       n_vector_calc=n_vector_calc, n_vector_orig=n_vector, $
                       time_vec=time_vec, time_win=time_win,sampling_time=sampling_time,$
                       blob_density_orig=blob_density_orig, blob_density_samp=blob_density_samp,$
                       decompose=decompose
 
                       
  
  
  default, shot, 14110
  default, spatial_pos, getcal_kstar_spat(shot)
  default, n_rad_meas, n_elements(spatial_pos[0,*,0])
  default, n_vert_meas, n_elements(spatial_pos[*,0,0])
  default, n_rad_calc, n_elements(spatial_pos[0,*,0])
  default, n_vert_calc, n_elements(spatial_pos[*,0,0])
  
  default, nocalib, 1.
  default, threshold, 0.001
  default, test, 1
  default, maxiter, 100.
  default, iterate, 1
  default, floating, 1
  default, decompose, 1
    
  default, fluct_amp, 1 ;*1e19
  default, blob_pos, [2200,40]
  default, blob_size, [10,10]
  default, noise_level, 0.01 ;% SNR=10
  default, electronic_noise_in, 0.003 ;V (sort of the KSTAR level)
  default, time_win, 50e-6
  default, hole_fluct_amp, 0.
  default, sampling_time, 0.5e-6    ;KSTAR level

  
  electronic_noise=electronic_noise_in/sqrt(20) ;1Mhz bandwidth --> 50kHz bandwidth for blobs, approx. power decrease
  ;Fill up of the error_vector

  s_vector=dblarr(n_vert_meas,n_rad_meas)
  n_vector=dblarr(n_vert_calc,n_rad_calc)
  n0=1. ;Characteristic density
  error_vector=dblarr(n_vert_meas,n_rad_meas) ;square of the standard error

  m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos, test=test) ;[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
  m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])
  n_vector=fluctuation_matrix_test(rad_pos=blob_pos[0],vert_pos=blob_pos[1], /moving, vert_size=blob_size[1], rad_size=blob_size[0], $
                                   vx=moving.vx, vy=moving.vy, v_size=moving.v_size, rot_freq=moving.rot_freq, time_vec=time_vec,$
                                   time_win=time_win,hole_fluct_amp=hole_fluct_amp, spatial_pos=spatial_pos, sampling_time=sampling_time,$
                                   blob_density=blob_density_orig, fluct_amp=fluct_amp)
  n_vector_calc=n_vector
  n_vector_calc[*,*,*]=0.
  loadct, 5       
  blob_density_samp=dblarr(n_elements(time_vec))
  for i=0,n_elements(time_vec)-1 do begin
      n_vector_ref = reform(transpose(reform(n_vector[*,*,i])), n_rad_meas*n_vert_meas)
      s_vector_no_noise=reform(m_matrix_ref ## n_vector_ref, n_rad_meas*n_vert_meas)
      
      ;randomn is mean=0 sigma=1 gaussian
      ;noise ref=
      sigma_vector_s_ref=abs(s_vector_no_noise*noise_level+electronic_noise) ;Ideally this should be estimated from data
      
      noise_ref=sigma_vector_s_ref*randomn(seed,n_rad_meas*n_vert_meas)
      
      blob_density_samp[i]=total(n_vector_ref)*(!pi*blob_size[0]*blob_size[1]) ;more like number of particles, not density
      s_vector_ref = abs(s_vector_no_noise + noise_ref) ;s_vector is a simulated measurement
      s_vector=reform(s_vector_ref, n_vert_meas, n_rad_meas)
      
      n_vector_calc[*,*,i]=calculate_decomposition(light_profile=s_vector, sigma_profile=sigma_vector_s_ref, m_matrix=m_matrix, spatial_pos=spatial_pos, $
                                                   floating=floating, iterate=iterate, k_factor=k_factor, $
                                                   contour=contour, ps=ps, nocalc=nocalc, save_file=save_file, $
                                                   n_vector_test=n_vector[*,*,i],$
                                                   decompose=decompose)
                                             
      if keyword_set(visual) then begin                                             
         if i eq 0 then range=range([n_vector,n_vector_calc])
         nr=n_elements(spatial_pos[0,*,0])
         nz=n_elements(spatial_pos[*,0,0])
         r_vec=reform((spatial_pos[nz/2-1,*,0]+spatial_pos[nz/2,*,0])/2)
         z_vec=reform((spatial_pos[*,nr/2-1,1]+spatial_pos[*,nr/2,1])/2)
         
         p_vector= reform(transpose(reform(n_vector_calc[*,*,i])), n_rad_calc*n_vert_calc)
         if i eq 0 then zmax_light=max([max(s_vector_ref),max(m_matrix_ref ## p_vector)])
         if i eq 0 then zmin_light=min([min(s_vector_ref),min(m_matrix_ref ## p_vector)])
         erase
         contour, transpose(n_vector[*,*,i]), r_vec, z_vec, /fill, nlevel=21, position=[0.05,0.8,0.45,0.95], $
           title="Original density", zrange=range, /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso
         contour, transpose(n_vector_calc[*,*,i]), r_vec,z_vec, /fill, nlevel=21, position=[0.05,0.55,0.45,0.7], /noerase, $
           title="Reconstructed density", zrange=range, /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso
         contour, reform(s_vector_ref, n_rad_meas,n_vert_meas), r_vec,z_vec,/fill, nlevel=21, position=[0.5,0.8,0.95,0.95], /noerase,$
           title="Original light with noise", /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso, zrange=[zmin_light,zmax_light]
         contour, reform(m_matrix_ref ## p_vector,n_rad_meas,n_vert_meas), r_vec,z_vec, /fill, nlevel=21, position=[0.5,0.55,0.95,0.7], /noerase, $
           title="Reconstructed light", /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso, zrange=[zmin_light,zmax_light]
         plot, r_vec, n_vector[3,*,i], /xstyle, /ystyle, xtitle="R [mm]", thick=3, position=[0.05,0.3,0.95,0.45], /noerase,$
           title="Original density profile"
         plot, r_vec, n_vector_calc[3,*,i], /xstyle, /ystyle, xtitle="R [mm]", thick=3, position=[0.05,0.05,0.95,0.2], /noerase,$
           title="Reconstructed density profile"
      endif
  endfor
end

