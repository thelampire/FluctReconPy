pro test_fluct_decomp, nocalibrate=nocalib, k_factor=k_factor, $
                       test=test, iterate=iterate, floating=floating,$
                       noise_level=noise_level, electronic_noise=electronic_noise_in, $
                       blob_size=blob_size, blob_pos=blob_pos, visual=visual, contour=contour, $
                       spatial_pos=spatial_pos, ps=ps, nocalc=nocalc, save_file=save_file,$
                       maxiter=maxiter, moving=moving, hole_fluct_amp=hole_fluct_amp,$
                       n_vector_calc=n_vector_calc, n_vector_orig=n_vector, $
                       time_vec=time_vec, time_win=time_win,sampling_time=sampling_time,$
                       blob_density_orig=blob_density_orig, blob_density_samp=blob_density_samp
 
                       
  
  
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
  default, blob_pos, [2200,40]
  default, blob_size, [10,10]
  default, noise_level, 0.01 ;%
  default, electronic_noise_in, 0.003 ;mV
  default, time_win, 50e-6
  default, hole_fluct_amp, 0.
  default, sampling_time, 0.5e-6
  
  electronic_noise=electronic_noise_in/sqrt(20) ;1Mhz bandwidth --> 50kHz bandwidth for blobs, approx. power decrease
  ;Fill up of the error_vector

  s_vector=dblarr(n_vert_meas,n_rad_meas)
  n_vector=dblarr(n_vert_calc,n_rad_calc)
  n0=1. ;Characteristic density
  error_vector=dblarr(n_vert_meas,n_rad_meas) ;square of the standard error
  if not keyword_set(moving) then begin
    m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos, test=test) ;[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
    m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])
    n_vector=fluctuation_matrix_test(spatial_pos=spatial_pos, rad_pos=blob_pos[0], vert_pos=blob_pos[1], rad_size=blob_size[0],vert_size=blob_size[1], blob_density=blob_density_orig)
    n_vector_ref = reform(transpose(n_vector), n_rad_meas*n_vert_meas)
    blob_density_samp=total(n_vector_ref)/(!pi*blob_size[0]*blob_size[1])
    s_vector_ref = m_matrix_ref ## n_vector_ref + (max(abs(m_matrix_ref ## n_vector_ref))*noise_level+electronic_noise)*randomn(seed,64) ;s_vector is a simulated measurement
    error_vector_s_ref=abs((m_matrix_ref ## n_vector_ref)*noise_level+electronic_noise)^2
    s_vector=reform(s_vector_ref,n_vert_meas, n_rad_meas)
    error_vector=transpose(reform(error_vector_s_ref,n_vert_meas, n_rad_meas))
    n_vector_calc=calculate_decomposition(light_profile=s_vector, error_profile=error_vector, m_matrix=m_matrix, spatial_pos=spatial_pos, $
                                          floating=floating, iterate=iterate, k_factor=k_factor, visual=visual, $
                                          contour=contour, ps=ps, nocalc=nocalc, save_file=save_file, n_vector_test=n_vector)
    
  endif else begin
    m_matrix=get_fluct_resp_matrix(spatial_pos=spatial_pos, test=test) ;[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
    m_matrix_ref = reform_matrix(m_matrix,[n_rad_calc*n_vert_calc,n_rad_meas*n_vert_meas])
    n_vector=fluctuation_matrix_test(rad_pos=blob_pos[0],vert_pos=blob_pos[1], /moving, vert_size=blob_size[1], rad_size=blob_size[0], $
                                     vx=moving.vx, vy=moving.vy, v_size=moving.v_size, rot_freq=moving.rot_freq, time_vec=time_vec,$
                                     time_win=time_win,hole_fluct_amp=hole_fluct_amp, spatial_pos=spatial_pos, sampling_time=sampling_time,$
                                     blob_density=blob_density_orig)
    n_vector_calc=n_vector
    loadct, 5
    blob_density_samp=dblarr(n_elements(time_vec))
    for i=0,n_elements(time_vec)-1 do begin
      value=max(abs(m_matrix_ref ## reform(transpose(reform(n_vector[*,*,0])), n_rad_meas*n_vert_meas)))
      noise=(value*noise_level+electronic_noise)*randomn(seed,n_rad_meas*n_vert_meas)
      n_vector_ref = reform(transpose(reform(n_vector[*,*,i])), n_rad_meas*n_vert_meas)
      blob_density_samp[i]=total(n_vector_ref)/(!pi*blob_size[0]*blob_size[1])
      s_vector_ref = m_matrix_ref ## n_vector_ref + noise ;s_vector is a simulated measurement
      error_vector_s_ref=abs((m_matrix_ref ## n_vector_ref)*noise_level+electronic_noise)^2
      s_vector=reform(s_vector_ref,n_vert_meas, n_rad_meas)
      error_vector=transpose(reform(error_vector_s_ref,n_vert_meas, n_rad_meas))
      n_vector_calc[*,*,i]=calculate_decomposition(light_profile=s_vector, error_profile=error_vector, m_matrix=m_matrix, spatial_pos=spatial_pos, $
                                                   floating=floating, iterate=iterate, k_factor=k_factor, $
                                                   contour=contour, ps=ps, nocalc=nocalc, save_file=save_file, n_vector_test=n_vector[*,*,i])
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
         contour, transpose(n_vector[*,*,i]), r_vec,z_vec, /fill, nlevel=21, position=[0.05,0.8,0.45,0.95], $
           title="Original density", zrange=range, /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso
         contour, transpose(n_vector_calc[*,*,i]), r_vec,z_vec, /fill, nlevel=21, position=[0.05,0.55,0.45,0.7], /noerase, $
           title="Reconstructed density", zrange=range, /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso
         contour, reform(s_vector_ref, n_rad_meas,n_vert_meas), r_vec,z_vec,/fill, nlevel=21, position=[0.5,0.8,0.95,0.95], /noerase,$
           title="Original light with noise", /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso, zrange=[zmin_light,zmax_light]
         contour, reform(m_matrix_ref ## p_vector,n_rad_meas,n_vert_meas), r_vec,z_vec, /fill, nlevel=21, position=[0.5,0.55,0.95,0.7], /noerase, $
           title="Reconstructed light", /xstyle, /ystyle, xtitle="R [mm]", ytitle="z [mm]", thick=3, /iso, zrange=[zmin_light,zmax_light]
         plot, r_vec, n_vector[1,*,i], /xstyle, /ystyle, xtitle="R [mm]", thick=3, position=[0.05,0.3,0.95,0.45], /noerase
         plot, r_vec, n_vector_calc[1,*,i], /xstyle, /ystyle, xtitle="R [mm]", thick=3, position=[0.05,0.05,0.95,0.2], /noerase
         ;print, total((n_vector_calc[*,*,i]-n_vector[*,*,i])^2)/n_vert_calc/n_rad_calc

      endif
    endfor
  endelse
end
