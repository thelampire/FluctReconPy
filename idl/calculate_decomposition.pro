function calculate_decomposition, light_profile=s_vector, error_profile=error_vector, m_matrix=m_matrix, $
                                  spatial_pos=spatial_pos, floating=floating, iterate=iterate, k_factor=k_factor, $
                                  visual=visual, contour=contour, ps=ps, nocalc=nocalc, save_file=save_file, $
                                  n_vector_test=n_vector, ksi=ksi

;*******************************************************
;***   calculate_decomposition       by M. Lampert   ***
;*******************************************************
;
; This procedure calculates the density matrix from 
; the measured BES light signal based on the arbitrary
; expansion method. The steps followed in the procedure
; are the following:
; 1: Arbitrary initial K value is considered.
; 2: Matrix N and L is calculated
; 3: Minimum point of X(P) is calculated
; 4: X is calculated
; 5: If X is not close enough to unity, then go back to
;    the second step
;
; INPUTs:
;   shot: shotnumber to be analyzed
;   time: time to be analyzed
;   t_win: time window around time, timerange: [time-t_win/2,time+t_win/2]
;   nocalibrate: No relative amplitude calibration
;   test: Test the code with dummy atomic physics
;   iterate: solve the equation with iteration
;   floating: the perimeter values of the undulation matrix are floating
;   For testing:
;   noise_level: the relative noise level of the signal
;   electronic_noise: electronic noise in mV for the calculated spectral range
;   blog_size: size of the simulated blob in mm
;   blob_pos: position of the blob in mm [R,z]
;   /visual: visualize the calculation
;   /contour: contour plot instead of simple plotting
; OUTPUTs:
;   return: calculated density fluctuation matrix
;   k_factor: smoothing factor for further calculation
;   n_vector: the original density distribution
;   
;

default, k_factor_first, 1e-10
default, k_factor_last, 1e3
default, nocalib, 1.
default, threshold, 0.001
default, maxiter, 100.
default, iterate, 0
default, floating, 1

increment=1.25
k_factor=k_factor_first
default, save_file, 'calculate_decomposition_save.sav'

if not keyword_set(nocalc) then begin
  if keyword_set(iterate) then begin
    k_factor_vector=[k_factor_first, k_factor_last]
    ksi_vector=dblarr(2)
    ksi_prev=0.
  endif else begin
    ksi_vector=dblarr(maxiter)
    k_factor_vector=dblarr(maxiter)
  endelse
 
  size_m_mat=size(m_matrix)
  
  n_vert_calc=double(size_m_mat[1])
  n_rad_calc=double(size_m_mat[2])
  n_vert_meas=double(size_m_mat[3])
  n_rad_meas=double(size_m_mat[4])
 
  s_vector_ref=reform((s_vector),n_vert_meas*n_rad_meas)
  error_vector_s_ref=reform((error_vector),n_vert_meas*n_rad_meas)
  m_matrix_ref = reform_matrix(m_matrix,[n_vert_calc*n_rad_calc,n_vert_meas*n_rad_meas])
  h_matrix=undulation_matrix(spatial_pos=spatial_pos, floating=floating) ;[n_meas_vert*n_meas_rad, n_calc_vert*n_calc_rad]
  h_matrix_ref = reform_matrix(h_matrix, [n_vert_calc*n_rad_calc, n_vert_meas*n_rad_meas])
  iter=1  
  repeat begin

    m_matrix_prime_ref=m_matrix_ref
    for i=0, n_rad_calc*n_vert_calc-1 do m_matrix_prime_ref[*,i]=m_matrix_ref[*,i]/error_vector_s_ref[*]
    
    p_vector = k_factor * invert(transpose(h_matrix_ref) + k_factor * transpose(m_matrix_prime_ref) ## m_matrix_ref, status) ## (transpose(m_matrix_prime_ref) ## s_vector_ref)
    
    if status ne 0 then begin
      print, 'The inversion is invalid. Matrix is singular.'
    endif
    
    ;Calculate ksi
    ksi=n_vert_meas^(-1)*n_rad_meas^(-1)*total((reform(s_vector_ref) - m_matrix_ref ## reform(p_vector))^2/reform(error_vector_s_ref))
    if ksi lt 0 then begin
      print, 'Ksi is lower than 0, iteration is aborted.'
      break
    endif
    if keyword_set(iterate) then begin
      if iter eq 1 then begin
        ksi_vector[0]=ksi
        k_factor=k_factor_vector[1]
      endif
      if iter eq 2 then begin
        ksi_vector[1]=ksi
        k_factor=(k_factor_vector[0]+k_factor_vector[1])/2.
      endif
      if iter ge 3 then begin
        if ksi lt 1 then begin
          ksi_vector[1]=ksi
          k_factor_vector[1]=k_factor
        endif else begin
          ksi_vector[0]=ksi
          k_factor_vector[0]=k_factor
        endelse
        k_factor=(k_factor_vector[0]+k_factor_vector[1])/2.
      endif
    endif else begin
      k_factor_vector[iter-1]=k_factor
      k_factor*=increment
      ksi_vector[iter-1]=ksi
    endelse
    iter+=1
  
    if iter eq maxiter+1 then break
  endrep until (abs(ksi - 1) lt threshold)
  
  save, /variables, filename=save_file
endif else begin
  if file_test(save_file) then begin
    restore, save_file 
  endif else begin
    print, 'The save file doesnt exist, please run the procedure without /nocalc'
    return, -1
  endelse
endelse
p_vector= transpose(reform(p_vector,n_rad_calc,n_vert_calc))
return, p_vector

end