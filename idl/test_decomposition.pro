pro test_decomposition, n_noise=n_noise, n_error=n_error, test_blob_position=test_blob_position,$
                        test_blob_size=test_blob_size

default, n_noise, 100.
noise_vec=dindgen(n_noise)/n_noise

if keyword_set(test_density) then begin

  default, n_error, 20.
  error_mat=dblarr(n_noise,n_error)
  for i=0,n_noise-1 do begin
    for j=0,n_error-1 do begin
      print, (i*n_error+j)/(n_error*n_noise)
      n_vec_calc=calculate_decomposition(14110, time=2.5, /nocalib, floating=0, blob_pos=[2200.,20.], blob_size=[20.,30.],$
                                         electronic_noise=0.003/25., noise_level=noise_vec[i], /iterate, n_vector=n_vec_sim,$
                                         spatial_pos=spatial_pos)
                                         
      error_mat[i,j]=sqrt(mean(1/double(n_elements(n_vec_calc))*((n_vec_calc-n_vec_sim))^2))/mean(n_vec_sim)
    endfor
  endfor
  plot, noise_vec,total(error_mat,2)/n_error
  error_vector=dblarr(n_noise)
  for i=0,n_noise-1 do error_vector[i]=sqrt(variance(error_mat[i,*]))
  errplot, noise_vec,total(error_mat,2)/n_error-error_vector, total(error_mat,2)/n_error+error_vector
endif

if keyword_set(test_blob_position) then begin
  default, n_error, 20.
  default, blob_pos, [2200.,20.]
  default, blob_size, [20.,30.]
  position_mat=dblarr(n_noise,n_error,2)
  size_mat=dblarr(n_noise,n_error,2)
  for i=0,n_noise-1 do begin
    for j=0,n_error-1 do begin
      print, (i*n_error+j)/(n_error*n_noise)
      n_vec_calc=calculate_decomposition(14110, time=2.5, /nocalib, floating=0, blob_pos=blob_pos, blob_size=blob_size,$
                                         electronic_noise=0.003/25., noise_level=noise_vec[i], /iterate, n_vector=n_vec_sim,$
                                          spatial_pos=spatial_pos)
        a=max(n_vec_calc, ind_max)
        ind=array_indices(n_vec_calc,ind_max)
      ;position_mat[i,j,*]=sqrt((spatial_pos[ind[0],ind[1],0]-blob_pos[0])^2+(spatial_pos[ind[0],ind[1],1]-blob_pos[1])^2) ;maximum intensity calculation
      position_mat[i,j,*]=abs([total(n_vec_calc*reform(spatial_pos[*,*,0]))/total(n_vec_calc),total(n_vec_calc*reform(spatial_pos[*,*,1]))/total(n_vec_calc)]-blob_pos) ;mass centre calculation
      size_mat[i,j,0]=sqrt(total(n_vec_calc*reform(spatial_pos[*,*,0]-position_mat[i,j,0]-blob_pos[0])^2)/total(n_vec_calc))-blob_size[0]
      size_mat[i,j,1]=sqrt(total(n_vec_calc*reform(spatial_pos[*,*,1]-position_mat[i,j,1]-blob_pos[1])^2)/total(n_vec_calc))-blob_size[1]
    endfor
  endfor
  plot, noise_vec,total(position_mat,2)/n_error
  error_vector=dblarr(n_noise)
  for i=0,n_noise-1 do error_vector[i]=sqrt(variance(position_mat[i,*,0]+variance(position_mat[i,*,1])))
  errplot, noise_vec,total(position_mat,2)/n_error-error_vector, total(position_mat,2)/n_error+error_vector
endif
end