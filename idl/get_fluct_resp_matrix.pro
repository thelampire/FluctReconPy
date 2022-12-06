function get_fluct_resp_matrix, timerange=timerange, test=test, reform=reform, spatial_pos=spatial_pos

default, test, 1
if keyword_set(test) then begin
  if not defined(spatial_pos) then begin
    errormess='The spatial coordinates need to be defined. Returning...'
    if not keyword_set(errormess) then print, errormess
    return, -1
  endif
  ;The matrix itself is a Lower triangular matrix
  nr=n_elements(spatial_pos[0,*,0])
  nz=n_elements(spatial_pos[*,0,0])
  tauvec=dblarr(nr,nz)
  matrix=dblarr(nz,nr,nz,nr)
  energy=90e3 ;eV
  e=1.6e-19
  mh=1.6727e-27
  v=sqrt(2*energy*e/mh)
  ;tauvec=[10.3,10.2,10.1,10.,9.,8.,7.,6.,5.,4.,3.,2.6,2.3,2.,1.5,1.]*1e-9 ;ns
  for i_rad=0,nr-1 do begin
    for j_vert=0,nz-1 do begin
      tauvec[i_rad,j_vert]=(tanh((spatial_pos[j_vert,nr-1,0]-spatial_pos[j_vert,i_rad,0])/20.-2.)+1.)*5.*1e-9
    endfor
  endfor

  for k=0,nz-1 do begin
    for j=0,nr-1 do begin
      for i=j,nr-1 do begin
        matrix[k,i,k,j]=1e-9/tauvec[i,k]*exp(-((i-j)*0.01)/(v*tauvec[i,k]))
      endfor
    endfor
  endfor
  return, matrix
endif else begin
  
endelse

end