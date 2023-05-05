pro test_undulation, nvec=nvec
  default, nvec, 100
  default, shot, 14110
  nr=16.
  nz=4.
  rres=10.
  zres=10.
  r_vec=findgen(nr)*rres ;radial coordinate in mm
  z_vec=findgen(nz)*zres ;horizontal coordinate in mm
  
  r_vec_interp=findgen(nr*nvec)*rres/nvec
  z_vec_interp=findgen(nz*nvec)*zres/nvec
  undulation=dblarr(16,4,16,4)
  diffmat=dblarr(nr*nvec,nz*nvec,4)
  if keyword_set(test) then begin
    for i1=1,nr-2 do begin
      for j1=1,nz-2 do begin
        for i2=1,nr-2 do begin
          for j2=1,nz-2 do begin
            diffmat[i1*nvec:(i1+1)*nvec-1,j1*nvec:(j1+1)*nvec-1,0] = base_function_diff(r_vec[i1-1:i1+1],z_vec[j1-1:j1+1], nvec=nvec, /rdiff)
            diffmat[i1*nvec:(i1+1)*nvec-1,j1*nvec:(j1+1)*nvec-1,1] = base_function_diff(r_vec[i1-1:i1+1],z_vec[j1-1:j1+1], nvec=nvec, /zdiff)
            diffmat[i2*nvec:(i2+1)*nvec-1,j2*nvec:(j2+1)*nvec-1,2] = base_function_diff(r_vec[i2-1:i2+1],z_vec[j2-1:j2+1], nvec=nvec, /rdiff)
            diffmat[i2*nvec:(i2+1)*nvec-1,j2*nvec:(j2+1)*nvec-1,3] = base_function_diff(r_vec[i2-1:i2+1],z_vec[j2-1:j2+1], nvec=nvec, /zdiff)
            undulation[i1,j1,i2,j2]=total(diffmat[*,*,0]*diffmat[*,*,2]+diffmat[*,*,1]*diffmat[*,*,3])*rres/nvec*zres/nvec
            if undulation[i1,j1,i2,j2] gt 0 then print, i1,j1,i2,j2, undulation[i1,j1,i2,j2]
            diffmat[*,*,*]=0
          endfor
        endfor
      endfor
    endfor
  endif else begin
    for i=1,nr-2 do begin
      for j=1,nz-2 do begin
        
        diffmat[i*nvec:(i+1)*nvec-1,j*nvec:(j+1)*nvec-1,0] = base_function_diff(r_vec[i-1:i+1],z_vec[j-1:j+1], nvec=nvec, /rdiff)
        diffmat[i*nvec:(i+1)*nvec-1,j*nvec:(j+1)*nvec-1,1] = base_function_diff(r_vec[i-1:i+1],z_vec[j-1:j+1], nvec=nvec, /zdiff)
        undulation[i1,j1,i2,j2]=total((diffmat[*,*,0]^2+diffmat[*,*,1]^2)*((r_vec[0:nr-3]-r_vec[2:nr-1])/2 # z_vec[0:2]-))
        
        if undulation[i1,j1,i2,j2] gt 0 then print, i1,j1,i2,j2, undulation[i1,j1,i2,j2]
        diffmat[*,*,*]=0
      endfor
    endfor
  endelse
stop
end

function base_function_diff, r, z, nvec=nvec, rdiff=rdiff, zdiff=zdiff
  default, nvec, 100
  if n_elements(r) ne 3 or n_elements(z) ne 3 then begin
    print, 'The boundary of the pyramid needs to be defined by 3-3 r z points.'
    return, 0
  endif

  r_vec=findgen(nvec)/nvec*(r[2]-r[0])+r[0]
  z_vec=findgen(nvec)/nvec*(z[2]-z[0])+z[0]

  base_diff=fltarr(nvec,nvec)
  for i_r=0,nvec-1 do begin
    for j_z=0,nvec-1 do begin
      
    if keyword_set(rdiff) then begin
      
      if i_r lt nvec/2 then begin
        if j_z lt nvec/2 then begin
          base_diff[i_r,j_z]=1./(r[1]-r[0])*(z_vec[j_z]-z[0])/(z[1]-z[0])
        endif else begin
          base_diff[i_r,j_z]=1./(r[1]-r[0])*(z_vec[j_z]-z[2])/(z[1]-z[2])
        endelse
      endif else begin
        if j_z lt nvec/2 then begin
          base_diff[i_r,j_z]=1./(r[1]-r[2])*(z_vec[j_z]-z[0])/(z[1]-z[0])
        endif else begin
          base_diff[i_r,j_z]=1./(r[1]-r[2])*(z_vec[j_z]-z[2])/(z[1]-z[2])
        endelse
      endelse
    endif
    
    if keyword_set(zdiff) then begin
      if i_r lt nvec/2 then begin
        if j_z lt nvec/2 then begin
          base_diff[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])/(z[1]-z[0])
        endif else begin
          base_diff[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])/(z[1]-z[2])
        endelse
      endif else begin
        if j_z lt nvec/2 then begin
          base_diff[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])/(z[1]-z[0])
        endif else begin
          base_diff[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])/(z[1]-z[2])
        endelse
      endelse
    endif
    
    endfor
  endfor
  return, base_diff
end