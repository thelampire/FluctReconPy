function base_function, r, z, nvec=nvec, r_interp=r_interp, z_interp=z_interp

if n_elements(r) ne 3 or n_elements(z) ne 3 then begin
  print, 'The boundary of the pyramid needs to be defined by 3-3 r z points.'
  return, 0
endif

if not defined(r_interp) and not defined(z_interp) then begin
  default, nvec, 100
  r_vec=findgen(nvec)/nvec*(r[2]-r[0])+r[0]
  z_vec=findgen(nvec)/nvec*(z[2]-z[0])+z[0]
  base=fltarr(nvec,nvec)
  for i_r=0,nvec-1 do begin
    for j_z=0,nvec-1 do begin
      if i_r lt nvec/2 then begin
        if j_z lt nvec/2 then begin
          base[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])*(z_vec[j_z]-z[0])/(z[1]-z[0])
        endif else begin
          base[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])*(z_vec[j_z]-z[2])/(z[1]-z[2])
        endelse
      endif else begin
        if j_z lt nvec/2 then begin
          base[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])*(z_vec[j_z]-z[0])/(z[1]-z[0])
        endif else begin
          base[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])*(z_vec[j_z]-z[2])/(z[1]-z[2])
        endelse
      endelse
    endfor
  endfor
  return, base
endif else begin
  if r_interp lt r[0] or r_interp gt r[2] or z_interp lt z[0] or z_interp gt z[2] then begin
    print, 'r_interp and z_interp should be between the given r,z boundaries. Returning -1...'
    return, -1.
  endif
  if r_interp lt r[1] then begin
    if z_interp lt z[1] then begin
      ret_val=(r_interp-r[0])/(r[1]-r[0])*(z_interp-z[0])/(z[1]-z[0])
    endif else begin
      ret_val=(r_interp-r[0])/(r[1]-r[0])*(z_interp-z[2])/(z[1]-z[2])
    endelse
  endif else begin
    if z_interp lt z[1] then begin
      ret_val=(r_interp-r[2])/(r[1]-r[2])*(z_interp-z[0])/(z[1]-z[0])
    endif else begin
      ret_val=(r_interp-r[2])/(r[1]-r[2])*(z_interp-z[2])/(z[1]-z[2])
    endelse
  endelse
  return, ret_val
endelse

end