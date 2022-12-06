function fluctuation_matrix_test, rad_size=rad_size, vert_size=vert_size, rad_pos=rad_pos, vert_pos=vert_pos,$
                                  fluct_amp=fluct_amp, spatial_pos=spatial_pos, plot=plot, $
                                  moving=moving, vx=vx, vy=vy, v_size=v_size, rot_freq=rot_freq,$
                                  hole_size=hole_size, hole_rad_pos=hole_rad_pos, hole_vert_pos=hole_vert_pos,$ 
                                  hole_fluct_amp=hole_fluct_amp,time_win=time,time_vec=time_vec,$
                                  sampling_time=sampling_time, blob_density=blob_density
                                  
  ;************************************************
  ;*** fluctuation_matrix_test    by M. Lampert    
  ;************************************************
  ;
  ;Returns a simulated Gaussian blob with a given
  ;radial and vertical size and a given radial and
  ;vertical position. The amplitude of the blob
  ;can be also set.
  ;
  ;************************************************
  ;
  ;INPUTs:
  ;  shot: shotnumber for spatial calibration
  ;  rad_size: radial size of the blob
  ;  vert_size: vertical size of the blog
  ;  rad_pos: radial position of the blob
  ;  vert_pos: vertical position of the blob
  ;  fluct_amp: amplitude of the blob
  ;  
  ;OUTPUTs:
  ;  Matrix of the blob with the same coordinates
  ;  as the BES measurement.
                                  
                                  
  default, rad_pos, 2200.
  default, vert_pos, 15.
  default, rad_size, 10.
  default, vert_size, 10.
  default, hole_rad_pos, 2200.
  default, hole_vert_pos, 15.
  default, hole_size, 10.
  default, fluct_amp, 0.05
  default, hole_fluct_amp, -0.02
  default, spatial_pos, getcal_kstar_spat(14110) 
  default, rot_freq, 100.
  default, vx, 100 ;m/s
  default, vy, 50 ;m/s
  default, hole_vx, -50.
  default, hole_vy, -100.
  default, sampling_time, 0.5e-6
  default, time, 1e-3
  default, v_size, 50.
  default, plot, 0
  
  
  
  n_rad=n_elements(spatial_pos[0,*,0])
  n_vert=n_elements(spatial_pos[*,0,0])
  
  if keyword_set(moving) then begin
    n_time=time/sampling_time
    return_matrix=dblarr(n_vert,n_rad,n_time)
    blob_density=dblarr(n_time) ; Gaussian integral divided by the ellipse defined by the sigma_r, sigma_z of the Gaussian (surface density)
    time_vec=dindgen(n_time)/n_time*time
    for i=0,n_time-1 do begin
      t=i*sampling_time
      arg=2*!pi*rot_freq*t
      arg2=0.
      for j=0,n_vert-1 do begin
        for k=0,n_rad-1 do begin
          ;blob
          x=spatial_pos[j,k,0]-(vx*1e3*t+rad_pos)
          y=spatial_pos[j,k,1]-(vy*1e3*t+vert_pos)
          a1=(cos(arg)/(rad_size+v_size*1e3*t))^2+(sin(arg)/(vert_size))^2
          b1=-0.5*sin(2*arg)/(rad_size+v_size*1e3*t)^2+0.5*sin(2*arg)/(vert_size)^2
          c1=(sin(arg)/(rad_size+v_size*1e3*t))^2+(cos(arg)/(vert_size))^2
          return_matrix[j,k,i]=fluct_amp*exp(-0.5*(a1*x^2+2*b1*x*y+c1*y^2))
          ;hole
          x=spatial_pos[j,k,0]-(hole_vx*1e3*t+hole_rad_pos)
          y=spatial_pos[j,k,1]-(hole_vy*1e3*t+hole_vert_pos)
          a2=(cos(arg2)/(hole_size))^2+(sin(arg2)/(hole_size))^2
          b2=-0.5*sin(2*arg2)/(hole_size)^2+0.5*sin(2*arg2)/(vert_size)^2
          c2=(sin(arg2)/(hole_size+v_size*1e3*t))^2+(cos(arg2)/(hole_size))^2
          return_matrix[j,k,i]+=hole_fluct_amp*exp(-0.5*(a2*x^2+2*b2*x*y+c2*y^2))
        endfor
     endfor
      min_axis=((a1+c1)-sqrt((a1+c1)^2 + 4*(b1^2/4 - a1*c1)))/2
      maj_axis=((a1+c1)+sqrt((a1+c1)^2 + 4*(b1^2/4 - a1*c1)))/2
      blob_density[i]=fluct_amp*(4*!pi/sqrt(4*a1*c1-b1^2))/(!pi*maj_axis*min_axis)
      if not i then zrange=[hole_fluct_amp, fluct_amp]
      if keyword_set(plot) then begin
        contour, return_matrix[*,*,i], reform(spatial_pos[*,*,0]),reform(spatial_pos[*,*,1]),nlevel=21, /fill, /irreg, /iso, zrange=zrange
        plots, reform(spatial_pos[*,*,0]),reform(spatial_pos[*,*,1]), psym=4
      endif
    endfor
  endif else begin
    return_matrix=fluct_amp*exp(-0.5*(((spatial_pos[*,*,0]-rad_pos)/rad_size)^2+((spatial_pos[*,*,1]-vert_pos)/vert_size)^2))
    blob_density=fluct_amp*2/sqrt(rad_size*vert_size)
    if keyword_set(plot) then begin
      contour, return_matrix, reform(spatial_pos[*,*,0]),reform(spatial_pos[*,*,1]),nlevel=21, /fill, /irreg, /iso
      plots, reform(spatial_pos[*,*,0]),reform(spatial_pos[*,*,1]), psym=4
    endif
  endelse
  
  return, return_matrix
end
