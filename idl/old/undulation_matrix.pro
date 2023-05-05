function undulation_matrix, spatial_pos=spatial_pos, r_vec=r_vec, z_vec=z_vec, floating=floating, onedim=onedim, errormess=errormess, silent=silent

;*********************************************************
;*** undulation_matrix.pro    by M. Lampert 2016-11-02 ***
;*********************************************************
;
; The procedure calculates the undulation matrix for a
; given geometry. The return matrix is a 4D matrix with
; dimensions of [n_r,n_z,n_r,n_z]. The current calculation
; is using bilinear base functions, which are placed on an
; equidistant grid.
; 
;  INPUTs:
;      edges instead of zero edges. (default: 0)
;    spatial_pos: [z_vec, r_vec, 2] matrix for the spatial
;      coordinates. The center line is taken for the 
;      calculation for both radial and vertical directions.
;    r_vec: 1D radial coordinate vector
;    v_vec: 1D vertical coordinate vector
;    /floating: the calculation is done with floating
;    /onedim: one dimensional coordinate vectors
;    errormess: error message
;    /silent: no printing (default: 0) 
;    
;  OUTPUT:
;    Returns the undulation matrix.
;    
;  Caveat:
;    The code only works for equidistant grid. Irregular
;    grids are averaged and the middle axis is taken as
;    radial and vertical coordinates.
; 

default, reform, 1
default, floating, 0
default, onedim, 1
default, silent, 0

if not defined(r_vec) and not defined(z_vec) then begin
  if not defined(spatial_pos) then begin
    errormess='No spatial coordinates are defined. Returning...'
    if not keyword_set(silent) then print, errormess
    return, 0
  endif
  
  nr=n_elements(spatial_pos[0,*,0])
  nz=n_elements(spatial_pos[*,0,0])
  if keyword_set(onedim) then begin
    r=reform(spatial_pos[nz/2-1,*,0]+spatial_pos[nz/2,*,0])/2.
    z=reform(spatial_pos[*,nr/2-1,1]+spatial_pos[*,nr/2,1])/2.
  endif else begin
    r=spatial_pos[*,*,0]
    z=spatial_pos[*,*,1]
  endelse
endif

if defined(r_vec) and defined(z_vec) then begin
    nr=n_elements(r_vec)
    nz=n_elements(z_vec)
    r=r_vec
    z=z_vec
endif 
 
  
undulation=fltarr(nz,nr,nz,nr)

for i=1,nr-2 do begin
  for j=1,nz-2 do begin
    
                            
    undulation[j,i,j,i-1]  =(z[j+1]-z[j-1])/3.*(1./(r[i-1]-r[i])) + $
                            (r[i]-r[i-1])  /6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            
    undulation[j,i,j,i]    =(z[j+1]-z[j-1])/3.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + $
                            (r[i+1]-r[i-1])/3.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            
    undulation[j,i,j,i+1]  =(z[j+1]-z[j-1])/3.*(1./(r[i]-r[i+1])) + $
                            (r[i+1]-r[i])  /6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
                            
    undulation[j,i,j-1,i-1]=(z[j]-z[j-1])  /6.*(1./(r[i-1]-r[i])) + $
                            (r[i]-r[i-1])  /6.*(1./(z[j-1]-z[j]))
                            
    undulation[j,i,j-1,i]  =(z[j]-z[j-1])  /6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + $
                            (r[i+1]-r[i-1])/3.*(1./(z[j-1]-z[j]))
                                                
    undulation[j,i,j-1,i+1]=(z[j]-z[j-1])  /6.*(1./(r[i]-r[i+1])) + $
                            (r[i+1]-r[i])  /6.*(1./(z[j-1]-z[j]))                            
                                  
    undulation[j,i,j+1,i-1]=(z[j+1]-z[j])  /6.*(1./(r[i-1]-r[i])) + $
                            (r[i]-r[i-1])  /6.*(1./(z[j]-z[j+1]))
    
    undulation[j,i,j+1,i]  =(z[j+1]-z[j])  /6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + $
                            (r[i+1]-r[i-1])/3.*(1./(z[j]-z[j+1]))
                            
    undulation[j,i,j+1,i+1]=(z[j+1]-z[j])  /6.*(1./(r[i]-r[i+1])) + $
                            (r[i+1]-r[i])  /6.*(1./(z[j]-z[j+1]))
                            
                            
    
  endfor
endfor

if keyword_set(floating) then begin
  for i=0,nr-1 do begin
    for j=0,nz-1 do begin

      if j eq 0 then begin
        if i eq 0 then begin
          undulation[j,i,j,i]=    (z[j+1]-z[j])/6.*(1./(r[i+1]-r[i])) + $
                                  (r[i+1]-r[i])/6.*(1./(z[j+1]-z[j]))
          undulation[j,i,j,i+1]=  (z[j+1]-z[j])/6.*(1./(r[i+1]-r[i])) + $
                                  (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))
          undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + $
                                  (r[i+1]-r[i])/6.*(1./(z[j+1]-z[j]))
          undulation[j,i,j+1,i+1]=(z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + $
                                  (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))
        endif else begin
          if i eq nr-1 then begin
            undulation[j,i,j,i-1]=  (z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j+1]-z[j]))
            undulation[j,i,j,i]=    (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j+1]-z[j])) ;corr
            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
          endif else begin
            undulation[j,i,j,i-1]=  (z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j+1]-z[j]))
            undulation[j,i,j,i]=    (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + $
                                    (r[i+1]-r[i-1])/3.*(1./(z[j+1]-z[j]))
            undulation[j,i,j,i+1]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j+1]-z[j]))
            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
            undulation[j,i,j+1,i+1]=(z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))
          endelse
        endelse
      endif else begin
        if j eq nz-1 then begin
          if i eq 0 then begin
            undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i+1]-r[i])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))
            undulation[j,i,j-1,i+1]=(z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))
            undulation[j,i,j,i]=    (z[j]-z[j-1])/6.*(1./(r[i+1]-r[i])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1]))
            undulation[j,i,j,i+1]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1]))
          endif else begin
            if i eq nr-1 then begin
              undulation[j,i,j-1,i-1]=(z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + $
                                      (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
              undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])) + $
                                      (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
              undulation[j,i,j,i-1]=  (z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + $
                                      (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1]))
              undulation[j,i,j,i]=    (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])) + $
                                      (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1])) ;corr
            endif else begin
              undulation[j,i,j-1,i-1]=(z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + $
                                      (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
              undulation[j,i,j-1,i]=  (z[j]-z[j-1])  /6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + $
                                      (r[i+1]-r[i-1])/3.*(1./(z[j-1]-z[j]))
                
              undulation[j,i,j-1,i+1]=(z[j]-z[j-1])  /6.*(1./(r[i]-r[i+1])) + $
                                      (r[i+1]-r[i])  /6.*(1./(z[j-1]-z[j]))
              undulation[j,i,j,i-1]=  (z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + $
                                      (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1]))
              undulation[j,i,j,i]=    (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])+1./(r[i+1]-r[i])) + $
                                      (r[i+1]-r[i-1])/3.*(1./(z[j]-z[j-1])) ;corr
              undulation[j,i,j,i+1]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + $
                                      (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1]))
            endelse
          endelse
        endif else begin
          if i eq 0 then begin
            undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i+1]-r[i])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j])) 
            undulation[j,i,j-1,i+1]=(z[j]-z[j-1])/6.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j-1]-z[j]))  
            undulation[j,i,j,i]=    (z[j+1]-z[j-1])/3.*(1./(r[i+1]-r[i])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
            undulation[j,i,j,i+1]=  (z[j+1]-z[j-1])/3.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i+1]-r[i])) + $
                                    (r[i+1]-r[i])/3.*(1./(z[j]-z[j+1]))
            undulation[j,i,j+1,i+1]=(z[j+1]-z[j])/6.*(1./(r[i]-r[i+1])) + $
                                    (r[i+1]-r[i])/6.*(1./(z[j]-z[j+1]))
          endif
          if i eq nr-1 then begin
            undulation[j,i,j-1,i-1]=(z[j]-z[j-1])/6.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
            undulation[j,i,j-1,i]=  (z[j]-z[j-1])/6.*(1./(r[i]-r[i-1])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j-1]-z[j]))
            undulation[j,i,j,i-1]=  (z[j+1]-z[j-1])/3.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
            undulation[j,i,j,i]=    (z[j+1]-z[j-1])/3.*(1./(r[i]-r[i-1])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j-1])+1./(z[j+1]-z[j]))
            undulation[j,i,j+1,i-1]=(z[j+1]-z[j])/6.*(1./(r[i-1]-r[i])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
            undulation[j,i,j+1,i]=  (z[j+1]-z[j])/6.*(1./(r[i]-r[i-1])) + $
                                    (r[i]-r[i-1])/6.*(1./(z[j]-z[j+1]))
          endif
        endelse
      endelse
    endfor
  endfor
endif
return, undulation
end