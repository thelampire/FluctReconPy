pro analyze_decomp_res, filename_anal=filename_anal
hardon, /color
loadct, 5
ymax=2.
imin=1.
number=20.
fnum=20.
default, filename_anal,'test_decomp_all_save_0-50noise'
restore, filename_anal+".sav"
!p.font=0
for j=0,n_elements(spatial_resolution)-1 do begin
  for i=imin,n_elements(blob_size)-1 do begin
    y=(abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,i,j,*,0]
    if i eq imin then begin
      plot, noise_vector, total(y,4)/n, $
        xtitle='Relative noise amplitude', ytitle='Radial FWHM error', title='Spat res:'+strtrim(spatial_resolution[j],2)+'mm',$
        charsize=2., thick=2, psym=-4, color=1, yrange=[0,ymax]
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=1
      xyouts, 0.8, 0.93, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=1
    endif else begin
      oplot, noise_vector, total(y,4)/n, color=fnum+i*number, psym=-4
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=fnum+i*number
      xyouts, 0.8, 0.93-(i-imin)*0.03, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=fnum+i*number
    endelse
  endfor
endfor
for j=0,n_elements(spatial_resolution)-1 do begin
  for i=imin,n_elements(blob_size)-1 do begin
    y=(abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,i,j,*,1]
    if i eq imin then begin
      plot, noise_vector, total(y,4)/n, $
        xtitle='Relative noise amplitude', ytitle='Vertical FWHM error', title='Spat res:'+strtrim(spatial_resolution[j],2)+'mm',$
        yrange=[0,ymax],charsize=2., thick=2, psym=-4, color=1
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=1
      xyouts, 0.8, 0.93, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=1
    endif else begin
      oplot, noise_vector, total(y,4)/n, color=fnum+i*number, psym=-4
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=fnum+i*number
      xyouts, 0.8, 0.93-(i-imin)*0.03, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=fnum+i*number
    endelse
  endfor
endfor

for j=0,n_elements(spatial_resolution)-1 do begin
  for i=imin,n_elements(blob_size)-1 do begin
    y=(abs(pos_calc-pos_orig)/(spatial_resolution[j]))[*,i,j,*,0]
    if i eq imin then begin
      plot, noise_vector, total(y,4)/n, $
        xtitle='Relative noise amplitude', ytitle='Radial position error', title='Spat res:'+strtrim(spatial_resolution[j],2)+'mm',$
        charsize=2., thick=2, psym=-4, color=1, yrange=[0,ymax]
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=1
      xyouts, 0.8, 0.93, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=1
    endif else begin
      oplot, noise_vector, total(y,4)/n, color=fnum+i*number, psym=-4
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=fnum+i*number
      xyouts, 0.8, 0.93-(i-imin)*0.03, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=fnum+i*number
    endelse
  endfor
endfor

for j=0,n_elements(spatial_resolution)-1 do begin
  for i=imin,n_elements(blob_size)-1 do begin
    y=(abs(pos_calc-pos_orig)/(spatial_resolution[j]))[*,i,j,*,1]
    if i eq imin then begin
      plot, noise_vector, total(y,4)/n, $
        xtitle='Relative noise amplitude', ytitle='Vertical position error', title='Spat res:'+strtrim(spatial_resolution[j],2)+'mm',$
        charsize=2., thick=2, psym=-4, color=1, yrange=[0,ymax]
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=1
      xyouts, 0.8, 0.93, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=1
    endif else begin
      oplot, noise_vector, total(y,4)/n, color=fnum+i*number, psym=-4
      errplot, noise_vector, total(y,4)/n-variance(y,dimension=4)/(n-1),total(y,4)/n+variance(y,dimension=4)/(n-1), color=fnum+i*number
      xyouts, 0.8, 0.93-(i-imin)*0.03, 'Blob size: '+strtrim(blob_size[i]*spatial_resolution[j],2)+'mm', /normal, color=fnum+i*number
    endelse
  endfor
endfor
hardfile, filename_anal+'.ps'
end
