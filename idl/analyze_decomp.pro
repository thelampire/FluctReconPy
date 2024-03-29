pro analyze_decomp, filename=filename

default, filename, '/Users/mlampert/IDLWorkspace83/test_decomp_all_save_0-50noise_nores_dens.sav'

restore, filename
loadct, 5
hardon, /color
fwhm_calc=reform(fwhm_calc)
fwhm_orig=reform(fwhm_orig)
n=n_elements(fwhm_calc[0,0,*,0])

for i=0,n_elements(chi2[0,*,0,0])-1 do begin
    y=total(chi2[*,i,*,*],4)
    if i eq 0 then begin
        plot, noise_vector, y, $
            xtitle='Relative noise amplitude', ytitle='Average density error [1e19m-3]', $
            yrange=range(total(chi2[*,*,*,*],4)),$
            charsize=2., thick=2
        xyouts, 0.8, 0.93, 'Blob size: 10mm', /normal
    endif else begin
        oplot, noise_vector, y, color=80+i*20
        xyouts, 0.8, 0.93-i*0.03, 'Blob size: '+strtrim((i+2)*5,2)+'mm', /normal, color=80+i*20
    endelse
endfor
for i=0,n_elements(fwhm_orig[0,*,0,0])-1 do begin
  if i eq 0 then begin
    plot, noise_vector, total((abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,i,*,0],3)/n, $
      xtitle='Relative noise amplitude', ytitle='Radial FWHM error', $
      yrange=range(reform(total((abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,*,*,0],3)/n)),$
      charsize=2., thick=2
    xyouts, 0.8, 0.93, 'Blob size: 10mm', /normal
  endif else begin
    oplot, noise_vector, total((abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,i,*,0],3)/n, color=80+i*20
    xyouts, 0.8, 0.93-i*0.03, 'Blob size: '+strtrim((i+2)*5,2)+'mm', /normal, color=80+i*20
  endelse
endfor

for i=0,n_elements(fwhm_orig[0,*,0,1])-1 do begin
  if i eq 0 then begin
    plot, noise_vector, total((abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,i,*,1],3)/n, $
      xtitle='Relative noise amplitude', ytitle='Vertical FWHM error', $
      yrange=[0,5],$
      charsize=2., thick=2
      ;yrange=range(reform(total((abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,*,*,1],3)/n))
    xyouts, 0.8, 0.93, 'Blob size: 10mm', /normal
  endif else begin
    oplot, noise_vector, total((abs(fwhm_calc-fwhm_orig)/(fwhm_orig))[*,i,*,1],3)/n, color=80+i*20
    xyouts, 0.8, 0.93-i*0.03, 'Blob size: '+strtrim((i+2)*5,2)+'mm', /normal, color=80+i*20
  endelse
endfor
pos_calc=reform(pos_calc)
pos_orig=reform(pos_orig)
for i=0,n_elements(fwhm_orig[0,*,0,1])-1 do begin
  if i eq 0 then begin
    plot, noise_vector, total((abs(pos_calc-pos_orig)/(pos_orig))[*,i,*,0],3)/n, $
      xtitle='Relative noise amplitude', ytitle='Radial position error', $
      yrange=range(reform(total((abs(pos_calc-pos_orig)/(pos_orig))[*,1:4,*,0],3)/n)),$
      charsize=2., thick=2
    xyouts, 0.8, 0.93, 'Blob size: 10mm', /normal
  endif else begin
    oplot, noise_vector, total((abs(pos_calc-pos_orig)/(pos_orig))[*,i,*,0],3)/n, color=80+i*20
    xyouts, 0.8, 0.93-i*0.03, 'Blob size: '+strtrim((i+2)*5,2)+'mm', /normal, color=80+i*20
  endelse
endfor

for i=0,n_elements(fwhm_orig[0,*,0,1])-1 do begin
  if i eq 0 then begin
    plot, noise_vector, total((abs(pos_calc-pos_orig)/(pos_orig))[*,i,*,1],3)/n, $
      xtitle='Relative noise amplitude', ytitle='Vertical position error', $
      yrange=range(reform(total((abs(pos_calc-pos_orig)/(pos_orig))[*,*,*,1],3)/n)),$
      charsize=2., thick=2
    xyouts, 0.8, 0.93, 'Blob size: 10mm', /normal
  endif else begin
    oplot, noise_vector, total((abs(pos_calc-pos_orig)/(pos_orig))[*,i,*,1],3)/n, color=80+i*20
    xyouts, 0.8, 0.93-i*0.03, 'Blob size: '+strtrim((i+2)*5,2)+'mm', /normal, color=80+i*20
  endelse
endfor

hardfile, filename+'.ps'
end