pro test_undulation_matrix
h_matrix_ref=reform_matrix(undulation_matrix(14110, floating=0),[64, 64])
while 1 do begin
  vector=(randomn(seed, 64)-0.5)*4.
  if (vector ## (h_matrix_ref ## vector)) lt 0 then begin
    print, vector ## (h_matrix_ref ## vector)
  endif
endwhile

end