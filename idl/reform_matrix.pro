function reform_matrix, matrix, dimensions, reverse=reverse, errormess=errormess

;**************************************************************
;*** reform_matrix.pro             by M. Lampert 08.31.2016 ***
;**************************************************************



if (size(matrix))[0] ne 4 and (size(matrix))[0] ne 2 then begin
  print, 'The software cannot be used for matrices other than 4 dimensions converting to 2 dimensions and backwards.'
  return, -1
endif

;Caveat: It is only working for matrices which need to be converted from 4D to 2D and backwards
return_matrix=dblarr(dimensions)
if not keyword_set(reverse) then begin
  
  n1=n_elements(matrix[*,0,0,0])
  n2=n_elements(matrix[0,*,0,0])
  n3=n_elements(matrix[0,0,*,0])
  n4=n_elements(matrix[0,0,0,*])
  
  if n_elements(dimensions) ne 2 then begin
    errormess= 'When forward reformation is used then the dimensions should be a two element vector.'
    print, errormess
    return, -1
  endif
  if n_elements(matrix[*,*,0,0]) ne dimensions[0] or n_elements(matrix[0,0,*,*]) ne dimensions[1] then begin
    errormess='The dimensions of the matrix must match the numbers in the dimensions vector.'
    print, errormess
    return, -1
  endif
  for i=0,n1-1 do begin
    for j=0,n2-1 do begin
      for k=0,n3-1 do begin
        for l=0,n4-1 do begin
          return_matrix[i*n2+j, k*n4+l]=matrix[i,j,k,l]
        endfor
      endfor
    endfor
  endfor
  return, return_matrix
endif else begin
  n1=dimensions[0]
  n2=dimensions[1]
  n3=dimensions[2]
  n4=dimensions[3]
  

  if n_elements(dimensions) ne 4 then begin
    errormess= 'When reverse reformation is used then the dimensions should be a two element vector.'
    print, errormess
    return, -1
  endif
  
  if n_elements(matrix[*,0]) ne dimensions[0]*dimensions[1] or n_elements(matrix[0,*]) ne dimensions[2]*dimensions[3] then begin
    errormess='The dimensions of the matrix must match the numbers in the dimensions vector.'
    print, errormess
    return, -1
  endif
  
  for i=0,n1-1 do begin
    for j=0,n2-1 do begin
      for k=0,n3-1 do begin
        for l=0,n4-1 do begin
          return_matrix[i,j,k,l]=matrix[i*n2+j, k*n4+l]
        endfor
      endfor
    endfor
  endfor
  return, return_matrix
endelse



end