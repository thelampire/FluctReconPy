pro test_fluct_resp_matrix


;The matrix itself is a Lower triangular matrix

tauvec=dblarr(16)
matrix=dblarr(16,16)
energy=90e3 ;eV
e=1.6e-19
mh=1.6727e-27
v=sqrt(2*energy*e/mh)
tauvec=[10.,10.,10.,10.,9.,8.,7.,6.,5.,4.,3.,2.,2.,2.,2.,2.]*1e-9 ;ns

for j=0,15 do begin
  for i=j,15 do begin
    matrix[i,j]=10e-9/tauvec[i]*exp(-((i-j)*0.01)/(v*tauvec[i]))
  endfor
endfor
matrix=transpose(matrix)

;density=(tanh(findgen(16)-4)+1)*0.95e20+1e19
;light = density ## matrix



stop
end