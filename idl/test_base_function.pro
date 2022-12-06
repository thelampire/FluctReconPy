pro test_base_function, shot, timerange=timerange

s=getcal_kstar_spat(shot)
nz=4
nr=16
z_inter=4
r_inter=4

r_vec=reform(s[1,*,0]+s[2,*,0])/2
z_vec=reform(s[*,7,1]+s[*,8,1])/2
si=dblarr((nz-1)*z_inter+1,(nr-1)*r_inter+1,2)
si[*,*,0]=interpolate(s[*,*,0], (dindgen((nz-1)*z_inter+1))/z_inter, (dindgen((nr-1)*r_inter+1))/r_inter,/grid)
si[*,*,1]=interpolate(s[*,*,1], (dindgen((nz-1)*z_inter+1))/z_inter, (dindgen((nr-1)*r_inter+1))/r_inter,/grid)

r_vec_inter=interpolate(r_vec, (dindgen((nr-1)*r_inter+1))/r_inter)
z_vec_inter=interpolate(z_vec, (dindgen((nz-1)*z_inter+1))/z_inter)

signal=dblarr(nz,nr)
signal_inter=dblarr((nz-1)*z_inter+1,(nr-1)*r_inter+1)
for i=0,nz-1 do begin
  for j=0,nr-1 do begin
    get_rawsignal, shot, 'BES-'+strtrim(i+1,2)+'-'+strtrim(j+1,2), t, d, timerange=timerange, /nocalib
    signal[i,j]=mean(d)
  endfor
endfor

for i=0,(nz-1)*z_inter-1 do begin
  for j=0,(nr-1)*r_inter-1 do begin
    ir=j/r_inter
    iz=i/z_inter
    p1=(r_vec_inter[j]-r_vec[ir+1])/(r_vec[ir]-r_vec[ir+1]) * (z_vec_inter[i]-z_vec[iz+1])/(z_vec[iz]-z_vec[iz+1])  
    p2=(r_vec_inter[j]-r_vec[ir])  /(r_vec[ir+1]-r_vec[ir]) * (z_vec_inter[i]-z_vec[iz])/  (z_vec[iz+1]-z_vec[iz])  
    p3=(r_vec_inter[j]-r_vec[ir+1])/(r_vec[ir]-r_vec[ir+1]) * (z_vec_inter[i]-z_vec[iz])/  (z_vec[iz+1]-z_vec[iz])
    p4=(r_vec_inter[j]-r_vec[ir])  /(r_vec[ir+1]-r_vec[ir]) * (z_vec_inter[i]-z_vec[iz+1])/(z_vec[iz]-z_vec[iz+1])
    signal_inter[i,j]=p1 * signal[iz,  ir]+$
                      p2 * signal[iz+1,ir+1]+$
                      p3 * signal[iz+1,ir]+$
                      p4 * signal[iz,  ir+1]
  endfor
endfor

for i=0,(nz-1)*z_inter-1 do begin
  j=(nr-1)*r_inter
  ir=j/r_inter
  iz=i/z_inter
  p1=(z_vec_inter[i]-z_vec[iz+1])/(z_vec[iz]-z_vec[iz+1])
  p3=(z_vec_inter[i]-z_vec[iz])/  (z_vec[iz+1]-z_vec[iz])
  signal_inter[i,j]=p1 * signal[iz,  ir] + p3 * signal[iz+1,ir]
endfor

for j=0,(nr-1)*r_inter-1 do begin
  i=(nz-1)*z_inter
  ir=j/r_inter
  iz=i/z_inter
  p1=(r_vec_inter[j]-r_vec[ir+1])  /(r_vec[ir]-r_vec[ir+1])
  p4=(r_vec_inter[j]-r_vec[ir])  /(r_vec[ir+1]-r_vec[ir])
  signal_inter[i,j]=p1 * signal[iz,  ir] + p4 * signal[iz,ir+1]
endfor

signal_inter[(nz-1)*z_inter,(nr-1)*r_inter]=signal[(nz-1),(nr-1)]


erase
loadct, 5
nlevel=21
min=min([min(signal),min(signal_inter)])
max=max([max(signal),max(signal_inter)])
levels=findgen(nlevel)/nlevel*(max-min)+min
contour, transpose(signal), transpose(reform(s[*,*,0])),transpose(reform(s[*,*,1])),        /iso, /fill, nlevel=nlevel, position=[0.05,0.40,0.95,0.70], /noerase, /xstyle, /ystyle, /irreg, levels=levels
plots, transpose(reform(s[*,*,0])),transpose(reform(s[*,*,1])), psym=4, color=120
contour, transpose(signal_inter), transpose(reform(si[*,*,0])),transpose(reform(si[*,*,1])),/iso, /fill, nlevel=nlevel, position=[0.05,0.05,0.95,0.35], /noerase, /xstyle, /ystyle, /irreg, levels=levels
plots, transpose(reform(s[*,*,0])),transpose(reform(s[*,*,1])), psym=4, color=120
stop
end