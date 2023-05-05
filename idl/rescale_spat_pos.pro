function rescale_spat_pos, shot=shot, nwin_x=nwin_x, nwin_y=nwin_y, rad_res=rad_res, vert_res=vert_res

default, shot, 14110
default, nwin_x, 8
default, nwin_y, 32
default, rad_res, 10.
default, vert_res, 10.

spatial_pos_init=getcal_kstar_spat(shot)
spat_pos_mid=[mean(spatial_pos_init[*,*,0]),mean(spatial_pos_init[*,*,1])]
spatial_pos=dblarr(nwin_x,nwin_y,2)

for i=0,nwin_x-1 do begin
    spatial_pos[i,*,1]=spat_pos_mid[1]-nwin_x/2*vert_res+i*vert_res
endfor
for j=0,nwin_y-1 do begin
    spatial_pos[*,j,0]=spat_pos_mid[0]-nwin_y/2*rad_res+j*rad_res
endfor

return, spatial_pos
end