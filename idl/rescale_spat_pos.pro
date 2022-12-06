function rescale_spat_pos, shot=shot, nwin=nwin, rad_res=rad_res, vert_res=vert_res

default, shot, 14110
default, nwin,16.
default, rad_res, 10.
default, vert_res, 10.

spatial_pos_init=getcal_kstar_spat(shot)
spat_pos_mid=[mean(spatial_pos_init[*,*,0]),mean(spatial_pos_init[*,*,1])]
spatial_pos=dblarr(nwin,nwin,2)


for i=0,nwin-1 do begin
    spatial_pos[*,i,0]=spat_pos_mid[0]-nwin/2*rad_res+i*rad_res
    spatial_pos[i,*,1]=spat_pos_mid[1]-nwin/2*vert_res+i*vert_res
  endfor
return, spatial_pos
end