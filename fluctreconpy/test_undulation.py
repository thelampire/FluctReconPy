import numpy
pro test_undulation, nvec=nvec
    default, nvec, 100
    default, shot, 14110
    nr=16.
    nz=4.
    rres=10.
    zres=10.
    r_vec=numpy.arange(nr)*rres #radial coordinate in mm
    z_vec=numpy.arange(nz)*zres #horizontal coordinate in mm
  
    r_vec_interp=numpy.arange(nr*nvec)*rres/nvec
    z_vec_interp=numpy.arange(nz*nvec)*zres/nvec
    undulation=numpy.zeroes(16,4,16,4)
    diffmat=numpy.zeroes(nr*nvec,nz*nvec,4)
    if (test is None) :
        for i1 in range(1, nr-2+1):
            for j1 in range(1, nz-2+1):
                for i2 in range(1, nr-2+1):
                    for j2 in range(1, nz-2+1):
                        diffmat[i1*nvec:(i1+1)*nvec-1,j1*nvec:(j1+1)*nvec-1,0] = base_function_diff(r_vec[i1-1:i1+1],z_vec[j1-1:j1+1], nvec=nvec, /rdiff)
                        diffmat[i1*nvec:(i1+1)*nvec-1,j1*nvec:(j1+1)*nvec-1,1] = base_function_diff(r_vec[i1-1:i1+1],z_vec[j1-1:j1+1], nvec=nvec, /zdiff)
                        diffmat[i2*nvec:(i2+1)*nvec-1,j2*nvec:(j2+1)*nvec-1,2] = base_function_diff(r_vec[i2-1:i2+1],z_vec[j2-1:j2+1], nvec=nvec, /rdiff)
                        diffmat[i2*nvec:(i2+1)*nvec-1,j2*nvec:(j2+1)*nvec-1,3] = base_function_diff(r_vec[i2-1:i2+1],z_vec[j2-1:j2+1], nvec=nvec, /zdiff)
                        undulation[i1,j1,i2,j2]=total(diffmat[*,*,0]*diffmat[*,*,2]+diffmat[*,*,1]*diffmat[*,*,3])*rres/nvec*zres/nvec
                        if undulation[i1,j1,i2,j2] > 0: print, i1,j1,i2,j2, undulation[i1,j1,i2,j2]
                        diffmat[*,*,*]=0
                            
                                
                        
                
         else begin
            for i in range(1, nr-2+1):
                for j in range(1, nz-2+1):
        
                    diffmat[i*nvec:(i+1)*nvec-1,j*nvec:(j+1)*nvec-1,0] = base_function_diff(r_vec[i-1:i+1],z_vec[j-1:j+1], nvec=nvec, /rdiff)
                    diffmat[i*nvec:(i+1)*nvec-1,j*nvec:(j+1)*nvec-1,1] = base_function_diff(r_vec[i-1:i+1],z_vec[j-1:j+1], nvec=nvec, /zdiff)
                    undulation[i1,j1,i2,j2]=total((diffmat[*,*,0]**2+diffmat[*,*,1]**2)*((r_vec[0:nr-3]-r_vec[2:nr-1])/2 # z_vec[0:2]-))
        
                    if undulation[i1,j1,i2,j2] > 0: print, i1,j1,i2,j2, undulation[i1,j1,i2,j2]
                    diffmat[*,*,*]=0
                        
                
        endelse
stop
end

function base_function_diff, r, z, nvec=nvec, rdiff=rdiff, zdiff=zdiff
    default, nvec, 100
    if len(r) != 3 | len(z) != 3 :
        print, 'The boundary of the pyramid needs to be defined by 3-3 r z points.'
        return, 0
        

        r_vec=numpy.arange(nvec)/nvec*(r[2]-r[0])+r[0]
        z_vec=numpy.arange(nvec)/nvec*(z[2]-z[0])+z[0]

        base_diff=numpy.zeroes(nvec,nvec)
        for i_r in range(0, nvec):
            for j_z in range(0, nvec):
      
            if (rdiff is None) :
      
                if i_r < nvec/2 :
                    if j_z < nvec/2 :
                        base_diff[i_r,j_z]=1./(r[1]-r[0])*(z_vec[j_z]-z[0])/(z[1]-z[0])
                             else begin
                                base_diff[i_r,j_z]=1./(r[1]-r[0])*(z_vec[j_z]-z[2])/(z[1]-z[2])
                                endelse
                         else begin
                            if j_z < nvec/2 :
                                base_diff[i_r,j_z]=1./(r[1]-r[2])*(z_vec[j_z]-z[0])/(z[1]-z[0])
                                 else begin
                                    base_diff[i_r,j_z]=1./(r[1]-r[2])*(z_vec[j_z]-z[2])/(z[1]-z[2])
                                endelse
                        endelse
                
    
                if (zdiff is None) :
                    if i_r < nvec/2 :
                        if j_z < nvec/2 :
                            base_diff[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])/(z[1]-z[0])
                                 else begin
                                    base_diff[i_r,j_z]=(r_vec[i_r]-r[0])/(r[1]-r[0])/(z[1]-z[2])
                                endelse
                         else begin
                            if j_z < nvec/2 :
                                base_diff[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])/(z[1]-z[0])
                                 else begin
                                    base_diff[i_r,j_z]=(r_vec[i_r]-r[2])/(r[1]-r[2])/(z[1]-z[2])
                                endelse
                        endelse
                
    
                
        
        return, base_diff
end