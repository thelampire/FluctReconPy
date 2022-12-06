import numpy
pro get_data_for_fluct_matrix, shot, timerange=timerange, thomson=thomson

default, shot, 14110
default, timerange, [3.5,4.5]
if (thomson is None):#Thomson data is useless most of the time
    #Te profile read from Thomson
    n_edge=15
    index_edge=indgen(15)+13
    n_core=12
    index_core=indgen(12)
    for i in range(0, n_edge):
        get_rawsignal, shot, '\TS_EDGE'+strtrim(index_edge[i],2)+':EDGE'+strtrim(index_edge[i],2)+'_TE', t, d, /store
        if i == 0: edge_te_profile=numpy.zeroes(len(d),n_edge)
        edge_te_profile[*,i]=d
        
  
        for i in range(0, n_core):
            get_rawsignal, shot,  '\TS_CORE'+strtrim(index_core[i],2)+':CORE'+strtrim(index_core[i],2)+'_TE', t, d, /store
            if i == 0: core_te_profile=numpy.zeroes(len(d),n_core)
            core_te_profile[*,i]=d
        
  
        #ne profile read from Thomson
  
        for i in range(0, n_edge):
            get_rawsignal, shot, '\TS_EDGE'+strtrim(index_edge[i],2)+':EDGE'+strtrim(index_edge[i],2)+'_NE', t, d, /store
            if i == 0: edge_ne_profile=numpy.zeroes(len(d),n_edge)
            edge_ne_profile[*,i]=d
        
  
        for i in range(0, n_core):
            get_rawsignal, shot,  '\TS_CORE'+strtrim(index_core[i],2)+':CORE'+strtrim(index_core[i],2)+'_NE', t, d, /store
            if i == 0: core_ne_profile=numpy.zeroes(len(d),n_core)
            core_ne_profile[*,i]=d
        
  
        thomson_calibration_edge=[2169.3, 2179.6, 2189.9, 2201.0, 2208.4, 2218.9, 2227.1, 2234.6, 2245.1, 2259.5, 2270.8, 2280.7, 2292.2, 2302.9, 2313.7]
        thomson_calibration_core=[1801.0,1835.9,1865.0,1898.0,1956.7,1986.3,2016.6,2048.1,2076.7,2108.6,2138.7,2173.0]
end

#EFIT read
efit=get_kstar_efit_mds(shot,/store,errormess=errormess)

#Read ECE data for temperature

ECE_channels=[string(indgen(10)+2, format='(i2.2)'),$
    string(indgen(15)+13, format='(i2.2)'),$
    string(indgen(20)+29, format='(i2.2)'),$
    string(indgen(6)+70, format='(i2.2)')]

n_ece=len(ece_channels)
ece_spatial_coordinates=numpy.zeroes(n_ece)
ece_channel_names=numpy.zeroes(n_ece)

for i in range(0, len(n_ece)):
    ece_channel_names[i]='\ECE'+ECE_channels[i]
    get_rawsignal, shot, ece_channel_names[i], timerange=timerange, t, d, /store

#Read ECE spatial calibration
if file_test(ece_spat_file) :
    restore, ece_spat_file
 else begin
    for i in range(0, len(n_ece)):
        read_kstar_mdsplus, 14110, ece_channel_names[i]+':RPOS2ND', time, data
        ece_spatial_coordinate[i]=data
        
        save, ece_channel_names, ece_spatial_coordinates, filename=strtrim(shot,2)+'_ECE_channels_spatial_coord.sav'
endelse

#Detector corner coordinates
det_corn_pos=getcal_kstar_spat(shot, /detcorn)

# Getting electron density profile

end