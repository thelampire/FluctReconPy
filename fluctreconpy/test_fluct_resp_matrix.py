import numpy
pro test_fluct_resp_matrix


#The matrix itself is a Lower triangular matrix

tauvec=numpy.zeroes(16)
matrix=numpy.zeroes(16,16)
energy=90e3 #eV
e=1.6e-19
mh=1.6727e-27
v=sqrt(2*energy*e/mh)
tauvec=[10.,10.,10.,10.,9.,8.,7.,6.,5.,4.,3.,2.,2.,2.,2.,2.]*1e-9 #ns

for j in range(0, 15+1):
    for i in range(j, 15+1):
        matrix[i,j]=10e-9/tauvec[i]*numpy.exp(-((i-j)*0.01)/(v*tauvec[i]))
        

matrix=transpose(matrix)

#density=(tanh(numpy.arange(16)-4)+1)*0.95e20+1e19
#light = density ## matrix



stop
end