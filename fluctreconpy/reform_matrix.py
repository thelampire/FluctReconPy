import numpy as np

def reform_matrix(matrix, dimensions, reverse=False):

#**************************************************************
#*** reform_matrix.pro --< python  by M. Lampert 08.31.2016 ***
#**************************************************************



    if (len(matrix.shape)) != 4 and (len(matrix.shape)) != 2 :
        print(f"The shape of the matrix is {matrix.shape}")
        raise ValueError('The software cannot be used for matrices other than 4 dimensions converting to 2 dimensions & backwards.')


    #Caveat: It is only working for matrices which need to be converted from 4D to 2D & backwards
    return_matrix=np.zeros(dimensions)
    if not reverse:

        n1=len(matrix[:,0,0,0])
        n2=len(matrix[0,:,0,0])
        n3=len(matrix[0,0,:,0])
        n4=len(matrix[0,0,0,:])

        if len(dimensions) != 2 :
            errormess= 'When forward reformation is used: the dimensions should be a two element vector.'
            print, errormess
            return -1

        if ((matrix[:,:,0,0].shape[0]*matrix[:,:,0,0].shape[1]) != dimensions[0] or
            (matrix[0,0,:,:].shape[0]*matrix[0,0,:,:].shape[1]) != dimensions[1]):
            raise ValueError('The dimensions of the matrix must match the numbers in the dimensions vector.')

        for i in range(0, n1):
            for j in range(0, n2):
                for k in range(0, n3):
                    for l in range(0, n4):
                        return_matrix[i*n2+j, k*n4+l]=matrix[i,j,k,l]

    else:
        n1=dimensions[0]
        n2=dimensions[1]
        n3=dimensions[2]
        n4=dimensions[3]

        if len(dimensions) != 4:
            raise ValueError('When reverse reformation is used: the dimensions should be a two element vector.')


        if (len(matrix[:,0]) != dimensions[0]*dimensions[1]) or (len(matrix[0,:]) != dimensions[2]*dimensions[3]):
            raise ValueError('The dimensions of the matrix must match the numbers in the dimensions vector.')


        for i in range(0, n1):
            for j in range(0, n2):
                for k in range(0, n3):
                    for l in range(0, n4):
                        return_matrix[i,j,k,l]=matrix[i*n2+j, k*n4+l]

    return return_matrix