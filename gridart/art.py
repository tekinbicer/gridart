# -*- coding: utf-8 -*-
import numpy as np
import ctypes
import os
import warnings
from skimage import io as skimage_io 

# Get the path and import the C-shared library.
libpath = os.path.abspath(os.path.join(os.path.dirname(__file__), './libgridart.so'))
libgridart = ctypes.CDLL(libpath)


def run(data, theta, num_grid=None, center=None):
    """
    Python wrapper for the art.c function.
    
    Parameters
    ----------
    data : raw data as 3-D ndarray
        1st dimension: Projections
        2nd dimension: Rows in projections (slices)
        3rd dimension: Columns in projections
        
    theta : projection angles 1-D ndarray
        The dimension of theta should match with
        the 1st dimension of data.
        
    Returns
    -------
    recon : Reconstructed 3-D volume
    """
    num_projs = np.array(data.shape[0], dtype='int32')
    num_rows = np.array(data.shape[1], dtype='int32')
    num_cols = np.array(data.shape[2], dtype='int32')
    
    if center is None:
        center = num_cols/2.
        
    # Make sure that inputs datatypes are correct
    data = np.array(data, dtype='float32')
    theta = np.array(theta, dtype='float32')
    center = np.array(center, dtype='float32')

    # Iterations count
    num_iter = 1

    # Init recon matrix.
    if num_grid is None:
        num_grid = num_cols
    recon = np.zeros((num_rows, num_grid, num_grid), dtype='float32')

    # Call C function to reconstruct recon matrix.
    c_float_p = ctypes.POINTER(ctypes.c_float)
    libgridart.myart.restype = ctypes.POINTER(ctypes.c_void_p)
    libgridart.myart(data.ctypes.data_as(c_float_p),
                 theta.ctypes.data_as(c_float_p),
                 ctypes.c_float(center),
                 ctypes.c_int(num_projs),
                 ctypes.c_int(num_rows),
                 ctypes.c_int(num_cols),
                 ctypes.c_int(num_grid),
                 ctypes.c_int(num_iter),
                 recon.ctypes.data_as(c_float_p))
    return recon
    

def data2tif(data, output_file=None, x_start=0,
             digits=5, axis=0, overwrite=True, 
             dtype='float32', data_min=None, data_max=None):
    """ 
    Write 3-D data to a stack of tif files.

    Parameters
    -----------
    output_file : str, optional
        Name of the output file.

    x_start : scalar, optional
        First index of the data on first dimension
        of the array.

    digits : scalar, optional
        Number of digits used for file indexing.
        For example if 4: test_XXXX.tiff
        
    axis : scalar, optional
        Imaages is read along that axis.
        
    overwrite: bool, optional
        if overwrite=True the existing files in the
        reconstruction folder will be overwritten
        with the new ones.
        
    dtype : bool, optional
        Export data type precision.
        
    data_min, data_max : scalar, optional
        User defined minimum and maximum values
        in the data that will be used to scale 
        the dataset when saving.
    
    Notes
    -----
    If file exists, saves it with a modified name.
    
    If output location is not specified, the data is
    saved inside ``recon`` folder where the input data
    resides. The name of the reconstructed files will
    be initialized with ``recon``
    """
    if output_file == None:
        output_file = "tmp/img_" 
    output_file =  os.path.abspath(output_file)
    dir_path = os.path.dirname(output_file)
        
    # Find max min of data for scaling
    if data_max is None:
        data_max = np.max(data)
    if data_min is None:
        data_min = np.min(data)
        
    if data_max < np.max(data):
        data[data>data_max] = data_max
    if data_min > np.min(data):
        data[data<data_min] = data_min
    
    # Remove TIFF extension if there is.
    if (output_file.endswith('tif') or
        output_file.endswith('tiff')) :
            output_file = output_file.split(".")[-2]
            
    # Create new folders.
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    # Select desired x from whole data.
    num_x, num_y, num_z = data.shape
    if axis == 0:
        x_end = x_start+num_x
    elif axis == 1:
        x_end = x_start+num_y
    elif axis == 2:
        x_end = x_start+num_z

    # Write data.
    file_index = ["" for x in range(digits)]
    for m in range(digits):
        file_index[m] = '0' * (digits - m - 1)
    ind = range(x_start, x_end)
    for m in range(len(ind)):
        for n in range(digits):
            if ind[m] < np.power(10, n + 1):
                file_body = output_file + file_index[n] + str(ind[m])
                file_name = file_body + '.tif'
                break

        # check if file exists.
        if not overwrite:
            if os.path.isfile(file_name):
                # generate new file unique name.
                indq = 1
                FLAG_SAVE = False
                while not FLAG_SAVE:
                    new_file_body = file_body + '-' + str(indq)
                    new_file_name = new_file_body + '.tif'
                    if not os.path.isfile(new_file_name):
                        FLAG_SAVE = True
                        file_name = new_file_name
                    else:
                        indq += 1

        if axis == 0:
            arr = data[m, :, :]
        elif axis == 1:
            arr = data[:, m, :]
        elif axis == 2:
            arr = data[:, :, m]

        if dtype is 'uint8':
            arr = ((arr*1.0 - data_min)/(data_max-data_min)*255).astype('uint8')
        elif dtype is 'uint16':
            arr = ((arr*1.0 - data_min)/(data_max-data_min)*65535).astype('uint16')
        elif dtype is 'float32':
            arr = arr.astype('float32')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            skimage_io.imsave(file_name, arr)
