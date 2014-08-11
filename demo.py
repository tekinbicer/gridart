#!/usr/bin/env python
# -*- coding: utf-8 -*-
from gridart import art
import numpy as np
import matplotlib.pylab as plt

# Import data
#data = np.load('/Users/bicer/Projects/data/gridart/2-ID-E-data.npy')
#theta = np.load('/Users/bicer/Projects/data/gridart/2-ID-E-theta.npy')
data = np.load('./data/data.npy')
theta = np.load('./data/theta.npy')

# Perform reconstruction
recon = art.run(data, theta)


# print sizes
print 'data shape  :' + str(data.shape)
print 'recon shape :' + str(recon.shape)

#np.save("orig.2-ID-E.i1.npy",recon);
np.save("recon.npy",recon);

#plt.figure()
#plt.imshow(data[:,3,:])
#plt.show()

#plt.figure()
#plt.imshow(recon[3,:,:])
#plt.show()

# Export data as TIFs.
#art.data2tif(recon, 'tmp/test_')
