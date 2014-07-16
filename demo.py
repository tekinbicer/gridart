#!/usr/bin/env python
# -*- coding: utf-8 -*-
from gridart import art
import numpy as np

# Import data
data = np.load('data/data.npy')
theta = np.load('data/theta.npy')

# Perform reconstruction
recon = art.run(data, theta)

# Export data as TIFs.
art.data2tif(recon, 'tmp/test_')
