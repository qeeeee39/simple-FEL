import matplotlib.pyplot as plt
import numpy as np

f  = 1.3e9
pi = 3.14
c  = 3e8
k  = 2*pi*f/c

phi = 5.5

Ef = 6804.928
de = 6804.928 - 542.697

R56 = Ef / (k*de*np.tan(phi))

print(R56)