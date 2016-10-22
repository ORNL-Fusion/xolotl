import numpy as np


# Reading the data to be normalized
x = np.loadtxt('xdata.txt', usecols = (0,), unpack=True)
y = np.loadtxt('ydata.txt', usecols = (0,), unpack=True)

# Applying Simpson's 1/3 Rule to determine the integration
n = len(x)
Sum1 = 0.0
Sum2 = 0.0

for i in range(1,n-2,2):
    Sum1 = Sum1 + y[i]

for i in range(2,n-2,2):
    Sum2 = Sum2 + y[i]

I = (x[n-1] - x[0])/(3*n) * (y[0] + 4 * Sum1 + 2 * Sum2 + y[n-1])

print I

# Normalizing the data
y = y / I
np.savetxt('Normalized_ydata.txt',y)


# Checking if the normalization succeeded
Sum1 = 0.0
Sum2 = 0.0

for i in range(1,n-2,2):
    Sum1 = Sum1 + y[i]

for i in range(2,n-2,2):
    Sum2 = Sum2 + y[i]

I = (x[n-1] - x[0])/(3*n) * (y[0] + 4 * Sum1 + 2 * Sum2 + y[n-1])

print 1.0-I