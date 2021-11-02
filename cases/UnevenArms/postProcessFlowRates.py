from postProcessPlotSingleFunction import hemeCollateFunction, hemeFinalProfileBuilder

import numpy as np
import matplotlib.pyplot as plt

stepSize = 500
DX = 0.00002

ForkIn = np.array([[23,23,203]])

ForkOut = np.array([[10.1,13.0,244.5],[187.0,13.0,138.4]])

ForkOutP, ForkOutV = hemeCollateFunction("results/outlet.txt", stepSize, DX, ForkOut)

##plt.figure(1)
##plt.subplot(2,1,1)
##for io in range(np.shape(ForkIn)[0]):
##    plt.plot(ForkInV[io][1::3],'-', label='IOlet #'+str(io))
##plt.legend()
##
##plt.subplot(2,1,2)
##for io in range(np.shape(ForkIn)[0]):
##    plt.plot(ForkInP[io][1::3],'-', label='IOlet #'+str(io))
##plt.legend()

plt.figure(2)
plt.suptitle("Pipe Outlet")
plt.subplot(3,1,1)
for io in range(np.shape(ForkOut)[0]):
    plt.plot(ForkOutV[io][1::3],'-', label='IOlet #'+str(io))
#plt.title("Velocity")
#plt.xlabel("Simulation Steps")
plt.ylabel("Velocity [Units]")
plt.legend()

plt.subplot(3,1,2)
for io in range(np.shape(ForkOut)[0]):
    plt.plot(ForkOutP[io][1::3],'-', label='IOlet #'+str(io))
#plt.title("Pressure")
plt.xlabel("Simulation Steps")
plt.ylabel("Pressure [Units]")
plt.legend()

plt.subplot(3,1,3)
plt.plot(np.divide(ForkOutV[0][1::3],ForkOutV[1][1::3]),'-')
#plt.title("Velocity")
#plt.xlabel("Simulation Steps")
plt.ylabel("Velocity [Units]")
plt.legend()
plt.show()
