from postProcessPlotSingleFunction import hemeCollateFunction, hemeFinalProfileBuilder

import numpy as np
import matplotlib.pyplot as plt

stepSize = 500
DX = 0.000025

ForkIn = np.array([[23,23,203]])

ForkOut = np.array([[8.0,13.0,73.0],[51.5,13.0,148.5],[51.5,13.0,3.0],[95.0,13.0,123.5],[94.5,13.0,73.0]])

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
plt.plot(np.divide(ForkOutV[1][1::3],ForkOutV[0][1::3]),'-',label='1/0')
plt.plot(np.divide(ForkOutV[2][1::3],ForkOutV[0][1::3]),'-',label='2/0')
plt.plot(np.divide(ForkOutV[3][1::3],ForkOutV[0][1::3]),'-',label='3/0')
plt.plot(np.divide(ForkOutV[4][1::3],ForkOutV[0][1::3]),'-',label='4/0')
#plt.title("Velocity")
#plt.xlabel("Simulation Steps")
plt.ylabel("Velocity [Units]")
plt.legend()
plt.show()
