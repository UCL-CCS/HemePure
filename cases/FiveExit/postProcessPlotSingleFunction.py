import numpy as np
import matplotlib.pyplot as plt

def hemeCollateFunction(path, stepSize, DX, lets):
    #stepSize = 100

    #DX = 0.05
    IOlets = DX*lets
    #IOlets = DX*np.array([[302.459,62.9992,4.78602],[352.459,62.9992,4.78601],[402.459,62.9992,4.786],[452.459,62.9992,4.786],[502.459,62.9992,4.786]])

    IOP = [[] for x in range(np.shape(IOlets)[0])]
    IOV = [[] for x in range(np.shape(IOlets)[0])]

    # steps gridX gridY gridZ velX velY velZ pressure

    print("starting to process", path)
    f = open(path, "r")
    #f = open("FluteFork/inletSecond.txt", "r")
    for i in range(22):
        next(f)
        
    stepCounter = 0

    outs=[]
    ioV = [[] for x in range(np.shape(IOlets)[0])]

    notStart = False

    for i in f:  
        i = np.float_(i.split())
        i = np.concatenate((i, np.linalg.norm(i[4:7])), axis=None)

        if int(i[0]) != stepCounter*stepSize:
            #print(stepCounter)
            stepCounter = stepCounter+1
            
            if notStart:       
                outs = np.array(outs)
                ioV = np.array(ioV)


                for io in range(np.shape(IOlets)[0]):
                    IOP[io] = IOP[io]+[np.max(ioV[io], axis=0)[7], np.average(ioV[io], axis=0)[7], np.min(ioV[io], axis=0)[7]]
                    IOV[io] = IOV[io]+[np.max(ioV[io], axis=0)[-1], np.average(ioV[io], axis=0)[-1], np.min(ioV[io], axis=0)[-1]]

            notStart=True
            
            outs=[]
            ioV = [[] for x in range(np.shape(IOlets)[0])]

        outs.append(i)

        dist = 1e10
        minIO = -1
        for io in range(np.shape(IOlets)[0]):
            if np.linalg.norm(i[1:4] - IOlets[io])<dist:
                dist = np.linalg.norm(i[1:4] - IOlets[io])
                minIO = io

        ioV[minIO].append(i)
          
    f.close()

    return IOP, IOV

def hemeFinalProfileBuilder(path, FinalStep, direction1,direction2, centre):
    # steps gridX gridY gridZ velX velY velZ pressure

    print("starting to process", path)
    f = open(path, "r")
    for i in range(22):
        next(f)
        
    stepCounter = 0
    
    ioV = []
    ioX = []

    for i in f:  
        i = np.float_(i.split())

        if int(i[0]) == FinalStep:
            if i[direction1+1] == centre:
                ioV.append(i[6])
                ioX.append(i[direction2+1])
        
    f.close()

    return ioV, ioX

##plt.figure(1)
##plt.subplot(1,2,1)
##plt.plot(P[0::3],'b.',P[1::3],'k-',P[2::3],'b.')
##
##plt.subplot(1,2,2)
##plt.plot(V[0::3],'r.',V[1::3],'k-',V[2::3],'r.')
##
##
##plt.figure(2)
##for io in range(np.shape(IOlets)[0]):
##    plt.plot(IOV[io][1::3],'-', label='IOlet #'+str(io))
##plt.legend()
##plt.show()

# normalised against inlet value.

# construct V profile??
