import numpy as np

#Since the system is isotropic, shear moduli can be computed from various components of the tensor. 
borntensor = np.loadtxt('borntensor.log',skiprows=1)

#First component is the xyxy component
shearmod1 = np.mean(borntensor[:,6])

#Second component comes from the xxxx and xxyyy component
shearmod2 = 0.5*np.mean(borntensor[:,1]-borntensor[:,2])

#Check if they're the same up some significant figures
print(shearmod1,shearmod2)
