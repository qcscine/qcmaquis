import TDDMRGResultFile
import matplotlib.pyplot as plt
import numpy as np

fileObject = TDDMRGResultFile.TDDMRGMeasurement("measure.27.1.resfile.h5")
listOfExcitations = fileObject.extractModeExcitationDegree(0)
numberOfModes = fileObject.getNumberOfModes()
numberOfSweeps = fileObject.getNumberOfSweeps()

fit, ax = plt.subplots()

xRange = np.arange(0, numberOfSweeps, 1.)

for i in range(numberOfModes):
    ax.plot(xRange, listOfExcitations[i], label="Mode "+str(i))

plt.legend()
plt.tight_layout()
plt.show()
