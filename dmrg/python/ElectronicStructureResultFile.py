import h5py
import argparse
import numpy as np

class ElectronicStructureResultFile:

    def __init__(self, resultFileName):
        self.__fileName = resultFileName
        self.__hdf5File = h5py.File(self.__fileName, 'r')
        self.__nOrbitals = self.__hdf5File['parameters']['L'][()]

    def getOneParticleRdm(self):
        oneParticleRdmMatrix = np.zeros([self.__nOrbitals, self.__nOrbitals])
        mapContainer = self.__hdf5File['spectrum']['results']['oneptdm']["mean"]["value"]
        for idx, i in enumerate(self.__hdf5File['spectrum']['results']['oneptdm']["labels_num"]):
            oneParticleRdmMatrix[i[0], i[1]] = mapContainer[0][idx]
            oneParticleRdmMatrix[i[1], i[0]] = mapContainer[0][idx]
        return oneParticleRdmMatrix

    def getNaturalOrbitalOccupationNumber(self, rdm = None):
        if rdm is None:
            rdm = self.getOneParticleRdm()
        eigenvalues, eigenvectors = np.linalg.eig(rdm)
        eigenvalues[::-1].sort()
	return eigenvalues

if __name__ == "__main__":
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--result", type=str, help="Name of the hdf5 file")
    args = parser.parse_args()
    resFileName = args.result
    resultFileObject = ElectronicStructureResultFile(resFileName)
    rdm = resultFileObject.getOneParticleRdm()
    occ = resultFileObject.getNaturalOrbitalOccupationNumber(rdm = rdm)
