from _dmrg import SRCASReal

from scine_qcmaquis import MaquisDmrg


class Srcas:
    def __init__(self, interface: MaquisDmrg):
        self._srcas = SRCASReal(interface._parameters.get_parameters(), interface._dmrg.get_dmrg())

    def run(self):
        self._srcas.run()

    def print(self):
        self._srcas.printSRCASSettings()
        self._srcas.printResults()


if __name__ == "__main__":
    dmrg = MaquisDmrg()
    dmrg.set_parameter("symmetry", "2u1pg")
    dmrg.init_dmrg("checkpoint_n2_triplet.2.2.h5", 6, 6, 2)

    print("hihi")
    srcas = Srcas(dmrg)
    srcas.run()
    srcas.print()
