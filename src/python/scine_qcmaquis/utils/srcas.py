from typing import Any

from _dmrg import SRCASReal

from scine_qcmaquis import QCMaquis


class Srcas:
    """SRCAS wrapper"""

    def __init__(self, interface: QCMaquis):
        """Create SRCAS wrapper

        Parameters
        ----------
        interface : QCMaquis
            wrapper around the qcmaquis interface
        """
        self._srcas = SRCASReal(interface._parameters.get_parameters(), interface._dmrg.get_dmrg())
        self._srcas.printSRCASSettings()

    def run(self) -> None:
        """Run SRCAS"""
        self._srcas.run()

    def print_results(self) -> None:
        """Print SRCAS results"""
        self._srcas.printResults()

    def get_srcas_obj(self) -> Any:
        """Get srcas object

        Returns
        -------
        _srcas : Any
            The srcas object
        """
        return self._srcas

    # TODO: srcas class has to be modified to make this work
    # def results(self):
    #     return self._srcas.getDetTable()
