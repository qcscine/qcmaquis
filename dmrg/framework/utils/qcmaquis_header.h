#include <iostream>
namespace maquis
{
void qcmaquis_header(bool parallel, int nthreads_mpi){
        std::cout << "  SCINE QCMaquis \n"
                  << "  Quantum Chemical Density Matrix Renormalization group\n"
                  << "  available from https://scine.ethz.ch/download/qcmaquis\n"
                  << "  based on the ALPS MPS codes from http://alps.comp-phys.org/\n"
                  << "  copyright (c) 2015-2018 Laboratory of Physical Chemistry, ETH Zurich\n"
                  << "  copyright (c) 2012-2016 by Sebastian Keller\n"
                  << "  copyright (c) 2016-2018 by Alberto Baiardi, Leon Freitag, \n"
                  << "                             Stefan Knecht, Yingjin Ma \n"
                  << "  copyright (c) 2018-     by Alberto Baiardi, Leon Freitag, \n"
                  << "                             Stefan Knecht \n"
                  << "  for details see the publication: \n"
                  << "  S. Keller et al., J. Chem. Phys. 143, 244118 (2015)\n"
                  << std::endl;
                  if(parallel){
                    std::cout << "  MPI parallel run with "<<  nthreads_mpi << " threads \n"
                    << std::endl;
                  }
    }
}
