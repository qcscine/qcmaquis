/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
*               2021 Alberto Baiardi <abaiardi@ethz.ch>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#define BOOST_TEST_MODULE MPS_INITIALIZER_ELECTRONIC

#include "Fixtures/BenzeneFixture.h"
#include "Fixtures/H2Fixture.h"
#include "Fixtures/LiHFixture.h"

#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps_initializers_helper.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/generate_mpo.hpp"

#include <iostream>
#include <boost/test/included/unit_test.hpp>


typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
#ifdef HAVE_TwoU1
, TwoU1
#endif
#ifdef HAVE_SU2U1PG
, SU2U1PG
#endif
#ifdef HAVE_SU2U1
, SU2U1
#endif
> symmetries;

BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_MPS_Initializers_Electronic_H2, S, symmetries, H2Fixture)
{
  parametersH2.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  auto lattice = Lattice(parametersH2);
  int latticeSize = lattice.size();
  auto model = Model<matrix, S>(lattice, parametersH2);
  auto mpo = make_mpo(lattice, model);
  parametersH2.set("init_type", "hf");
  parametersH2.set("hf_occ", "4,1");
  auto mpsTHF = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energyTHF = expval(mpsTHF, mpo)/norm(mpsTHF);
  parametersH2.set("hf_occ", "3,2");
  auto mpsSHF = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energySHF = expval(mpsSHF, mpo)/norm(mpsSHF);
  parametersH2.set("init_type", "basis_state_generic");
  parametersH2.set("init_basis_state", "4,1");
  auto mpsTHF_g = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energyTHF_g = expval(mpsTHF_g, mpo)/norm(mpsTHF_g);
  BOOST_CHECK_CLOSE(energyTHF, energyTHF_g, 1.0E-10);
  parametersH2.set("init_basis_state", "3,2");
  auto mpsSHF_g = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energySHF_g = expval(mpsSHF_g, mpo)/norm(mpsSHF_g);
  BOOST_CHECK_CLOSE(energySHF, energySHF_g, 1.0E-10);
  parametersH2.set("init_type", "const");
  auto mpsC = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energyC = expval(mpsC, mpo)/norm(mpsC);
  parametersH2.set("init_type", "basis_state_generic_const");
  parametersH2.set("init_space", "5,5");
  auto mps_gc = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energy_gc = expval(mps_gc, mpo)/norm(mps_gc);
  BOOST_CHECK_CLOSE(energyC, energy_gc, 1.0E-10);
  parametersH2.set("init_space", "3,5");
  auto mpsSHF_gc = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersH2)));
  auto energySHF_gc = expval(mpsSHF_gc, mpo)/norm(mpsSHF_gc);
  BOOST_CHECK_CLOSE(energySHF, energySHF_gc, 1.0E-10);
}
  
BOOST_FIXTURE_TEST_CASE_TEMPLATE(Test_MPS_Initializers_Electronic_Benzene, S, symmetries, BenzeneFixture)
{
  parametersBenzene.set("symmetry", symm_traits::SymmetryNameTrait<S>::symmName());
  auto lattice = Lattice(parametersBenzene);
  int latticeSize = lattice.size();
  auto model = Model<matrix, S>(lattice, parametersBenzene);
  auto mpo = make_mpo(lattice, model);
  parametersBenzene.set("init_type", "hf");
  parametersBenzene.set("hf_occ", "4,4,4,1,1,1");
  auto mpsTHF = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energyTHF = expval(mpsTHF, mpo)/norm(mpsTHF);
  parametersBenzene.set("hf_occ", "3,3,3,2,2,2");
  auto mpsSHF = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energySHF = expval(mpsSHF, mpo)/norm(mpsSHF);
  parametersBenzene.set("init_type", "basis_state_generic");
  parametersBenzene.set("init_basis_state", "4,4,4,1,1,1");
  auto mpsTHF_g = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energyTHF_g = expval(mpsTHF_g, mpo)/norm(mpsTHF_g);
  BOOST_CHECK_CLOSE(energyTHF, energyTHF_g, 1.0E-10);
  parametersBenzene.set("init_basis_state", "3,3,3,2,2,2");
  auto mpsSHF_g = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energySHF_g = expval(mpsSHF_g, mpo)/norm(mpsSHF_g);
  if (symm_traits::SymmetryNameTrait<S>::symmName() == "2u1PG" || symm_traits::SymmetryNameTrait<S>::symmName() == "2u1") // as this is not true for SU2
    BOOST_CHECK_CLOSE(energySHF, energySHF_g, 1.0E-10);
  parametersBenzene.set("init_bond_dimension", 100); // so that we don't truncate anything --> otherwise reshuffling the order has an effect on the energy!
  parametersBenzene.set("init_type", "const");
  auto mpsC = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energyC = expval(mpsC, mpo)/norm(mpsC);
  parametersBenzene.set("init_type", "basis_state_generic_const");
  parametersBenzene.set("init_space", "5,5,5,5,5,5");
  auto mps_gc = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energy_gc = expval(mps_gc, mpo)/norm(mps_gc);
  BOOST_CHECK_CLOSE(energyC, energy_gc, 1.0E-10);
  parametersBenzene.set("init_space", "3,3,3,2,2,5");
  auto mpsSHF_gc = MPS<matrix, S>(latticeSize, *(model.initializer(lattice, parametersBenzene)));
  auto energySHF_gc = expval(mpsSHF_gc, mpo)/norm(mpsSHF_gc);
  if (symm_traits::SymmetryNameTrait<S>::symmName() == "2u1PG" || symm_traits::SymmetryNameTrait<S>::symmName() == "2u1") // as this is not true for SU2
    BOOST_CHECK_CLOSE(energySHF_gc, energySHF_g, 1.0E-10);  
  parametersBenzene.set("init_space", "3,3,3,2,2,2");
  parametersBenzene.set("orbital_order", "6,5,4,3,2,1"); // reshuffling the orbitals
  auto lattice_r = Lattice(parametersBenzene);
  auto model_r = Model<matrix, S>(lattice_r, parametersBenzene);
  auto mpo_r = make_mpo(lattice_r, model_r);
  auto mpsSHFr_gc = MPS<matrix, S>(latticeSize, *(model_r.initializer(lattice_r, parametersBenzene)));
  auto energySHFr_gc = expval(mpsSHFr_gc, mpo_r)/norm(mpsSHFr_gc);
  BOOST_CHECK_CLOSE(energySHFr_gc, energySHF_gc, 1.0E-10);
  parametersBenzene.set("init_space", "4,4,4,5,5,5");
  auto mpsTHFr_gc = MPS<matrix, S>(latticeSize, *(model_r.initializer(lattice_r, parametersBenzene)));
  auto energyTHFr_gc = expval(mpsTHFr_gc, mpo_r)/norm(mpsTHFr_gc);
  BOOST_CHECK_CLOSE(energyTHFr_gc, energyTHF, 1.0E-10);
  parametersBenzene.set("init_space", "5,5,5,1,1,1");
  auto mpsTHFr_gc2 = MPS<matrix, S>(latticeSize, *(model_r.initializer(lattice_r, parametersBenzene)));
  auto energyTHFr_gc2 = expval(mpsTHFr_gc2, mpo_r)/norm(mpsTHFr_gc2);
  BOOST_CHECK_CLOSE(energyTHFr_gc2, energyTHF, 1.0E-10);
}