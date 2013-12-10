/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASURE_ON_MPS_H
#define MEASURE_ON_MPS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/models/lattice.h"
#include "dmrg/models/measurements.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/meas_eval.hpp"

#include <boost/shared_ptr.hpp>
#include <stdexcept>

template<class Matrix, class SymmGroup>
void measure_on_mps(MPS<Matrix, SymmGroup> const& mps, Lattice const& lat,
                    Measurements<Matrix, SymmGroup> const& meas,
                    std::string const& h5name, std::string basepath = std::string("/spectrum/results/"))
{
	if (basepath[basepath.size()-1] != '/')
        basepath += "/";
    
    if (meas.n_terms() > 0) {
        bool super_meas=meas.is_super_meas();
        
        boost::scoped_ptr<meas_eval::LocalMPSMeasurement<Matrix, SymmGroup> > local_measurement;
        if (!super_meas)
            local_measurement.reset( new meas_eval::LocalMPSMeasurement<Matrix, SymmGroup>(mps, lat) );
        
        // + storage::encode(meas[i].name)
        for (int i = 0; i < meas.n_terms(); ++i)
        {
            maquis::cout << "Calculating " << meas[i].name << std::endl;
            switch (meas[i].type)
            {
                case Measurement_Term<Matrix, SymmGroup>::Local:
                    assert(meas[i].operators.size() == 1  || meas[i].operators.size() == 2);
                    if (!super_meas && meas[i].operators.size() == 1) // Local measurements are fast and efficient!
                        local_measurement->site_term(meas[i].operators[0],
                                                     h5name, basepath, meas[i].name);
                    else
                        meas_eval::measure_local(mps, lat,
                                                 meas.identity_matrices(), meas.filling_matrices(),
                                                 meas[i], h5name, basepath, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::MPSBonds:
                    assert(meas[i].operators.size() == 2);
                    local_measurement->bond_term(meas[i].operators,
                                                 h5name, basepath, meas[i].name);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Average:
                    assert(meas[i].operators.size() == 1  || meas[i].operators.size() == 2);
                    meas_eval::measure_average(mps, lat,
                                               meas.identity_matrices(), meas.filling_matrices(),
                                               meas[i], h5name, basepath, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Overlap:
                    meas_eval::measure_overlap(mps, dynamic_cast<OverlapMeasurement<Matrix, SymmGroup> const & >(meas[i]).bra_ckp,
                                               h5name, basepath, meas[i].name);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Correlation:
                    meas_eval::measure_correlation(mps, lat, meas.identity_matrices(), meas.filling_matrices(),
                                                   meas[i], h5name, basepath, false, false, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::HalfCorrelation:
                    meas_eval::measure_correlation(mps, lat, meas.identity_matrices(), meas.filling_matrices(),
                                                   meas[i], h5name, basepath, true, false, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::CorrelationNN:
                    if (meas[i].operators.size() % 2 != 0)
                        throw std::runtime_error("Nearest neighbors correlators have to have even number of operators");
                    meas_eval::measure_correlation(mps, lat, meas.identity_matrices(), meas.filling_matrices(),
                                                   meas[i], h5name, basepath, false, true, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::HalfCorrelationNN:
                    if (meas[i].operators.size() % 2 != 0)
                        throw std::runtime_error("Nearest neighbors correlators have to have even number of operators");
                    meas_eval::measure_correlation(mps, lat, meas.identity_matrices(), meas.filling_matrices(),
                                                   meas[i], h5name, basepath, true, true, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Custom:
                    meas_eval::measure_custom(mps, lat,
                                              meas.identity_matrices(), meas.filling_matrices(), meas[i].custom_ops,
                                              h5name, basepath, meas[i].name);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::LocalAt:
                    meas_eval::measure_local_at(mps, lat, meas.identity_matrices(), meas.filling_matrices(),
                                                meas[i], h5name, basepath);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::DMOverlap:
                case Measurement_Term<Matrix, SymmGroup>::DMMultioverlap:
                    DMOverlapMeasurement<Matrix, SymmGroup> const & cast_meas = dynamic_cast<DMOverlapMeasurement<Matrix, SymmGroup> const & >(meas[i]);
                    bool is_multi_overlap = (meas[i].type == Measurement_Term<Matrix, SymmGroup>::DMMultioverlap);
                    meas_eval::dm_overlap(mps, cast_meas.mps_ident, cast_meas.overlaps_mps, cast_meas.labels, is_multi_overlap,
                                          h5name, basepath, meas[i].name);
                    break;
            }
        }
    }
}

#endif
