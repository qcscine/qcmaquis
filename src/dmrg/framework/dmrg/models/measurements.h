/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "dmrg/models/measurements/average.h"
#include "dmrg/models/measurements/local.h"
#include "dmrg/models/measurements/local_at.h"
#include "dmrg/models/measurements/correlations.h"
#include "dmrg/models/measurements/custom.h"
#include "dmrg/models/measurements/overlap.h"
#include "dmrg/models/measurements/entanglement.h"


template <class Matrix, class SymmGroup>
class measure_and_save {
public:
    measure_and_save(std::string const& rfile_, std::string const& archive_path_,
                     MPS<Matrix, SymmGroup> const& mps_, int eigenstate_=0)
    : rfile(rfile_)
    , archive_path(archive_path_)
    , eigenstate(eigenstate_)
    , mps(mps_)
    , rmps(mps)
    { }
    
    void operator()(measurement<Matrix, SymmGroup> & meas) const
    {
        maquis::cout << "Measuring " << meas.name() << std::endl;
        meas.eigenstate_index() = eigenstate;
        meas.evaluate(mps, rmps);
        storage::archive ar(rfile, "w");
        ar[archive_path] << meas;
    }
    
private:
    std::string rfile, archive_path;
    int eigenstate;
    MPS<Matrix, SymmGroup> const& mps;
    reduced_mps<Matrix, SymmGroup> rmps;
};


namespace detail {
    class name_not_in_list {
    public:
        name_not_in_list(std::vector<std::string> const& list_)
        : list(list_)
        { }
        
        template <class Matrix, class SymmGroup>
        bool operator() (measurement<Matrix, SymmGroup> const& term) const
        {
            return std::find(list.begin(), list.end(), term.name()) == list.end();
        }
        
    private:
        std::vector<std::string> const& list;
    };
}

template <class Matrix, class SymmGroup>
boost::ptr_vector<measurement<Matrix, SymmGroup> > &
operator<< (boost::ptr_vector<measurement<Matrix, SymmGroup> > & lhs,
            boost::ptr_vector<measurement<Matrix, SymmGroup> > const& rhs)
{
    lhs.insert(lhs.end(), rhs.begin(), rhs.end());
    return lhs;
}

template <class Matrix, class SymmGroup>
boost::ptr_vector<measurement<Matrix, SymmGroup> >
meas_sublist(boost::ptr_vector<measurement<Matrix, SymmGroup> > const& m,
             std::vector<std::string> const& meas_list)
{
    boost::ptr_vector<measurement<Matrix, SymmGroup> > sublist(m.clone());
    sublist.erase_if( detail::name_not_in_list(meas_list) );
    return sublist;
}

//template <class Matrix, class SymmGroup>
//class DMOverlapMeasurement : public Measurement_Term<Matrix, SymmGroup> {
//public:
//    typedef Measurement_Term<Matrix, SymmGroup> base;
//    
//    MPS<Matrix, SymmGroup> mps_ident;
//    std::vector<MPS<Matrix, SymmGroup> > overlaps_mps;
//    std::vector<std::string> labels;
//    
//protected:
//    virtual Measurement_Term<Matrix, SymmGroup> * do_clone() const
//    {
//        return new DMOverlapMeasurement<Matrix, SymmGroup>(*this);
//    }
//};


template <class Matrix, class SymmGroup>
boost::ptr_vector<measurement<Matrix, SymmGroup> >
overlap_measurements(BaseParameters const & parms, boost::optional<size_t> sweep = boost::none)
{
    /* Syntax for MEASURE_OVERLAP:
     *  (1) MEASURE_OVERLAP[obsname] = "/path/to/ckp.h5"
     *  (2) MEASURE_OVERLAP[obsname(sweep)] = "/path/to/ckp.h5"
     *
     * `obsname` is the name that will be given in archive output.
     * if `sweep` is prensent, the overlap will only be computed when the sweep number
     * matches the given one.
     */
    boost::ptr_vector<measurement<Matrix, SymmGroup> > meas;
    boost::regex expression("^MEASURE_OVERLAP\\[([a-zA-Z]+)(\\(([0-9]+)\\))?\\]$");
    boost::smatch what;
    for (BaseParameters::const_iterator it=parms.begin();it != parms.end();++it) {
        std::string lhs = it->key();
        if (boost::regex_match(lhs, what, expression)) {
            if (sweep && !what[3].matched) continue;
            if (sweep && what[3].matched && boost::lexical_cast<long>(what.str(3)) != sweep.get()) continue;
            
            std::string name = what.str(1), bra_chkp = it->value();
            meas.push_back( new measurements::overlap<Matrix, SymmGroup>(name, bra_chkp) );
        }
    }
    return meas;
}

#endif
