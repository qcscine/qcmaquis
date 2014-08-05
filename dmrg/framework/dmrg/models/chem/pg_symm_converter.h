/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef PG_SYMM_CONVERTER_H
#define PG_SYMM_CONVERTER_H

#include <vector>
#include <set>
#include <map>

#include <boost/shared_ptr.hpp>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/models/op_handler.h"
#include "dmrg/models/chem/pg_util.h"

template <class Matrix, class SymmGroup>
class PGSymmetryConverter
{
public:
    PGSymmetryConverter(std::vector<int> const & site_irreps_) {}
    void convert_tags_to_symm_tags(MPO<Matrix, SymmGroup> & mpo_in) {} 
};

template <class Matrix>
class PGSymmetryConverter<Matrix, TwoU1PG>
{
    typedef TwoU1PG::subcharge subcharge;

    typedef tag_detail::tag_type tag_type;
    typedef typename OPTable<Matrix, TwoU1PG>::op_t op_t;

public:
    PGSymmetryConverter(std::vector<subcharge> const & site_irreps_)
        : site_irreps(site_irreps_)
    {
        // collect all irreducible representations in a set
        std::set<subcharge> irrep_set;
        for(typename std::vector<subcharge>::const_iterator it=site_irreps.begin(); it != site_irreps.end(); ++it)
            irrep_set.insert(*it); 

        // convert set to vector (faster)
        for(typename std::set<subcharge>::const_iterator it=irrep_set.begin(); it != irrep_set.end(); ++it)
            irrep_vector.push_back(*it);
    }

    void convert_tags_to_symm_tags(MPO<Matrix, TwoU1PG> & mpo_in)
    {

        alps::numeric::matrix<subcharge> translation_table;

        /*** Set up a new OPTable for the symmetry-enabled operators via TagHandler ****/
        // TODO: assert that op_table is valid for all tensors in the MPO!
        boost::shared_ptr<OPTable<Matrix, TwoU1PG> > op_table = mpo_in[0].get_operator_table();

        TagHandler<Matrix, TwoU1PG> symm_handler;

        subcharge highest_irrep = *std::max_element(irrep_vector.begin(), irrep_vector.end());
        translation_table.resize(op_table->size(), highest_irrep+1); 

        /*** Create the tag translation table ****/
        PGDecorator<TwoU1PG> set_symm;
        for (std::size_t i = 0; i < op_table->size(); ++i)
        {
            op_t base_op = (*op_table)[i];
            for(typename std::vector<subcharge>::const_iterator it=irrep_vector.begin(); it != irrep_vector.end(); ++it)
            {
                op_t modified(set_symm(base_op.basis(), *it));
                for (std::size_t p = 0; p < modified.n_blocks(); ++p)
                    modified[p] = base_op[p];

                std::pair<tag_type, typename Matrix::value_type> symm_tag = symm_handler.checked_register(modified, tag_detail::bosonic);
                translation_table(i, *it) = symm_tag.first;
            }
        }

        typedef typename MPOTensor<Matrix, TwoU1PG>::row_proxy row_proxy;
        
        // TODO: traverse by columns instead of rows
        /*** Perform tag translation to symmetry tags ****/
        for (std::size_t i = 0; i < mpo_in.length(); ++i)
        {
            mpo_in[i].operator_table = symm_handler.get_operator_table();
            for (typename MPOTensor<Matrix, TwoU1PG>::index_type b1 = 0; b1 < mpo_in[i].row_dim(); ++b1)
            {
                row_proxy cur_row = mpo_in[i].row(b1);
                for (typename row_proxy::const_iterator row_it = cur_row.begin(); row_it != cur_row.end(); ++row_it)
                {
                    typename MPOTensor<Matrix, TwoU1PG>::index_type b2 = row_it.index();

                    // access the mpo tensor at (b1, b2)
                    typename MPOTensor<Matrix, TwoU1PG>::internal_value_type access = mpo_in[i].col_tags(b1, b2);
                    tag_type base_tag = access.first;
                    tag_type symm_tag = translation_table(base_tag, site_irreps[i]);

                    // replace 'base_tag' with a tag corresponding to the same operator, but with the symmetry
                    // irreducible representation of site 'i', which is site_irreps[i]
                    mpo_in[i].col_tags(b1, b2) = typename MPOTensor<Matrix, TwoU1PG>::internal_value_type(symm_tag, access.second);
                }
            }
        }
    }

private:
    std::vector<subcharge> irrep_vector; // a list of all irreps present
    std::vector<subcharge> site_irreps;  // irreps of the sites
};

#endif
