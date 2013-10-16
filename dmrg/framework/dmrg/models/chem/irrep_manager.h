/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef IRREP_MANAGER_H
#define IRREP_MANAGER_H

#include <vector>
#include <set>
#include <map>

#include <boost/shared_ptr.hpp>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/models/op_handler.h"
#include "dmrg/models/chem/chem_detail.h"
#include "dmrg/models/chem/pg_util.h"


template <class Matrix, class SymmGroup>
class PGSymmetryConverter : public HamiltonianTraits
{
    typedef typename PGDecorator<SymmGroup>::irrep_t irrep_t;

    typedef typename tag_type<Matrix, SymmGroup>::type tag_type;
    typedef typename op_t<Matrix, SymmGroup>::type op_t;

public:
    PGSymmetryConverter(std::vector<irrep_t> const & site_irreps_,
                        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_base,
                        boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_symm)
        : site_irreps(site_irreps_)
    {
        // collect all irreducible representations in a set
        std::set<irrep_t> irrep_set;
        for(typename std::vector<irrep_t>::const_iterator it=site_irreps.begin(); it != site_irreps.end(); ++it)
            irrep_set.insert(*it); 

        // convert set to vector (faster)
        for(typename std::set<irrep_t>::const_iterator it=irrep_set.begin(); it != irrep_set.end(); ++it)
            irrep_vector.push_back(*it);

        translation_table.resize(tag_handler_base->size(), irrep_vector.size()); 

        PGDecorator<SymmGroup> set_symm;
        for (std::size_t i=0; i < tag_handler_base->size(); ++i)
        {
            op_t base_op = tag_handler_base->get_op(i);
            for(typename std::vector<irrep_t>::const_iterator it=irrep_vector.begin(); it != irrep_vector.end(); ++it)
            {
                op_t modified = base_op;
                set_symm(modified.left_basis(), *it);
                set_symm(modified.right_basis(), *it);

                std::pair<tag_type, typename Matrix::value_type> symm_tag = tag_handler_symm->checked_register(modified);
                translation_table(i, *it) = symm_tag.first;

            }
        }

    }

    /*
    tag_type operator()(tag_type base_op, std::size_t site)
    {
        return modified_ops[base_op][irreps[site]];
    }

    tag_type register_op(typename op_t<Matrix, SymmGroup>::type op, tag_detail::operator_kind kind)
    {
        PGDecorator<SymmGroup> set_symm;

        // register point-group-modified operators in the tag_handler
        modified_ops.push_back( std::map<irrep_t, tag_type>() ); 
        tag_type ret = modified_ops.size() - 1;
        for (typename std::set<irrep_t>::iterator it = irrep_set.begin(); it != irrep_set.end(); ++it)
        {
            typename op_t<Matrix, SymmGroup>::type modified = op; 
            set_symm(modified.left_basis(), *it);
            set_symm(modified.right_basis(), *it);
            tag_type mod_id = tag_handler->register_op(modified, kind);

            assert (ret == modified_ops.size()-1);
            modified_ops[ret][*it] = mod_id;
        }
        return ret;
    }
    */

private:
    alps::numeric::matrix<irrep_t> translation_table;

    std::set<irrep_t> irrep_vector;
    std::vector<irrep_t> site_irreps;
};

#endif
