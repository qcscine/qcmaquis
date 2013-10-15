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

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/models/op_handler.h"
#include "dmrg/models/chem/chem_detail.h"
#include "dmrg/models/chem/pg_util.h"


template <class Matrix, class SymmGroup>
class IrrepManager : public HamiltonianTraits
{
    typedef typename PGDecorator<SymmGroup>::irrep_t irrep_t;

    typedef typename tag_type<Matrix, SymmGroup>::type tag_type;

public:
    IrrepManager(std::vector<irrep_t> const & irreps, boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_)
        : tag_handler(tag_handler_)
    {
        // collect all irreducible representations in a set
        for(typename std::vector<irrep_t>::const_iterator it=irreps.begin(); it != irreps.end(); ++it)
            irrep_set.insert(*it); 
    }

    tag_type operator()(tag_type base_op, irrep_t irrep)
    {
        return modified_ops[base_op][irrep];
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

private:
    std::vector<typename op_t<Matrix, SymmGroup>::type> base_ops;
    std::vector< std::map<irrep_t, tag_type> > modified_ops;
    std::set<irrep_t> irrep_set;
    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
};

#endif
