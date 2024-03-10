/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "dmrg/mp_tensors/reduced_mps.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/site_operator_algorithms.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/models/lattice/lattice.h"
#include <alps/parser/xmlstream.h>

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include <stdexcept>

#include <regex>
#include <boost/static_assert.hpp>

/// TODO: 1) move data to new object measurement_result.
///       2) store only measurement description, i.e. no matrix, and pass model+lattice in evaluate method.

template<class Matrix, class SymmGroup>
class measurement {
public:
    typedef typename Matrix::value_type value_type;
    typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;

    measurement(std::string const& n="")
    : cast_to_real(true), is_super_meas(false), name_(n), eigenstate(0)
    {}

    virtual ~measurement() { }

    virtual void evaluate(MPS<Matrix, SymmGroup> const&, boost::optional<reduced_mps<Matrix, SymmGroup> const&> = boost::none) =0;
    template <class Archive>
    void save(Archive &) const;
    void write_xml(alps::oxstream &) const;
    virtual void print(std::ostream& os) const;

    std::string const& name() const { return name_; }
    int& eigenstate_index() { return eigenstate; }
    int eigenstate_index() const { return eigenstate; }

    void set_super_meas(Index<SymmGroup> const& phys_psi_);

    measurement* clone() const { return do_clone(); }

    std::vector<std::vector<int> > get_labels_num() { return labels_num; };
    std::vector<value_type> get_vec_results() { return vector_results; };

protected:
    virtual measurement* do_clone() const =0;

    bool cast_to_real, is_super_meas;

    // Labels in a string
    std::vector<std::string> labels;

    // Plain indices
    std::vector<std::vector<int> > labels_num;
    value_type result;
    std::vector<value_type> vector_results;

    Index<SymmGroup> phys_psi;

private:
    std::string name_;
    int eigenstate;
};

template<class Matrix, class SymmGroup>
inline measurement<Matrix, SymmGroup>* new_clone( const measurement<Matrix, SymmGroup>& m )
{
    return m.clone();
}

template <class Archive, typename T>
void save_val_at_index(Archive & ar, std::string const& archive_path, T const& val, std::size_t eig)
{
    std::vector<T> vals;
    if (ar.is_data(archive_path.c_str())) ar[archive_path] >> vals;
    vals.resize( std::max(vals.size(), eig+1) );
    vals[eig] = val;
    ar[archive_path] << vals;
}

template<class Matrix, class SymmGroup>
template <class Archive>
void measurement<Matrix, SymmGroup>::save(Archive & ar) const
{
    if (vector_results.size() > 0) {
        if (cast_to_real) {
            save_val_at_index(ar, storage::encode(name()) + std::string("/mean/value"),  maquis::real(vector_results), eigenstate_index());
        } else {
            save_val_at_index(ar, storage::encode(name()) + std::string("/mean/value"), vector_results, eigenstate_index());
        }
        if (labels.size() > 0)
            ar[storage::encode(name()) + std::string("/labels")] << labels;
        // Save also the numerical values of the labels for those who do not want to parse the strings on reading
        if (labels_num.size() > 0)
            ar[storage::encode(name()) + std::string("/labels_num")] << labels_num;
    } else {
        if (cast_to_real) {
            save_val_at_index(ar, storage::encode(name()) + std::string("/mean/value"), maquis::real(result), eigenstate_index());
        } else {
            save_val_at_index(ar, storage::encode(name()) + std::string("/mean/value"), result, eigenstate_index());
        }
    }
}

template<class Matrix, class SymmGroup>
void measurement<Matrix, SymmGroup>::write_xml(alps::oxstream & out) const
{
    if (labels.size() > 0) {
        out << alps::start_tag("VECTOR_AVERAGE") << alps::attribute("name", name());
        for (int i=0; i<labels.size(); ++i) {
            out << alps::start_tag("SCALAR_AVERAGE");
            out << alps::attribute("indexvalue",labels[i]) << alps::no_linebreak;

            if (cast_to_real) {
                out << alps::start_tag("MEAN") << alps::no_linebreak << maquis::real(vector_results[i]) << alps::end_tag("MEAN");
            } else {
                out << alps::start_tag("MEAN") << alps::no_linebreak << vector_results[i] << alps::end_tag("MEAN");
            }
            out << alps::end_tag("SCALAR_AVERAGE");

        }
        out << alps::end_tag("VECTOR_AVERAGE");
    } else {
        out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name", name()) << alps::no_linebreak;
        out << alps::start_tag("MEAN") <<  alps::no_linebreak << ((cast_to_real) ?  maquis::real(result) : result) << alps::end_tag("MEAN");
        out << alps::end_tag("SCALAR_AVERAGE");
    }
}


template<class Matrix, class SymmGroup>
void measurement<Matrix, SymmGroup>::set_super_meas(Index<SymmGroup> const& phys_psi_)
{
    phys_psi = phys_psi_;
    is_super_meas = true;
}

template<class Matrix, class SymmGroup>
void measurement<Matrix, SymmGroup>::print(std::ostream& os) const
{
    os << "MEASURE[" << name_ << "]";
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, measurement<Matrix, SymmGroup> const& m)
{
    m.print(os);
    return os;
}


/// UTILITIES

template<class BlockMatrix>
bool is_hermitian_meas(std::vector<BlockMatrix> const & ops)
{
    return all_true(ops.begin(), ops.end(),
                    boost::bind(static_cast<bool (*)(BlockMatrix const&)>(&is_hermitian), _1));
    return true;
}

template<class BlockMatrix>
bool is_hermitian_meas(std::vector<std::pair<std::vector<BlockMatrix>, bool> > const & ops)
{
    bool is_herm = true;
    for (int i=0; i<ops.size() && is_herm; ++i)
        is_herm = is_hermitian_meas(ops[i].first);
    return is_herm;
}
// Create label strings from indices including orbital reordering
inline std::string label_string(const Lattice& lat, const std::vector<Lattice::pos_t> & indices)
{
    std::string s;
    for (auto it = indices.begin(); it != indices.end(); it++)
        s += "( " + boost::lexical_cast<std::string>(lat.get_prop<std::string>("label", *it))
          + (it != std::prev(indices.end()) ? " ) -- " : " )") ;
    return s;
}

inline std::vector<std::string> label_strings (const Lattice& lat, const std::vector<std::vector<Lattice::pos_t> >& labels)
{
    std::vector<std::string> ret(labels.size());
    std::transform(labels.begin(), labels.end(), ret.begin(),
        [&](const std::vector<Lattice::pos_t> & l) { return label_string(lat, l); });
    return ret;
}

// Create label strings from indices WITHOUT reordering
inline std::string label_string(const std::vector<Lattice::pos_t> & indices)
{
    std::string s;
    for (auto it = indices.begin(); it != indices.end(); it++)
        s += "( " + boost::lexical_cast<std::string>(*it) + (it != std::prev(indices.end()) ? " ) -- " : " )") ;
    return s;
}
inline std::vector<std::string> label_strings(const std::vector<std::vector<Lattice::pos_t> >& labels)
{
    std::vector<std::string> ret(labels.size());
    std::transform(labels.begin(), labels.end(), ret.begin(),
        [&](const std::vector<Lattice::pos_t> & l) { return label_string(l);});
    return ret;
}

// Order orbital labels according to their ordering in the lattice
template <class T>
inline std::vector<T> order_labels(const Lattice& lat, const std::vector<T> & labels)
{
    // std::vector<T> ret(labels.size());
    // std::transform(labels.begin(), labels.end(), ret.begin(), [&](T i){return lat.get_prop<T>("label_int", i);});
    // return ret;
    return labels;
}


#endif
