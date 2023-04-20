/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPOTENSOR_DETAIL_H
#define MPOTENSOR_DETAIL_H

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>

template<class Matrix, class SymmGroup>
class MPOTensor;

namespace MPOTensor_detail
{
    template <class T, bool C>
    struct const_type { typedef T type; };

    template <class T>
    struct const_type<T, true> { typedef const T type; };

    template <class Matrix, class SymmGroup, bool Const>
    class term_descriptor
    {
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::internal_value_type internal_value_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::op_table_ptr op_table_ptr;

    public:
        term_descriptor() {}
        term_descriptor(typename const_type<internal_value_type, Const>::type & term_descs,
                        op_table_ptr op_tbl_)
            : operator_table(op_tbl_), term_descriptors(term_descs) {}

        std::size_t size() const { return term_descriptors.size(); }
        typename const_type<op_t, Const>::type & op(std::size_t i=0) { return (*operator_table)[term_descriptors[i].first]; }
        typename const_type<value_type, Const>::type & scale(std::size_t i=0) { return term_descriptors[i].second; }

    private:
        typename const_type<internal_value_type, Const>::type & term_descriptors;
        op_table_ptr operator_table;
    };

    template <class ConstIterator>
    class IteratorWrapper : public std::iterator<std::forward_iterator_tag, typename std::iterator_traits<ConstIterator>::value_type>
    {
        typedef ConstIterator internal_iterator;

    public:
        typedef IteratorWrapper<ConstIterator> self_type;
        typedef typename std::iterator_traits<internal_iterator>::value_type value_type;

        IteratorWrapper(internal_iterator i) : it_(i) { }

        void operator++() { ++it_; }
        void operator++(int) {it_++; }
        bool operator!=(self_type const & rhs) { return it_ != rhs.it_; }

        value_type index() const { return *it_; }
        value_type operator*() const {
            throw std::runtime_error("direct MPOTensor access via row iterators currently not implemented\n");
            return *it_;
        }

    private:
       internal_iterator it_;
    };

    template <class ConstIterator>
    class row_proxy : public std::pair<ConstIterator, ConstIterator>
    {
        typedef ConstIterator internal_iterator;
        typedef std::pair<internal_iterator, internal_iterator> base;

    public:
        typedef IteratorWrapper<ConstIterator> const_iterator;
        row_proxy(internal_iterator b, internal_iterator e) : base(b, e) { }

        const_iterator begin() const { return const_iterator(base::first); }
        const_iterator end() const { return const_iterator(base::second); }
    };

    using namespace boost::tuples;

    template<class Tuple>
    struct row_cmp
    {
        bool operator() (Tuple const & i, Tuple const & j) const
        {
            if (get<0>(i) < get<0>(j))
                return true;
            else if (get<0>(i) > get<0>(j))
                return false;
            else
                return get<1>(i) < get<1>(j);
        }
    };

    template<class Tuple>
    struct col_cmp
    {
        bool operator() (Tuple const & i, Tuple const & j) const
        {
            if (get<1>(i) < get<1>(j))
                return true;
            else if (get<1>(i) > get<1>(j))
                return false;
            else
                return get<0>(i) < get<0>(j);
        }
    };

    class Hermitian
    {
        typedef std::size_t index_type;

        friend Hermitian operator * (Hermitian const &, Hermitian const &);

    public:
        Hermitian(index_type ld, index_type rd)
        {
            LeftHerm.resize(ld);
            RightHerm.resize(rd);
            LeftPhase = std::vector<int>(ld, 1);
            RightPhase = std::vector<int>(rd, 1);

            index_type z=0;
            std::generate(LeftHerm.begin(), LeftHerm.end(), boost::lambda::var(z)++);
            z=0;
            std::generate(RightHerm.begin(), RightHerm.end(), boost::lambda::var(z)++);
        }

        Hermitian(std::vector<index_type> const & lh,
                  std::vector<index_type> const & rh,
                  std::vector<int> const & lp,
                  std::vector<int> const & rp)
        : LeftHerm(lh), RightHerm(rh), LeftPhase(lp), RightPhase(rp)
        {}

        bool left_skip(index_type b1) const { return LeftHerm[b1] < b1; }
        bool right_skip(index_type b2) const { return RightHerm[b2] < b2; }

        index_type  left_conj(index_type b1) const { return  LeftHerm[b1]; }
        index_type right_conj(index_type b2) const { return RightHerm[b2]; }

        std::size_t left_size() const { return LeftHerm.size(); }
        std::size_t right_size() const { return RightHerm.size(); }

        int left_phase(std::size_t i) const { return LeftPhase[i]; }
        int right_phase(std::size_t i) const { return RightPhase[i]; }

    private:
        std::vector<index_type> LeftHerm;
        std::vector<index_type> RightHerm;

        std::vector<int> LeftPhase;
        std::vector<int> RightPhase;
    };

    inline Hermitian operator * (Hermitian const & a, Hermitian const & b)
    {
        return Hermitian(a.LeftHerm, b.RightHerm, a.LeftPhase, b.RightPhase);
    }

    template <class Matrix, class SymmGroup>
    symm_traits::disable_if_su2_t<SymmGroup, int> get_spin(MPOTensor<Matrix, SymmGroup> const & mpo,
                                                                                   typename MPOTensor<Matrix, SymmGroup>::index_type k, bool left)
    {
        return 0;
    }

    template <class Matrix, class SymmGroup>
    symm_traits::enable_if_su2_t<SymmGroup, int> get_spin(MPOTensor<Matrix, SymmGroup> const & mpo,
                                                                                  typename MPOTensor<Matrix, SymmGroup>::index_type k, bool left)
    {
        if (left)
        return mpo.left_spin(k).get();
        else
        return mpo.right_spin(k).get();
    }
}

#endif
