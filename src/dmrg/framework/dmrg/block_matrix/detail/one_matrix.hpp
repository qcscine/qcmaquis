/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
 *
 *****************************************************************************/

#ifndef MAQUIS_ONE_MATRIX_HPP
#define MAQUIS_ONE_MATRIX_HPP

namespace maquis {
    namespace dmrg {

    // Dummy Matrix of constant size 1 for MPO construction and compression
    template <typename T>
    class one_matrix {
    public:
        typedef T                       value_type;
        typedef T&                      reference;
        typedef T const&                const_reference;
        typedef std::size_t             size_type;
        typedef std::ptrdiff_t          difference_type;  
        // TODO: Introduce iterator classes that support *, !=, ++, ...
        typedef reference               element_iterator;
        typedef const_reference         const_element_iterator;

        explicit one_matrix(size_type rows = 1, size_type cols = 1, T init_value = T()) {
            assert(cols==1 && rows==1);
            val_ = init_value;
        }

        void swap(one_matrix & r) { std::swap((*this)(0,0), r(0,0)); }

        friend void swap(one_matrix & x, one_matrix & y)
        {
            x.swap(y);
        }

        inline value_type& operator()(const size_type i, const size_type j) {
            assert(i==0 && j==0);
            return val_;
        }
        inline value_type const& operator()(const size_type i, const size_type j) const {
            assert(i==0 && j==0);
            return val_;
        }

        inline bool operator == (one_matrix const& rhs) const { return this->val_ == rhs(0,0); }

        inline one_matrix<T>& operator += (one_matrix const& rhs) { this->val_ += rhs(0,0); return *this; }
        inline one_matrix<T>& operator -= (one_matrix const& rhs) { this->val_ -= rhs(0,0); return *this; }

        template <typename T2>
        inline one_matrix<T>& operator *= (T2 const& t) { this->val_ *= t; return *this; }

        template <typename T2>
        inline one_matrix<T>& operator /= (T2 const& t) { this->val_ /= t; return *this; }

        inline bool empty() const { return false; }

        inline size_type num_rows() const { return 1; }
        inline size_type num_cols() const { return 1; }

        std::pair<element_iterator, element_iterator> elements() { return std::make_pair(val_, NULL); }
        std::pair<const_element_iterator, const_element_iterator> elements() const { return std::make_pair(val_, NULL); }

        //MemoryBlock const& get_values() const;
        //MemoryBlock & get_values();

    private:

        T val_;
    };

    }
}

template <typename T>
inline std::size_t num_rows(maquis::dmrg::one_matrix<T> const & m) { return 1; }

template <typename T>
inline std::size_t num_cols(maquis::dmrg::one_matrix<T> const & m) { return 1; }

template <typename T>
inline void gemm(maquis::dmrg::one_matrix<T> const & a, maquis::dmrg::one_matrix<T> const & b, maquis::dmrg::one_matrix<T> & c)
{
    c(0,0) = a(0,0) * b(0,0);
} 

namespace maquis {
    namespace dmrg {

    template <typename T>
    const one_matrix<T> operator + (one_matrix<T> m1, one_matrix<T> const& m2)
        { return one_matrix<T>(1,1, m1(0,0) + m2(0,0)); }

    template <typename T>
    const one_matrix<T> operator - (one_matrix<T> m1, one_matrix<T> const& m2)
        { return one_matrix<T>(1,1, m1(0,0) - m2(0,0)); }

    template <typename T>
    const one_matrix<T> operator - (one_matrix<T> a)
        { return one_matrix<T>(1,1, -a(0,0)); }

    template<typename T>
    const one_matrix<T> operator * (one_matrix<T> const& m1, one_matrix<T> const& m2)
        { return one_matrix<T>(1,1, m1(0,0) * m2(0,0)); }

    template<typename T>
    std::size_t size_of(one_matrix<T> const & m) { return 1; }

    template <typename T>
    std::ostream& operator << (std::ostream& o, one_matrix<T> const& m) { o << m(0,0); return o; }

    }
}

#endif //MAQUIS_ONE_MATRIX_HPP
