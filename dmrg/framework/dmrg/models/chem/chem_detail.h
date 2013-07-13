/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_CHEM_DETAIL_H
#define QC_CHEM_DETAIL_H

struct HamiltonianTraits {

    template <class M>
    struct hamterm_t { typedef typename Hamiltonian<M, TwoU1>::hamterm_t type; };

    template <class M>
    struct hamtagterm_t { typedef typename Hamiltonian<M, TwoU1>::hamtagterm_t type; };

    template <class M>
    struct op_t { typedef typename Hamiltonian<M, TwoU1>::op_t type; };

    template <class M>
    struct op_tag_t { typedef typename OPTagTable<M, TwoU1>::op_tag_t type; };

    template <class M>
    struct op_pair_t { typedef typename generate_mpo::Operator_Term<M, TwoU1>::op_pair_t type; };

    template <class M>
    struct op_tag_pair_t { typedef typename generate_mpo::Operator_Tag_Term<M, TwoU1>::op_pair_t type; };
};

namespace chem_detail {

    class IndexTuple : public NU1Charge<4>
    {
    public:
        IndexTuple() {}
        IndexTuple(int i, int j, int k, int l) {
            (*this)[0] = i; (*this)[1] = j; (*this)[2] = k; (*this)[3] = l;
        }
    };

    IndexTuple align(int i, int j, int k, int l) {
        if (i<j) std::swap(i,j);
        if (k<l) std::swap(k,l);
        if (i<k) { std::swap(i,k); std::swap(j,l); }
        return IndexTuple(i,j,k,l);
    }
    
    IndexTuple align(IndexTuple const & rhs) {
        return align(rhs[0], rhs[1], rhs[2], rhs[3]);
    }

    int sign(IndexTuple const & idx)
    {
        int inv_count=0, n=4;
        for(int c1 = 0; c1 < n - 1; c1++)
            for(int c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;  

        return 1 - 2 * (inv_count % 2);
    }

    std::ostream& operator<<(std::ostream & os, IndexTuple const & c) {
        os << "<";
        for (int i = 0; i < 4; ++i) {
            os << c[i];
            if (i+1 < 4)
                os << ",";
        }
        os << ">";
        return os;
    }

    class TermTuple : public NU1Charge<8>
    {
    public:
        TermTuple() {}
        TermTuple(IndexTuple const & a, IndexTuple const & b) {
            for (int i=0; i<4; i++) { (*this)[i] = a[i]; (*this)[i+4] = b[i]; }
        }
    };

    template <typename M>
    class ChemHelper : HamiltonianTraits
    {
    public:
        typedef typename M::value_type value_type;
        typedef typename op_tag_t<M>::type tag_t;

        ChemHelper(BaseParameters & parms, Lattice const & lat,
                   tag_t ident_, tag_t fill_, boost::shared_ptr<OPTagTable<M, TwoU1> > op_table_) 
            : ident(ident_), fill(fill_), op_table(op_table_)
        {
            this->parse_integrals(parms, lat);

            for (std::size_t m=0; m < matrix_elements.size(); ++m) {
                IndexTuple pos;
                std::copy(idx_[m].begin(), idx_[m].end(), pos.begin());
                coefficients[pos] = matrix_elements[m];
            }
        }

        std::vector<value_type> & getMatrixElements() { return matrix_elements; }
        
        int idx(int m, int pos) { return idx_[m][pos]; }

        void commit_three_terms(std::vector<typename hamtagterm_t<M>::type> & tagterms) {
            for (typename std::map<TermTuple, typename hamtagterm_t<M>::type>::const_iterator it = three_terms.begin();
                    it != three_terms.end(); ++it)
                tagterms.push_back(it->second);
        }

        void add_term(std::vector<typename hamtagterm_t<M>::type> & tagterms,
                           value_type scale, int s, int p1, int p2, tag_t op_i, tag_t op_k, tag_t op_l, tag_t op_j) {

            typename hamtagterm_t<M>::type
            term = TermMaker<M>::three_term(ident, fill, scale, s, p1, p2, op_i, op_k, op_l, op_j, op_table);
            TermTuple id(IndexTuple(s,s,p1,p2),IndexTuple(op_i,op_k,op_l,op_j));
            if (three_terms.count(id) == 0) {
                three_terms[id] = term;
            }
            else
                three_terms[id].scale += term.scale;
    
        }

        void add_term(std::vector<typename hamtagterm_t<M>::type> & tagterms,
                           int i, int k, int l, int j, tag_t op_i, tag_t op_k, tag_t op_l, tag_t op_j)
        {
            // Collapse terms with identical operators and different scales into one term
            if (op_i == op_k && op_j == op_l) {

                // if i>j, we switch l,j to get the related term
                // if j<i, we have to switch i,k, otherwise we get a forbidden permutation
                IndexTuple self(i,j,k,l), twin(i,l,k,j);
                if (i<j) twin = IndexTuple(k,j,i,l);

                if (self > twin) {
                
                    typename hamtagterm_t<M>::type
                    term = TermMaker<M>::four_term(ident, fill, coefficients[align(i,j,k,l)], i,k,l,j,
                                                   op_i, op_k, op_l, op_j, op_table);

                    term.scale += value_type(sign(twin)) * coefficients[align(twin)];

                    tagterms.push_back(term);
                }
                //else: we already have the term
            }
            else {
                tagterms.push_back( TermMaker<M>::four_term(ident, fill, coefficients[align(i,j,k,l)], i,k,l,j,
                                   op_i, op_k, op_l, op_j, op_table) );
            }
        }
    
    private:
        void parse_integrals(BaseParameters & parms, Lattice const & lat) {

            // ********************************************************************
            // *** Parse orbital data *********************************************
            // ********************************************************************

            std::ifstream orb_file;
            orb_file.open(parms.get<std::string>("integral_file").c_str());
            for (int i = 0; i < 4; ++i){
                orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
            }

            std::vector<double> raw;
            std::copy(std::istream_iterator<double>(orb_file), std::istream_iterator<double>(),
                        std::back_inserter(raw));

            std::vector<double>::iterator it = raw.begin();
            while (it != raw.end()) {
                matrix_elements.push_back(*it++);
                std::vector<int> tmp;
                std::transform(it, it+4, std::back_inserter(tmp), boost::lambda::_1-1);
                idx_.push_back(tmp);
                it += 4;
            }

            #ifndef NDEBUG
            for (std::size_t m = 0; m < matrix_elements.size(); ++m)
            {
                assert( *max_element(idx[m].begin(), idx[m].end()) <= lat.size() );
            }
            #endif
            
        }

        tag_t ident, fill;
        boost::shared_ptr<OPTagTable<M, TwoU1> > op_table;

        std::vector<value_type> matrix_elements;
        std::vector<std::vector<int> > idx_;

        std::map<IndexTuple, value_type> coefficients;

        std::map<TermTuple, typename hamtagterm_t<M>::type> three_terms;

    };
}

#endif
