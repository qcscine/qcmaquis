/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
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

#include <iostream>
#include <algorithm>
#include <boost/lexical_cast.hpp>

#include "dmrg/sim/matrix_types.h"
#include "dmrg/models/model.h"
//#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/lattice.h"
#include "alps/numeric/matrix.hpp"
#include "dmrg/models/chem/util.h"


#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> matrix;
#endif


#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif

namespace deas_detail
{
    template <class SymmGroup, class = void>
    class charge_from_int
    {
        typedef Lattice::pos_t pos_t;
        typedef typename SymmGroup::charge charge;
        typedef std::vector<Index<SymmGroup> > index_vec;
        typedef std::vector<typename SymmGroup::subcharge> site_vec;
    public:
        std::vector<charge> operator()(int sc_input, pos_t p, index_vec const & phys_dims, site_vec const & site_types)
        {
            std::vector<charge> site_charges;
            switch(sc_input) 
            {
                case 4:
                    site_charges.push_back(phys_dims[site_types[p]][0].first); // updown                   
                    break;
                case 3:
                    site_charges.push_back(phys_dims[site_types[p]][1].first); // up
                    break;
                case 2:
                    site_charges.push_back(phys_dims[site_types[p]][2].first); // down
                    break;
                case 1:
                    site_charges.push_back(phys_dims[site_types[p]][3].first); // empty
                    break;
            }   
            return site_charges;
         }
    };
   
    template <class SymmGroup>
    class charge_from_int<SymmGroup, typename boost::enable_if< symm_traits::HasSU2<SymmGroup> >::type>
    {
        typedef Lattice::pos_t pos_t;
        typedef typename SymmGroup::charge charge;
        typedef std::vector<Index<SymmGroup> > index_vec;
        typedef std::vector<typename SymmGroup::subcharge> site_vec;
    public:
        std::vector<charge> operator()(int sc_input, pos_t p, index_vec const & phys_dims, site_vec const & site_types)
        {
            std::vector<charge> site_charges;
            switch(sc_input) {
                case 4:
                    site_charges.push_back(phys_dims[site_types[p]][0].first); // doubly-occ
                    break;
                case 3:
                    site_charges.push_back(phys_dims[site_types[p]][1].first); // singly-occ
                    site_charges.push_back(phys_dims[site_types[p]][2].first); // singly-occ
                    break;
                case 2:
                    site_charges.push_back(phys_dims[site_types[p]][1].first); // singly-occ
                    site_charges.push_back(phys_dims[site_types[p]][2].first); // singly-occ
                    break;
                case 1:
                    site_charges.push_back(phys_dims[site_types[p]][3].first); // empty
                    break;
            }
            return site_charges;
        }
    };
}

template<class Matrix, class SymmGroup, class=void>
struct deas_mps_init : public mps_initializer<Matrix,SymmGroup>
{
    deas_mps_init(BaseParameters parms_,
                std::vector<Index<SymmGroup> > const& phys_dims_,
                typename SymmGroup::charge right_end_,
                std::vector<int> const& site_type,
                std::vector<std::vector<std::size_t> > const& det_list_)
    : parms(parms_)
    , phys_dims(phys_dims_)
    , site_types(site_type)
    , di(parms, phys_dims_, right_end_, site_type)
    , right_end(right_end_)
    , det_list(det_list_)
    , determinants(det_list.size())
    {
        // convert det_list to vec<vec<charge>>
        for (size_t i = 0; i < det_list.size(); ++i)
            for (size_t j = 0; j < det_list[i].size(); ++j)
                determinants[i].push_back(deas_detail::charge_from_int<SymmGroup>()(det_list[i][j], j, phys_dims, site_types));
    }

    typedef Lattice::pos_t pos_t;
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        pos_t L = mps.length();
        if (determinants[0].size() != L)
            throw std::runtime_error("HF occupation vector length != MPS length\n");

        charge doubly_occ = phys_dims[0].begin()->first, empty = phys_dims[0].rbegin()->first;

        // initialize objects required 
        std::vector<std::vector<int > > rows_to_fill(determinants.size(), std::vector <int> (L));
        std::vector<std::map<charge, std::map<std::string, int> > > str_to_col_map(L);

        //main loop
        for(int d = 0; d < determinants.size(); ++d)
        {
            std::vector<charge> bond_charge, accumulated_charge, site_charge;
            bond_charge.push_back(right_end);
            for(int s = L - 1; s > 0; --s)
            {
                site_charge = determinants[d][s];
                for (typename std::vector<charge>::const_iterator it = bond_charge.begin(); it != bond_charge.end(); ++it)
                {
                    for (typename std::vector<charge>::const_iterator it2 = site_charge.begin(); it2 != site_charge.end(); ++it2)
                    {
                        charge current_charge = SymmGroup::fuse(*it, -*it2);
                        if (charge_detail::physical<SymmGroup>(current_charge))
                        {
                            accumulated_charge.push_back(current_charge);

                            std::string str = det_string(s, det_list[d]);
  
                            std::map<std::string, int> & str_map = str_to_col_map[s-1][current_charge];

                            if (str_map[str])
                                rows_to_fill[d][s] = str_map[str] - 1;

                            else
                            {
                                //get largest element in map
                                int max_value = str_map.size();
                                str_map[str] = max_value;
                                rows_to_fill[d][s] = max_value - 1;
                            }
                        }
                    }
                }
                bond_charge = accumulated_charge;
                accumulated_charge.clear();
            }
            bond_charge.clear();
        }

        //this now calls a function which is part of this structure
        init_sect(mps, str_to_col_map, true, 0); 

        //this here is absolutely necessary
        for(pos_t i = 0; i < L; ++i)
            mps[i].multiply_by_scalar(0.0);

        //fill loop
        for(int d = 0; d < determinants.size(); ++d)
            {
            std::vector<charge> bond_charge, accumulated_charge, site_charge;
            bond_charge.push_back(right_end);
            int prev_row = 0;
            for(int s = L - 1; s > 0; --s)
            {
                site_charge = determinants[d][s];
                for (typename std::vector<charge>::const_iterator it = bond_charge.begin(); it != bond_charge.end(); ++it)
                {
                    for (typename std::vector<charge>::const_iterator it2 = site_charge.begin(); it2 != site_charge.end(); ++it2)
                    {
                        charge current_charge = SymmGroup::fuse(*it, -*it2);
                        if (charge_detail::physical<SymmGroup>(current_charge) && mps[s].row_dim().has(current_charge))
                        {
                           int nrows_fill = mps[s].row_dim().size_of_block(current_charge);
                            //get current matrix
                            size_t max_pos = mps[s].data().left_basis().position(*it);
                            Matrix & m_insert = mps[s].data()[max_pos];
                     
                            int nrows = m_insert.num_rows(), off = 0;
                     
                            //get additional offsets for subsectors
                            if(*it2 == phys_dims[site_types[s]][1].first){
                     
                                if(mps[s].row_dim().has(SymmGroup::fuse(*it, -doubly_occ)))
                                    off =  mps[s].row_dim().size_of_block(SymmGroup::fuse(*it, -doubly_occ));
                          
                            }
                            else if(*it2 == phys_dims[site_types[s]][2].first)
                            {
                                if(mps[s].row_dim().has(*it)){
                                    off = nrows - nrows_fill - mps[s].row_dim().size_of_block(*it);
                                }
                                else {
                                    off = nrows - nrows_fill;
                                } 
                            }  
                            else if(*it2 == phys_dims[site_types[s]][3].first){
                                off = nrows - nrows_fill;
                            }
                            //actual insertion
                            m_insert(off+rows_to_fill[d][s],prev_row) = 1;
                            prev_row = rows_to_fill[d][s];//hier wird es kompliziert!!! Er geht ja mehrere Kombinationen durch
                            accumulated_charge.push_back(current_charge); 
                        }
                    }
                }
                bond_charge = accumulated_charge;
                accumulated_charge.clear();
            }
            bond_charge.clear();
       }

       //first site needs to be filled as well
       for(int d = 0; d < determinants.size(); d++){
          charge first_charge = determinants[d][0][0];
          size_t first_pos = mps[0].data().left_basis().position(first_charge);
          Matrix & m_first = mps[0].data()[first_pos];
          m_first(0,0) = 1;
       }

    }//end of main initialization function 


    //function to get string of left or right part from det
    std::string det_string(int s, std::vector<size_t> det){
       std::string str;
       char c;
       int L = det.size();
       if(s > L/2){
          for(int i = s; i<det.size(); i++){
             c = boost::lexical_cast<char>(det[i]);
             str.push_back(c);
          }
       }else{
          for(int i = 0; i<s; i++){
             c = boost::lexical_cast<char>(det[i]);
             str.push_back(c);
          }
       }
    return str;
    }

    //function to initalize sectors -> copied from mps-initializers

    void init_sect(MPS<Matrix, SymmGroup> & mps,
                   const std::vector<std::map<charge, std::map<std::string, int> > > & str_to_col_map, 
                   bool fillrand = true, 
                   typename Matrix::value_type val = 0)
    {
        parallel::scheduler_balanced scheduler(mps.length());
        std::size_t L = mps.length();

        std::vector<Index<SymmGroup> > allowed = allowed_sectors(site_types, phys_dims, right_end, 5);
        maquis::cout << "allowed sectors created" << std::endl;
        allowed = adapt_allowed(allowed, str_to_col_map);
        maquis::cout << "size of sectors succesfully adapted" << std::endl;

        omp_for(size_t i, parallel::range<size_t>(0,L), {
            parallel::guard proc(scheduler(i));
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_types[i]], allowed[i], allowed[i+1], fillrand, val);
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        });
    }

    std::vector<Index<SymmGroup> > adapt_allowed(std::vector<Index<SymmGroup> > allowed,
                                                     const std::vector<std::map<charge, std::map<std::string, int > > > &str_to_col_map)
    {
        std::vector<Index<SymmGroup> > adapted = allowed;
        pos_t L = str_to_col_map.size();
        for(pos_t i = 1; i < L+1; ++i)
        {
            for(typename Index<SymmGroup>::iterator it = adapted[i].begin();
                     it != adapted[i].end(); ++it)
            {
                if (str_to_col_map[i-1].count(it->first) != 0)
                    it->second = str_to_col_map[i-1].at(it->first).size();
            }
        }
        return adapted;
    }


    BaseParameters parms;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<typename SymmGroup::subcharge> site_types;
    default_mps_init<Matrix, SymmGroup> di;
    std::vector<std::vector<size_t> > det_list;
    std::vector<std::vector<std::vector<charge> > > determinants;
    charge right_end;
};

//parse determinants
std::vector<std::vector<std::size_t> > dets_from_file(std::string file){
    std::ifstream config_file;
    config_file.open(file.c_str());

    std::vector<std::vector<std::size_t> > configs;

    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_coeff;
        boost::split(det_coeff, line, boost::is_any_of(" "));

        std::string det = det_coeff[0];

        std::vector<size_t> tmp;
        for (std::size_t i = 0; i < det.size(); ++i) {
            int occ = boost::lexical_cast<size_t>(det[i]);
            switch(occ) {
                case 4:
                    tmp.push_back(4); // doubly occ
                    break;
                case 3:
                    tmp.push_back(3); // up
                    break;
                case 2:
                    tmp.push_back(2); // down 
                    break;
                case 1:
                    tmp.push_back(1); // empty
                    break;
            }
        }
        configs.push_back(tmp);
    }
    return configs;
}


int main(int argc, char ** argv){
    try {
        if (argc != 3) {
            std::cout << "Usage: " << argv[0] << " <result.h5>" << "determinants.txt" << std::endl;
            return 1;
        }

        std::string rfile(argv[1]);
        std::string det_file(argv[2]);

        storage::archive ar(rfile, "r");
        BaseParameters parms;
        ar["/parameters"] >> parms;
        Lattice lat(parms);

        /***Create TEST environment not needed in actual implementation***/
        size_t L = parms["L"];

        //create symmetry vector -> site_types
        int max_site_type = 0;
        std::vector<int> site_types(lat.size());
        for (int i = 0; i<lat.size(); i++){
           site_types[i] = lat.get_prop<int>("type", i);
           max_site_type =std::max(site_types[i],max_site_type);
        }
        std::cout <<"site_types: ";
        for(int i = 0; i<site_types.size(); i++){
           std::cout << site_types[i] <<", "; 
        }
        std::cout <<std::endl;
        maquis::cout << "maximal site type: " << max_site_type << std::endl;

    //create vector of indices -> phys_dims
    grp::charge a(2), b(1), c(1), d(0);
    a[1]= 0;
    c[1] = -1;

   std::vector<Index<grp> > phys_dims(4);
   for(int i = 0; i< max_site_type+1; i++){
      b[2]=c[2]=i;
      phys_dims[i].insert(std::make_pair(a,1));
      phys_dims[i].insert(std::make_pair(b,1));
      phys_dims[i].insert(std::make_pair(c,1));
      phys_dims[i].insert(std::make_pair(d,1));
   }
   for(int i = 0; i<max_site_type+1; i++){
      std::cout << "phys_dims["<<i<<"] = " <<phys_dims[i] <<std::endl; 
   }




   //get right end charge
   grp::charge right_end = chem_detail::qn_helper<grp>().total_qn(parms);
   maquis::cout << "Right end: " << right_end <<std::endl;

    //det_list as function parameter
    std::vector<std::vector<std::size_t> >  det_list = dets_from_file(det_file);

   std::cout <<"hf_determinant = 1st determinant in list: ";
   for(int i = 0; i<det_list[0].size(); i++){
      std::cout << det_list[0][i] <<", "; 
   }
   std::cout <<std::endl;
   


   //create MPS
   MPS<matrix,grp> hf_mps(L);
   deas_mps_init<matrix,grp> hf(parms,phys_dims,right_end,site_types,det_list);
   hf(hf_mps); 
   maquis::cout << "MPS created" << std::endl;
   save("test",hf_mps);
   maquis::cout << "MPS saved" << std::endl; 


    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
