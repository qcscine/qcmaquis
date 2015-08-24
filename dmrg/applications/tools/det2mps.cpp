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
//#include <boost/bimap.hpp>
//#include <boost/bimap/multiset_of.hpp>
#include <boost/lexical_cast.hpp>
#include "dmrg/sim/matrix_types.h"
//#include "dmrg/mp_tensors/compression.h"
#include "dmrg/models/model.h"
//#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/models/lattice.h"
#include "alps/numeric/matrix.hpp"
#include "dmrg/models/chem/util.h"

namespace deas_detail
{
    // free functions

    /*
        template <class SymmGroup, class = void>
        class PGCharge
        {
        public:
            typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, int irr)
            {
                return rhs;
            }
        };

        template <class SymmGroup>
        class  PGCharge<SymmGroup, typename boost::enable_if<symm_traits::HasPG<SymmGroup> >::type>
        {
        public:
            typedef typename SymmGroup::subcharge subcharge;
            typename SymmGroup::charge operator()(typename SymmGroup::charge rhs, subcharge irr)
            {
                SymmGroup::irrep(rhs) = irr;
                return rhs;
            }
        };
    */
}

template<class Matrix, class SymmGroup, class=void>
struct deas_mps_init : public mps_initializer<Matrix,SymmGroup>
{
    deas_mps_init(BaseParameters parms_,
                std::vector<Index<SymmGroup> > const& phys_dims_,
                typename SymmGroup::charge right_end_,
                std::vector<int> const& site_type,
                std::vector<std::vector<std::size_t> > const& det_list_ )
    : parms(parms_)
    , phys_dims(phys_dims_)
    , site_types(site_type)
    , di(parms, phys_dims_, right_end_, site_type)
    , det_list(det_list_)
    , right_end(right_end_)
    {}

    typedef Lattice::pos_t pos_t;
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {

     //get hf determinant, should be extended to a vector of determinants
        std::vector<std::vector<std::size_t> > dets=det_list;

     //order will be treated while creating the determinants!!!
      //check for correct size
        if (dets[0].size() != mps.length())
            throw std::runtime_error("HF occupation vector length != MPS length\n");

      //actual initialization; check for sector and fill it with ones
      //max_charge is accumulated charge of the sites
        typename SymmGroup::charge max_charge, search_charge, first_charge = SymmGroup::IdentityCharge;

        // NO! -- access mps[i].site_dim() gives the Index with correct charges
        //initialize updown and empty inidices that will be needed later fo determination of the target size
        typename SymmGroup::charge ud(1), em(0);
        typename SymmGroup::charge up(0), down(0);
        up[0] = down[1] = 1;

        //go through determinants and get possible charges
        std::map<typename SymmGroup::charge, int> charge_to_int;
        typename std::map<typename SymmGroup::charge, int>::iterator it;
        int sc = 1;
        int count = 0;
        for(int d = 0; d < dets.size(); d++)
        {
            accumulated_charge = SymmGroup::IdentityCharge;           
            for(int i = 0; i < dets[0].size(); i++)
            {
                //search if accumulated_charge is already in map
                it = charge_to_int.find(accumulated_charge);
                if(it == charge_to_int.end())
                {
                  charge_to_int[accumulated_charge] = count;
                  count++;
                }
                charge site_charge = charge_from_int(dets[d][i]);
                
                if (SymmGroup::particleNumber(site_charge) %2 == 1)
                    PGCharge<SymmGroup>()(site_charge, site_types[i]);

                accumulated_charge = SymmGroup::fuse(accumulated_charge, site_charge);
            }
        }

        //New Approach including subsectors

        // initialize objects required 
        std::vector<std::vector<int > > rows_to_fill(dets.size(), std::vector <int> (dets[0].size()));
        std::vector<std::vector<std::map<std::string, int> > > str_to_col_map(dets[0].size(),
                                                                              std::vector<std::map<std::string, int> > (charge_to_int.size())
                                                                              );
        std::string str;
        int ifc, max_value, prev_row, nrows, nrows_fill, off = 0;

        // needed?
        int Mmax = 0;

         //main loop
         for(int d= 0; d<dets.size(); d++)
         {
             max_charge = right_end;
             for(int s = dets[0].size()-1; s > 0; s--){
                 charge site_charge = charge_from_int(dets[d][s]);

                if (SymmGroup::particleNumber(site_charge) %2 == 1)
                    PGCharge<SymmGroup>()(site_charge, site_types[i]);

                 max_charge = SymmGroup::fuse(max_charge, -site_charge);
                 str = det_string(s, dets[d]);
                 ifc = charge_to_int[max_charge];

                 if (str_to_col_map[s-1][ifc][str])
                    rows_to_fill[d][s] = str_to_col_map[s-1][ifc][str]-1;

                 else
                 {
                    //get largest element in map
                    max_value = str_to_col_map[s-1][ifc].size();
                    str_to_col_map[s-1][ifc][str] = max_value;
                    rows_to_fill[d][s] = max_value - 1;
                    //get size of largest sector
                    Mmax = std::max(Mmax, max_value);
                 }
             }
         } 

          //this now calls a function which is part of this structure
          init_sect(mps, str_to_col_map, charge_to_int, true, 0); 

          //this here is absolutely necessary
          for(pos_t i = 0; i < mps.length(); ++i)
             mps[i].multiply_by_scalar(0.0);


       //fill loop
       for(int d= 0; d<dets.size();d++){
          max_charge = right_end;
          prev_row = 0;
          for(int s = dets[0].size()-1;s > 0; s--){
              charge site_charge = charge_from_int(dets[d][s]);
              if(dets[d][s]==2||dets[d][s]==3){site_charge[2]=site_types[s];}
              search_charge = SymmGroup::fuse(max_charge,-site_charge);
              nrows_fill = mps[s].row_dim().size_of_block(search_charge);
             //get current matrix
             size_t max_pos = mps[s].data().left_basis().position(max_charge);
             Matrix & m_insert = mps[s].data()[max_pos];
             nrows = m_insert.num_rows();
             //get additional offsets for subsectors
             off = 0;
             if(dets[d][s] == 3){
                if(mps[s].row_dim().has(SymmGroup::fuse(max_charge,-ud))){
                  off =  mps[s].row_dim().size_of_block(SymmGroup::fuse(max_charge,-ud));
                }
             }
             else if(dets[d][s] == 2){
               if(mps[s].row_dim().has(SymmGroup::fuse(max_charge,-em))){
                 off = nrows - nrows_fill - mps[s].row_dim().size_of_block(SymmGroup::fuse(max_charge,-em));
               }else{off = nrows-nrows_fill;}
             }
             else if(dets[d][s] == 1){off = nrows- nrows_fill;}
          //actual insertion
             m_insert(off+rows_to_fill[d][s],prev_row) = 1;
             prev_row = rows_to_fill[d][s];
             max_charge = SymmGroup::fuse(max_charge,-site_charge);
          }
       }

         //first site needs to be filled as well
       int fill = 0;
       for(int d = 0; d<dets.size(); d++){
          fill = dets[d][0];
          first_charge = charge_from_int(fill);
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


    //function to get charge from int
    typename SymmGroup::charge charge_from_int (int sc_input){
      typename SymmGroup::charge site_charge(0);
       switch(sc_input) {
           case 4:
           site_charge = phys_dims[site_types[0]][0].first; // updown                   
           break;
           case 3:
           site_charge = phys_dims[site_types[0]][1].first; // up
           break;
           case 2:
           site_charge = phys_dims[site_types[0]][2].first; // down
           break;
           case 1:
           site_charge = phys_dims[site_types[0]][3].first; // empty
           break;
       }   
       return site_charge;
    }

    //function to initalize sectors -> copied from mps-initializers

    void init_sect(MPS<Matrix, SymmGroup> & mps,
                   const std::vector<std::vector<std::map<std::string, int > > > & str_to_col_map, 
                   std::map <typename SymmGroup::charge, int> &charge_to_int,
                   bool fillrand = true, 
                   typename Matrix::value_type val = 0)
    {
        parallel::scheduler_balanced scheduler(mps.length());
        std::size_t L = mps.length();

        std::vector<Index<SymmGroup> > allowed = allowed_sectors(site_types, phys_dims, right_end, 5);
        maquis::cout << "allowed sectors created" <<std::endl;
        allowed = adapt_allowed(allowed, str_to_col_map, charge_to_int);

        maquis::cout << "size of sectors succesfully adapted" << std::endl;

        omp_for(size_t i, parallel::range<size_t>(0,L), {
            parallel::guard proc(scheduler(i));
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_types[i]], allowed[i], allowed[i+1], fillrand, val);
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        });
    }

    std::vector<Index<SymmGroup> > adapt_allowed(std::vector<Index<SymmGroup> > allowed,
                                                     const std::vector<std::vector<std::map<std::string, int > > > &str_to_col_map,
                                                     std::map <typename SymmGroup::charge, int> &charge_to_int)
    {
       std::vector<Index<SymmGroup> > adapted = allowed;
       int ifc = 0, Mmax = 0;
       int L = str_to_col_map.size();
       for(int i = 1; i<L+1; ++i){
          for(typename Index<SymmGroup>::iterator it = adapted[i].begin();
                     it != adapted[i].end(); ++it){
             ifc = charge_to_int[it->first];
             Mmax = str_to_col_map[i-1][ifc].size();
             if(Mmax>0){
               it->second = str_to_col_map[i-1][ifc].size();
             }
          }
       }
       return adapted;
    }


    BaseParameters parms;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<int> site_types;
    default_mps_init<Matrix, SymmGroup> di;
    std::vector<std::vector<std::size_t> > det_list;
    typename SymmGroup::charge right_end;
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

        typedef alps::numeric::matrix<double> Matrix;


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
TwoU1PG::charge a(0), b(0), c(0), d(1);
b[0]=c[1] = 1;

   std::vector<Index<TwoU1PG> > phys_dims(4);
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
 TwoU1PG::charge right_end = chem_detail::qn_helper<TwoU1PG>().total_qn(parms);
 maquis::cout << "Right end: " << right_end <<std::endl;

//det_list as function parameter
  std::vector<std::vector<std::size_t> >  det_list = dets_from_file(det_file);

   std::cout <<"hf_determinant = 1st determinant in list: ";
   for(int i = 0; i<det_list[0].size(); i++){
      std::cout << det_list[0][i] <<", "; 
   }
   std::cout <<std::endl;
   


//create MPS
   MPS<Matrix,TwoU1PG> hf_mps(L);;
   deas_mps_init<Matrix,TwoU1PG> hf(parms,phys_dims,right_end,site_types,det_list);
   hf(hf_mps); 
   maquis::cout << "MPS created" << std::endl;
   save("test",hf_mps);
   maquis::cout << "MPS saved" << std::endl; 


    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
