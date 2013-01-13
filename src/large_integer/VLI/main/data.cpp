#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include "boost/lexical_cast.hpp"

typedef boost::tuple<int, int, int, double, double, double> mytuple;

bool mycompare (const mytuple &lhs, const mytuple &rhs){
  return (boost::get<0>(lhs) < boost::get<0>(rhs) );
}
// I could boost mpl later ... 
template<int numvar>
class data{
   public:
   data(std::string const& name):name_(name),maxorder_(14){
       datanumbits_.resize(3); // 3 because 128, 192 et 256 input vli
   }; 

   void operator()(){
       load();
       sort();
       save();
   }

   void load(){
       for( int i=1; i<=maxorder_; ++i){
           std::string file(name_);
           file +=  boost::lexical_cast<std::string>(i);
           std::ifstream is(file.c_str(),std::ios::in);
           int a,b,c;
           double d,e,f;
           int pos;
           while (is.good()){ // this duplicate the last element of the opening file
               is >> a >> b >> c >> d >> e >> f; 
               pos = b/64-2;
               datanumbits_[pos].push_back(boost::make_tuple<int, int, int, double, double, double>(a,b,c,d,e,f));
           }
           std::vector< boost::tuple<int, int, int, double, double, double> >::iterator it = datanumbits_[pos].end(); //delete the duplicate element 
           datanumbits_[pos].erase(it);
           is.close();
       }
   } 

   void sort(){
       std::vector< std::vector< boost::tuple<int, int, int, double, double, double> > >::iterator it = datanumbits_.begin(); //delete the duplicate element 
       for(it; it != datanumbits_.end(); ++it)
           std::sort((*it).begin(), (*it).end(), mycompare);
   }
 
   void save(){
       std::vector< std::vector< boost::tuple<int, int, int, double, double, double> > >::iterator it1 = datanumbits_.begin(); 
       for(; it1 != datanumbits_.end(); ++it1){
           std::string name(name_);
           std::vector< boost::tuple<int, int, int, double, double, double> >::const_iterator it2 = (*it1).begin();
           name += boost::lexical_cast<std::string>(boost::get<1>((*it2))) + "bits.dat";
           std::ofstream os(name.c_str(),std::ios::app);
           os.setf(std::ios::fixed,std::ios::floatfield);
           for(; it2 < (*it1).end(); it2+=numvar){ // MaxOrderEach 3 var, MaxOrderCombined 4 var
               os << boost::get<0>((*it2))<< " " << boost::get<1>((*it2)) ; /* give info on the poly + vli */
               int i =0;
               int nv = 1;
               while(i<numvar){
                    if(nv == boost::get<2>(*(it2+i))){
                        os << " " << boost::get<2>(*(it2+i)) << " " << boost::get<3>(*(it2+i)) << " " << boost::get<4>(*(it2+i)) << " " << boost::get<5>(*(it2+i)); /* give #var, tgmp, tvlicpu, tvligpu */
                        nv++;i=0;
                    }else
                        i++;
               };  
               os << std::endl;
               
           }
      }
   }
   private:
   std::vector< std::vector< boost::tuple<int, int, int, double, double, double> > > datanumbits_;
   
   std::string name_;
   int maxorder_;
};

int main(){
    data<4>  MaxOrderEach("MaxOrderEachTime");
    data<4>  MaxOrderCombined("MaxOrderCombinedTime"); 
    MaxOrderEach();
    MaxOrderCombined();
}
