
#define BOOST_TEST_MODULE ambient_c_kernels
//#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include"utils/gnuplot_i.hpp"


#define SLEEP_LGTH 2  // sleep time in seconds
#define NPOINTS    50 // length of array

void wait_for_key(); // Programm halts until keypress

using std::cout;
using std::endl;

typedef boost::mpl::list< bool> plot;

template<typename T>
void Load(std::string name, std::vector<T> & size, std::vector<T> &gflop){
  std::string line;
  std::ifstream myfile(name.c_str());
  T a,b,c,d,e;
  if (myfile.is_open()){
    while (!myfile.eof()) {
      myfile >> a >> b >> c >> d >> e;
      gflop.push_back(b);
      size.push_back(c);
    }   
    myfile.close();
  } else cout << "Unable to open file";
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test, T, plot  )
{
    std::vector<float> SizeBlas,GflopBlas;
    std::vector<float> SizeAmbient,GflopAmbient;
    Load<float>("TimeGemmBlas.txt", SizeBlas, GflopBlas);
    Load<float>("TimeGemmAmbient.txt", SizeAmbient, GflopAmbient);
    try
    {
        Gnuplot g1("lines");
        g1.savetops("test_output");
        g1.set_style("lines").plot_xy(SizeBlas,GflopBlas,"user-defined points 2d");
        g1.set_style("lines").plot_xy(SizeAmbient,GflopAmbient,"user-defined points 2d");
        g1.showonscreen(); // window output
    } catch (GnuplotException ge) {
        cout << ge.what() << endl;
    }
    BOOST_CHECK(true==true );
}
