#include <iostream>
#include <fstream>
#include <string>

#include <alps/model/hamiltonian_matrix.hpp>

// for the moment it cannot be used
//#include <alps/numeric/matrix.hpp>
//typedef alps::numeric::matrix<double> Matrix;

// we have to switch back to ublas
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
typedef boost::numeric::ublas::matrix<double> Matrix;

int main(int argc, char ** argv)
{
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <model>" << std::endl;
        return 0;
    }
    try {
        std::ifstream pfile(argv[1]);
        alps::Parameters parms(pfile);
        parms["N_total"] = 1; // force the system to one particle!
        
        Matrix H = alps::hamiltonian_matrix<Matrix>(parms).matrix();
        std::cout << H << std::endl;
    }
    catch (std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return -1;
    }
    catch (...) {
        std::cerr << "Fatal Error: Unknown Exception!\n";
        return -2;
    }
    return 0;
}
