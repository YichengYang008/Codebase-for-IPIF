
#include <string>
#include <vector>


class Rfn_test { 
public: 
//- Constructors and destructor. 
// 
	Rfn_test(double* x, int* r, int* nrow_x, int* ncol_x, double* k, 
	         double* d, int* M);

	


//- A function for deallocating memory is made publicly available because 
// I want to not only run test this code in a pure C++ test 
// program outside of R, testRMat.cc, but I also want to run this code 
// from within R. In the latter case, I will let R handle memory 
// management. This means that I can’t automatically deallocate memory 
// in the destructor, which is the usual thing one might do. 
// For symmetry, the function for memory allocation is also made publicly 
// available, although it didn’t have to be. (Maybe it should be private!) 
// 


private: 
};

