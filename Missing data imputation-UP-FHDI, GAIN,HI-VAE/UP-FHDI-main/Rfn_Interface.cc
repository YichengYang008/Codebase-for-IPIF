//local functions 
//#include "Rfn_test.h" //this is called by Rfn_test.cc before 2017_0112

#include "Rfn_test.cc" 

//R built-in functions
//#include "R.h"     // R functions 
//#include "Rmath.h" // Rmath

//----------------------------------------------------------------------------
//- This file defines TWO functions. 
// 
// The first one is a C++ function that directly accesses 
// the C++ class. 
// 
// The second one is a C function named ’CWrapper’, that INdirectly accesses 
// the C++ class, via a call to the function.
//----------------------------------------------------------------------------
void Rfn_test_call(double* x, int* r, int * nrow_x, int * ncol_x, 
                   double* k, double* d, int * M)

{ 

	//Rfn_test Rfn_test(x, r, nrow_x, ncol_x, k, d, M); //before 2017_0112
	Rfn_test(x, r, nrow_x, ncol_x, k, d, M); 
	
	return ;

}

//----------------------------------------------------------------------------
//- CWrapper 
// 
// C function that in turn invokes the above C++ function . 
// R can access C code but can’t access C++ code directly. 
// This C function provides a C interface to the C++ code that R can access. 
// See: http://www.parashift.com/c++-faq-lite/mixing-c-and-cpp.html 
// In this C function, you must NOT include class declarations, 
// instantiate any C++ objects, or do any oher obvious C++ 
// things. 
// The EXTERN statement tells the C++ compiler that the 
// enclosed function ’CWrapper’ is a C function. 
// Although apparently we can insert C++ style comments and 
// we can even declare variables in the middle of the function, 
// which I thought you can’t do in regular C. 
// 
extern "C" 
{ 
	void CWrapper(double* x, int* r, int * nrow_x, int * ncol_x, double* k, 
	              double* d, int* M) 
	{ 
		//- Invoke second function which internally can do C++ things. 
		// 
		
		Rfn_test_call(x, r, nrow_x, ncol_x, k, d, M); 

	} 
}

