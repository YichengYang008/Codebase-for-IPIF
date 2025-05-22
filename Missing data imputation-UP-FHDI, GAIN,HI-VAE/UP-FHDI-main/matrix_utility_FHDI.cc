#ifndef MATRIX_UTILITY_FHDI
#define MATRIX_UTILITY_FHDI
//---------------------
//Collection of basic matrix and vector utilities
//for FHDI program 
//October 5, 2016
//
//Developed by I. Cho
//All rights reserved
//----------------------
//#include <string>
//#define R_NO_REMAP //so that R does not define length(x) which may cause many complie error with fstream

#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>

using std::cout; 
using std::cerr; 
using std::endl; 


void Copy_dVector(double Source[], int n, double Target[])
{
   for(int i=0; i<n; i++)
   {
      Target[i] = Source[i];
   }
   return;
}
//==============================================================================
//==============================================================================
void Copy_iVector(int Source[], int n, int Target[])
{
   for(int i=0; i<n; i++)
   {
      Target[i] = Source[i];
   }
   return;
}
//==============================================================================
//==============================================================================



void Copy_dMatrix(double** Source, int n_row, int n_col, double** Target)
{
   for(int i_row=0; i_row<n_row; i_row++)
   {
      for(int i_col=0; i_col<n_col; i_col++)
      {
         Target[i_row][i_col] = Source[i_row][i_col];
      }

   }
   return;
}
//==============================================================================
//==============================================================================



void Copy_iMatrix(int** Source, int n_row, int n_col, int** Target)
{
   for(int i_row=0; i_row<n_row; i_row++)
   {
      for(int i_col=0; i_col<n_col; i_col++)
      {
         Target[i_row][i_col] = Source[i_row][i_col];
      }

   }
   return;
}
//==============================================================================
//==============================================================================



double ** New_dMatrix(int n_row, int n_col)
//Description================================
// make new double MATRIX
//============================================
{

   //===============================================
   // NON contiguous dynamic multidimensional array.
   // for MPI version ineffective from June 15 09
   //===============================================
   /*double ** matrix;

   matrix = new double*[n_row];
   for(int i=0; i<n_row;i++)
   {
      matrix[i] = new double[n_col];

      for(int j=0; j<n_col; j++)
      {
         matrix[i][j] = 0.0;
      }
   }*/


   
   //==========================================
   //Contiguous dynamic multidimensional array
   //which is essential for MPI usage
   //for MPI version effective from June 15 09
   //==========================================
   double ** matrix;

   matrix = new double*[n_row];
   matrix[0] = new double[n_row*n_col];//allocate the total storage as a contiguous block
   for(int i=1; i<n_row; i++)
   {
      matrix[i] = matrix[0] + i*n_col; //matrix[i] points to the entire block of ith row
   }

   for(int i=0; i<n_row;i++)
   {
      for(int j=0; j<n_col; j++)
      {
         matrix[i][j] = 0.0;
      }
   }
   


   return matrix;
}

char ** New_cMatrix(int n_row, int n_col)
//Description================================
// make new character MATRIX
//============================================
{

	//===============================================
	// NON contiguous dynamic multidimensional array.
	// for MPI version ineffective from June 15 09
	//===============================================
	/*double ** matrix;

	matrix = new double*[n_row];
	for(int i=0; i<n_row;i++)
	{
	matrix[i] = new double[n_col];

	for(int j=0; j<n_col; j++)
	{
	matrix[i][j] = 0.0;
	}
	}*/



	//==========================================
	//Contiguous dynamic multidimensional array
	//which is essential for MPI usage
	//for MPI version effective from June 15 09
	//==========================================
	char ** matrix;

	matrix = new char*[n_row];
	matrix[0] = new char[n_row*n_col];//allocate the total storage as a contiguous block
	for (int i = 1; i<n_row; i++)
	{
		matrix[i] = matrix[0] + i*n_col; //matrix[i] points to the entire block of ith row
	}

	for (int i = 0; i<n_row;i++)
	{
		for (int j = 0; j<n_col; j++)
		{
			matrix[i][j] = 0.0;
		}
	}



	return matrix;
}
//==============================================================================
//==============================================================================


void Del_dMatrix(double ** matrix, int n_row, int n_col)
//Description================================
// delete the double MATRIX
//============================================
{
   /*for(int i=0; i<n_row; i++)
   {
      delete[] matrix[i];
   }
   delete[] matrix;*/

   //==========================================
   //Contiguous dynamic multidimensional array
   //which is essential for MPI usage
   // effective from June 15 09
   //==========================================
   delete[] matrix[0];
   delete[] matrix;
}
//==============================================================================
//==============================================================================

void Del_cMatrix(char ** matrix, int n_row, int n_col)
//Description================================
// delete the double MATRIX
//============================================
{
	/*for(int i=0; i<n_row; i++)
	{
	delete[] matrix[i];
	}
	delete[] matrix;*/

	//==========================================
	//Contiguous dynamic multidimensional array
	//which is essential for MPI usage
	// effective from June 15 09
	//==========================================
	delete[] matrix[0];
	delete[] matrix;
}
//==============================================================================
//==============================================================================

int** New_iMatrix(int n_row, int n_col)
//Description================================
// make new integer MATRIX
//============================================
{

   //===============================================
   // NON contiguous dynamic multidimensional array.
   // for MPI version ineffective from June 15 09
   //===============================================
   /*int ** matrix;
   matrix = new int*[n_row];
   for(int i=0; i<n_row;i++)
   {
      matrix[i] = new int[n_col];

      for(int j=0; j<n_col; j++)
      {
         matrix[i][j] = 0;
      }
   }*/


   
   //==========================================
   //Contiguous dynamic multidimensional array
   //which is essential for MPI usage
   //for MPI version effective from June 15 09
   //==========================================
   int** matrix ;
   matrix = new int*[n_row];
   matrix[0] = new int[n_row*n_col]; //allocate the total storage of entire block

   for(int i=1; i<n_row; i++)
   {
      matrix[i] = matrix[0] + i*n_col; //allocate sequential addr starting from [0]
   }

   for(int i=0; i<n_row; i++)
   {
      for(int j=0; j<n_col; j++)
      {
         matrix[i][j] = 0 ;
      }
   }
   


   return matrix;
}
//==============================================================================
//==============================================================================



void Del_iMatrix(int ** matrix, int n_row, int n_col)
//Description================================
// delete the integer MATRIX
//============================================

{
   /*for(int i=0; i<n_row; i++)
   {
      delete[] matrix[i];
   }
   delete[] matrix;*/

   //==========================================
   //Contiguous dynamic multidimensional array
   //which is essential for MPI usage
   // effective from June 15 09
   //==========================================
	//matrix[0] = NULL;
	//matrix = NULL;

   delete[] matrix[0];
   delete[] matrix;
}
//==============================================================================
//==============================================================================




int Find_iValue(int** i_matrix, int n_row, int n_col,
                char s_rowcol, int n_rowcol, int i_value)
// Description======================
//IN   : int** i_matrix[n_row][n_col]
//OUT  : (i, n_rowcol) = position of the 'i_value'
//     or(n_rowcol, i)
//
//   by searching i_matrix[:][n_rowcol] when s_rowcol="r"
//             or i_matrix[n_rowcol][:] when s_rowcol="c"
//
//Note: -1 is returned when there is no matched value
//===================================

{
   if(s_rowcol == 'r')
   {
      for(int i=0; i<n_row; i++)
      {
         if(i_matrix[i][n_rowcol] == i_value)
         {
            return i;
         }
      }
   }
   else if(s_rowcol == 'c')
   {
      for(int i=0; i<n_col; i++)
      {
         if(i_matrix[n_rowcol][i] == i_value)
         {
            return i;
         }
      }

   }
   return -1; //when no matched value

}
//==============================================================================
//==============================================================================


int Find_dValue(double** d_matrix, int n_row, int n_col,
                char s_rowcol, int n_rowcol, double d_value)
// Description======================
//IN   : double** d_matrix[n_row][n_col]
//OUT  : (i, n_rowcol) = position of the 'd_value'
//     or(n_rowcol, i)
//   by searching d_matrix[:][n_rowcol] when s_rowcol="row"
//             or d_matrix[n_rowcol][:] when s_rowcol="col"
//===================================

{
   if(s_rowcol == 'r')
   {
      for(int i=0; i<n_row; i++)
      {
         double d_temp=0.0;
         d_temp = fabs(d_matrix[i][n_rowcol] - d_value);
         if(d_temp < 10.E-10)
         {
            return i;
         }
      }
   }
   else if(s_rowcol == 'c')
   {
      for(int i=0; i<n_col; i++)
      {
         double d_temp=0.0;
         d_temp = fabs(d_matrix[n_rowcol][i] - d_value);
         if(d_temp < 10.E-10)
         {
            return i;
         }
      }

   }
   return -1; //when no matched value

}
//==============================================================================
//==============================================================================


int iMaxValue(int** i_matrix, int n_row, int n_col, char s_where,
              int n_begin, int n_end, int n_at)
//Description=================
//
// find maximum int value in the int matrix
// (1) s_where ="row" ; search through row(within n_begin~n_end) at n_at col
// (2) s_where ="col" ; search through col(within n_begin~n_end) at n_at row
// (3) s_where ="all" ; search all matrix
//
//IN   : int** i_matrix[n_row][n_col]
//OUT  : int maximum value
//============================

{
   int i_maximum =0;
   int i_temp=0;

   if(s_where =='r')
   {
      for(int i=n_begin;i<=n_end; i++)
      {
         if(i_temp < i_matrix[i][n_at]) i_temp = i_matrix[i][n_at] ;
      }
   }
   else if(s_where =='c')
   {
      for(int i=n_begin;i<=n_end; i++)
      {
         if(i_temp < i_matrix[n_at][i]) i_temp = i_matrix[n_at][i] ;
      }
   }
   else if(s_where =='a')
   {
      for(int i=0;i<n_row; i++)
      {
         for(int j=0; j<n_col; j++)
         {
            if(i_temp < i_matrix[i][j]) i_temp = i_matrix[i][j] ;
         }
      }

   }

   i_maximum = i_temp;
   return i_maximum;
}
//==============================================================================
//==============================================================================


int iMinValue(int** i_matrix, int n_row, int n_col, char s_where,
              int n_begin, int n_end, int n_at)
//Description=================
//
// find minimum int value in the int matrix
// note: do search among only positive values
//
// (1) s_where ="row" ; search through row(within n_begin~n_end) at n_at col
// (2) s_where ="col" ; search through col(within n_begin~n_end) at n_at row
// (3) s_where ="all" ; search all matrix
//
//IN   : int** i_matrix[n_row][n_col]
//OUT  : int minimum value
//============================

{
   int i_minimum =0;
   int i_temp=0;//max_previous;

   if(s_where =='r')
   {
      //initiallize i_temp value
      for(int i=n_begin;i<=n_end; i++)
      {
         if(i_matrix[i][n_at]>0)
         {
            i_temp = i_matrix[i][n_at] ;
            break; //exit this loop
         }
      }

      for(int i=n_begin;i<=n_end; i++)
      {

         if(i_temp > i_matrix[i][n_at]&& i_matrix[i][n_at]>0)
            i_temp = i_matrix[i][n_at] ;
      }
   }
   else if(s_where =='c')
   {
      //initialize i_temp value
      for(int i=n_begin;i<=n_end; i++)
      {
         if(i_matrix[n_at][i]>0)
         {
            i_temp = i_matrix[n_at][i] ;
            break; //exit this loop
         }
      }

      for(int i=n_begin;i<=n_end; i++)
      {
         if(i_temp > i_matrix[n_at][i]&& i_matrix[n_at][i]>0)
            i_temp = i_matrix[n_at][i] ;
      }
   }
   else if(s_where =='a')
   {
      //initialize i_temp
      for(int i=0;i<n_row; i++)
      {
         for(int j=0; j<n_col; j++)
         {
            if(i_matrix[i][n_at]>0)
            {
              i_temp = i_matrix[i][j] ;
              break; //exit this loop
            }
         }
      }

      for(int i=0;i<n_row; i++)
      {
         for(int j=0; j<n_col; j++)
         {
            if(i_temp > i_matrix[i][j]&& i_matrix[i][n_at]>0)
               i_temp = i_matrix[i][j] ;
         }
      }
   }

   i_minimum = i_temp;
   return i_minimum;
}
//==============================================================================
//==============================================================================


void Fill_dVector(double *d_vector, const int n_size, const double value)
{
   for(int i=0; i<n_size; i++)
   {
      d_vector[i] = value;
   }

   return;
}
//==============================================================================
//==============================================================================
void Fill_iVector(int i_vector[], const int n_size, const int value)
{
   for(int i=0; i<n_size; i++)
   {
      i_vector[i] = value;
   }

   return;
}

void Fill_iVector(std::vector<int> i_vector, const int value)
{
	int n_size = i_vector.size();

	for (int i = 0; i<n_size; i++)
	{
		i_vector[i] = value;
	}

	return;
}
//==============================================================================
//==============================================================================



void Fill_dMatrix(double** d_Matrix, int n_row, int n_col, double value)
{
   for(int i=0; i<n_row; i++)
   {
      for(int j=0; j<n_col; j++) d_Matrix[i][j] = value;
   }
   return;
}
//==============================================================================
//==============================================================================



void Fill_iMatrix(int** i_Matrix, int n_row, int n_col, int value)
{
   for(int i=0; i<n_row; i++)
   {
      for(int j=0; j<n_col; j++) i_Matrix[i][j] = value;
   }
   return;
}
//==============================================================================
//==============================================================================



void Inverse_dMatrix(double** d_Mat, const int n, double** d_Inv)
//Description===================
//return inverse matrix of the n*n matrix
//using Gauss-Jordan elimination
//
// if diagonal term is zero and too small
// perform pivoting with largest value on the columns
//
//IN    : double** d_Mat[n][n]
//OUT   : double** d_Inv[n][n]
//==============================
{
   const double eps=1.e-15 ;


   //make d_Inv unity matrix
   for(int i=0; i<n; i++)
   {
      for(int j=0; j<n; j++) d_Inv[i][j] =0.0;

      d_Inv[i][i]=1.0 ;
   }



   double c=0.0 ; //coeff.

   for(int i_diag=0; i_diag<n; i_diag++)
   {
      c = d_Mat[i_diag][i_diag];

      //================
      //when diagonal term is too small or zero, then needs pivoting
      //================
      if(fabs(c)<eps )
      {
         //find max on current column ======
         double d_temp= c;  int i_loc=i_diag;

         for(int i=(n-1); i>i_diag; i--)//from the nth ~ (i_diag+1)
         {
            if(fabs(d_temp)< fabs(d_Mat[i][i_diag]) )
            {
               i_loc = i ;
               d_temp = d_Mat[i][i_diag] ;
            }
         }

         //When Pivot is necessary!
         if(i_loc != i_diag)
         {
            for(int i=0; i<n; i++)
            {
               d_temp    = d_Mat[i_diag][i] ; //store temporarily
               d_Mat[i_diag][i] = d_Mat[i_loc][i] ;//exchange with the max.
               d_Mat[i_loc][i]  = d_temp ;

               d_temp    = d_Inv[i_diag][i] ; //store temporarily
               d_Inv[i_diag][i] = d_Inv[i_loc][i] ;//exchange with the max.
               d_Inv[i_loc][i]  = d_temp ;
            }
         }
         else if(i_loc == i_diag) //can't find max value than current diagonal term
         {
            cerr<<"Error! no pivoting is possible with current mat. in invers matrix";
            return;
         }
      }


      c = d_Mat[i_diag][i_diag];  //get original or exchanged one


      //make current diag. term 1.0
      for(int i=0; i<n; i++)
      {
         d_Mat[i_diag][i]= d_Mat[i_diag][i] /c; //divide [i_diag]th row terms with the diagonal term
         d_Inv[i_diag][i]= d_Inv[i_diag][i] /c;
      }

      //cout<<'\n'<< "i_diag = "<<i_diag <<endl;
      //cout<< "   d_Mat[i_diag][]   d_Inv[i_diag][] "<<endl;
      //for(int i11=0; i11<n; i11++) cout<<d_Mat[i_diag][i11]<<'\t';
      //for(int i11=0; i11<n; i11++) cout<<d_Inv[i_diag][i11]<<'\t';
      //cout << endl;

      //Lower term elimination================
      //eliminate below terms of [i_diag]th col
      if(i_diag == n-1) continue; // don't need below for the last diagonal

      for(int j=i_diag+1; j<n; j++)
      {
         c=d_Mat[j][i_diag] ; //jth row first term

         for(int i=0; i<n; i++)   //all jth row terms
         {
            d_Mat[j][i] = d_Mat[j][i] -c* d_Mat[i_diag][i];
            d_Inv[j][i] = d_Inv[j][i] -c* d_Inv[i_diag][i];
         }
      }

      //cout<< "afterLOW d_Mat[][]"<<endl;
      //for(int i22=0; i22<n; i22++)
      //{
      //   for(int i11=0; i11<n; i11++)
      //   {
      //      cout<<d_Mat[i22][i11]<<'\t';
      //   }
      //   cout<<endl;
      //}

      //cout<< "afterLOW   d_Inv[][] "<<endl;
      //for(int i22=0; i22<n; i22++)
      //{
      //   for(int i11=0; i11<n; i11++)
      //   {
      //      cout<<d_Inv[i22][i11]<<'\t';
      //   }
      //   cout<<endl;
      //}

   }

   //Upper term elimination================
   //subtract upper terms of [i_diag]th col
   for(int i_diag=1; i_diag<n; i_diag++) //note: begin from the second row
   {
      for(int j=0; j<i_diag; j++)
      {
         c=d_Mat[j][i_diag] ; //jth row first term

         for(int i=0; i<n; i++)   //all jth row terms
         {
            d_Mat[j][i] = d_Mat[j][i] -c* d_Mat[i_diag][i];
            d_Inv[j][i] = d_Inv[j][i] -c* d_Inv[i_diag][i];
         }
      }
   }

   return; 
}
//==============================================================================
//==============================================================================



void dMatrix_Mul_AB(double** A, int n_row, int n_col1,
                double** B, int n_col2,
                double** AB)
//Description================================
//  matrix multiplication
//  C = A*B
//IN   :double** A(n_row, n_col1)
//     :double** B(n_col1, n_col2)
//OUT  :double** AB(n_row, n_col2)
//===========================================
{
   const double tolerance = 10.E-15;
   double d_temp=0.0;

   for(int ic=0; ic<n_col2; ic++)
   {
      for(int ir=0; ir<n_row; ir++)
      {
         d_temp=0.0;

         for(int i=0; i<n_col1; i++)
         {
            d_temp = d_temp + A[ir][i]*B[i][ic] ;
         }
         if(fabs(d_temp) < tolerance ) d_temp =0.0 ; //delete numerical error

         AB[ir][ic] = d_temp ;
      }
   }

   return;
}
//==============================================================================
//==============================================================================


void dMatrix_Mul_AtB(double** A, int n_row, int n_col1,
                double** B, int n_col2,
                double** AtB)
//Description================================
//  matrix multiplication
//  AtB = transpose(A)*B
//IN   :double** A(n_row, n_col1)
//     :double** B(n_row, n_col2)
//OUT  :double** AtB(n_col1, n_col2)
//===========================================
{
   const double tolerance = 10.E-15;
   double d_temp=0.0;

   for(int ic=0; ic<n_col2; ic++)
   {
      for(int ir=0; ir<n_col1; ir++)
      {
         d_temp=0.0;

         for(int i=0; i<n_row; i++)
         {
            d_temp = d_temp + A[i][ir]*B[i][ic] ;
         }

         if(fabs(d_temp) < tolerance ) d_temp =0.0 ; //delete numerical error
         AtB[ir][ic] = d_temp ;
      }
   }

   return;
}

//==============================================================================
//==============================================================================

void dMatrix_Mul_AtBA(double** A, const int n_row, const int n_col,
                      double** B,
                      double** AtBA)
//Description================================
//  matrix multiplication
//  AtBA = transpose(A)*B*A
//IN   :double** A(n_row, n_col)
//     :double** B(n_row, n_row)
//OUT  :double** AtBA(n_col, n_col)
//===========================================
{

   const double tolerance = 10.E-15;
   double d_temp=0.0;

   //double AtB[n_col][n_row] ;
   double** AtB = New_dMatrix(n_col,n_row) ;

   for(int i=0; i<n_col; i++) //initialize
   {
      for(int j=0; j<n_row; j++)
      {
         AtB[i][j] = 0.0 ;
      }
   }

   //AtB
   for(int ic=0; ic<n_row; ic++)
   {
      for(int ir=0; ir<n_col; ir++)
      {
         d_temp=0.0;

         for(int i=0; i<n_row; i++)
         {
            d_temp = d_temp + A[i][ir]*B[i][ic] ;
         }

         if(fabs(d_temp) < tolerance ) d_temp =0.0 ; //delete numerical error
         AtB[ir][ic] = d_temp ;
      }
   }

   //AtBA
   d_temp=0.0;

   for(int ic=0; ic<n_col; ic++)
   {
      for(int ir=0; ir<n_col; ir++)
      {
         d_temp=0.0;

         for(int i=0; i<n_row; i++)
         {
            d_temp = d_temp + AtB[ir][i]*A[i][ic] ;
         }
         if(fabs(d_temp) < tolerance ) d_temp =0.0 ; //delete numerical error

         AtBA[ir][ic] = d_temp ;
      }
   }
   
   Del_dMatrix(AtB, n_col,n_row) ;

   return;

}

void dMatrix_Mul_AtBA_Yicheng(double** A, const int n_row, const int n_col,
	double** B,
	double* AtBA)
	//Description================================
	//  matrix multiplication
	//  AtBA = transpose(A)*B*A
	//IN   :double** A(n_row, n_col)
	//     :double** B(n_row, n_row)
	//OUT  :double** AtBA(n_col, n_col)
	//===========================================
{

	const double tolerance = 10.E-15;
	double d_temp = 0.0;

	//double AtB[n_col][n_row] ;
	double** AtB = New_dMatrix(n_col, n_row);

	for (int i = 0; i<n_col; i++) //initialize
	{
		for (int j = 0; j<n_row; j++)
		{
			AtB[i][j] = 0.0;
		}
	}

	//AtB
	for (int ic = 0; ic<n_row; ic++)
	{
		for (int ir = 0; ir<n_col; ir++)
		{
			d_temp = 0.0;

			for (int i = 0; i<n_row; i++)
			{
				d_temp = d_temp + A[i][ir] * B[i][ic];
			}

			if (fabs(d_temp) < tolerance) d_temp = 0.0; //delete numerical error
			AtB[ir][ic] = d_temp;
		}
	}

	//AtBA
	d_temp = 0.0;

	for (int ic = 0; ic<n_col; ic++)
	{

		d_temp = 0.0;

		for (int i = 0; i<n_row; i++)
		{
			d_temp = d_temp + AtB[ic][i] * A[i][ic];
		}
		if (fabs(d_temp) < tolerance) d_temp = 0.0; //delete numerical error

		AtBA[ic] = d_temp;

	}

	Del_dMatrix(AtB, n_col, n_row);

	return;

}

//==============================================================================
//==============================================================================


void dMatrix_dVector_Mul_Av(double** A, int n_row, int n_col,
                            double   v[],
                            double  Av[])
//Description================================
//  matrix & Vector multiplication
//  Av = A*v
//IN   :double** A(n_row, n_col)
//     :double   v(n_col)
//OUT  :double   Av(n_row)
//===========================================
{
   const double tolerance = 10.E-15;
   double d_temp=0.0;

   for(int i_r=0; i_r<n_row; i_r++)
   {
      d_temp =0.0;
      for(int i_c=0; i_c<n_col; i_c++)
      {
         d_temp = d_temp + A[i_r][i_c]*v[i_c] ;
      }
      if(fabs(d_temp) < tolerance) d_temp =0.0;

      Av[i_r] = d_temp ;
   }
   return;
}


//==============================================================================
//==============================================================================

void dMatrix_dVector_Mul_Atv(double** A, int n_row, int n_col,
                            double   v[],
                            double  Atv[])
//Description================================
//  Transpose(matrix) & Vector multiplication
//  Atv = Transpose(A)*v
//IN   :double** A(n_row, n_col)
//     :double   v(n_row)
//OUT  :double   Atv(n_col)
//===========================================
{
   const double tolerance = 10.E-15;
   double d_temp=0.0;

   for(int i_c=0; i_c<n_col; i_c++)
   {
      d_temp =0.0;
      for(int i_r=0; i_r<n_row; i_r++)
      {
         d_temp = d_temp + A[i_r][i_c]*v[i_r] ;
      }
      if(fabs(d_temp) < tolerance) d_temp =0.0;

      Atv[i_c] = d_temp ;
   }
   return;
}

//==============================================================================
//==============================================================================

double dMaxValue(double** d_matrix, int n_row, int n_col, char s_where,
              int n_begin, int n_end, int n_at)
//Description=================
//
// find maximum double value in the double matrix
// (1) s_where ="row" ; search through row(within n_begin~n_end) at n_at col
// (2) s_where ="col" ; search through col(within n_begin~n_end) at n_at row
// (3) s_where ="all" ; search all matrix
//
//IN   : double** d_matrix[n_row][n_col]
//OUT  : double maximum value
//============================

{
   double d_maximum =0.0;
   double d_temp=0.0;

   if(s_where =='r')
   {
      for(int i=n_begin;i<=n_end; i++)
      {
         if(d_temp < d_matrix[i][n_at]) d_temp = d_matrix[i][n_at] ;
      }
   }
   else if(s_where =='c')
   {
      for(int i=n_begin;i<=n_end; i++)
      {
         if(d_temp < d_matrix[n_at][i]) d_temp = d_matrix[n_at][i] ;
      }
   }
   else if(s_where =='a')
   {
      for(int i=0;i<n_row; i++)
      {
         for(int j=0; j<n_col; j++)
         {
            if(d_temp < d_matrix[i][j]) d_temp = d_matrix[i][j] ;
         }
      }

   }

   d_maximum = d_temp;
   return d_maximum;
}

//==============================================================================
//==============================================================================

void Compare_Two_dMatrix(double** A, double** B, int n_row, int n_col)
{
   for(int i=0; i<n_row; i++)
   {
      for(int j=0; j<n_col; j++)
      {
         if(fabs(A[i][j] - B[i][j]) != 0.0)
         {
            cout << "difference component at (i,j)=("<<i+1<<","<<j+1<<")" <<endl;
            cout << "Aij = " <<A[i][j] <<"   Bij ="<<B[i][j]<<endl;
         }
      }
   }
   system("PAUSE") ;
}


//=============================================================================
//=============================================================================
void c1A_p_c2B(const double c1, double** A, const int n_row, const int n_col,
               const double c2, double** B,
               double ** M)
//Description=========================================
//  perform
//     M = c1*[A] + c2*[B]
//====================================================
{
   for(int i=0; i<n_row; i++)
   {
      for(int j=0; j<n_col; j++)
      {
         M[i][j] = c1*A[i][j] + c2*B[i][j] ;
      }
   }
}

//=============================================================================
//=============================================================================
double my_dot(const int n, const double* u, const double * v)
//Description===============
//  inner product of two vectors 
//  double = {u}.{v}
//==========================
{
   double d=0.0;

   for(int i=0; i<n; i++)
      d += u[i]*v[i];
   return d;
}


bool Inverse_dMatrix_FHDI(double** d_Mat, const int n, double** d_Inv)
//Description===================
//return inverse matrix of the n*n matrix
//using Gauss-Jordan elimination
//
// if diagonal term is zero and too small
// perform pivoting with largest value on the columns
//
// Note: for FHDI, n = 1 and n =2 cases are separately handled. 
//
//IN    : double** d_Mat[n][n]
//OUT   : double** d_Inv[n][n]
//OUT   : bool b_success = 0 when abrupt exit due to zero digonal term 
//==============================
{
   bool b_success = true; 
   const double eps=1.e-15 ;

   //------------------
   //if n = 1
   //------------------
    if(fabs(d_Mat[0][0]) > eps) 
    {
	   d_Inv[0][0] = 1.0/d_Mat[0][0];
	   return b_success; 
	}
    if(fabs(d_Mat[0][0]) <= eps) 
    {
	   d_Inv[0][0] = 1.0;
	   b_success = false; 
	   return b_success; 
	}	

   //------------------
   //if n = 2
   //------------------
   const double det2 = d_Mat[0][0]*d_Mat[1][1] - d_Mat[0][1]*d_Mat[1][0];
    if(fabs(det2) > eps) 
    {
	   d_Inv[0][0] = d_Mat[1][1]/det2;
	   d_Inv[0][1] = -1.0*d_Mat[0][1]/det2;
	   d_Inv[1][0] = -1.0*d_Mat[1][0]/det2;
	   d_Inv[1][1] = d_Mat[0][0]/det2;
	   
	   return b_success; 
	}
    if(fabs(det2) <= eps) 
    {
	   Fill_dMatrix(d_Inv, 2, 2, 1.0);
	   b_success = false; 
	   return b_success; 
	}	
	
	//------------------
	// below is for n > 2
	//------------------
   //make d_Inv unity matrix
   for(int i=0; i<n; i++)
   {
      for(int j=0; j<n; j++) d_Inv[i][j] =0.0;

      d_Inv[i][i]=1.0 ;
   }



   double c=0.0 ; //coeff.

   for(int i_diag=0; i_diag<n; i_diag++)
   {
      c = d_Mat[i_diag][i_diag];

      //================
      //when diagonal term is too small or zero, then needs pivoting
      //================
      if(fabs(c)<eps )
      {
         //find max on current column ======
         double d_temp= c;  int i_loc=i_diag;

         for(int i=(n-1); i>i_diag; i--)//from the nth ~ (i_diag+1)
         {
            if(fabs(d_temp)< fabs(d_Mat[i][i_diag]) )
            {
               i_loc = i ;
               d_temp = d_Mat[i][i_diag] ;
            }
         }

         //When Pivot is necessary!
         if(i_loc != i_diag)
         {
            for(int i=0; i<n; i++)
            {
               d_temp    = d_Mat[i_diag][i] ; //store temporarily
               d_Mat[i_diag][i] = d_Mat[i_loc][i] ;//exchange with the max.
               d_Mat[i_loc][i]  = d_temp ;

               d_temp    = d_Inv[i_diag][i] ; //store temporarily
               d_Inv[i_diag][i] = d_Inv[i_loc][i] ;//exchange with the max.
               d_Inv[i_loc][i]  = d_temp ;
            }
         }
         else if(i_loc == i_diag) //can't find max value than current diagonal term
         {
			//----
			//below condition is added for FHDI
			//----
			if(fabs(c) < eps)
			{ 
				cerr<<"Error! no pivoting is possible with current mat. in invers matrix";
				b_success = false; 
				return b_success;
			}
			if(fabs(c) >= eps)
			{
				//keep going with current Non-zero diagonal value 
			}
			
         }
      }


      c = d_Mat[i_diag][i_diag];  //get original or exchanged one


      //make current diag. term 1.0
      for(int i=0; i<n; i++)
      {
         d_Mat[i_diag][i]= d_Mat[i_diag][i] /c; //divide [i_diag]th row terms with the diagonal term
         d_Inv[i_diag][i]= d_Inv[i_diag][i] /c;
      }

      //cout<<'\n'<< "i_diag = "<<i_diag <<endl;
      //cout<< "   d_Mat[i_diag][]   d_Inv[i_diag][] "<<endl;
      //for(int i11=0; i11<n; i11++) cout<<d_Mat[i_diag][i11]<<'\t';
      //for(int i11=0; i11<n; i11++) cout<<d_Inv[i_diag][i11]<<'\t';
      //cout << endl;

      //Lower term elimination================
      //eliminate below terms of [i_diag]th col
      if(i_diag == n-1) continue; // don't need below for the last diagonal

      for(int j=i_diag+1; j<n; j++)
      {
         c=d_Mat[j][i_diag] ; //jth row first term

         for(int i=0; i<n; i++)   //all jth row terms
         {
            d_Mat[j][i] = d_Mat[j][i] -c* d_Mat[i_diag][i];
            d_Inv[j][i] = d_Inv[j][i] -c* d_Inv[i_diag][i];
         }
      }

      //cout<< "afterLOW d_Mat[][]"<<endl;
      //for(int i22=0; i22<n; i22++)
      //{
      //   for(int i11=0; i11<n; i11++)
      //   {
      //      cout<<d_Mat[i22][i11]<<'\t';
      //   }
      //   cout<<endl;
      //}

      //cout<< "afterLOW   d_Inv[][] "<<endl;
      //for(int i22=0; i22<n; i22++)
      //{
      //   for(int i11=0; i11<n; i11++)
      //   {
      //      cout<<d_Inv[i22][i11]<<'\t';
      //   }
      //   cout<<endl;
      //}

   }

   //Upper term elimination================
   //subtract upper terms of [i_diag]th col
   for(int i_diag=1; i_diag<n; i_diag++) //note: begin from the second row
   {
      for(int j=0; j<i_diag; j++)
      {
         c=d_Mat[j][i_diag] ; //jth row first term

         for(int i=0; i<n; i++)   //all jth row terms
         {
            d_Mat[j][i] = d_Mat[j][i] -c* d_Mat[i_diag][i];
            d_Inv[j][i] = d_Inv[j][i] -c* d_Inv[i_diag][i];
         }
      }
   }

   return b_success; 
}
#endif /* MATRIX_UTILITY_FHDI */