#include "binary.cpp"
void order_FHDI_binary(int* i_original, const int n)

//Description ================================

// Order the positive integer array in ascending order

//

//INOUT   : int i_original_0(n) returned with the ordered (Actual) cell numbers     

//          i_original > 0

//=============================================

{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//cout<<"Yeah2 how do you feel"<<endl;

	int* i_source = new int[n];

	int* i_order = new int[n];



	for (int i = 0; i<n; i++)

	{

		i_source[i] = i_original[i]; //backup

		i_order[i] = i + 1; //default

	}



	//-----------

	//leverage sorting library

	//-----------

	std::sort(i_source, i_source + n);

	int i_now = 0;

	//for (int k = 0;k < n;k++) {
	//	cout<<"i_source["<<k<<"]: "<< i_source [k] <<endl;
	//}


	i_order[0] = 1; //first cell location as default

	//for (int i = 0; i<n; i++)

	//{

	//	i_now = i_source[i];

	//	//----------------

	//	//comparisons from the first entiry to now 

	//	//----------------

	//	for (int j = 0; j<n; j++)

	//	{

	//		if (fabs_FHDI(i_now - i_original[j])<1e-3)

	//		{

	//			i_order[i] = j + 1; //Actual location

	//			i_original[j] = -1; //dummy value

	//			break;

	//		}

	//	}

	//}
	//for (int k = 0;k < n;k++) {
	//	cout << "i_original[" << k << "]: " << i_original[k] << endl;
	//}


	for (int i = 0;i < n;i++) {
		i_now = i_original[i];

		if (i == 0) {
			int index = i;
			int left = 0;
			int right = 0;
			binary(left, right, i_now, n, i_source);
			for (int k = left;k < (right + 1);k++) {
				i_order[k] = index + 1;
				index++;
			}
		}

		if (i != 0) {
			if (i_now != i_original[i - 1]) {
				int index = i;
				int left = 0;
				int right = 0;
				binary(left, right, i_now, n, i_source);
				for (int k = left;k < (right + 1);k++) {
					i_order[k] = index + 1;
					index++;
				}
			}
		}

	}


	//---prep return

	for (int i = 0; i<n; i++)

	{

		i_original[i] = i_order[i]; //backup

	}



	delete[] i_source;

	delete[] i_order;



	return;

}