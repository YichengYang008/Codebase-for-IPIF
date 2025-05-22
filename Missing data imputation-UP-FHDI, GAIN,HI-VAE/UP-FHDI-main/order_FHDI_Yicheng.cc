#include "tune_points.cc"
void order_FHDI_Yicheng(int* i_original, const int n)
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

	double order_begin = MPI_Wtime();
	//cout << "Order_FHDI here" << endl;
	int* i_source = new int[n];
	int* i_order = new int[n];

	for (int i = 0; i < n; i++)
	{
		i_source[i] = i_original[i]; //backup
		i_order[i] = i + 1; //default
	}
	//if (mynode == 1) {
	//	for (int i = 0; i < n; i++) {
	//		cout << "i_source[" << i << "]: " << i_source[i] << endl;
	//	}
	//}
	//-----------
	//leverage sorting library
	//-----------
	std::sort(i_source, i_source + n);

	//if (mynode == 1) {
	//	for (int i = 0; i < n; i++) {
	//		cout << "i_source[" << i << "]: " << i_source[i] << endl;
	//	}
	//}

	int i_now = 0;
	//cout << "Order_FHDI1 time is " << MPI_Wtime() - order_begin << endl;
	//double order_begin2 = MPI_Wtime();
	i_order[0] = 1; //first cell location as default
	//cout << "n value is: " << n << endl;
	int numWorkPerProc = (int)floor(1.0*n / (1.0*totalnodes - 1));
	int numWorkLocalLast = n - numWorkPerProc * (totalnodes - 2);
	//cout << "numWorkPerProc on node " << mynode << "is " << numWorkPerProc << endl;
	//cout << "numWorkLocalLast on node " << mynode << "is " << numWorkLocalLast << endl;
	int startpoint = 0;
	int endpoint = 0;

	if (mynode != 0 && mynode != (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = mynode*numWorkPerProc;
		//cout << "Strating point and ending point of order_FHDI before on node " << mynode << " are: " << startpoint << ", " << endpoint << endl;
	}

	if (mynode == (totalnodes - 1)) {
		startpoint = (mynode - 1)*numWorkPerProc;
		endpoint = (mynode - 1)*numWorkPerProc + numWorkLocalLast;
		//cout << "Strating point and ending point of order_FHDI before on node " << mynode << " are: " << startpoint << ", " << endpoint << endl;
	}

	tune_points(startpoint, endpoint, i_source, mynode, totalnodes, n);

	//if (mynode == 1) {
	//	cout << "Successful out tunepoints" << endl;
	//}
	//if (mynode != 0 && mynode != (totalnodes - 1)) {
	//	cout << "Strating point and ending point of order_FHDI after on node " << mynode << " are: " << startpoint << ", " << endpoint << endl;
	//}

	//if (mynode == (totalnodes - 1)) {
	//	cout << "Strating point and ending point of order_FHDI after on node " << mynode << " are: " << startpoint << ", " << endpoint << endl;
	//}


	//--------------------------------------------------------------------------------------------
	int* i_order_temp = new int[endpoint - startpoint];
	int clock = 0;

	for (int i = startpoint; i < endpoint; i++)
	{
		i_now = i_source[i];
		//----------------
		//comparisons from the first entiry to now 
		//----------------
		for (int j = 0; j < n; j++)
		{
			if (fabs(i_now - i_original[j]) < 1e-3)
			{
				i_order_temp[clock] = j + 1; //Actual location
				i_original[j] = -1; //dummy value
				break;
			}
		}
		clock++;
	}


	//if (mynode == 1) {
	//	cout << "Order_check1 " << endl;
	//}

	if (mynode != 0) {
		MPI_Send(&startpoint, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&endpoint, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	}


	if (mynode != 0) {
		MPI_Send(i_order_temp, (endpoint - startpoint), MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	//int *i_order_recv;
	int * i_order_total = new int[n];
	if (mynode == 0) {
		for (int j = 1; j < totalnodes; j = j + 1) {
			MPI_Recv(&startpoint, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&endpoint, 1, MPI_INT, j, 2, MPI_COMM_WORLD, &status);

			int* i_order_recv = new int[endpoint - startpoint];
			MPI_Recv(i_order_recv, (endpoint - startpoint), MPI_INT, j, 1, MPI_COMM_WORLD, &status);

			int time_counter = 0;
			for (int k = startpoint; k < endpoint; k++) {
				i_order_total[k] = i_order_recv[time_counter];
				time_counter++;
			}
			delete[] i_order_recv;
		}


		for (int i = 0; i<n; i++)
		{
			i_original[i] = i_order_total[i]; //backup
			//cout<<"i_order_total["<<i<<"]: "<< i_original[i]<<endl;
		}

	}

	//if (mynode == 1) {
	//	cout << "Order_check2 " << endl;
	//}
	//-------------------------------------------------------------------------------------------------
	//---prep return
	//cout << "Order_FHDI2 time is " << MPI_Wtime() - order_begin2 << endl;
	//double order_begin3 = MPI_Wtime();

	MPI_Bcast(i_original,n,MPI_INT,0,MPI_COMM_WORLD);
	
	//for (int i = 0; i < n; i++)
	//{
	//	cout << "mynode: " << mynode << " i_original [" << i << "]: " << i_original[i] << endl;
	//}

	//if (mynode == 1) {
	//	cout << "Order_check3 " << endl;
	//}
	//cout << "Order_FHDI3 time is " << MPI_Wtime() - order_begin3 << endl;
	delete[] i_source;
	delete[] i_order;
	delete[] i_order_temp;
	delete[] i_order_total;
	//delete[] i_order_recv;
	return;
}