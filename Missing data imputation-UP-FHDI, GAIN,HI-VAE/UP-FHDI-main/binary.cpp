#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
//#include "order.cpp"
using namespace std;
int binarySearch(int arr[], int l, int r, int x) {
	if (r >= l) {
		int mid = l + (r - l) / 2;

		// If the element is present at the middle 
		// itself 
		if (arr[mid] == x)
			return mid;

		// If element is smaller than mid, then 
		// it can only be present in left subarray 
		if (arr[mid] > x)
			return binarySearch(arr, l, mid - 1, x);

		// Else the element can only be present 
		// in right subarray 
		return binarySearch(arr, mid + 1, r, x);
	}

	// We reach here when element is not 
	// present in array 
	return -1;
}


void binary(int& left, int& right, int target, int n,int* arr)
{
	int mynode, totalnodes;
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Status status;

	//cout<<"BINARY"<<endl;
	//int n = sizeof(arr) / sizeof(arr[0]);
	//cout << "n binary: " << n << endl;
	int result = binarySearch(arr, 0, n - 1, target);


	//cout << "result: " << result << endl;
	if (result == -1) {
		left = result;
		right = result;
	}

	if (result != -1) {
		int left_temp = result;
		while (arr[left_temp] == target) {
			left_temp--;
			//cout << "left: " << left << endl;
		}
		left = left_temp + 1;
		//cout << "left: " << left + 1 << endl;
		int right_temp = result;
		while (arr[right_temp] == target) {
			right_temp++;
		}
		right = right_temp - 1;
	}

	//cout << "left: "<<left<<"and right: " << right << endl;
}

