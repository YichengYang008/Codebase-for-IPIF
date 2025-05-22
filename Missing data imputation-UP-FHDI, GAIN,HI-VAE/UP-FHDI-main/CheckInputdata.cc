
void CheckInputdata(double **x_raw, int **r_raw, const int nrow_x, const int ncol_x, ofstream& TestOut) {

	//=====================================================================
	// Check if the input data has instances without observed values
	//=====================================================================

	for (int i = 0;i < nrow_x; i++) {
		int sum_temp = 0;
		for (int j = 0;j < ncol_x;j++) {
			sum_temp = sum_temp + r_raw[i][j];
		}

		if (sum_temp == 0) {
			TestOut << "Caution! The " << i + 1 << "th row of data set has no observed values" << endl;
		}
	}


	//========================================================================
	// Check if any variables with almost 0 variance
	// variance = sum(x-x_bar)^2/(n-1)
	//=======================================================================

	for (int j = 0; j < ncol_x; j++) {

		//----------------
		//Prepare fully observed variable
		//---------------------
		std::vector<int> ol;
		ol.clear();

		for (int i = 0; i < nrow_x; i++) {

			//if (i == 8) cout << "raw[" << i << "][" << j << "]: " << r_raw[j][i] << endl;
			if (fabs(r_raw[i][j]) > 1e-15) //this row has no missing cells
			{
				ol.push_back(i);
			} //actual number of the row having no missing cells
		}

		int ol_size = ol.size();

		double mean = 0.0;
		double sum_temp2 = 0.0;
		for (int i = 0; i < ol_size; i++) {
			sum_temp2 = sum_temp2 + x_raw[ol[i]][j];
		}
		mean = sum_temp2 / ol_size;

		double sum_temp3 = 0.0;
		for (int i = 0; i < ol_size; i++) {
			sum_temp3 = sum_temp3 + (x_raw[ol[i]][j] - mean)*(x_raw[ol[i]][j] - mean);
		}

		//cout<<"The variance at variable "<<j<<" is "<< sum_temp3 / (ol_size - 1) <<", where ol is "<< ol_size <<" and mean is "<<mean<<endl;
		double variance = sum_temp3 / (ol_size - 1);

		if (variance < 1e-3) {
			TestOut << "Error! The variance of " << j << "th variable is "<< variance <<" approaching 0. Please consider removing this variable from daty" << endl;
		}

	}


}