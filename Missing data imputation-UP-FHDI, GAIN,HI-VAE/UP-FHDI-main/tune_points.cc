void tune_points(int& startpoint, int& endpoint, int* i_source, int mynode, int totalnodes,int n) {
	if (mynode != 0 && mynode != 1 && mynode != (totalnodes - 1)) {
		if (i_source[startpoint] == i_source[startpoint - 1]) {
			for (int i = startpoint + 1; i < n;i++) {
				if (i_source[startpoint] != i_source[i]) {
					startpoint = i;
					break;
				}
				if (i == (n - 1)) {
					startpoint = i + 1;
				}
			}
		}

		if (i_source[endpoint - 1] == i_source[endpoint]) {
			for (int i = endpoint;i < n;i++) {
				//cout << "i1:" << i << endl;
				if (i_source[endpoint - 1] != i_source[i]) {
					//cout << "i2:" << i << endl;
					endpoint = i;
					break;
				}
				if (i == (n - 1)) {
					endpoint = i + 1;
				}
			}
		}
	}

	if (mynode == 1) {
		//for (int i = 0; i < n; i++) {
		//	cout << "i_source[" << i << "]: " << i_source[i] << endl;
		//}

		//cout << "endpoint: " << endpoint <<endl;
		//cout << "I'm here" << endl;
		//cout <<"endpoint: "<< endpoint <<"    i_source[endpoint-1]: " << i_source[endpoint - 1] << ", i_source[endpoint]: " << i_source[endpoint] << endl;
		if (i_source[endpoint] == i_source[endpoint - 1]) {
			for (int i = endpoint;i < n;i++) {
				//cout << "i1:" << i << endl;
				if (i_source[i] != i_source[endpoint - 1]) {
					//cout << "i2:" << i << endl;
					endpoint = i;
					break;
				}
				if (i == (n - 1)) {
					endpoint = i + 1;
				}
			}
		}
	}

	if (mynode == (totalnodes - 1)) {
		if (i_source[startpoint] == i_source[startpoint - 1]) {
			for (int i = startpoint + 1; i < n;i++) {
				if (i_source[startpoint] != i_source[i]) {
					startpoint = i;
					break;
				}
				if (i == (n - 1)) {
					startpoint = i + 1;
				}
			}
		}
	}
	//cout<<"I'm out"<<endl;
}