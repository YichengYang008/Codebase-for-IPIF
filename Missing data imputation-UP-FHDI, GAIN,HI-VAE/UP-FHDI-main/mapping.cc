#include <vector>
#include <map>

// voteCand: key, C++ map value is actually a vector that allows multiple values for a single key
void mapping(map< vector<double>, int >& freqMap, const vector<double>& voteCand) {
	// If the key does not exist create new key, value pair and insert it to map
	if (freqMap.find(voteCand) == freqMap.end()) {
		freqMap.insert(pair<vector<double>, int >(voteCand, 1));
	}
	// Otherwise, update value (count)
	else {
		freqMap.find(voteCand)->second++;
	}
	return;
}

//void mapping(int index, map< vector<double>, int > *arr, vector<double> voteCand) {
//	for (vector<double>::iterator it = voteCand.begin(); it != voteCand.end(); ++it) {
//		map< vector<double>, int >::iterator mapIter;
//		mapIter = arr[index].find(*it);
//		// If the key does not exist create new key, value pair and insert it to map
//		if (mapIter == arr[index].end()) { // *it ---> vector<double>
//			arr[index].insert(pair<vector<double>, int >(*it, 1));
//		}
//		// Otherwise, update value (count)
//		else {
//			mapIter->second++;
//		}
//	}
//	return;
//}