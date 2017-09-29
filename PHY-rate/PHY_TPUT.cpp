//============================================================================
// Name        : PHY_TPUT.cpp
// Author      : Ali Abedi
// Version     :
// Copyright   : Your copyright notice
// Description : 802.11n PHY layer Throughput in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstdlib>
#include "Ratetable.h"

using namespace std;

int main(int argc, char* argv[]) {

	if (argc == 1){
		cout << "Usage: Rate[0..95] or Rate[0..7] sgi[0..1] ht40[0..1]" << endl;
		return -1;
	}
	Rate_table rt = Rate_table();

	if (argc == 2){
		cout << rt.get_tput(atoi(argv[1])) << endl;
	}
	else if (argc == 4){
		cout << rt.get_tput(atoi(argv[1]), atoi(argv[2]), atoi(argv[3])) << endl;
	}

	return 0;
}
