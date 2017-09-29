/*
 * Ratetable.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: ali
 */

#include "Ratetable.h"
#include <iostream>

using namespace std;
// 1SS rates need to be hard-coded 2SS and 3SS can be calculated from 1SS
float Rate_table::tput[24][2][2] = { { {6.5 , 7.2}, {13.5, 15} },
							{ {13, 14.4}, {27, 30} },
							{ {19.5, 21.7}, {40.5, 45} },
							{ {26, 28.9}, {54, 60} },
							{ {39, 43.3}, {81, 90} },
							{ {52, 57.8}, {108, 120} },
							{ {58.5, 65}, {121.5, 135} },
							{ {65, 72.2}, {135, 150} } };

Rate_table::Rate_table() {
	// TODO Auto-generated constructor stub
	//initializing the table
	//2SS rates
	for (int i=8; i<16; i++){
		for (int sgi=0; sgi<=1; sgi++){
			for (int ht40=0; ht40<=1; ht40++){
				tput[i][sgi][ht40]= 2 * tput[i-8][sgi][ht40];
			}
		}
	}

	//3SS rates
	for (int i=16; i<24; i++){
		for (int sgi=0; sgi<=1; sgi++){
			for (int ht40=0; ht40<=1; ht40++){
				tput[i][sgi][ht40]= 3 * tput[i-16][sgi][ht40];
			}
		}
	}

/*	for (int i=0; i<24; i++){
		for (int sgi=0; sgi<=1; sgi++){
			for (int ht40=0; ht40<=1; ht40++){
				cout << tput[i][sgi][ht40] << "\t";
			}
		}
		cout << endl;
	}*/
}

Rate_table::~Rate_table() {
	// TODO Auto-generated destructor stub
}

float Rate_table::get_tput (int idx, int sgi, int ht40)
{
	return tput[idx][ht40][sgi];
}

float Rate_table::get_tput (int idx)
{
	int ht40 , sgi;

	if (idx < 48){
		ht40 = 0;
	}
	else{
		ht40 = 1;
	}

	if ( (idx % 48) < 24 ){
		sgi = 0;
	}
	else{
		sgi = 1;
	}

	return Rate_table::get_tput(idx % 24, sgi, ht40);
}




