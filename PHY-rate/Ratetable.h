/*
 * Ratetable.h
 *
 *  Created on: Jul 29, 2015
 *      Author: ali
 */

#ifndef RATETABLE_H_
#define RATETABLE_H_

class Rate_table {

public:
	Rate_table();
	virtual ~Rate_table();
	float get_tput (int idx, int sgi, int ht40);
	float get_tput (int idx);

private:
	static float tput[24][2][2];

};



#endif /* RATETABLE_H_ */
