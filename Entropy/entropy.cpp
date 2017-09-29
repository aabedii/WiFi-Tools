#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

#include <boost/math/distributions/students_t.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace std;
using namespace boost::math;
using namespace boost::accumulators;

#define NUM_ARGS 8
#define NUM_LEVELS 101
#define FILES_PATH "results/"
#define NUM_RATES 24

void calculate_pdf_joint(vector<int> &x, vector<int> &y, float pdf_xy[][NUM_LEVELS])
{
  for (int i=0; i<x.size(); i++){
    pdf_xy[ x[i] ][ y[i] ] += 1.0 / (float)x.size();
  }
}

void calculate_pdf(vector<int> &x, float  pdf_x[])
{

  for (int i=0; i<x.size(); i++){
    pdf_x[ x[i] ]+= 1.0 / (float)x.size();
  }


}


// Convert char * with format hh:mm:ss to a unique integer
long int time_to_value(char* time){
	  char * pch;
	  long int result = 0;

	  try{
		  pch = strtok (time,":");
		  for (int i=2; i>=0 ; i--){
			  result += atoi(pch) * (pow(60,i));
			  pch = strtok (NULL,":");

		  }
	  }
	  catch (int e){
		  return -1;
	  }

	  return result;
}

void read_file(ifstream &input, vector<float> &x, int start_time, int end_time)
{
  if (!input.is_open()){
    cout << "Could not open the file!" << endl;
    return;
  }

  float tmp;
  string junk, time;
  int relative_time;


  input >> time >> tmp  >> junk;
  long int trace_start_time = time_to_value((char*)time.c_str());

  //cout << "trace start = " << trace_start_time<< "\t" << (char*)time.c_str()<<endl;

  while (!input.eof()) {
    input >> time >> tmp  >> junk;
    relative_time = time_to_value((char*)time.c_str()) - trace_start_time;
   //cout << relative_time << endl;
    if ( relative_time < start_time){
    	continue;
    }
    if ( relative_time > end_time){
    	//cout << "break time = " << relative_time << endl;
    	break;
    }

    x.push_back ( tmp*100 );
    //multiples of 5
 // x.push_back ( (((int)round(tmp*100)+2)/10)*10  );
  }

}

void read_file(ifstream &input, vector<float> &x)
{
  if (!input.is_open()){
    cout << "Could not open the file!" << endl;
    return;
  }

  float tmp;
  string junk, time;


  while (!input.eof()) {
    input >> time >> tmp  >> junk;
    x.push_back ( tmp*100 );
    //multiples of 5
 // x.push_back ( (((int)round(tmp*100)+2)/10)*10  );
  }

}
void read_file(ifstream &input, vector<int> &x)
{
  if (!input.is_open()){
    cout << "Could not open the file!" << endl;
    return;
  }

  float tmp;
  string junk, time;

  while (!input.eof()) {
    input >> time >> tmp  >> junk;
    x.push_back ( round(tmp*100)  );
    //multiples of 5
 // x.push_back ( (((int)round(tmp*100)+2)/10)*10  );
  }

}

float calculate_entropy (float pdf_x[])
{
  float entropy = 0;
  for (int i=0; i<NUM_LEVELS; i++){
    if (pdf_x[i] == 0)
      continue;
    entropy += pdf_x[i] * log2 (pdf_x[i]) ;
  }

  return -entropy;

}

float calculate_conditional_entropy( float pdf_x[], float pdf_xy[][NUM_LEVELS] )
{
  float entropy = 0;

  for (int i=0; i<NUM_LEVELS; i++){
    for (int j=0; j<NUM_LEVELS; j++){
      if (pdf_xy[i][j] == 0)
	continue;

      entropy += pdf_xy[i][j] * log2 (pdf_x[i] / pdf_xy[i][j]);
    }
  }

  return entropy;
}

float get_mean (vector<float> x){
	float sum =0;
	for (int i=0; i< x.size(); i++){
		sum += x[i];
	}
	return sum / (float)x.size();
}

float get_mean (vector<int> x){
	float sum =0;
	for (int i=0; i< x.size(); i++){
		sum += x[i];
	}
	return sum / (float)x.size();
}

float get_max_dispersion(int r1_mcs, int r1_sgi, int r1_ht40, int r2_mcs, int r2_sgi, int r2_ht40 , char* path, int start_time, int end_time)
{
	vector<float>  X,Y;


	char file1[200];
	char file2[200];

	sprintf(file1, "%s/%d-%d-%d.dat", path, r1_mcs, r1_sgi, r1_ht40 );
	sprintf(file2, "%s/%d-%d-%d.dat", path, r2_mcs, r2_sgi, r2_ht40 );

	//Reading the data files and initialize the X, Y, and XY arrays
	ifstream r1, r2;

	r1.open(file1);
	r2.open(file2);

	read_file(r1, X, start_time, end_time);
	read_file(r2, Y,  start_time, end_time);

	if (X.size() <= 2 || Y.size() <= 2 )
		return 0;

	float max_dispersion = 0;

	const float step = 1;

	float min, max;

	for (float e1=0; e1<=100; e1+=step){
		min = 100;
		max = 0;
		for (int i = 0; i < X.size() && i < Y.size(); i++){
			if (e1-(step/2.0) < X[i] && X[i] < e1+(step/2.0)){
				if (Y[i] < min){
					min = Y[i];
				}
				if (Y[i] > max){
					max = Y[i];
				}
			}
		}
		if (max_dispersion < max-min){
			max_dispersion = max-min;
		}
	}

	return max_dispersion;

}


float get_correlation_coeffcient(int r1_mcs, int r1_sgi, int r1_ht40, int r2_mcs, int r2_sgi, int r2_ht40 , char* path)
{


  vector<float>  X,Y;


  char file1[200];
  char file2[200];

  sprintf(file1, "%s/%d-%d-%d.dat", path, r1_mcs, r1_sgi, r1_ht40 );
  sprintf(file2, "%s/%d-%d-%d.dat", path, r2_mcs, r2_sgi, r2_ht40 );

  //Reading the data files and initialize the X, Y, and XY arrays
  ifstream r1, r2;

  r1.open(file1);
  r2.open(file2);

  read_file(r1, X);
  read_file(r2, Y);

  if (X.size() <= 2 || Y.size() <= 2 )
	  	  return 0;

  float mean_X = get_mean(X);
  float mean_Y = get_mean(Y);

  float S_x = 0;
  float S_y = 0;

  for (int i=0; i< X.size(); i++){
	  S_x += pow( X[i]-mean_X, 2);
	  S_y += pow( Y[i]-mean_Y, 2);
  }

 // S_x *= 1.0/(float)(X.size()-1);
  S_x = sqrt (S_x);

 // S_y *= 1.0/(float)(Y.size()-1);
  S_y = sqrt (S_y);

  if (S_x == 0 && S_y == 0 )
	  return 1;

  if (S_x == 0 || S_y == 0 )
  	  return 0;


  float r=0;

  for (int i=0; i< X.size(); i++){
	  r+= (X[i]-mean_X) * (Y[i]-mean_Y);
  }

  r = r / (S_x * S_y);

 // r = r /(float)(X.size()-1);

  //this is not corect, abs should be removed
  return abs(r);
}


float get_conditional_entropy(int r1_mcs, int r1_sgi, int r1_ht40, int r2_mcs, int r2_sgi, int r2_ht40 , char* path)
{


  vector<int>  X,Y;
  float PDF_X [NUM_LEVELS]={0};
  float PDF_Y [NUM_LEVELS]={0};













  float PDF_XY[NUM_LEVELS][NUM_LEVELS]={0};
  vector< vector<int> > XY;

  char file1[200];
  char file2[200];

  sprintf(file1, "%s/%d-%d-%d.dat", path, r1_mcs, r1_sgi, r1_ht40 );
  sprintf(file2, "%s/%d-%d-%d.dat", path, r2_mcs, r2_sgi, r2_ht40 );

  //Reading the data files and initialize the X, Y, and XY arrays
  ifstream r1, r2;

  r1.open(file1);
  r2.open(file2);

  read_file(r1, X);
  read_file(r2, Y);

  if (X.size() < 2 || Y.size() < 2 )
	  	  return 0;

  calculate_pdf(X,PDF_X);
  calculate_pdf(Y,PDF_Y);
  calculate_pdf_joint(X, Y, PDF_XY);

  return calculate_conditional_entropy(PDF_X, PDF_XY);

}

float get_entropy (int r1_mcs, int r1_sgi, int r1_ht40, char* path)
{

	  vector<int>  X;
	  float PDF_X [NUM_LEVELS]={0};

	  char file1[200];


	  sprintf(file1, "%s/%d-%d-%d.dat", path, r1_mcs, r1_sgi, r1_ht40 );

	  //Reading the data files and initialize the X, Y, and XY arrays
	  ifstream r1, r2;

	  r1.open(file1);

	  read_file(r1, X);

	  calculate_pdf(X,PDF_X);

	  return calculate_entropy (PDF_X);

}


void print_joint_data_heat (int r1_mcs, int r1_sgi, int r1_ht40, int r2_mcs, int r2_sgi, int r2_ht40 , char* path)
{

	vector<int>  X,Y;
	char file1[200];
	char file2[200];
	const int num_levels = 101;
	int heat [num_levels][num_levels]= {0};

	sprintf(file1, "%s/%d-%d-%d.dat", path, r1_mcs, r1_sgi, r1_ht40 );
	sprintf(file2, "%s/%d-%d-%d.dat", path, r2_mcs, r2_sgi, r2_ht40 );

	//Reading the data files and initialize the X, Y, and XY arrays
	ifstream r1, r2;

	r1.open(file1);
	r2.open(file2);

	read_file(r1, X);
	read_file(r2, Y);



	for (int i=0; i< X.size(); i++){
		for (int j=0; j<Y.size(); j++){
			heat[X[i]][Y[j]]++;
		}

	}

	int num_pairs = X.size() * Y.size();
	for (int i=0; i< num_levels; i++){
		for (int j=0; j< num_levels; j++){
			cout << 100.0 * heat[i][j] / (float)num_pairs << " ";
		}
		cout << endl;
	}


}


void print_joint_data (int r1_mcs, int r1_sgi, int r1_ht40, int r2_mcs, int r2_sgi, int r2_ht40 , char* path)
{

	vector<int>  X,Y;
	char file1[200];
	char file2[200];

	sprintf(file1, "%s/%d-%d-%d.dat", path, r1_mcs, r1_sgi, r1_ht40 );
	sprintf(file2, "%s/%d-%d-%d.dat", path, r2_mcs, r2_sgi, r2_ht40 );

	//Reading the data files and initialize the X, Y, and XY arrays
	ifstream r1, r2;

	r1.open(file1);
	r2.open(file2);

	read_file(r1, X);
	read_file(r2, Y);

	//cout << "H[X] = " << get_entropy(r1_mcs,r1_sgi,r1_ht40, path)<< endl;
	//cout << "H[Y] = " << get_entropy(r2_mcs,r2_sgi,r2_ht40, path)<< endl;

	for (int i=0; i<X.size() && i<Y.size(); i++){
		cout << X[i] << "\t" << Y[i]<< endl;
	}


}

int main(int argc, char* argv[])
{

  if (argc <= NUM_ARGS) {
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr,
//	    "\entropy [R1 file] [R2 file] \n");
	    "entropy [path to data files] [MODE] [OUTPUT_MODE] [MCS_1] [SGI_1] [Ht40_1] [MCS_2] [SGI_2] [Ht40_2] ( [start_time] [end_time] ) \n");
        exit(EXIT_FAILURE);
  }


  char* path = argv[1];
  char* mode = argv[2];
  char* output_mode = argv[3];
  int mcs_1 = atoi(argv[4]);
  int sgi_1 = atoi(argv[5]);
  int ht40_1 = atoi(argv[6]);
  int mcs_2 = atoi(argv[7]);
  int sgi_2 = atoi(argv[8]);
  int ht40_2 = atoi(argv[9]);
  int start_time = atoi(argv[10]);
  int end_time = atoi(argv[11]);




 // int rate = atoi(argv[2]);
  //int ht_40 = 1;


  if (strcmp(mode,"joint") == 0){
	  print_joint_data(mcs_1, sgi_1, ht40_1, mcs_2, sgi_2, ht40_2, path);
  }
  else if (strcmp(mode,"conditional-entropy") == 0){

	  // Enable this to print the conditional entropy of all combinations
	  if (strcmp(output_mode,"all") == 0){
		  for (int r1_ht_40=0; r1_ht_40<2; r1_ht_40++){
			  for (int r1_mcs=0; r1_mcs<NUM_RATES; r1_mcs++){
				  for (int r1_sgi=0; r1_sgi<2; r1_sgi++){
					  for(int r2_ht_40=0; r2_ht_40<2; r2_ht_40 ++){
						  for (int r2_mcs=0; r2_mcs<NUM_RATES; r2_mcs++){
							  for (int r2_sgi=0; r2_sgi<2; r2_sgi++){
								  if (r1_sgi != sgi_1 || r2_sgi != sgi_2 || r1_ht_40 != ht40_1 || r2_ht_40 != ht40_2){
									  continue;
								  }
								  cout << get_conditional_entropy(r1_mcs, r1_sgi, r1_ht_40, r2_mcs, r2_sgi, r2_ht_40 , path ) << "\t";
							  }
						  }
					  }
					  cout << endl;
				  }
			  }
		  }
	  }
	  else{
		  for (int r1_mcs=0; r1_mcs<NUM_RATES; r1_mcs++){
			  for (int r2_mcs=0; r2_mcs<NUM_RATES; r2_mcs++){
				  cout << get_conditional_entropy(r1_mcs, sgi_1, ht40_1, r2_mcs, sgi_2, ht40_2 , path ) << "\t";
			  }
			  cout << endl;
		  }

	  }
  }

  else if (strcmp(mode,"correlation-coefficient") == 0){
	  // Enable this to print the conditional entropy of all combinations
	  if (strcmp(output_mode,"all") == 0){
		  for (int r1_ht_40=0; r1_ht_40<2; r1_ht_40++){
			  for (int r1_mcs=0; r1_mcs<NUM_RATES; r1_mcs++){
				  for (int r1_sgi=0; r1_sgi<2; r1_sgi++){
					  for(int r2_ht_40=0; r2_ht_40<2; r2_ht_40 ++){
						  for (int r2_mcs=0; r2_mcs<NUM_RATES; r2_mcs++){
							  for (int r2_sgi=0; r2_sgi<2; r2_sgi++){
								  cout << get_correlation_coeffcient(r1_mcs, r1_sgi, r1_ht_40, r2_mcs, r2_sgi, r2_ht_40 , path ) << "\t";
							  }
						  }
					  }
					  cout << endl;
				  }
			  }
		  }
	  }
	  else{
		  for (int r1_mcs=0; r1_mcs<NUM_RATES; r1_mcs++){
			  for (int r2_mcs=0; r2_mcs<NUM_RATES; r2_mcs++){
				  cout << get_correlation_coeffcient(r1_mcs, sgi_1, ht40_1, r2_mcs, sgi_2, ht40_2 , path ) << "\t";
			  }
			  cout << endl;
		  }
	  }

  }
  else if (strcmp(mode,"dispersion") == 0){
	  //Enable to output dispersion of all combinations of rates (specified by the input)
	  if (strcmp(output_mode,"all") == 0){
		  for (int r1_mcs=0; r1_mcs<NUM_RATES; r1_mcs++){
			  for (int r2_mcs=0; r2_mcs<NUM_RATES; r2_mcs++){
				  cout << get_max_dispersion(r1_mcs, sgi_1, ht40_1, r2_mcs, sgi_2, ht40_2 , path , start_time, end_time) << "\t";
			  }
			  cout << endl;
		  }
	  }
	  else{

		  for (int r1_mcs=0; r1_mcs<NUM_RATES; r1_mcs++){
			  cout <<  r1_mcs << "\t" << get_max_dispersion(r1_mcs, sgi_1, ht40_1, r1_mcs, sgi_2, ht40_2, path, start_time, end_time) << "\t"
					  << get_max_dispersion(r1_mcs, sgi_2, ht40_1, r1_mcs, sgi_1, ht40_2, path, start_time, end_time )  << endl;
		  }
	  }



  }
  else{
	  cout << "Unkonwn MODE" << endl;
	  return EXIT_FAILURE;
  }






  cout.precision(3);
  cout << fixed;

/*
  cout <<"#The following output shows the conditional entropy of two rates that are different only in the GI setting" <<endl;
  cout << "#X = LGI, Y = SGI" <<endl;
  cout << "#MCS\tH(X)\tH(Y)\tH(Y|X)\tH(X|Y)\tH(X)\tH(Y)\tH(Y|X)\tH(X|Y)"<<endl;
  for (int mcs=0; mcs<NUM_RATES; mcs++){
	  cout << mcs << "\t";
	  for (int ht_40=0; ht_40<2; ht_40++){
		  cout <<get_entropy(mcs,0,ht_40, path) << "\t" << get_entropy(mcs,1,ht_40, path)<< "\t"
				  <<get_conditional_entropy(mcs, 0, ht_40, mcs, 1, ht_40, path) << "\t"
				  <<get_conditional_entropy(mcs, 1, ht_40, mcs, 0, ht_40, path) << "\t";

	  }
	  cout << endl;

  }
*/

 /* cout <<endl<<"#The following output shows the conditional entropy of two rates that are different only in the Channel Width setting" <<endl;
  cout << "#X = 20, Y = 40" <<endl;
  cout << "#MCS\tH(X)\tH(Y)\tH(Y|X)\tH(X|Y)\tH(X)\tH(Y)\tH(Y|X)\tH(X|Y)"<<endl;
  for (int mcs=0; mcs<NUM_RATES; mcs++){
	  cout << mcs << "\t";
	  for (int sgi=0; sgi<2; sgi++){
		  cout <<get_entropy(mcs,sgi, 0, path) << "\t" << get_entropy(mcs,sgi,1, path)<< "\t"
				  <<get_conditional_entropy(mcs, sgi, 0, mcs, sgi, 1, path) << "\t"
				  <<get_conditional_entropy(mcs, sgi, 1, mcs, sgi, 0, path) << "\t";
	  }
	  cout << endl;

  }*/


  //#printing the entropy of all rates
  /*for (int ht_40=0; ht_40<2; ht_40++){
	  for (int sgi=0; sgi<2; sgi++){
		  for (int mcs=0 ; mcs < NUM_RATES ; mcs++){
		  cout << mcs << " " << sgi << " " << ht_40 << "\t" <<  get_entropy(mcs, sgi, ht_40, path) << endl;
		  }
	  }
  }*/


  /* for (int i=0; i<NUM_LEVELS; i++){
    cout << PDF_X[i]<< "\t" << PDF_Y[i]<< endl;
  }

  cout.precision(2);
  cout << fixed;
  //  cout << "H(X) = "  << calculate_entropy(PDF_X)<< endl;
  //  cout << "H(Y) = "  << calculate_entropy(PDF_Y)<< endl;
  //  cout <<"H(Y|X) = " << calculate_conditional_entropy(PDF_X, PDF_XY) << endl;
   cout << calculate_conditional_entropy(PDF_X, PDF_XY) << endl;
  */

   // print_joint_data(19,0,0,19,1,0, path);
  //print_joint_data_heat(rate,0,0,rate,1,0, path);
    return EXIT_SUCCESS;
}
