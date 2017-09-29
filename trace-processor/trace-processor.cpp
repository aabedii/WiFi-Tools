/*
 * Copyright 2015 Ali Abedi, Andrew Heard, University of Waterloo.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <boost/math/distributions/students_t.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <ctime>
#include <thread>

#include <map>


using ::boost::accumulators::accumulator_set;
using ::boost::accumulators::features;
using ::boost::accumulators::stats;
using ::boost::math::students_t;
using ::std::string;
using ::std::abs;
using ::std::cout;
using ::std::cerr;
using ::std::endl;
using ::std::ifstream;
using ::std::vector;
using ::std::setw;
using ::std::setprecision;
using ::std::fixed;
using ::std::map;

#define PACKET_LENGTH 1536.0
#define CONFIDENCE_LEVEL 0.95
#define NUM_ARGS 9
//#define DOWNSAMPLE_WINDOW 0.15  // The time gap betweeen downsampled samples
#define ROUND_DURATION 0.350       // A complete round-robin round time
#define AVERAGE_WINDOW 1 // used to find correlation between two rates


// Pseudo random string for implementing Randomized multiple interleaved trials
// This works for 10 or 20 trials case
const int rand_stream[40] = {1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0,
                             1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
                             0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1};
/*
0	0	0	1	0	0	0	0	1
0	1	0	0	1	1	1	1	1
0	0
1	0	0	1	1	0	1	1	1
0	0	0	0	1	1	1	0	1
0	0
1	0	0	0	1	1	1	1	1
0	1	1	0	1	0	0	1	1
0	0
1	0	1	0	1	0	0	1	0
1	0	1	1	0	1	0	1	1
1	1
1	1	0	0	1	0	1	1	1
1	0	1	0	0	1	0	1	1
1	0
 */

struct line {
  float m;  // slope
  float b;
  float R_squared;
};

struct average_result {
  int rssi;
  float error;
};

struct sample_metric_measured {
  string readable_time;
  string time;
  double metric_measure;
};

struct sample_rssi_received {
  string readable_time;
  double time;
  int mcs;
  int sgi;
  int ht40;
  int acked;
  int rssi1;
  int rssi2;
  int rssi3;
  int rssi_effective;
  int length;
  int seqno;
  string type;
  string src_mac;
  string dst_mac;
};

struct sample_aggr_sent {
  string readable_time;
  double time;
  int mcs;
  int sgi;
  int ht40;
  int acked;
  int total_sent;
  int total_failed;
  int ack_rssi;

  sample_aggr_sent(sample_aggr_sent* s)
     {
	  readable_time = s->readable_time;
	   time = s->time;
	   mcs = s->mcs;
	   sgi = s->sgi;
	   ht40 = s->ht40;
	   acked = s->acked;
	   total_sent = s->total_sent;
	   total_failed = s->total_failed;
	   ack_rssi = s->ack_rssi;
     }
  sample_aggr_sent (){}
};

struct output_throughput_value {
  double time;
  double throughput;
  double confidence_interval;
};

map<double,sample_aggr_sent*> cache;
const double epsilon = 0.000001;

void print_packet(sample_aggr_sent packet) {
  cout << packet.time << "\t" << packet.mcs << "\t" << packet.total_sent << "\t"
       << packet.total_failed << endl;
}

/*
 void print_data()
 {
 cout.precision(10);
 cout << "time \t mcs \t sgi \t ht40 \t acked \t ack rssi "<<endl;
 for (int i=0; i< trace_sent.size(); i++)
 {
 print_packet(trace_sent[i]);
 }

 }
 */

double confidence_interval(
    const accumulator_set<double,
                          features<boost::accumulators::tag::mean,
                                   boost::accumulators::tag::lazy_variance,
                                   boost::accumulators::tag::count>> &stats,
    double confidence_level) {
  double probability = 1.0 - confidence_level;

  double count = boost::accumulators::extract::count(stats);
  if (count > 1) {
    students_t t_distribution(count - 1);

    double t_constant =
        boost::math::quantile(complement(t_distribution, probability / 2));
    double variance = boost::accumulators::extract::variance(stats);
    double stdev = sqrt(variance);

    return t_constant * stdev / sqrt(count);
  } else {
    return 0;
  }
}

bool trace_get_metric(ifstream *input, struct sample_metric_measured *sm) {
  string junk, metric;
  if (!input->eof()) {
    // input >> junk >> sm.readable_time >> junk;
    // input >> sm.metric_measure;
    // input >> junk >> junk >> junk >> junk >> junk ;
    // input >> sm.time;
    // while (input.peek() != '\n' && !input.eof()) {
    //     input >> sm.time;
    // }
    *input >> sm->readable_time >> sm->metric_measure >> sm->time;

    return true;
  }

  return false;
}

bool trace_get_packet(ifstream *input, struct sample_aggr_sent *sp,
                      double *start_time) {
  string junk, type, time_s;

  if (!input->eof()) {
    *input >> junk >> junk;
    *input >> sp->readable_time;
    *input >> junk >> junk;
    *input >> time_s;

    if (time_s.length() < 2) {
      *input >> time_s;
    }

    // removing [ and ] from the string
    for (int i = 0; i < time_s.length(); i++) {
      if (time_s[i] == ']' || time_s[i] == '[') {
        // cout << time_s[i] << endl;
        time_s.erase(i, 1);
        i--;
      }
    }

    *input >> type;

    if (type.compare("[AGGR]") == 0) {
      sp->time = atof(time_s.c_str()) - *start_time;
      if (*start_time <= 0) {
        *start_time = sp->time;
        sp->time = 0;
      }

      *input >> junk;  // HT rate
      *input >> sp->mcs;
      *input >> sp->sgi;
      *input >> sp->ht40;
      *input >> junk;  // RTS
      *input >> junk;  // rates[0].count
      *input >> sp->total_failed;
      *input >> sp->total_sent;

      if (input->peek() != '\n') {
        *input >> sp->acked;
        *input >> sp->ack_rssi;
      }

      while (input->peek() != '\n') {
        *input >> junk;
      }

      return true;
    } else if (type.compare("[TRACE]") == 0) {
      // Do stuff for round robin traces
      sp->time = atof(time_s.c_str()) - *start_time;
      if (*start_time <= 0) {
        *start_time = sp->time;
        sp->time = 0;
      }
      *input >> junk;  // HT rate
      *input >> sp->mcs;
      *input >> sp->acked;
      *input >> sp->sgi;
      *input >> sp->ht40;
      *input >> junk;  // RTS
      *input >> junk;  // Retries
      *input >> junk;  // rates[1] index
      *input >> junk;  // rates[1] retries
      *input >> junk;  // TX
      *input >> junk;  // RX
      *input >> junk;  // Busy
      *input >> junk;  // Cycles
      *input >> junk;  // Ack Signal
      *input >> junk;  // Sequence number
      *input >> junk;  // Source MAC address
      *input >> junk;  // Destination MAC address

      sp->total_sent = 1;

      if (sp->acked) {
        sp->total_failed = 0;
      } else {
        sp->total_failed = 1;
      }

      return true;
    }
  }

  return false;
}

void add_to_cache(struct sample_aggr_sent *sp){
	sample_aggr_sent *temp = new sample_aggr_sent(sp);
	cache[sp->time] = temp;
}

sample_aggr_sent* is_in_cache(double time){
	//is there a packet with grather (not equal) time in cache?
	sample_aggr_sent* result = NULL;
	if (cache.empty()){
		return result;
	}
	std::map<double,sample_aggr_sent*>::iterator ithigh;
	ithigh=cache.upper_bound (time);
	if (ithigh != cache.end() ){
		result = ithigh->second;
	}
	return result;
}
void clear_cache (double time){
	//delete all elements with a time less than (not equal) the specified time
	std::map<double,sample_aggr_sent*>::iterator itlow;
	itlow=cache.lower_bound (time);
	cache.erase(cache.begin(),itlow);        // erases [itlow,itup)
}


bool trace_get_packet(ifstream *input, struct sample_aggr_sent *sp,
		double *start_time, double *current_computation_time) {
	string junk, type, time_s;

	//check cache, if a packet with higher time than current_computation_time
	//exists return that rather than reading a new packet from the trace
	//This helps to implement the sliding window mode
	sample_aggr_sent* cached_packet = is_in_cache(*current_computation_time + epsilon);
	if (cached_packet != NULL){
		*sp = new sample_aggr_sent(cached_packet);
		//cout << endl << "Time before " << *current_computation_time ;
		*current_computation_time = sp->time;
		//cout << "\t Time after " << *current_computation_time << endl;
		return true;
	}


	*input >> junk >> junk;
	*input >> sp->readable_time;
	*input >> junk >> junk;
	*input >> time_s;

	if (input->eof()) {
		return false;
	}

	if (time_s.length() < 2) {
		*input >> time_s;
	}

	// removing [ and ] from the string
	for (int i = 0; i < time_s.length(); i++) {
		if (time_s[i] == ']' || time_s[i] == '[') {
			// cout << time_s[i] << endl;
			time_s.erase(i, 1);
			i--;
		}
	}

	*input >> type;

	if (type.compare("[AGGR]") == 0) {
		sp->time = atof(time_s.c_str()) - *start_time;
		*current_computation_time = sp->time;
		if (*start_time <= 0) {
			*start_time = sp->time;
			sp->time = 0;
			*current_computation_time = 0;
		}

		*input >> junk;  // HT rate
		*input >> sp->mcs;
		*input >> sp->sgi;
		*input >> sp->ht40;
		*input >> junk;  // RTS
		*input >> junk;  // rates[0].count
		*input >> sp->total_failed;
		*input >> sp->total_sent;

		if (input->peek() != '\n') {
			*input >> sp->acked;
			*input >> sp->ack_rssi;
		}

		while (input->peek() != '\n') {
			*input >> junk;
		}

		//Add the newly read packet to the cache for future reuse
		*current_computation_time = sp->time;
		add_to_cache(sp);


		return true;
	} else if (type.compare("[TRACE]") == 0) {
		// Do stuff for round robin traces
		sp->time = atof(time_s.c_str()) - *start_time;
		if (*start_time <= 0) {
			*start_time = sp->time;
			sp->time = 0;
		}
		*input >> junk;  // HT rate
		*input >> sp->mcs;
		*input >> sp->acked;
		*input >> sp->sgi;
		*input >> sp->ht40;
		*input >> junk;  // RTS
		*input >> junk;  // Retries
		*input >> junk;  // rates[1] index
		*input >> junk;  // rates[1] retries
		*input >> junk;  // TX
		*input >> junk;  // RX
		*input >> junk;  // Busy
		*input >> junk;  // Cycles
		*input >> junk;  // Ack Signal
		*input >> junk;  // Sequence number
		*input >> junk;  // Source MAC address
		*input >> junk;  // Destination MAC address

		sp->total_sent = 1;

		if (sp->acked) {
			sp->total_failed = 0;
		} else {
			sp->total_failed = 1;
		}

		return true;
	}


}

bool trace_get_rssi_log(ifstream *input, struct sample_rssi_received *rp,
                        double *start_time) {
  string junk, type, time_s;
  int protocol;

  while (!input->eof()) {
    *input >> junk >> junk;
    *input >> rp->readable_time;
    *input >> junk >> junk;
    *input >> time_s;
    if (time_s.length() < 2) {
      *input >> time_s;
    }

    // removing [ and ] from the string
    for (int i = 0; i < time_s.length(); i++) {
      if (time_s[i] == ']' || time_s[i] == '[') {
        // cout << time_s[i] << endl;
        time_s.erase(i, 1);
        i--;
      }
    }

    *input >> type >> protocol;

    if (type.compare("[TRACE]") == 0 ||
        type.find("[AGGR]") != string::npos) {  // SENT PACKET
      getline(*input, junk);
      continue;
    } else {  // RECEIVED PACKET
      rp->time = atof(time_s.c_str()) - *start_time;
      if (*start_time <= 0) {
        *start_time = rp->time;
        rp->time = 0;
      }

      *input >> rp->mcs >> rp->sgi >> rp->ht40;
      *input >> rp->rssi1 >> rp->rssi2 >> rp->rssi3;
      *input >> junk;  // RSSI 4
      *input >> rp->rssi_effective;
      *input >> rp->length;
      *input >> rp->type;
      *input >> rp->seqno;
      *input >> rp->src_mac;
      *input >> rp->dst_mac;

      return true;
    }
  }
  return false;
}

bool trace_get_rssi_tcpdump(ifstream *input, struct sample_rssi_received *rp,
                            double *start_time) {
  string junk, type;

  if (!input->eof()) {
    *input >> junk >> junk >> junk;
    *input >> rp->readable_time >> rp->time >> rp->rssi_effective;
    rp->readable_time =
        rp->readable_time.substr(0, rp->readable_time.find(".", 0));

    // cout <<  rp.readable_time << endl;
    if (*start_time <= 0) {
      *start_time = rp->time;
      rp->time = 0;
    }

    return true;
  }

  return false;
}

bool trace_get_rssi(ifstream *input, struct sample_rssi_received *rp,
                    double *start_time, bool use_tcpdump) {
  if (use_tcpdump) {
    return trace_get_rssi_tcpdump(input, rp, start_time);
  } else {
    return trace_get_rssi_log(input, rp, start_time);
  }
}

void print_stats(
    const vector<string> &readable_time_1,
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>> *stats_1,
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>> *error_1,
    const vector<string> &readable_time_2,
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>> *stats_2,
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>> *error_2,
    const vector<string> &readable_time_3,
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>> *stats_3,
    vector<accumulator_set<
        double, features<boost::accumulators::tag::mean,
                         boost::accumulators::tag::lazy_variance,
                         boost::accumulators::tag::count>>> *error_3,
    bool relative_timestamps) {
  bool aabb = stats_2 != NULL;
  bool has_error_rate = (error_1 != NULL);
  bool has_random =
      (stats_3 != NULL);  // has_random can be removed once all functions
  // implement randomized

  vector<double> mean_1;
  vector<double> ci_1;
  vector<double> mean_2;
  vector<double> ci_2;
  vector<double> mean_3;
  vector<double> ci_3;

  // Extract 1_A
  mean_1.push_back(boost::accumulators::extract::mean(stats_1->front()));
  ci_1.push_back(confidence_interval(stats_1->front(), CONFIDENCE_LEVEL));

  if (aabb) {
    // Extract 2_A
    mean_2.push_back(boost::accumulators::extract::mean(stats_2->front()));
    ci_2.push_back(confidence_interval(stats_2->front(), CONFIDENCE_LEVEL));
    if (has_random) {
      mean_3.push_back(boost::accumulators::extract::mean(stats_3->front()));
      ci_3.push_back(confidence_interval(stats_3->front(), CONFIDENCE_LEVEL));
    }

    // Extract 1 and 2 B..X
    for (int i = 1; i < stats_1->size(); i++) {
      mean_1.push_back(boost::accumulators::extract::mean(stats_1->at(i)));
      ci_1.push_back(confidence_interval(stats_1->at(i), CONFIDENCE_LEVEL));
      mean_2.push_back(boost::accumulators::extract::mean(stats_2->at(i)));
      ci_2.push_back(confidence_interval(stats_2->at(i), CONFIDENCE_LEVEL));
      if (has_random) {
        mean_3.push_back(boost::accumulators::extract::mean(stats_3->at(i)));
        ci_3.push_back(confidence_interval(stats_3->at(i), CONFIDENCE_LEVEL));
      }
    }
  }

  vector<double> mean_error_1;
  vector<double> ci_error_1;
  vector<double> mean_error_2;
  vector<double> ci_error_2;
  vector<double> mean_error_3;
  vector<double> ci_error_3;

  if (has_error_rate) {
    // Extract 1_A
    mean_error_1.push_back(
        boost::accumulators::extract::mean(error_1->front()));
    ci_error_1.push_back(
        confidence_interval(error_1->front(), CONFIDENCE_LEVEL));

    if (aabb) {
      // Extract 2_A
      mean_error_2.push_back(
          boost::accumulators::extract::mean(error_2->front()));
      ci_error_2.push_back(
          confidence_interval(error_2->front(), CONFIDENCE_LEVEL));

      // Extract 1 and 2 B..X
      for (int i = 1; i < error_1->size(); i++) {
        mean_error_1.push_back(
            boost::accumulators::extract::mean(error_1->at(i)));
        ci_error_1.push_back(
            confidence_interval(error_1->at(i), CONFIDENCE_LEVEL));
        mean_error_2.push_back(
            boost::accumulators::extract::mean(error_2->at(i)));
        ci_error_2.push_back(
            confidence_interval(error_2->at(i), CONFIDENCE_LEVEL));
        if (has_random) {
          mean_error_3.push_back(
              boost::accumulators::extract::mean(error_3->at(i)));
          ci_error_3.push_back(
              confidence_interval(error_3->at(i), CONFIDENCE_LEVEL));
        }
      }
    }
  }

  // Print 1_A
  cout << setw(10) << readable_time_1.front() << " " << setw(10)
       << setprecision(4) << fixed << mean_1.front() << " " << setw(10)
       << ci_1.front();

  if (has_error_rate) {
    cout << " " << setw(10) << mean_error_1.front() << " " << setw(10)
         << ci_error_1.front();

  } else {
    cout << " " << setw(10) << "nan"
         << " " << setw(10) << 0;
  }

  if (aabb) {
    // Print 1_B .. 1_X
    for (int i = 1; i < stats_1->size(); i++) {
      cout << setw(10) << readable_time_1.at(i) << " " << setw(10)
           << mean_1.at(i) << " " << setw(10) << ci_1.at(i);
      if (has_error_rate) {
        cout << " " << setw(10) << mean_error_1.at(i) << " " << setw(10)
             << ci_error_1.at(i);
      } else {
        cout << " " << setw(10) << "nan"
             << " " << setw(10) << 0;
      }
    }

    // Print 2_A .. 2_X
    for (int i = 0; i < stats_2->size(); i++) {
      cout << setw(10) << readable_time_2.at(i) << " " << setw(10)
           << mean_2.at(i) << " " << setw(10) << ci_2.at(i);
      if (has_error_rate) {
        cout << " " << setw(10) << mean_error_2.at(i) << " " << setw(10)
             << ci_error_2.at(i);
      } else {
        cout << " " << setw(10) << "nan"
             << " " << setw(10) << 0;
      }
    }

    // Print 3_A .. 3_X
    if (has_random) {
      for (int i = 0; i < stats_3->size(); i++) {
        cout << setw(10) << readable_time_3.at(i) << " " << setw(10)
             << mean_3.at(i) << " " << setw(10) << ci_3.at(i);
        if (has_error_rate) {
          cout << " " << setw(10) << mean_error_3.at(i) << " " << setw(10)
               << ci_error_3.at(i);
        } else {
          cout << " " << setw(10) << "nan"
               << " " << setw(10) << 0;
        }
      }
    }
  }

  cout << endl;
}

float get_mean(vector<float> x) {
  float sum = 0;
  for (int i = 0; i < x.size(); i++) {
    sum += x[i];
  }
  return sum / static_cast<float>(x.size());
}

float get_standard_deviation(vector<float> X) {
  float S_x = 0;
  float mean_X = get_mean(X);
  for (int i = 0; i < X.size(); i++) {
    S_x += pow(X[i] - mean_X, 2);
  }

  S_x *= 1.0 / static_cast<float>(X.size() - 1);
  S_x = sqrt(S_x);

  return S_x;
}

// For the case of a linear model with a single independent variable,
// the coefficient of determination (R squared) is the square of r,
// Pearson's product-moment coefficient.
// Pearson's product-moment coefficient: returns r_xy
float get_correlation_coeffcient(vector<float> X, vector<float> Y) {
  float mean_X = get_mean(X);
  float mean_Y = get_mean(Y);

  float S_x = get_standard_deviation(X);
  float S_y = get_standard_deviation(Y);

  // These conditions should be checked
  if (S_x == 0 && S_y == 0) {
    return 1;
  }

  if (S_x == 0 || S_y == 0) {
    return 0;
  }

  float r = 0;

  for (int i = 0; i < X.size(); i++) {
    r += ((X[i] - mean_X) / S_x) * ((Y[i] - mean_Y) / S_y);
  }

  r *= 1.0 / (static_cast<float>(X.size()) - 1);

  return r;
}

line fit(vector<float> X, vector<float> Y) {
  // y = m.x +b
  struct line fitted_line;

  // m = r_xy (S_y/S_x)
  float r_xy = get_correlation_coeffcient(X, Y);
  fitted_line.m =
      r_xy * (get_standard_deviation(Y) / get_standard_deviation(X));

  // b = mean_y - m * mean_x
  fitted_line.b = get_mean(Y) - fitted_line.m * get_mean(X);
  fitted_line.R_squared = pow(r_xy, 2);

  return fitted_line;
}

float linearInterpolation(float x1, float y1, float x2, float y2, float x){
	//y-y1 = m (x-x1)
	//m = y2-y1/x2-x1

	return ((y2-y1)/(x2-x1))* (x-x1) + y1;

}

vector<float> findFunction(vector<float> X, vector<float> Y){

	const float dispersion_limit = 0.4;
	const float dispersion_window = 0.01;
	const float step = 0.01;


	vector<float> result;
/*
	if (X.size() != Y.size()){
		return result;
	}
 */

	float sum;
	int count;
	float min, max;

	for (float e1=0; e1<=1; e1+=step){
		sum = count = 0;
		min = 1;
		max = 0;
		for (int i = 0; i < X.size() && Y.size() ; i++){
			if (e1-(dispersion_window) <= Y[i] && Y[i] <= e1+(dispersion_window)){
				if (X[i] < min){
					min = X[i];
				}
				if (X[i] > max){
					max = X[i];
				}
				count++;
				sum += X[i];
			}
		}

		if (max - min < dispersion_limit){

			if (count > 0){
				result.push_back(sum/count);
			}
			else{
				result.push_back(-1);
			}
		}
		else{
			result.push_back(-2);
		}
	}

	// error rate of 0 or 1 are not reliable for estimation as they don't indiacte
	// how good or bad the channel can be beyond the 0 or 1 limit.

	//result[0] = -3;
	//result[result.size()-1] = -3;


	/*
	if (result[0] == -1){
		for (int i=1; i< result.size(); i++){
			if (result[i] >= 0 ){
				result[0] = result[i];
				break;
			}
		}
	}

	if (result[result.size()-1] == -1){
		for (int i=result[result.size()-2]; i>0; i--){
			if (result[i] >= 0 ){
				result[result.size()-1] = result[i];
				break;
			}
		}
	}


	for (int i=0; i< result.size(); i++){
		if (result[i] == -1){
			for (int j=i+1; j< result.size(); j++){
				if (result[j] >= 0){
					result[i] = linearInterpolation ( i-1, result[i-1], j, result[j], i);
					break;
				}
			}
		}
	}
*/
	return result;
}


// finds the packet with the nearst time to input "time"
int find_index(vector<sample_aggr_sent> data, double time) {
	if (time < data[0].time ){
		return 0;
	}

  if ( time > data[data.size() - 1].time){
	  return data.size() - 1;
  }
  const float epsilon = 0.001;
  int start = 0;
  int end = static_cast<int>(data.size() - 1);
  int index = end / 2;

  while (true) {
    if (abs(data[index].time - time) <= epsilon) break;

    if ((end - start) <= 2) {
      if (abs(data[start].time - time) < abs(data[end].time - time)) {
        index = start;
      } else {
        index = end;
      }
      break;
    }

    if (data[index].time > time) {
      end = (start + end) / 2;
    } else {
      start = (start + end) / 2;
    }
    index = (start + end) / 2;
  }  // while

  return index;
}

// |..window/2..| ...offset... (time) ...offset... | ..window/2..|
average_result get_estimation(vector<sample_aggr_sent> data, double time,
                           double window, double offset, bool distro_based) {

  average_result result;
  int rssi = 0;
  int rssi_count = 0;
  int sent = 0;
  int failed = 0;

  double start_time = time + offset;
  int time_idx = find_index(data, start_time);

  const int NUM_BINS = 101;
  float bin [NUM_BINS] = {0};
  int i;
  float error_rate;

  for (i = time_idx; (abs(data[i].time - start_time) < (window / 2)) &&
                             (i < data.size()) && (time_idx != -1);
       i++) {

      error_rate = static_cast<float>(data[i].total_failed) / static_cast<float>(data[i].total_sent);
      bin [static_cast<int> (error_rate * (NUM_BINS - 1))]+= data[i].total_failed;
      failed += data[i].total_failed;
      sent+= data[i].total_sent;

      if (data[i].ack_rssi != -128) {  //  -128 indicates no ack received
      rssi += data[i].ack_rssi;
      rssi_count++;
    }
  }

 //cout << time << "\t" << sent << "\t" << data[i-1].time << endl;

  start_time = time - offset;
  time_idx = find_index(data, start_time);

  for (i = time_idx - 1; (abs(data[i].time - start_time) < (window / 2)) &&
                                 (i >= 0) && (time_idx != -1);
       i--) {

    error_rate = static_cast<float>(data[i].total_failed) / static_cast<float>(data[i].total_sent);
    bin [static_cast<int>(error_rate * (NUM_BINS - 1))]+= data[i].total_failed;
    failed += data[i].total_failed;
    sent+= data[i].total_sent;


    if (data[i].ack_rssi != -128) {  //  -128 indicates no ack received
      rssi += data[i].ack_rssi;
      rssi_count++;
    }
  }
  //cout << time << "\t" << sent << "\t" << data[i+1].time << endl << endl;
  /*
  for (int i = 0; i < NUM_BINS; i++) {
    bin[i] = bin[i] / static_cast<float>(sent);
    }*/

  if (sent > 0) {
	  if (distro_based){
		  int random = rand() % sent + 1;

		  int index = 0;
		  int threshold = 0;
		  do{
			  threshold += bin[index];
			  index++;
		  }while (threshold < random);

		  result.error = (index - 1) / 100.0;
	  }
	  else{
		  result.error = static_cast<float>(failed) / static_cast<float>(sent);
	  }
  } else {
	  result.error = -1;
  }

  if (rssi_count == 0) {
    result.rssi = 0;
  } else {
    result.rssi = rssi / rssi_count;
  }

  if ( result.error < -1 )
	  cout << time << "\t" << result.error << endl;

  return result;
}
/*
// |..window/2..| ...offset... (time) ...offset... | ..window/2..|
average_result get_average(vector<sample_aggr_sent> data, double time,
                           double window, double offset) {
  average_result result;
  int failed = 0;
  int sent = 0;
  int rssi = 0;
  int rssi_count = 0;

  double start_time = time + offset;
  int time_idx = find_index(data, start_time);

  for (int i = time_idx; (i < data.size()) && (abs(data[i].time - start_time) < (window / 2)) && (time_idx != -1); i++) {
	  failed += data[i].total_failed;
	  sent += data[i].total_sent;
	  if (data[i].ack_rssi != -128) {  //  -128 indicates no ack received
		  rssi += data[i].ack_rssi;
		  rssi_count++;
    }
  }
  start_time = time - offset;
  time_idx = find_index(data, start_time);

  for (int i = time_idx - 1; (i >= 0) && (abs(data[i].time - start_time) < (window / 2)) && (time_idx != -1);i--) {
    failed += data[i].total_failed;
    sent += data[i].total_sent;
    if (data[i].ack_rssi != -128) {  //  -128 indicates no ack received
      rssi += data[i].ack_rssi;
      rssi_count++;
    }
  }


  if (sent > 0) {
    result.error = static_cast<float>(failed) / static_cast<float>(sent);
  } else {
    result.error = -1;
  }

  if (rssi_count == 0) {
    result.rssi = 0;
  } else {
    result.rssi = rssi / rssi_count;
  }

  return result;
}
*/

/*average_result get_average (vector<sample_aggr_sent> data, double time, double
window){

        average_result result;
        int failed = 0;
        int sent = 0;
        int rssi = 0;
        int rssi_count = 0;

        int time_idx = find_index(data, time);

        for (int i= time_idx ; abs(data[i].time - time) < (window/2) &&
i<data.size(); i++ ){
                failed += data[i].total_failed;
                sent += data[i].total_sent;
                if (data[i].ack_rssi != -128){ //  -128 indicates no ack
received
                        rssi += data[i].ack_rssi;
                        rssi_count ++;
                }
        }

        for (int i= time_idx-1 ; abs(data[i].time - time) < (window/2) && i>=0;
i-- ){
                        failed += data[i].total_failed;
                        sent += data[i].total_sent;
                        if (data[i].ack_rssi != -128){ //  -128 indicates no ack
received
                                rssi += data[i].ack_rssi;
                                rssi_count ++;
                        }
        }

        if (sent > 0){
                result.error = (float)failed/(float)sent;
        }
        else{
                result.error = -1;
        }

        if (rssi_count == 0){
                result.rssi = 0;
        }
        else{
                result.rssi = rssi / rssi_count;
        }

        return result;

}*/

void read_trace(char *input_files[], int num_files,
                vector<sample_aggr_sent> *trace) {
  int current_file = 0;
  double start_time = 0;

  ifstream input;
  input.open(input_files[current_file]);

  if (!input.is_open()) {
    cout << "Unable to open " << input_files[0] << endl;
    exit(EXIT_FAILURE);
  }

  // Reading the input file(s) data into a vector

  while (true) {
    sample_aggr_sent sample;

    bool got_packet = trace_get_packet(&input, &sample, &start_time);

    if (!got_packet) {  // no more packets in the current file
      input.close();
      if (++current_file < num_files) {
        input.close();
        input.open(input_files[current_file]);
        got_packet = trace_get_packet(&input, &sample, &start_time);
      }
    }

    if (!got_packet) {  // no more packet in any file
      break;
    }

    trace->push_back(sample);
  }  // while (true)
  input.close();
}

void coherence(char *input_files[], int num_files) {
  vector<sample_aggr_sent> original_trace;
  read_trace(input_files, num_files, &original_trace);
  // Reading file(s) finished files. All the data is now in original_trace

  if (original_trace.size() == 0) {
    cout << "No data read" << endl;
    exit(EXIT_FAILURE);
  }

  // We now calculate how the average error rate diviates as we move farther
  // from an arbitary center time

  float current_time = 0;
  const int levels = 60;
  double sum[levels] = {0};
  int count[levels] = {0};
  double E0, error;

  while (current_time < original_trace[original_trace.size() - 1].time) {
    current_time += 1;  // sliding window steps

    E0 = get_estimation(original_trace, current_time, 0.1, 0, 0).error;

    for (int i = 0; i < levels; i++) {
      error = get_estimation(original_trace, current_time, 0.1, 0.05 * i, 0).error;

      // cout << i<< "\t" << error << "\t" << E0 << endl;
      if (error == -1) {
        continue;
      }
      sum[i] += abs(error - E0);
      count[i]++;
    }
    // cin >> x;
  }  // while

  for (int i = 0; i < levels; i++) {
    cout << i * 0.05 << "\t" << sum[i] / static_cast<double>(count[i]) << endl;
  }
}

void variation(char *input_files[], int num_files, float window) {
  vector<sample_aggr_sent> original_trace;
  read_trace(input_files, num_files, &original_trace);
  // Reading file(s) finished files. All the data is now in original_trace

  if (original_trace.size() == 0) {
    cout << "No data read" << endl;
    exit(EXIT_FAILURE);
  }

  // We now calculate how the max variation over "window".

  const float sub_window = 0.030;

  double current_time = original_trace[0].time;
  double error;
  double max = 0;
  double min = 1;
  while (current_time <
         original_trace[original_trace.size() - 1].time - window) {
    current_time += 1;  // sliding window steps

    for (double i = current_time; i < (current_time + window);
         i += sub_window) {
      error = get_estimation(original_trace, i, sub_window, 0, 0).error;

      if (error == -1) {
        continue;
      }

      if (error < min) {
        min = error;
      }

      if (error > max) {
        max = error;
      }
    }

    if (min > max) {
      continue;
    }

    cout << current_time << "\t" << max - min << endl;

    max = 0;
    min = 1;
  }  // while
}


// Input round-robin trace
// Output number of rates used in the round-robin trace
int get_num_rates(vector<sample_aggr_sent> data)
{
	// There are zero round-robin rates in en empty trace!
	if (data.size() == 0 )
	{
		return 0;
	}

	// It sould be enough to check the few first seconds
	// Check the first 5 seconds

	int i=0;
	vector<sample_aggr_sent> rates;

	while ( (data[i].time - data[0].time < 5) && i <data.size() )
	{
		bool found = false;
		for (int j=0; j< rates.size(); j++ )
		{
			if (data[i].mcs == rates[j].mcs && data[i].sgi == rates[j].sgi && data[i].ht40 == rates[j].ht40)
			{
				found = true;
				break;
			}
		}
		if (!found)
		{
			rates.push_back(data[i]);
		}

		i++;

	}

	return rates.size();
}

void print_time(string s){
	  time_t rawtime;
	  struct tm * timeinfo;
	  char buffer[80];

	  time (&rawtime);
	  timeinfo = localtime(&rawtime);

	  strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
	  std::string str(buffer);

	  std::cout << "# " <<  s << " : " << str << endl;


}


void find_correlation(const vector<vector<sample_aggr_sent>> &trace , int conf, vector< vector<float> > *correlations, vector<line> *correlations_line){
	// Separating/averaging downsampled to X (sgi) and Y (lgi) vectors for
	// finding the correlation.
	vector<float> X;  // conf - base
	vector<float> Y;  // other rates

	float sum = 0;
	int count = 0;
	double last_timestamp = trace[conf][0].time;

	for (int i = 0; i < trace[conf].size(); i++) {
		if ( (trace[conf][i].time - last_timestamp) > AVERAGE_WINDOW )  {
			if (count == 0) {
				X.push_back(0);
				cout << "Missing an packet at " << trace[conf][i].time
						<< endl;
			} else {
				X.push_back(sum / static_cast<float>(count));
			}

			sum = count = 0;
			last_timestamp = trace[conf][i].time;
		}  // if

		count += trace[conf][i].total_sent;
		sum += trace[conf][i].total_failed;


	}  // for

	// Now let's compute Y for each config and compare it agains X
	for (int j=0; j<trace.size(); j++)
	{
		Y.clear();
		sum = count = 0;
		last_timestamp = trace[j][0].time ;

		for (int i = 0; i < trace[j].size(); i++) {
			if ((trace[j][i].time - last_timestamp) > AVERAGE_WINDOW) {
				if (count == 0) {
					Y.push_back(0);
					cout << "Missing an packet at " << trace[j][i].time
							<< endl;
				} else {
					Y.push_back(sum / static_cast<float>(count));
				}

				sum = count = 0;
				last_timestamp = trace[j][i].time;
			}  // if

			count += trace[j][i].total_sent;
			sum += trace[j][i].total_failed;

		}  // for

		// Y computation is done
		// fit a line to X and Y
		 struct line fitted_line = fit (Y, X);
		 //cout << j << "\t" << fitted_line.R_squared << "\t" << fitted_line.m << " x + " << fitted_line.b << endl;
		 correlations_line->push_back(fitted_line); // it might not be a good fitted line should be checked latter if can be used
		 correlations->push_back(findFunction(X,Y));

/*
	 	if (j == 1 && conf == 0){
			for (int k=0; k< correlations[correlations.size()-1].size(); k++){
				cout << k << "\t" << correlations[correlations.size()-1][k] << endl;
			}
	 	}
*/
		 //cout << conf << " and " << j << endl;
		//cout << fitted_line.m <<" x + "<< fitted_line.b << endl;
		//cout << "Coefficient of determination: "<< fitted_line.R_squared << endl;
		 //cout << conf << "\t" << j << endl;

	} //for



}

int round (float num){
	int result = floor(num + 0.5);
	if (result < 0 || result > 100)
		cout << result << endl;
	return result;
}

void generate_cor(vector<sample_aggr_sent> *generated, const vector<sample_aggr_sent> &original,
		const vector<vector<sample_aggr_sent>> &downsampled, int conf, double W_gen, double start_time, double end_time){

	vector < vector<float> > correlations;
	vector < line > correlations_line;

	find_correlation(downsampled, conf, &correlations, &correlations_line);
	average_result estimated;
	float sum;
	int count;
	double current_time;
	/*
	for (int i=0; i< correlations[0].size(); i++){
		cout << i << "\t" << correlations[0][i] << "\t" << correlations[1][i] << "\t" << correlations[2][i] << "\t" << correlations[3][i] << "\t" <<endl;
	}*/
	for (int j=0; j<original.size(); j++){
		if ( original[j].time < start_time || original[j].time > end_time){
				continue;
			}

		current_time = original[j].time;
		sum=count=0;
		double window;

		for (int i=0; i<downsampled.size(); i++){

			estimated.error = -1;
			window = W_gen;

			while (estimated.error < 0) {
				estimated = get_estimation(downsampled[i], current_time, window, 0, 0);
				if (window > 8 * W_gen){
					cout << conf << "\t" <<current_time << endl;
				}
				window = window * 2;
			}

/*
			if (estimated.error < 0){
				continue;
			}
*/
			//estimated = get_estimation(downsampled[i], current_time, W_gen, 0, 0);
			/*if (conf==1){
						cout << "raw error "<< i  << " = " << estimated.error << endl;
			}*/

			float predicted_error;
			if (conf != i){
				predicted_error = correlations_line[i].m * estimated.error + correlations_line[i].b;
			}
			else{ // don't estimate for conf
				predicted_error = estimated.error;
			}

/*
			if (estimated.error == 1 && predicted_error < 0.95)
				continue;

			if (estimated.error == 0 && predicted_error > 0.05)
				continue;
*/

			if (estimated.error == 1 && predicted_error != 1)
				predicted_error = (predicted_error + 1) / 2.0;

			if (estimated.error == 0 && predicted_error != 0)
				predicted_error = (predicted_error) / 2.0;


			if (predicted_error < 0)
				predicted_error = 0;
			else if (predicted_error > 1)
				predicted_error = 1;


			//predicted_error = 0;
			//sum += correlations[i][round(estimated.error*100)];
			sum += predicted_error;
			count++;
		}

		sample_aggr_sent s;
		float error;

		if (count != 0){
			error = sum / count;
		}
		else{ // is this good to do?
			error = 0;
		}

		/*
		if (conf==1){
			cout << error << endl;
		}*/
		s.time = current_time;
		s.total_failed = error * 100;
		s.total_sent = 100;

		generated->push_back(s);
	}


}

void generate(vector<sample_aggr_sent> *generated, const vector<sample_aggr_sent> &original,
		const vector<sample_aggr_sent> &downsampled, double W_gen, double start_time, double end_time){
	double current_time;
	average_result estimated;

	double window;

	for (int j=0; j<original.size(); j++){
		if ( original[j].time < start_time || original[j].time > end_time){
			continue;
		}

		current_time = original[j].time;
		estimated.error = -1;
		//get_estimation(downsampled, current_time, W_gen, 0, 0);
		window = W_gen;
		// increase the window until we find a simple
		while (estimated.error < 0) {
			estimated = get_estimation(downsampled, current_time, window, 0, 0);
			if (window > 8 * W_gen){
				cout << window << "\t" <<current_time << endl;
			}
			window = window * 2;
		}

		sample_aggr_sent s;
		s.time = current_time;
		s.total_failed = estimated.error * 100;
		s.total_sent = 100;

		generated->push_back(s);
	}

}

void compare(vector<float> *results_original, vector<float> *results_estimated,
		vector<float> *results_estimated_cor, const vector<sample_aggr_sent> &original_trace,
		const vector<sample_aggr_sent> &downsampled_trace,const vector<sample_aggr_sent> &downsampled_trace_cor,
		double time_step, double start_time, double end_time, double W_rep){
	double current_time = start_time;
	float ground , gen_vanilla, gen_cor;
	double window;

	while ( current_time < end_time){

		ground = gen_vanilla = gen_cor = -1;

		window = W_rep;
		while (ground < 0){
			ground = get_estimation(original_trace, current_time, window, 0, 0).error;
			window = window * 2;
		}

		window = W_rep;
		while (gen_vanilla < 0){
			gen_vanilla = get_estimation(downsampled_trace , current_time, window, 0, 0).error;
			window = window * 2;
		}

		window = W_rep;
		while (gen_cor < 0){
			gen_cor = get_estimation(downsampled_trace_cor , current_time, window, 0, 0).error;
			window = window * 2;
		}


		assert (gen_vanilla >= 0);
		assert (gen_cor >= 0);
		assert (ground >= 0);

		///// if no data assume the error rate is 0
/*
		if (gen_vanilla == -1){
			gen_vanilla = 0;
		}

		if (gen_cor < 0){
			gen_cor = 0;
		}

		if (ground == -1){
			ground = 0;
		}
*/
		results_original->push_back(ground);
		results_estimated->push_back(gen_vanilla);
		results_estimated_cor->push_back(gen_cor);
		current_time += time_step;
	}
}


void upsample(const vector<vector<sample_aggr_sent>> &trace, int conf, vector<sample_aggr_sent> *generated_trace)
{

	//Find the correlation with other confs and use them when estimating the error rate

	//vector<line> correlations;
	vector < vector<float> > correlations;
	vector < line > correlations_line;
	// find corerlations betwean rate "conf" and all other rates
	find_correlation(trace, conf, &correlations, &correlations_line);

	// for a given rate (conf) we can compute the merged downsampled trace which merges are correlated
	// rates using the linear correlation just found


	vector<long int> cur_idx (trace.size());
	double residue = 0;

	while (true)
	{
		int min_trace = 0;
		bool finished = true;
		for (int i=0; i< trace.size(); i++)
		{
			if (cur_idx[i] < trace[i].size())
			{
				finished = false;
				if (trace[min_trace].size() == cur_idx[min_trace] ||
						trace[i][cur_idx[i]].time < trace[min_trace][cur_idx[min_trace]].time)
				{
					min_trace = i;
				}
			}
		}

		if (finished)
		{
			break;
		}

		sample_aggr_sent temp;
		temp.time = trace[min_trace][cur_idx[min_trace]].time;
		int sent = trace[min_trace][cur_idx[min_trace]].total_sent;
		int failed = trace[min_trace][cur_idx[min_trace]].total_failed;
		float error = static_cast<float>(failed)/sent;

		if (min_trace == conf){
				temp.total_failed = failed;
				temp.total_sent = sent;
				generated_trace->push_back(temp);
				cur_idx[min_trace]++;
				continue;
		}

		//float error_mapped = error * correlations[min_trace].m + correlations[min_trace].b + residue;
		float error_mapped = correlations[min_trace][(int)(error*100)];

		if (error_mapped < 0){ // dispersion is high so not a good correlation to use
			cur_idx[min_trace]++;
			continue;
		}


		if (error_mapped < 0){
			error_mapped = 0;
		}
		else if (error_mapped > 1){
			error_mapped = 1;
		}


		temp.total_failed = 100 * error_mapped + 0.5; // +0.5 to round after casting to int
		temp.total_sent = 100;

//		if (min_trace == conf){
			generated_trace->push_back(temp);
	//	}
		cur_idx[min_trace]++;

	}
}

void downsample_mean_based(char *input_files[], int num_files, float W_rep, float W_gen) {

	const double start_time = 180;
	const double end_time = 250;
	vector<sample_aggr_sent> original_trace;
	vector<vector<sample_aggr_sent>> original_trace_sub; //original per rate
	vector<vector<sample_aggr_sent>> downsampled_trace_sub; //downsampled per rate

	int rates_num; // number of rates found in the collected trace

	read_trace(input_files, num_files, &original_trace);
	// Reading file(s) finished files. All the data is now in original_trace

	if (original_trace.size() == 0) {
		cout << "No data read" << endl;
		exit(EXIT_FAILURE);
	}


	// Creating original sub traces
	bool found;
	for (int i=0; i<original_trace.size(); i++){
		found = false;
		//search if same type packet already has a sub trace
		for (int j=0; j< original_trace_sub.size(); j++ )
		{
			if (original_trace[i].mcs == original_trace_sub[j][0].mcs &&
					original_trace[i].sgi == original_trace_sub[j][0].sgi &&
					original_trace[i].ht40 == original_trace_sub[j][0].ht40)
			{
				original_trace_sub[j].push_back(original_trace[i]);
				found = true;
				break;
			}
		}
		//this type of packet has not seen yet so create a new sub trace for it
		if (!found)
		{
			vector<sample_aggr_sent> temp;
			original_trace_sub.push_back(temp);
			original_trace_sub[original_trace_sub.size()-1].push_back(original_trace[i]);
		}

	}

	// spliting rates into different sub vectors finished

	rates_num = original_trace_sub.size();

	// Creating the downsampled trace
	print_time("Creating the downsampled trace");
	double last_sample_time ;

	for (int conf=0; conf < rates_num; conf++)
	{
		vector<sample_aggr_sent> temp;
		downsampled_trace_sub.push_back(temp);

		//last_sample_time = -INFINITY;
		last_sample_time = original_trace_sub[0][0].time  + conf * (ROUND_DURATION/rates_num);
		cout << "#last_sample_time = " << last_sample_time  << "\t conf = " << conf << endl;

		for (int i = 0; i < original_trace_sub[conf].size(); i++) {
			if ((original_trace_sub[conf][i].time - last_sample_time) > ROUND_DURATION ) {
				downsampled_trace_sub[conf].push_back(original_trace_sub[conf][i]);
				last_sample_time = original_trace_sub[conf][i].time;
			}
		}
		cout << "#DS Conf " << conf << " trace size = " << downsampled_trace_sub[conf].size() << endl;
		cout << "#OG Conf " << conf << " trace size = " << original_trace_sub[conf].size() << endl;
	}

	print_time("Creating the downsampled trace done.");

	vector<std::thread> threads;

	// Generating the final trace using the downsampled traces
	print_time("Generating the final trace using the downsampled traces");
	vector<vector<sample_aggr_sent>> generated_trace_sub;
	vector<vector<sample_aggr_sent>> generated_trace_sub_cor;

	for (int i=0; i< rates_num; i++)
	{
		vector<sample_aggr_sent> temp;
		generated_trace_sub.push_back(temp);
		vector<sample_aggr_sent> temp2;
		generated_trace_sub_cor.push_back(temp2);
	}

	threads.clear();
	double current_time;
	// For every rate used for trace collection do:
	for (int i=0; i< rates_num; i++)
	{
		threads.push_back(std::thread (generate, &generated_trace_sub[i], original_trace_sub[i], downsampled_trace_sub[i], W_gen, start_time, end_time));
		threads.push_back(std::thread(generate_cor, &generated_trace_sub_cor[i], original_trace_sub[i], downsampled_trace_sub, i, W_gen, start_time, end_time));
		//threads.push_back(std::thread(generate_cor, &generated_trace_sub_cor[i], original_trace_sub[i], original_trace_sub, i, W_gen));
	}
	for (int i=0; i< threads.size(); i++){
		threads[i].join();
	}

	print_time("Generating the final trace using the downsampled traces done.");

	// Comparision and printing out the results
	print_time("Comparision");


	const float time_step = 1; // This is the gap between two consecutive comparision points

	float ground, gen_vanilla, gen_cor;
	current_time = start_time;//original_trace[0].time;
	last_sample_time = end_time; //original_trace[original_trace.size()-1].time;

	vector<vector<float>> results_original;
	vector<vector<float>> results_estimated;
	vector<vector<float>> results_estimated_cor;

	for (int i=0; i< downsampled_trace_sub.size(); i++)
	{
		vector<float> temp;
		results_original.push_back(temp);
		vector<float> temp2;
		results_estimated.push_back(temp2);
		vector<float> temp3;
		results_estimated_cor.push_back(temp3);
	}


	threads.clear();


	for (int i=0; i< downsampled_trace_sub.size(); i++){

		threads.push_back(std::thread(compare, &results_original[i],&results_estimated[i], &results_estimated_cor[i],
				original_trace_sub[i], generated_trace_sub[i], generated_trace_sub_cor[i],
				time_step, start_time, end_time , W_rep));

		}
	for (int i=0; i< threads.size(); i++){
		threads[i].join();
	}
	print_time("Comparision done");

	for (int i=0; i< results_original[0].size(); i++){
		cout << i* time_step << "\t";
		for (int j=0; j< results_original.size(); j++){
			cout << results_original[j][i] << "\t" << results_estimated[j][i] << "\t" << results_estimated_cor[j][i] << "\t";
		}
		cout << endl;
	}


	// Calculating and printing out the overal statistics

	cout << "# Time \t";
	for (int i=0; i< downsampled_trace_sub.size(); i++){
		cout << "Original-"<< downsampled_trace_sub[i][0].mcs << "-" << downsampled_trace_sub[i][0].sgi << "-" << downsampled_trace_sub[i][0].ht40 << "\t"
			 << "Gen-vanilla \t Gen-Cor \t";
	}
	cout << endl;

	float actual_error, estimated_error, estimated_cor_error;

	for (float threshold = 0; threshold <= 0.3 ; threshold+=0.1){
		cout << "# " << threshold << " ";
		for (int i=0; i< rates_num; i++){
			double sum = 0;
			double sum_cor = 0;
			int count= 0;

			for (int j=0; j< results_original[0].size(); j++ ){
				if (results_original[i][j] >= threshold){
					sum += abs ( results_original[i][j] - results_estimated[i][j]);
					sum_cor += abs (results_original[i][j] - results_estimated_cor[i][j] );
					count ++;
				}
			}
			cout << sum / static_cast<float>(count) << "\t" << sum_cor / static_cast<float>(count) << "\t";
		}
		cout << endl;
	}
}
void downsample(char *input_files[], int num_files, float W_rep, float W_gen) {
	// Reading the input file(s) data into a vector

	vector<sample_aggr_sent> original_trace;
	vector<vector<sample_aggr_sent>> original_trace_sub; //original per rate
	vector<vector<sample_aggr_sent>> downsampled_trace_sub; //downsampled per rate
	vector<vector<sample_aggr_sent>> cor_upsampled_trace_sub; // downsampled and then upsampled using correlation per rate (not generated yet just the number of samples increased)

	int rates_num; // number of rates found in the collected trace

	read_trace(input_files, num_files, &original_trace);
	// Reading file(s) finished files. All the data is now in original_trace

	if (original_trace.size() == 0) {
		cout << "No data read" << endl;
		exit(EXIT_FAILURE);
	}


	// Creating original sub traces
	bool found;
	for (int i=0; i<original_trace.size(); i++){
		found = false;
		//search if same type packet already has a sub trace
		for (int j=0; j< original_trace_sub.size(); j++ )
		{
			if (original_trace[i].mcs == original_trace_sub[j][0].mcs &&
					original_trace[i].sgi == original_trace_sub[j][0].sgi &&
					original_trace[i].ht40 == original_trace_sub[j][0].ht40)
			{
				original_trace_sub[j].push_back(original_trace[i]);
				found = true;
				break;
			}
		}
		//this type of packet has not seen yet so create a new sub trace for it
		if (!found)
		{
			vector<sample_aggr_sent> temp;
			original_trace_sub.push_back(temp);
			original_trace_sub[original_trace_sub.size()-1].push_back(original_trace[i]);
		}

	}

	// spliting rates into different sub vectors finished

	rates_num = original_trace_sub.size();

	// Creating the downsampled trace
	print_time("Creating the downsampled trace");
	double last_sample_time ;

	for (int conf=0; conf < rates_num; conf++)
	{
		vector<sample_aggr_sent> temp;
		downsampled_trace_sub.push_back(temp);

		last_sample_time = -INFINITY;

		for (int i = 0; i < original_trace_sub[conf].size(); i++) {
			if ((original_trace_sub[conf][i].time - last_sample_time) > ROUND_DURATION ) {
				downsampled_trace_sub[conf].push_back(original_trace_sub[conf][i]);
				last_sample_time = original_trace_sub[conf][i].time;
			}
		}
	}

	// Upsample the downsampled traces using correlation if existed
	print_time("Upsample the downsampled traces using correlation if existed");

	vector<std::thread> threads;

	for (int i=0; i< rates_num; i++)
	{
		vector<sample_aggr_sent> temp;
		cor_upsampled_trace_sub.push_back(temp);
		upsample(downsampled_trace_sub, i, &cor_upsampled_trace_sub[i]);
		//threads.push_back(std::thread (upsample,downsampled_trace_sub, i, &cor_upsampled_trace_sub[i]));
	}
/*
	for (int i=0; i< threads.size(); i++){
		threads[i].join();
	}
*/

	for (int i=0; i< rates_num; i++)
	{
		std::stringstream filename;
		filename << "ds-data/trace-downsampled-" << i << ".log";
		std::ofstream fs(filename.str().c_str());
		if(fs)
		{
			for (int j=0; j< downsampled_trace_sub[i].size(); j++ ){
				fs << downsampled_trace_sub[i][j].time << "\t" <<  downsampled_trace_sub[i][j].total_failed /  static_cast<float>(downsampled_trace_sub[i][j].total_sent) << endl;
			}
		}
		fs.close();

	}


	for (int i=0; i< rates_num; i++)
	{
		std::stringstream filename;
		filename << "ds-data/trace-upsampled-" << i << ".log";
		std::ofstream fs(filename.str().c_str());
		if(fs)
		{
			for (int j=0; j< cor_upsampled_trace_sub[i].size(); j++ ){
				fs << cor_upsampled_trace_sub[i][j].time << "\t" <<  cor_upsampled_trace_sub[i][j].total_failed /  static_cast<float>(cor_upsampled_trace_sub[i][j].total_sent) << endl;
			}
		}
		fs.close();

	}
	//return;



	// Generating the final trace using the downsampled traces
	print_time("Generating the final trace using the downsampled traces");
	vector<vector<sample_aggr_sent>> generated_trace_sub;
	vector<vector<sample_aggr_sent>> generated_trace_sub_cor;

	for (int i=0; i< rates_num; i++)
	{
		vector<sample_aggr_sent> temp;
		generated_trace_sub.push_back(temp);
		vector<sample_aggr_sent> temp2;
		generated_trace_sub_cor.push_back(temp2);
	}

	threads.clear();
	double current_time;
	// For every rate used for trace collection do:
	for (int i=0; i< rates_num; i++)
	{
		//generate(&generated_trace_sub[i], original_trace_sub[i], downsampled_trace_sub[i], W_gen);
		//std::thread t (generate, &generated_trace_sub[i], original_trace_sub[i], downsampled_trace_sub[i], W_gen);
		threads.push_back(std::thread (generate, &generated_trace_sub[i], original_trace_sub[i], downsampled_trace_sub[i], W_gen, 0 , 600));
		//Generating without coordination for comparision purposes
		/*for (int j=0; j<original_trace_sub[i].size(); j++){

			sample_aggr_sent s;
			current_time = original_trace_sub[i][j].time;

			average_result estimated =
					get_estimation(downsampled_trace_sub[i], current_time, W_gen, 0);

			if (estimated.error == -1) {
				continue;
			}

			s.time = current_time;
			s.total_failed = estimated.error * 25;
			s.total_sent = 25;

			generated_trace_sub[i].push_back(s);
		}
		 */

		//generate(&generated_trace_sub_cor[i], original_trace_sub[i], cor_upsampled_trace_sub[i], W_gen);
		//std::thread t_cor (generate, &generated_trace_sub_cor[i], original_trace_sub[i], cor_upsampled_trace_sub[i], W_gen);
		threads.push_back(std::thread(generate, &generated_trace_sub_cor[i], original_trace_sub[i], cor_upsampled_trace_sub[i], W_gen, 0 , 600));

		//Generating using coordination
		/*for (int j=0; j<original_trace_sub[i].size(); j++){

				sample_aggr_sent s;
				current_time = original_trace_sub[i][j].time;

				average_result estimated =
						get_estimation(cor_upsampled_trace_sub[i], current_time, W_gen, 0);

				if (estimated.error == -1) {
					continue;
				}

				s.time = current_time;
				s.total_failed = estimated.error * 25;
				s.total_sent = 25;

				generated_trace_sub_cor[i].push_back(s);

			}*/
	}
	for (int i=0; i< threads.size(); i++){
		threads[i].join();
	}

	// Comparision and printing out the results
	print_time("Comparision and printing out the results");

	cout << "# Time \t";
	for (int i=0; i< downsampled_trace_sub.size(); i++){
		cout << "Original-"<< downsampled_trace_sub[i][0].mcs << "-" << downsampled_trace_sub[i][0].sgi << "-" << downsampled_trace_sub[i][0].ht40 << "\t"
			 << "Gen-vanilla \t Gen-Cor \t";
	}
	cout << endl;

	const float time_step = 1; // This is the gap between two consecutive comparision points
	current_time = original_trace[0].time;
	float ground, gen_vanilla, gen_cor;
	last_sample_time = original_trace[original_trace.size()-1].time;

	vector<vector<float> >results_original;
	vector<vector<float>> results_estimated;
	vector<vector<float>> results_estimated_cor;

	for (int i=0; i< downsampled_trace_sub.size(); i++)
	{
		vector<float> temp;
		results_original.push_back(temp);
		vector<float> temp2;
		results_estimated.push_back(temp2);
		vector<float> temp3;
		results_estimated_cor.push_back(temp3);
	}

/*	while ( current_time < last_sample_time){
		current_time += time_step;
		//cout << current_time << "\t";
		for (int i=0; i< downsampled_trace_sub.size(); i++){
			ground = get_average(original_trace_sub[i], current_time, W_rep, 0).error;
			gen_vanilla = get_average(downsampled_trace_sub[i] , current_time, W_rep, 0).error;
			gen_cor = get_average(cor_upsampled_trace_sub[i] , current_time, W_rep, 0).error;

			///// if no data assume the error rate is 0
			if (gen_vanilla == -1){
				gen_vanilla = 0;
			}

			if (gen_cor == -1){
				gen_cor = 0;
			}

			if (ground == -1){
				ground = 0;
			}
			results_original[i].push_back(ground);
			results_estimated[i].push_back(gen_vanilla);
			results_estimated_cor[i].push_back(gen_cor);

		}
		//cout << endl;
	}
	*/

	threads.clear();

	for (int i=0; i< downsampled_trace_sub.size(); i++){

		threads.push_back(std::thread(compare, &results_original[i],&results_estimated[i], &results_estimated_cor[i],
				original_trace_sub[i], generated_trace_sub[i], generated_trace_sub_cor[i],
				time_step, current_time, last_sample_time , W_rep));

		}
	for (int i=0; i< threads.size(); i++){
		threads[i].join();
	}

	for (int i=0; i< results_original[0].size(); i++){
		cout << i* time_step << "\t";
		for (int j=0; j< results_original.size(); j++){
			cout << results_original[j][i] << "\t" << results_estimated[j][i] << "\t" << results_estimated_cor[j][i] << "\t";
		}
		cout << endl;
	}

	// Calculating and printing out the overal statistics
	float actual_error, estimated_error, estimated_cor_error;

	for (float threshold = 0; threshold <= 0.2 ; threshold+=0.1){
		cout << "# " << threshold << " ";
		for (int i=0; i< rates_num; i++){
			double sum = 0;
			double sum_cor = 0;
			int count= 0;

			for (int j=0; j< results_original[0].size(); j++ ){
				if (results_original[i][j] >= threshold){
					sum += abs ( results_original[i][j] - results_estimated[i][j]);
					sum_cor += abs (results_original[i][j] - results_estimated_cor[i][j] );
					count ++;
				}
			}
			cout << sum / static_cast<float>(count) << "\t" << sum_cor / static_cast<float>(count) << "\t";
		}
		cout << endl;
	}
}

void print_raw_average(char *input_files[], int num_files, float window) {
  int current_file = 0;
  double start_time = 0;

  ifstream input;
  input.open(input_files[current_file]);

  if (input.is_open()) {
  } else {
    cout << "Unable to open " << input_files[0] << endl;
    exit(EXIT_FAILURE);
  }

  // Reading the input file(s) data into a vector

  vector<sample_aggr_sent> original_trace;
  vector<sample_aggr_sent> downsampled_trace;

  while (true) {
    sample_aggr_sent sample;

    bool got_packet = trace_get_packet(&input, &sample, &start_time);

    if (!got_packet) {  // no more packets in the current file
      input.close();
      if (++current_file < num_files) {
        input.close();
        input.open(input_files[current_file]);
        got_packet = trace_get_packet(&input, &sample, &start_time);
      }
    }

    if (!got_packet) {  // no more packet in any file
      break;
    }

    original_trace.push_back(sample);
  }  // while (true)
  input.close();

  // Reading file(s) finished files. All the data is now in original_trace

  if (original_trace.size() == 0) {
    cout << "No data read" << endl;
    exit(EXIT_FAILURE);
  }

  int failed = 0;
  int count = 0;
  int rssi = 0;
  int rssi_count = 0;

  double time = original_trace[0].time;

  for (int i = 0; i < original_trace.size(); i++) {
    if (original_trace[i].time - time < window) {
      failed += original_trace[i].total_failed;
      count += original_trace[i].total_sent;
      if (original_trace[i].ack_rssi !=
          -128) {  // -128 represents no ack received
        rssi_count++;
        rssi += original_trace[i].ack_rssi;
      }
    } else {
      if (count == 0) {
        time = original_trace[i].time;
        continue;
      }
      if (rssi_count == 0) {
        cout << original_trace[i].time << "\t"
             << static_cast<float>(failed) / static_cast<float>(count) << "\t"
             << 0 << endl;
      } else {
        cout << original_trace[i].time << "\t"
             << static_cast<float>(failed) / static_cast<float>(count) << "\t"
             << rssi / rssi_count << endl;
      }
      count = failed = rssi = rssi_count = 0;
      time = original_trace[i].time;
    }
  }
}

void metric_trial_stats_continuous(char *input_files[], int num_files,
                                   int trial_count, int num_comparisons,
                                   bool az_technique) {
  int current_file = 0;

  ifstream input;
  input.open(input_files[current_file]);

  if (input.is_open()) {
    cout << "#\t\tAABB\t\t\t\t\t\tABABt\t\t\t\t\tABAB-Rand" << endl;
    cout << "# Time\tMean_A\tCI_A\tMean_B\tCI_B\tMean_C\tCI_C" << endl;

  } else {
    cout << "Unable to open " << input_files[0] << endl;
    exit(EXIT_FAILURE);
  }

  vector<string> readable_time_1(num_comparisons);
  vector<string> readable_time_2(num_comparisons);
  vector<string> readable_time_3(num_comparisons);

  int experiment_index = 1;

  sample_metric_measured sample_metric;

  bool complete_experiment;

  while (!input.eof()) {
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_2(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_3(num_comparisons);

    int trial_index = 0;

    complete_experiment = true;

    while (trial_index < trial_count * 2) {
      bool got_packet = trace_get_metric(&input, &sample_metric);

      // This is to avoid problematic measurements in their dataset
      if (got_packet) {
        if (sample_metric.metric_measure == 0) {
          continue;
        }
      }

      if (!got_packet) {
        input.close();
        if (++current_file < num_files) {
          input.open(input_files[current_file]);
          got_packet = trace_get_metric(&input, &sample_metric);
        }
      }

      if (!got_packet) {
        complete_experiment = false;
        break;
      }

      // cout << sample_metric.readable_time << "\t" <<
      // sample_metric.metric_measure << endl;

      if (az_technique) {
        stats_1[0](sample_metric.metric_measure);
        readable_time_1[0] = sample_metric.time;
      } else {
        // AABB
        int comparison_index =
            (trial_index % (trial_count * num_comparisons)) / trial_count;
        stats_1[comparison_index](sample_metric.metric_measure);
        readable_time_1[comparison_index] = sample_metric.time;

        // ABAB
        comparison_index = trial_index % num_comparisons;
        stats_2[comparison_index](sample_metric.metric_measure);
        readable_time_2[comparison_index] = sample_metric.time;

        // RANDOMIZED ABAB
        comparison_index = rand_stream[trial_index - 1];
        stats_3[comparison_index](sample_metric.metric_measure);
        readable_time_3[comparison_index] = sample_metric.time;
      }

      trial_index++;
    }  // while ( trial_index < trial_count * 2 )

    if (complete_experiment) {
      if (az_technique) {
        print_stats(readable_time_1, &stats_1, NULL, readable_time_2, NULL,
                    NULL, readable_time_3, NULL, NULL, false);
      } else {
        print_stats(readable_time_1, &stats_1, NULL, readable_time_2, &stats_2,
                    NULL, readable_time_3, &stats_3, NULL, false);
      }
    }
    experiment_index++;
  }  // while (!input.eof())
}

void rssi_trial_stats_discontinuous(char *input_files[], int num_files,
                                    double trial_duration, int trial_count,
                                    int num_comparisons, bool use_tcpdump) {
  ifstream input;
  int sum, count;
  int failed;
  double start_time;
  vector<string> readable_time_1;
  vector<string> readable_time_2;
  vector<string> readable_time_3;

  sample_rssi_received sample_received;
  cout << "#\t\tAABB\t\t\t\t\t\tABAB" << endl;
  cout << "# Time\tMean_A\tCI_A\tMean_B\tCI_B" << endl;

  for (int i = 0; i < num_files - (num_comparisons * trial_count);
       i += trial_count * num_comparisons) {
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_2(num_comparisons);

    for (int trial_index = i; trial_index < i + (num_comparisons * trial_count);
         trial_index++) {
      sum = count = start_time = failed = 0;
      double rssi_avg;

      input.open(input_files[trial_index]);

      if (!input.is_open()) {
        cout << "Unable to open " << input_files[trial_index] << endl;
        exit(EXIT_FAILURE);
      }

      while (
          trace_get_rssi(&input, &sample_received, &start_time, use_tcpdump)) {
        sum += sample_received.rssi_effective;
        count++;
      }

      input.close();

      rssi_avg = sum / static_cast<double>(count);

      // AABB
      int comparison_index =
          (trial_index % (trial_count * num_comparisons)) / trial_count;
      stats_1[comparison_index](rssi_avg);
      readable_time_1[comparison_index] = sample_received.readable_time;

      // ABAB
      comparison_index = trial_index % num_comparisons;
      stats_2[comparison_index](rssi_avg);
      readable_time_2[comparison_index] = sample_received.readable_time;
    }

    print_stats(readable_time_1, &stats_1, NULL, readable_time_2, &stats_2,
                NULL, readable_time_3, NULL, NULL, false);
  }
}

void rssi_trial_stats_continuous(char *input_files[], int num_files,
                                 vector<double> trial_durations,
                                 int trial_count, int num_comparisons,
                                 bool az_technique, bool use_tcpdump) {
  int current_file = 0;

  ifstream input;
  input.open(input_files[current_file]);

  if (input.is_open()) {
    cout << "#\t\tAABB\t\t\t\t\t\tABAB" << endl;
    cout << "# Time\tMean_A\tCI_A\tMean_B\tCI_B" << endl;

  } else {
    cout << "Unable to open " << input_files[0] << endl;
    exit(EXIT_FAILURE);
  }

  double start_time = 0;
  double t_start = 0;

  vector<string> readable_time_1(num_comparisons);
  vector<string> readable_time_2(num_comparisons);
  vector<string> readable_time_3(num_comparisons);

  int experiment_index = 1;

  sample_rssi_received sample_received;

  int experiment_duration = trial_durations[0] * trial_count;
  if (!az_technique) {
    for (int i = 1; i < num_comparisons; i++) {
      experiment_duration += trial_durations[i] * trial_count;
    }
  }

  bool complete_experiment;

  while (!input.eof()) {
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_2(num_comparisons);

    int sum = 0;
    int count = 0;

    double rssi_avg;
    int trial_index = 0;

    complete_experiment = false;

    while (t_start < (experiment_duration * experiment_index)) {
      bool got_packet =
          trace_get_rssi(&input, &sample_received, &start_time, use_tcpdump);

      if (t_start == 0) {
        t_start = sample_received.time;
      }

      if (!got_packet) {
        input.close();
        if (++current_file < num_files) {
          input.open(input_files[current_file]);
          got_packet = trace_get_rssi(&input, &sample_received, &start_time,
                                      use_tcpdump);
        }
      }

      double trial_duration =
          trial_durations[trial_index % trial_durations.size()];

      if ((sample_received.time - t_start) > trial_duration) {
        t_start += trial_duration;

        complete_experiment = true;
        rssi_avg = sum / static_cast<double>(count);

        if (az_technique) {
          stats_1[0](rssi_avg);
          readable_time_1[0] = sample_received.readable_time;
        } else {
          // AABB
          int comparison_index =
              (trial_index % (trial_count * num_comparisons)) / trial_count;
          stats_1[comparison_index](rssi_avg);
          readable_time_1[comparison_index] = sample_received.readable_time;

          // ABAB
          comparison_index = trial_index % num_comparisons;
          stats_2[comparison_index](rssi_avg);
          readable_time_2[comparison_index] = sample_received.readable_time;
        }

        trial_index++;
        sum = count = 0;
      }

      if (got_packet) {
        sum += sample_received.rssi_effective;
        count++;
      } else {
        complete_experiment = false;
        break;
      }
    }

    if (complete_experiment) {
      if (az_technique) {
        print_stats(readable_time_1, &stats_1, NULL, readable_time_2, NULL,
                    NULL, readable_time_3, NULL, NULL, false);
      } else {
        print_stats(readable_time_1, &stats_1, NULL, readable_time_2, &stats_2,
                    NULL, readable_time_3, NULL, NULL, false);
      }
    }
    experiment_index++;
  }
}

void trial_stats_discontinuous(char *input_files[], int num_files,
                               double trial_duration, int trial_count,
                               int num_comparisons) {
  ifstream input;
  int sum;
  int failed;
  double start_time;
  vector<string> readable_time_1(num_comparisons);
  vector<string> readable_time_2(num_comparisons);
  vector<string> readable_time_3(num_comparisons);

  sample_aggr_sent sample_sent;

  for (int i = 0; i < num_files - (num_comparisons * trial_count);
       i += trial_count * num_comparisons) {
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_2(num_comparisons);

    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        error_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        error_2(num_comparisons);

    cout << "#\t\tAABB\t\t\t\t\t\tABAB" << endl;
    cout << "# Time\tMean_A\tCI_A\tMean_B\tCI_B" << endl;

    for (int trial_index = i; trial_index < i + (num_comparisons * trial_count);
         trial_index++) {
      sum = start_time = failed = 0;
      double throughput, error;

      input.open(input_files[trial_index]);

      if (!input.is_open()) {
        cout << "Unable to open " << input_files[trial_index] << endl;
        exit(EXIT_FAILURE);
      }

      while (trace_get_packet(&input, &sample_sent, &start_time)) {
        sum += sample_sent.total_sent - sample_sent.total_failed;
        failed += sample_sent.total_failed;
      }

      input.close();

      throughput = sum * (PACKET_LENGTH * 8) / (trial_duration * 1024 * 1024);
      error = (static_cast<double>(failed)) / (failed + sum);

      // AABB
      int comparison_index =
          (trial_index % (trial_count * num_comparisons)) / trial_count;
      stats_1[comparison_index](throughput);
      error_1[comparison_index](error);
      readable_time_1[comparison_index] = sample_sent.readable_time;

      // ABAB
      comparison_index = trial_index % num_comparisons;
      stats_2[comparison_index](throughput);
      error_2[comparison_index](error);
      readable_time_2[comparison_index] = sample_sent.readable_time;
    }

    print_stats(readable_time_1, &stats_1, &error_1, readable_time_2, &stats_2,
                &error_2, readable_time_3, NULL, NULL, false);
  }
}

void print_rate_stat(int rate[][3][2][2], int64_t total_packets) {
  int sum = 0;

  cout << "MCS \t LGI-20 \t LGI-40 \t SGI-20 \t SGI-40" << endl;
  for (int ss = 0; ss <= 2; ss++) {
    for (int mcs = 0; mcs < 8; mcs++) {
      cout << mcs << "\t";
      for (int sgi = 0; sgi <= 1; sgi++) {
        for (int ht40 = 0; ht40 <= 1; ht40++) {
          sum += rate[mcs][ss][sgi][ht40];
          cout << rate[mcs][ss][sgi][ht40] /
                      static_cast<double>(total_packets) * 100 << "\t";
        }
      }
      cout << endl;
    }
  }

  cout << "Total packets  " << total_packets << "\t Calculate Total packets "
       << sum << endl;
}



void trial_stats_continuous(char *input_files[], int num_files,
                            vector<double> trial_durations, int trial_count,
                            int num_comparisons, bool az_technique,
                            bool relative_timestamps, double sliding_window) {
  int current_file = 0;
  double current_sliding_window_start_time = 0;
  double current_computation_time = 0;


  // rate [MCS][SS][SGI][HT40]
  //int rate[8][3][2][2] = {0};
  int64_t total_packets = 0;

  if (az_technique) {
    num_comparisons = 1;
  }

  ifstream input;
  input.open(input_files[current_file]);

  if (input.is_open()) {
    if (num_comparisons > 1) {
      cout << "#\t\tAABB\t\t\t\t\t\tABAB" << endl;
    }

    for (int i = 0; i < num_comparisons; i++) {
      string comparison_identifier = "";
      comparison_identifier += ('A' + i);
/*
      cout << "#" << setw(9) << "Time_" + comparison_identifier << " ";
      cout << setw(10) << "Mean_" + comparison_identifier << " ";
      cout << setw(10) << "CI_" + comparison_identifier << " ";
      cout << setw(10) << "Error_" + comparison_identifier << " ";
      cout << setw(10) << "CI_" + comparison_identifier;
 */
    }
    // cout << endl;
    // cout << "# Time\tMean_A\tCI_A\tMean_B\tCI_B" << endl;

  } else {
    cout << "Unable to open " << input_files[0] << endl;
    exit(EXIT_FAILURE);
  }

  double start_time = 0;
  double t_start = 0; // trial start

  vector<string> readable_time_1(num_comparisons);
  vector<string> readable_time_2(num_comparisons);
  vector<string> readable_time_3(num_comparisons);

  int experiment_index = 1;

  sample_aggr_sent sample_sent;

  int experiment_duration = trial_durations[0] * trial_count;
  if (!az_technique) {
    for (int i = 1; i < num_comparisons; i++) {
      experiment_duration += trial_durations[i] * trial_count;
    }
  }

  while (!input.eof()) {
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        stats_2(num_comparisons);

    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        error_1(num_comparisons);
    vector<accumulator_set<double,
                           features<boost::accumulators::tag::mean,
                                    boost::accumulators::tag::lazy_variance,
                                    boost::accumulators::tag::count>>>
        error_2(num_comparisons);

    int sum = 0;
    int failed = 0;

    double throughput, error;
    int trial_index = 0;

    bool complete_experiment = true;
    double window_advancement;

    if ( sliding_window <= 0 ){
    	window_advancement = experiment_duration;
    }
    else{
    	window_advancement = sliding_window;
    }

    while (t_start < (experiment_duration +  ((experiment_index -1) * window_advancement)) ) {

    	bool got_packet = trace_get_packet(&input, &sample_sent, &start_time, &current_computation_time);
    	//cout << "Retruning from  trace_get_packet " << sample_sent.time << endl;
    	//cout << sample_sent.time << "\t";

    	if (t_start == 0) {
    		t_start = sample_sent.time;
    	}

    	if (!got_packet) {
    		input.close();
    		if (++current_file < num_files) {
    			input.open(input_files[current_file]);
    			got_packet = trace_get_packet(&input, &sample_sent, &start_time, &current_computation_time);
    		}
    	}

    	double trial_duration =
    			trial_durations[trial_index % trial_durations.size()];

    	if ((sample_sent.time - t_start) > trial_duration) {

    		t_start += trial_duration;


    		throughput = sum * (PACKET_LENGTH * 8) / (trial_duration * 1024 * 1024);

    		if ((failed + sum) == 0) {
    			error = 1;
    		} else {
    			error = (static_cast<double>(failed) / (failed + sum));
    		}
    		// AABB
    		int comparison_index =
    				(trial_index % (trial_count * num_comparisons)) / trial_count;
    		stats_1[comparison_index](throughput);
    		error_1[comparison_index](error);
    		if (relative_timestamps) {
    			readable_time_1[comparison_index] =
    					std::to_string(static_cast<int>(sample_sent.time));
    		} else {
    			readable_time_1[comparison_index] = sample_sent.readable_time;
    		}


    		// ABAB
    		comparison_index = trial_index % num_comparisons;
    		stats_2[comparison_index](throughput);
    		error_2[comparison_index](error);
    		readable_time_2[comparison_index] = sample_sent.readable_time;

    		trial_index++;
    		sum = 0;
    		failed = 0;
    	}

    	if (got_packet) {
    		//rate[sample_sent.mcs % 8][sample_sent.mcs / 8][sample_sent.sgi]
    		//   [sample_sent.ht40] += sample_sent.total_sent;
    		total_packets += sample_sent.total_sent;
    		sum += sample_sent.total_sent - sample_sent.total_failed;
    		failed += sample_sent.total_failed;

    	} else {
    		complete_experiment = false;
    		break;
    	}
    } // end of  while (t_start < (experiment_duration * experiment_index))
    if (sliding_window > 0){
    	current_sliding_window_start_time += sliding_window;
    	current_computation_time = current_sliding_window_start_time;
    	t_start = current_sliding_window_start_time;
    	clear_cache(current_sliding_window_start_time);
    }
  // cout << endl << "#################################### " << current_sliding_window_start_time << endl;

    if (complete_experiment) {
    	if (az_technique) {
    		print_stats(readable_time_1, &stats_1, &error_1, readable_time_2, NULL,
    				NULL, readable_time_3, NULL, NULL, false);
    	} else {
    		print_stats(readable_time_1, &stats_1, &error_1, readable_time_2,
    				&stats_2, &error_2, readable_time_3, NULL, NULL, false);
    	}
    }

    experiment_index++;
  }// while (!input.eof())

  // print_rate_stat( rate, total_packets);
}

void write_output_value(const double duration,
                        const output_throughput_value &output) {
  cout << setprecision(8) << fixed << setw(12) << output.throughput << " ";
  cout << setprecision(8) << fixed << setw(12) << output.confidence_interval
       << " ";
}

void split_aabb(const vector<double> &measurements, int experiment_count,
                vector<double> *measurements_a,
                vector<double> *measurements_b) {
  for (int i = 0; i < measurements.size(); i++) {
    if (i % (experiment_count * 2) < experiment_count) {
      measurements_a->push_back(measurements[i]);
    } else {
      measurements_b->push_back(measurements[i]);
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc <= NUM_ARGS) {
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr,
            "\toutput-processor [data type] [window duration] [trial count]"
            "[continuous] [a..z technique] [num_comparisons] "
            "[use_tcpdump] [relative timestamp?] [trace file 1] [trace file 2] ... "
            "[trace file n] \n");
    exit(EXIT_FAILURE);
  }


  /* initialize random seed: */
  srand (time(NULL));

  char *bin_name = argv[0];
  char *data_type = argv[1];

  vector<double> trial_durations;
  char trial_delimiter = ',';
  char *trial_duration = NULL;
  char *saveptr;
  trial_duration = strtok_r(argv[2], &trial_delimiter, &saveptr);

  while (trial_duration) {
    trial_durations.push_back(atof(trial_duration));
    trial_duration = strtok_r(NULL, &trial_delimiter, &saveptr);
  }

  int trial_count = atoi(argv[3]);
  int continuous = atoi(argv[4]);
  int az_technique = atoi(argv[5]);
  int num_comparisons = atoi(argv[6]);
  int use_tcpdump = atoi(argv[7]);
  bool relative_timestamps = atoi(argv[8]);
  double sliding_window = atoi(argv[9]);
  char **trace_file = &argv[10];

  cout << "# " << bin_name << " " << data_type << " " << trial_durations[0];
  cout << " " << trial_count << " " << continuous << " " << az_technique;
  cout << " " << num_comparisons << " " << use_tcpdump << " " << relative_timestamps;
  cout << " " << sliding_window;

  for (int i = NUM_ARGS + 1; i < argc; i++) {
    cout << " " << argv[i];
  }
  cout << endl;

  if (num_comparisons % trial_durations.size() != 0) {
    cerr << "Warning: number of compared alternatives and number of "
         << "trials do not divide evenly." << endl;
  }

  int num_files = argc - NUM_ARGS;

  if (strcmp(data_type, "throughput") == 0) {
    if (continuous) {
    	trial_stats_continuous(trace_file, num_files, trial_durations,
                             trial_count, num_comparisons, az_technique,
                             relative_timestamps, sliding_window);

    } else {
      trial_stats_discontinuous(trace_file, num_files, trial_durations[0],
                                trial_count, num_comparisons);
    }
  } else if (strcmp(data_type, "rssi") == 0) {
    if (continuous) {
      rssi_trial_stats_continuous(trace_file, num_files, trial_durations,
                                  trial_count, num_comparisons, az_technique,
                                  use_tcpdump);
    } else {
      rssi_trial_stats_discontinuous(trace_file, num_files, trial_durations[0],
                                     trial_count, num_comparisons, use_tcpdump);
    }
  } else if (strcmp(data_type, "EC2") == 0) {
    metric_trial_stats_continuous(trace_file, num_files, trial_count,
                                  num_comparisons, az_technique);
  } else if (strcmp(data_type, "downsample") == 0) {
    float window_average_report = atof(argv[2]);  // results averaging window
    float window_generation = atof(argv[3]);  // Window used to find the estimated error when generating
    cout << "# Downsampling mode ... window_average_report = " << window_average_report << "  window_generation = " << window_generation << endl;
    downsample_mean_based(trace_file, num_files, window_average_report, window_generation);
  } else if (strcmp(data_type, "print") == 0) {
    float window = atof(argv[2]);  // The window over which the average error
                                   // rate should be calculated
    print_raw_average(trace_file, num_files, window);
  } else if (strcmp(data_type, "coherence") == 0) {
    // coherence(trace_file, num_files);
    double window = atof(argv[2]);
    variation(trace_file, num_files, window);
  } else {
    cout << "Data type must be: throughput or rssi or EC2" << endl;
  }

  return EXIT_SUCCESS;
}
