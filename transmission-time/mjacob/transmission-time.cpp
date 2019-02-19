/*
 * transmission-time.cpp
 *
 *  Created on: Jun 25, 2015
 *      Author: ali
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <string>
extern "C" {
#include "rates.h"
#include "rate-configs.h"

}

using namespace std;


/* Data structures */
enum Standard {
    g,
    n
};

//802.11n rate structure
struct rate_80211 {
    Standard standard;
    float rate_short; //rate in Mbps like 6.5 Mbps
    float rate_long; //rate in Mbps like 7.2 Mbps
    int MCS;
    int nSS; //number of streams 1 to 4
    int nDBPS; // number of data bits per symbol
    int HT40; // 1 40 Mhz channel - 0 20 Mhz channel
};


/*Progrm general constants */
#define ARG_NUM 1
#define MEGA 1000000
#define RATETABLE_FILE "/mnt/c/Users/MIDUL JACOB/Documents/Projects/set-cover/tools/WiFi-Tools/transmission-time/80211n_rates"


/* 802.11 general constants */
#define N_NUM_RATES 64
#define G_NUM_RATES 8
#define DATA_LEN  1470 //1560
#define ACK_LEN 14
#define BLOCK_ACK_LEN 20

#define SIFS 0.000010              // 10 micro-seconds
#define slot_short 0.000009        // 9 micro-seconds

/* PLCP constants */
#define L_STF 0.000008
#define L_LTF 0.000008
#define L_SIG 0.000004
#define HT_STF 0.000004
#define HT_LTF 0.000004
#define HT_SIG 0.000008

/* OFDM constants*/
#define service_bits 16
#define tail_bits  6                   //tail bits per spatial stream
#define OFDM_signal_ext 0.000006        // 6 micro-seconds

/* Aggregations constants*/
#define MPDU_Delimiter_bits 32


/* Global variables*/
rate_80211 RATE_TABLE_80211[N_NUM_RATES + G_NUM_RATES];


bool init_ratetable(string fn) {
    ifstream input(fn.c_str());
    int mcs, n_dbps;
    float rate_s, rate_l;
    string dummy;

    if (input.is_open()) {
        for (int i = 0; i < (N_NUM_RATES + G_NUM_RATES); i++) {
            //sample line: 0 BPSK 1/2 1 52 4 52 26 6.5 7.2
            // dbps - total databits per symbol
            // rate_l - rate of LGI
            // rate_s - rate of SGI
            input >> mcs >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
                  >> n_dbps >> rate_l >> rate_s;
            if (i < G_NUM_RATES)
                RATE_TABLE_80211[i].standard = g;
            else
                RATE_TABLE_80211[i].standard = n;

            RATE_TABLE_80211[i].rate_short = rate_s;
            RATE_TABLE_80211[i].rate_long = rate_l;
            RATE_TABLE_80211[i].MCS = mcs;
            RATE_TABLE_80211[i].nSS = mcs / 8 + 1;
            RATE_TABLE_80211[i].nDBPS = n_dbps;
            RATE_TABLE_80211[i].HT40 = (i >= (32 + 8)); //the first 32+8 are
            // 20 Mhz and the second 32+8 are 40 Mhz

        }
    } else {
        cout
                << "Could not initialize the rate table. Set the path to "
                   "802.11n_rate in #define RATETABLE_FILE"
                << endl;
        return false;
    }

    return true;

}

double get_PSDU_duration(int rate_idx, int SGI, int length, int aggr_num) {
    float symbol_length;

    if (SGI == 0)
        symbol_length = 4; // 4 microseconds
    else
        symbol_length = 3.6;

    //service and tail bits
    double PLCP_bits =
            service_bits + RATE_TABLE_80211[rate_idx].nSS * tail_bits;


    if (aggr_num == 1) { //no-Aggregation
        return ceil((length * 8.0 + PLCP_bits) /
                    static_cast<float> (RATE_TABLE_80211[rate_idx].nDBPS)) *
               symbol_length / static_cast<float>(MEGA);
    } else { //frame aggregation enable
        int mpdu_length = length * 8 + MPDU_Delimiter_bits;
        //adding pad. the length should be a multiple of 4 bytes
        mpdu_length = mpdu_length + (mpdu_length % 32);

        return ceil((aggr_num * mpdu_length + PLCP_bits) /
                    static_cast<float> (RATE_TABLE_80211[rate_idx].nDBPS)) *
               symbol_length / static_cast<float> (MEGA);
    }
}

double get_PSDU_data_duration(int rate_idx, int SGI, int length, int aggr_num) {
    float symbol_length;

    if (SGI == 0)
        symbol_length = 4; // 4 microseconds
    else
        symbol_length = 3.6;

    //service and tail bits
    double PLCP_bits = 0.0;
    // service_bits + RATE_TABLE_80211[rate_idx].nSS * tail_bits;


    if (aggr_num == 1) { //no-Aggregation
        return ceil((length * 8.0) /
                    static_cast<float> (RATE_TABLE_80211[rate_idx].nDBPS)) *
               symbol_length / static_cast<float>(MEGA);
    } else { //frame aggregation enable
        int mpdu_length = length * 8;
        //adding pad. the length should be a multiple of 4 bytes
        mpdu_length = mpdu_length + (mpdu_length % 32);

        return ceil((aggr_num * mpdu_length + PLCP_bits) /
                    static_cast<float> (RATE_TABLE_80211[rate_idx].nDBPS)) *
               symbol_length / static_cast<float> (MEGA);
    }
}


double get_PLCP_duration(int rate_idx, int SGI, int length, int aggr_num) {
    //cout << get_PSDU_duration(rate_idx, SGI, length, aggr_num) << endl;
    //non-HT mode (also used for 802.11g)


    if (rate_idx <= 7) {
        return L_STF + L_LTF + L_SIG +
               get_PSDU_duration(rate_idx, SGI, length, aggr_num) +
               OFDM_signal_ext;
    }
        //Mixed mode
    else {
        // cout << "MPDU transmision time = " << get_PSDU_duration(rate_idx,
        // SGI, length, aggr_num) / aggr_num * MEGA << endl;
        //1 stream 1 HT_LTF, 2 streams 2 HT_LTFs, 3 or 4 streams 4 HT_LTFs
        int HT_LTF_sum;
        if (RATE_TABLE_80211[rate_idx].nSS <= 2)
            HT_LTF_sum = RATE_TABLE_80211[rate_idx].nSS * HT_LTF;
        else
            HT_LTF_sum = 4 * HT_LTF;
        return L_STF + L_LTF + L_SIG + HT_SIG + HT_STF + HT_LTF_sum +
               get_PSDU_duration(rate_idx, SGI, length, aggr_num) +
               OFDM_signal_ext;
    }
}

//return the average backoff time
double get_backoff() {
    return 7.5 * slot_short;
}

int get_ack_rate(int rate_idx) {
    if (rate_idx >= 4)
        return 4;
    if (rate_idx >= 2)
        return 2;
    return 0;
}


double get_duration(int rate_idx, int SGI, int length,
                    int aggr_num)//int rate_idx, int rts, int acked)
{
    double DIFS = SIFS + 3 *
                         slot_short; // 37 micro second <--changed from 2 to 3
    // based on QoS docs //  verified: this is a 802.11e
    // feature for "best effort" class

    int ack_idx;
    double ack_time;
    // ACK rate is device dependent
    ack_idx = get_ack_rate(rate_idx);

    double backoff = get_backoff();

    double data_time = get_PLCP_duration(rate_idx, SGI, length, aggr_num);


    if (aggr_num > 1) {
        ack_time = get_PLCP_duration(ack_idx, 0, BLOCK_ACK_LEN,
                                     1); //needs to be checked
    } else {
        ack_time = get_PLCP_duration(ack_idx, 0, ACK_LEN,
                                     1); //needs to be checked
    }

    return DIFS + backoff + data_time + SIFS + ack_time;

}


void benchmark(int rate_idx, int SGI, int length) {
//	const int aggr_pattern_len = 26;
//	int aggr_pattern [2][aggr_pattern_len] =
//	{
//	        {32, -1, 1, 1, 1, 32, -1, 1, 1, 1, 32, 32, -1, 1, 1, 1, 32, -1,
//          1, 1, 1, 32, -1, 1, 1, 1}, {32, 32, 32, 32, -1, 1, 1, 1, 32, 32,
//                                      32, 32, 32, -1, 1, 1, 32, 32, 32, 32,
//                                      32, 32,-1, 1, 1, 1}
//	};


    const int aggr_pattern_len = 7;
    int aggr_pattern[2][aggr_pattern_len] =
            {
                    {32, -1, 1,  32, 32, -1, 1},
                    {32, 32, 32, 32, 32, -1, 1}
            };

    float time = 0;
    int count_packets_sent = 0;
    int inner_counter = 0;

    for (int p = 0; p < 2; p++) {
        while (time < MEGA) {
            int aggr_num = aggr_pattern[p][inner_counter++];

            if (inner_counter == aggr_pattern_len) {
                inner_counter = 0;
            }

            if (aggr_num == -1) {
                aggr_num = rand() % 31 + 1;
            }
            //cout << aggr_num << " ";
            if (aggr_num != 1 ||
                (aggr_num == 1 && aggr_pattern[p][inner_counter] == 1))
                count_packets_sent += aggr_num;
            time += get_duration(rate_idx, SGI, length, aggr_num) * MEGA;
        }
        cout << "Pattern " << p << " : "
             << count_packets_sent * length * 8 / static_cast<float>(MEGA)
             << endl;
        inner_counter = count_packets_sent = time = 0;
    }
}

int compute(int ht_rate, int index, int sgi, int ht40, int num_aggr) {

    if (index != 32) {

        float duration;
        float transmission_time, expected_throughput, data_duration, efficiency;
        if (ht_rate == 0) {
            duration = static_cast<float>(get_duration(index, 0, DATA_LEN, 1)) *
                       MEGA;
            data_duration = static_cast<float>(get_PSDU_data_duration(index, 0,
                                                                      DATA_LEN,
                                                                      1)) *
                            MEGA;
        } else {
            duration = static_cast<float>(get_duration(
                    (index + (ht40 * N_NUM_RATES / 2)) + G_NUM_RATES, sgi,
                    DATA_LEN,
                    num_aggr)) * MEGA;
            data_duration =
                    static_cast<float>(get_PSDU_data_duration(
                            (index + (ht40 * N_NUM_RATES / 2)) + G_NUM_RATES,
                            sgi,
                            DATA_LEN, num_aggr)) * MEGA;
        }
        transmission_time = duration;
        expected_throughput =
                static_cast<float>(MEGA * num_aggr) / duration * DATA_LEN * 8 /
                MEGA;

        efficiency = (data_duration / duration) * 100;

        if (ht_rate == 0) {
            cout << "|" << setw(8)<<"g | " << setw(6) <<  index;
        } else {
            cout << "|" << setw(8)<<"n | " << setw(6) <<  index;
        }
        if (sgi == 0) {
            cout << " |"<< setw(8) << "SGI";
        } else {
            cout << " |" << setw(8) << "LGI";
        }
        if (ht40 == 0) {
            cout << " |"<< setw(14) << " 20 MHz |";
        } else {
            cout << " |"<< setw(14) << " 40 MHz |";
        }

        cout << std::setprecision
                (4) << std::fixed << setw(24) << transmission_time << " |"<< setw(26)
                << expected_throughput << std::setprecision
                (4) << std::fixed << " |" << setw(10) <<  efficiency
             ;
    } else {

        return 0;
    }
}

int main(int argc, char *argv[]) {
    cout << "\033[2J\033[1;1H"; // This command clears the screen

    cout << "Input Parameters are: [NUM_AGGR]\n\n" << endl;

//    cout << "Input Parameters are: " << ARG_NUM
//         << " [IS_HT_RATE] [INDEX] [SGI] [HT_40] [NUM_AGGR]\n\n" << endl;

    if (argc < ARG_NUM) {
//        cout << "Provide " << ARG_NUM
//             << " inputs: [IS_HT_RATE] [INDEX] [SGI] [HT_40] [NUM_AGGR]"
//             << endl;
        cout << "Provide " << ARG_NUM
             << " input: [NUM_AGGR]"
             << endl;
        return -1;
    }

//    int ht_rate = atoi(argv[1]);
//    int index = atoi(argv[2]);
//    int sgi = atoi(argv[3]);
//    int ht40 = atoi(argv[4]);ls

    int num_aggr = atoi(argv[1]);




    cout << "|" << setw(8) << "n/g | " << setw(6) << "INDEX  | " <<
    setw(8) << "SGI/LGI |" << setw(14)
    << "20MHz/40MHz |" << setw(26) << "Transmission time(ms) |" << setw
    (28) << "Expected Throughput(Mbps) |" << setw(12) << std::setprecision
    (4) << std::fixed << "Efficiency |" << setw(16) << "Physical Rate |" <<
    setw(18) << "Rate Configuration |"
    <<
    setw(12) << "Filename |" << endl;
    cout
            << "-----------------------------------------------------"
               "--------------------------------------------------------------"
               "----------------------------------------\n";


    //Consistency checks
//    if (!ht_rate && index >= G_NUM_RATES) {
//        cout << "Inconsistent configuration: index < 8 for 802.11g" << endl;
//        return -1;
//    }
//
//    if (ht_rate && index >= N_NUM_RATES / 2) {
//        cout << "Inconsistent configuration: index < 32 for 802.11n" << endl;
//        return -1;
//    }

//    if (!(sgi != 1 || sgi != 0) || !(ht40 != 1 || ht40 != 0) || !(ht_rate !=
//                                                                  1 ||
//                                                                  ht_rate !=
//                                                                  0)) {
//        cout << "Inconsistent configuration: SGI,HT40, and IS_HT_RATE must be"
//                " 0 or 1" << endl;
//        return -1;
//    }

    if (num_aggr < 1) {
        cout << "Inconsistent configuration: NUM_AGGR must be at least 1"
             << endl;
        return -1;
    }

    //cout << "Initializing the rate table: " << RATETABLE_FILE << endl;
    //Initializing rate table
    if (!init_ratetable(RATETABLE_FILE)) {
        return -1;
    }

    // (int ht_rate, int index, int sgi, int ht40, int num_aggr)
    // 1S-I1-LG-40M=13.5

    int rate_select = 0;
    string rate_string;
    int test;
    char ratestr[80];
    int index = 0;
    int i = 0;
    char filename[80];

    // int test = compute(0, 0, 0, 0, num_aggr);
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= 31; j++) {
            for (int k = 0; k <= 1; k++) {
                for (int l = 0; l <= 1; l++) {
                    test = compute(i, j, k, l, num_aggr);
                    if( j >=0 && j < 8 ) {
                            rate_select = 1 + 0 + j % 8;
                    }
                    if( j >=8 && j < 16 ) {
                        rate_select = 1 + 8 + j % 8;
                    }
                    if( j >=16 && j < 24 ) {
                        rate_select = 1 + 16 + j % 8;
                    }
                    if( l == 1  && j < 24) {
                        rate_select += 48;
                    }
                    if( k == 0 && j < 24 ) {
                        rate_select += 24;
                    }
                    // cout << " |" << setw(14) << "4.4.44";
                    universal_index_to_str(rate_select, ratestr);
                    rate_string = ratestr;
                    cout << " |" << setw(14) << std::setprecision
                            (1) << std::fixed << strtof((rate_string.substr
                            (rate_string.find
                            ("=") + 1)).c_str(),0)
                    << " |" << setw(18)
                    << rate_string;
                    // printf("r = %2d %18s   ", rate_select, ratestr);

                    ratestr_to_filename(ratestr, filename);
                    cout << " |" << setw(10) << filename << " |";
                    // printf("filename = %s", filename);
                    cout << endl;

                }
            }
        }
    }
    cout
            << "-----------------------------------------------------"
               "--------------------------------------------------------------"
               "----------------------------------------\n";


//
//    printf("----------------------------------------------------------------------\n");
//    printf("Testing universal rates strings\n");
//    fflush(stdout);
//    for (i = 1; i <= MAX_RATES; i ++) {
//        universal_index_to_str(i, ratestr);
//        printf("r = %2d %18s   ", i, ratestr);
//
//        ratestr_to_filename(ratestr, filename);
//        printf("filename = %s\n", filename);
//
//        if ((i % 8) == 0) {
//            printf("\n");
//        }
//
//        index = ratestr_to_universal_index(ratestr);
//        assert(index == i);
//        if (index != i) {
//            printf("Index is not as expected index = %d i = %d\n", index, i);
//            exit(1);
//        }
//    }

//
//    benchmark((index + (ht40 * N_NUM_RATES / 2)) + G_NUM_RATES, sgi, DATA_LEN);
//
//    float duration;
//    float transmission_time, expected_throughput, data_duration, efficiency;
//    if (ht_rate == 0) {
//        duration =
//                static_cast<float>(get_duration(index, 0, DATA_LEN, 1)) * MEGA;
//        data_duration =
//                static_cast<float>(get_PSDU_data_duration(index, 0, DATA_LEN,
//                                                          1)) * MEGA;
//    } else {
//        duration = static_cast<float>(get_duration(
//                (index + (ht40 * N_NUM_RATES / 2)) + G_NUM_RATES, sgi, DATA_LEN,
//                num_aggr)) * MEGA;
//        data_duration = static_cast<float>(get_PSDU_data_duration(
//                (index + (ht40 * N_NUM_RATES / 2)) + G_NUM_RATES, sgi, DATA_LEN,
//                num_aggr)) * MEGA;
//    }
//    transmission_time = duration;
//    expected_throughput =
//            static_cast<float>(MEGA * num_aggr) / duration * DATA_LEN * 8 /
//            MEGA;
//
//    efficiency = (data_duration / duration) * 100;
//    cout << "Transmission Duration (ms) = " << transmission_time << endl;
//
//    //cout << "Single packet duration = " << duration << endl;
//    cout << "Expected throughput (Mbps) = " << expected_throughput << endl;
//    cout << "Efficiency = " << efficiency << "%";
//
//    cout
//            << "\n\n| 802.11 n/g \t\t | INDEX \t\t | SGI/LGI \t\t | "
//               "20MHz/40MHz \t\t | Transmission Duration(ms) \t\t"
//               " | Expected Throughput(Mbps) \t\t | Efficiency\n";
//    cout
//            << "-----------------------------------------------------------"
//               "-------------------------------------------------"
//               "-----------------------------------------------------\n";
//    if (ht_rate == 0) {
//        cout << "| g\t\t\t\t\t | " << index;
//    } else {
//        cout << "| n\t\t\t\t\t | " << index;
//    }
//    if (sgi == 0) {
//        cout << "\t\t\t | SGI\t\t ";
//    } else {
//        cout << "\t\t\t | LGI\t\t ";
//    }
//    if (ht40 == 0) {
//        cout << "\t | 20 MHz\t\t";
//    } else {
//        cout << "\t | 40 MHz\t\t";
//    }
//
//    cout << "\t | " << duration << "\t\t\t\t\t\t\t | " << expected_throughput
//         << "\t\t\t\t\t\t\t | " << efficiency;

    if (test == 0) {
        return 0;
    } else {
        return -1;
    }
}


