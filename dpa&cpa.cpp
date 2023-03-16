#include <iostream>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <chrono>
#include <math.h>
#include <vector>

using namespace std;
#define KEY_SIZE 16

#define NB_SAMPLES 3000
#define NB_TRACES 100 // max == 400

const char *path_trace  = "/home/bayo/CLionProjects/sca_tme5/traces.raw";
const char *path_key    = "/home/bayo/CLionProjects/sca_tme5/key.raw";
const char *path_textin = "/home/bayo/CLionProjects/sca_tme5/textin.raw";

double traces[NB_SAMPLES * NB_TRACES];
double keys[KEY_SIZE * NB_TRACES];
double textin[KEY_SIZE * NB_TRACES];

typedef struct inter {
    int x;
    int y;
} ;

int intermediate(int pt, int keyguess);
void DPA_attack(vector<double> trace, vector<int> textin, int nb_trace);
int attack_bytes(vector<double> tra, vector<int> text, int nb_traces, int bnum, inter intervalle);
void cpa_attack(vector<double> tr, vector<int> text, int nb_traces, inter * intervals);
double pearson(double *hypo, double *tr);

int avg_point[16] = {1841, 2057, 2273, 2489, 1889, 2105, 2305, 2537,
                     1937, 2153, 2369, 2585, 2009, 2225, 2441, 2620};

inter intervalles[KEY_SIZE] = { {1836, 2057}, {2052, 2062} , { 2268, 2278 },
                                {2484, 2494}, {1884, 1894}, {2100, 2110},
                                {2300, 2310},{2532, 2542}, {1932, 1942},
                                {2148, 2158}, {2364, 2374}, {2580, 2590},
                                {2004, 2014}, {2220, 2230},
                                {2436, 2446}, {2615, 2625}};

int sbox[256] = {0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
                0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
                0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
                0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
                0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
                0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
                0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
                0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
                0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
                0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
                0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
                0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
                0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
                0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
                0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
                0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
};

int HW[256] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3,
               3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3,
               4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4,
               3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3,
               4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6,
               5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3,
               4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
               4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2,
               3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
               4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6,
               7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};


int main(int argc, char const *argv[]) {
    int fd1, fd2, fd3;
    vector<double> tr;
    vector<int> tex;
    vector<int> ke;

    if ((fd1 = open(path_trace, O_RDONLY)) == -1) {
        perror("open");
        return -1;
    }

    if (read(fd1, traces, sizeof(double) * NB_SAMPLES * NB_TRACES) == -1) {
        perror("read");
        return -1;
    }

    for (int k = 0; k < NB_TRACES; ++k) {
        for (int i = 0; i < NB_SAMPLES; i += 1) {
            tr.push_back(traces[i + (NB_SAMPLES * k)]);
        }
    }

    if ((fd2 = open(path_key, O_RDONLY)) == -1) {
        perror("open");
        return -1;
    }

    if (read(fd2, keys, sizeof(double) * NB_TRACES * 16) == -1) {
        perror("read");
        return -1;
    }


    if ((fd3 = open(path_textin, O_RDONLY)) == -1) {
        perror("open");
        return -1;
    }

    if (read(fd3, textin, sizeof(double) * NB_TRACES * KEY_SIZE) == -1) {
        perror("read");
        return -1;
    }

    for (int k = 0; k < NB_TRACES; ++k) {
        for (int i = 0; i < KEY_SIZE; i += 1) {
            tex.push_back(textin[i + (KEY_SIZE * k)]);
        }
    }

    printf("valeur exacte de la clé : \n");
    printf(" 43 126 21 22 40 174 159 166 171 247 21 136 9 207 79 60\n");
    auto start = std::chrono::high_resolution_clock::now();
    DPA_attack(tr, tex, NB_TRACES);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout << "DPA attack " << float_ms.count() << " milliseconds" << std::endl;
    printf("################################\n");
    printf("valeur avec l'attaque cpa\n");
    auto start2 = std::chrono::high_resolution_clock::now();
    cpa_attack(tr, tex, NB_TRACES, intervalles);
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms2 = end2 - start2;
    std::cout << "CPA attack " << float_ms2.count() << " milliseconds" << std::endl;

    return EXIT_SUCCESS;
}


int intermediate(int pt, int keyguess){
    return sbox[pt ^ keyguess];
}

double pearson(double *hypo, double *tr){
    double hypo_mean = 0.0;
    double mean_tr = 0.0;
    for (int i = 0; i < NB_TRACES; ++i) {
        hypo_mean += hypo[i];
        mean_tr += tr[i];
    }
    hypo_mean /= NB_TRACES;
    mean_tr /= NB_TRACES;

    double cov_h = 0.0;
    double cov_t = 0.0;
    double cov = 0.0;
    for (int i = 0; i < NB_TRACES; ++i) {
        cov +=  (hypo[i] - hypo_mean) * (tr[i] - mean_tr);
        cov_h += (hypo[i] - hypo_mean) * (hypo[i] - hypo_mean);
        cov_t += (tr[i] - mean_tr) * (tr[i] - mean_tr);
    }
    return cov / sqrt(cov_h * cov_t);
}


int attack_bytes(vector<double> tra, vector<int> text, int nb_traces, int bnum, inter intervalle){
    double best_coef = 0.0;
    int best_key = 0;
    for (int guess = 0; guess < 256; ++guess) {
        for (int k =  intervalle.x; k < intervalle.y; ++k) {
            double hypo[nb_traces];
            double tr[nb_traces];
            for (int i = 0; i < nb_traces; ++i) {
                hypo[i] = HW[intermediate(text[bnum + i*KEY_SIZE ], guess)];
                tr[i] = tra[k + i*NB_SAMPLES];
            }
            double coef = fabs(pearson(hypo, tr));
            if( best_coef < coef){
                best_coef = coef;
                best_key = guess;
            }
        }
    }
    return best_key;
}


void cpa_attack(vector<double> tr, vector<int> text, int nb_traces, inter * intervals){
    int recovered_key[KEY_SIZE];
    for (int bnum = 0; bnum < KEY_SIZE; ++bnum) {
        recovered_key[bnum] = attack_bytes(tr, text , nb_traces, bnum, intervals[bnum]);
    }
    for (int i = 0; i < KEY_SIZE; ++i) {
        cout << recovered_key[i] << " ";
    }
    cout << endl;
}

void DPA_attack(vector<double> trace, vector<int> textin, int nb_trace){
    int result[16];
    for (int bnum = 0; bnum < 16; ++bnum) {
        double Max = -1.0;
        int key_found = -1;

        for (int i = 0; i < 256; ++i) {
            vector<double>  bit_one;
            vector<double>  bit_zero;

            for (int j = 0; j < nb_trace; ++j) {
                int sbox_out = HW[intermediate(textin[bnum + j*KEY_SIZE], i)];
                if (sbox_out > 4){
                    bit_one.push_back(trace[avg_point[bnum] + j*NB_SAMPLES]);
                } else {
                    bit_zero.push_back(trace[avg_point[bnum] + j*NB_SAMPLES] );
                }
            }

            double mean_one = 0.0;
            for (int j = 0; j < bit_one.size(); ++j) {
                mean_one += bit_one[j];
            }
            mean_one = mean_one / bit_one.size();
            double mean_zero = 0.0;
            for (int j = 0; j < bit_zero.size(); ++j) {
                mean_zero += bit_zero[j];
            }
            mean_zero = mean_zero / bit_zero.size();
            double diff = fabs(abs(mean_one) - abs(mean_zero));

            if (diff > Max){
                Max = diff;
                key_found = i;
            }
        }
        result[bnum] = key_found;
    }

    for (int i = 0; i < KEY_SIZE; ++i) {
        cout << result[i] << " ";
    }
    cout << endl;
}