#ifndef WAVELET_UNIT_H
#define WAVELET_UNIT_H

#include <armadillo>

using namespace std;
using namespace arma;

//analysis half lengths
#define HL_A53_LP 2
#define HL_A53_HP 1
#define HL_A97_LP 4
#define HL_A97_HP 3

//synthesis half lengths
#define HL_S53_LP 1
#define HL_S53_HP 2
#define HL_S97_LP 3
#define HL_S97_HP 4

int analysis_5_3(char* inputFile, char* outputFile, int levels);
void analysis_5_3(my_image_comp *in, my_image_comp *out, int levels);
int analysis_9_7(char* inputFile, char* outputFile, int levels);
void analysis_9_7(my_image_comp *in, my_image_comp *out, int levels);
int synthesis_5_3(char* inputFile, char* outputFile, int levels);
void synthesis_5_3(my_image_comp *in, my_image_comp *out, int levels);
int synthesis_9_7(char* inputFile, char* outputFile, int levels);
void synthesis_9_7(my_image_comp *in, my_image_comp *out, int levels);
void analysis(my_image_comp *in, my_image_comp *out, int spacing, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL);
void synthesis(my_image_comp *in, my_image_comp *out, int spacing, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL);
int increaseWaveletLevel(char* inputFile, char* outputFile, int levels);
void increaseWaveletLevel(my_image_comp *in, my_image_comp *out, int scale);
void copyBuffer(my_image_comp *copyFrom, my_image_comp *copyTo);

#endif //WAVELET_UNIT_H