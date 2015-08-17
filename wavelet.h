#ifndef WAVELET_H
#define WAVELET_H

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

//performs 5_3 wavelet analysis on an image to the given number of levels
int analysis_5_3(char* inputFile, char* outputFile, int levels);
int analysis_9_7(char* inputFile, char* outputFile, int levels);
int synthesis_5_3(char* inputFile, char* outputFile, int levels);
int synthesis_9_7(char* inputFile, char* outputFile, int levels);
int increaseWaveletLevel(char* inputFile, char* outputFile, int scale);

void analysis_5_3(my_image_comp *in, my_image_comp *out, int levels, Mat<float> offset);
void analysis_9_7(my_image_comp *in, my_image_comp *out, int levels, Mat<float> offset);
void synthesis_5_3(my_image_comp *in, my_image_comp *out, int levels, Mat<float> offset);
void synthesis_9_7(my_image_comp *in, my_image_comp *out, int levels, Mat<float> offset);
void increaseWaveletLevel(my_image_comp *in, my_image_comp *out, int scale);

#endif //WAVELET_H