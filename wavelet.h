#ifndef WAVELET_H
#define WAVELET_H

//performs 5_3 wavelet analysis on an image to the given number of levels
int analysis_5_3(char* inputFile, char* outputFile, int levels);
int analysis_9_7(char* inputFile, char* outputFile, int levels);
int synthesis_5_3(char* inputFile, char* outputFile, int levels);
int synthesis_9_7(char* inputFile, char* outputFile, int levels);

#endif //WAVELET_H