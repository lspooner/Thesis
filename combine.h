#ifndef COMBINE_IMAGES_H
#define COMBINE_IMAGES_H

//wavelet types
#define WAVELET_53 0
#define WAVELET_97 1

//wavelet combination methods
#define MAX_VALUE 0

int combineImages(char* LRinputFile, char* HRinputFile, char* outputFile, int waveletType, int combineType);
int combineWavelets(char* LRinputFile, char* HRinputFile, char* outputFile, int LRspacing, int method, int Roffset, int Coffset);

#endif //COMBINE_IMAGES_H