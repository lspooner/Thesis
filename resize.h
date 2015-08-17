#ifndef RESIZE_H
#define RESIZE_H

//reduces an image by an integer scale
int reduceImage(char* inputFile, char* outputFile, int H, int scale);

//enlarges an image by an integer scale
int enlargeImage(char* inputFile, char* outputFile, int H, int scale);

void enlargementFilter(my_image_comp *in, my_image_comp *out, int halfLength, int scale);
void reductionFilter(my_image_comp *in, my_image_comp *out, int halfLength, int scale);

//filter functions for print debugging
float* makeSincFilter(int halfLength, float centre, float scaling);
float* makeWindow(int halfLength);
float* makeSinc(int halfLenght, float centre, float scaling);
#endif //RESIZE_H