#ifndef AFFINE_UNIT_H
#define AFFINE_UNIT_H
#include <armadillo>

using namespace std;
using namespace arma;

Mat<float> matchImages(my_image_comp *LR, my_image_comp *HR, float *LRpoints, float *HRpoints);
Mat<float> computeAffineTransform(float *LRpoints, float *HRpoints);
Mat<float> scaleTransform(Mat<float> transform, int* scaleDiff);
void resampleImage(my_image_comp *in, my_image_comp *out, Mat<float> transform, Mat<float> offset);
double getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff);
float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11);
float bilinear_interpolation_1D(float x, float x0_val, float x1_val);

int matchImages(char *LRfilename, char *HRfilename, float* LRpoints, float *HRpoints);
int resampleImage(char* inputFile, char* outputFile, arma::Mat<float> transform);

#endif //AFFINE_UNIT_H