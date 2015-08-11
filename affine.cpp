#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include "image_comps.h"
#include "io_bmp.h"
#include "affine.h"

#define PI 3.141592

using namespace std;
using namespace arma;

Mat<float> matchImages(my_image_comp *LR, my_image_comp *HR, float *LRpoints, float *HRpoints);
Mat<float> computeAffineTransform(float *LRpoints, float *HRpoints);
Mat<float> scaleTransform(Mat<float> transform);
void resampleImage(my_image_comp *in, my_image_comp *out, Mat<float> transform);
float getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<float> offset);
float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11);
float bilinear_interpolation_1D(float x, float x0_val, float x1_val);

int matchImages(char *LRfilename, char *HRfilename, float* LRpoints, float *HRpoints){

    // Read the input image
    int err_code;
    my_image_comp *LR_input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(LRfilename, &LR_input_comps, 1, &num_comps)) != 0){
        return err_code;
    }

    my_image_comp *HR_input_comps = NULL;
    if ((err_code = readBMP(HRfilename, &HR_input_comps, 1, &num_comps)) != 0){
        return err_code;
    }    

    for(int n=0; n < num_comps; n++){
        Mat<float> transform = matchImages(LR_input_comps+n, HR_input_comps+n, LRpoints, HRpoints);
        cout << transform << endl;
    }

    delete[] LR_input_comps;
    delete[] HR_input_comps;
    return 0;
}

Mat<float> matchImages(my_image_comp *LR, my_image_comp *HR, float *LRpoints, float *HRpoints){
    Mat<float> transform_offset = computeAffineTransform(LRpoints, HRpoints);
    Mat<float> transform = transform_offset.submat(0, 0, 1, 1);
    Mat<float> offset = transform_offset.submat(0, 2, 1, 2);
    transform = scaleTransform(transform);

    my_image_comp *HRTransformed = new my_image_comp;
    HRTransformed->init(HR->height, HR->width, 0);
    resampleImage(HR, HRTransformed, transform);

    float bestMatchQuality = getMatchMetric(LR, HRTransformed, offset);

    for(int i = 0; i < 6; i++){
        for(float diff = -0.5; diff <= 0.5; diff+=1){
            HRpoints[i] += diff;
            transform_offset = computeAffineTransform(LRpoints, HRpoints);
            transform = transform_offset.submat(0, 0, 1, 1);
            offset = transform_offset.submat(0, 2, 1, 2);
            transform = scaleTransform(transform);

            resampleImage(HR, HRTransformed, transform);
            float matchQuality = getMatchMetric(LR, HRTransformed, offset);

            if(matchQuality < bestMatchQuality){
                //restart looking into coordinates
                i = -1;
                diff = 1;
            } else {
                HRpoints[i] -= diff;
            }
        }
    }
    return computeAffineTransform(LRpoints, HRpoints);
}

//arrays expected to contain 6 values in the form:
//x1, y1, x2, y2, x3, y3
Mat<float> computeAffineTransform(float *LRpoints, float *HRpoints){
    Mat<float> LRmatrix(2, 3);
    Mat<float> HRmatrix(3, 3);

    for(int c = 0; c < 3; c++){
        for(int r = 0; r < 2; r++){
            LRmatrix(r, c) = LRpoints[2*c+r];
            HRmatrix(r, c) = HRpoints[2*c+r];
        }
    }

    for(int c = 0; c < 3; c++){
        HRmatrix(2, c) = 1;
    }

    return LRmatrix*inv(HRmatrix);
}

Mat<float> scaleTransform(Mat<float> transform){
    assert(transform.n_rows == 2);
    assert(transform.n_cols == 2);

    Mat<float> tInv = inv(transform);

    Mat<float> *corners = new Mat<float>[4];
    for(int i = 0; i < 4; i++){
        corners[i].resize(2, 1);
    }
    corners[0](0,0) = PI;
    corners[0](1,0) = PI;
    corners[1](0,0) = PI;
    corners[1](1,0) = -PI;
    corners[2](0,0) = -PI;
    corners[2](1,0) = PI;
    corners[3](0,0) = -PI;
    corners[3](1,0) = -PI;

    float maxCorner = PI;

    for(int i = 0; i < 4; i++){
        Mat<float> outCorner = tInv*corners[i];
        maxCorner = max(maxCorner, abs(outCorner(0,0)));
        maxCorner = max(maxCorner, abs(outCorner(1,0)));
    }

    int scale = pow(2, (int)(log2(maxCorner/PI)+0.5));

    return scale*transform;
}

int resampleImage(char* inputFile, char* outputFile, arma::Mat<float> transform){

    assert(transform.n_rows == 2);
    assert(transform.n_cols == 2);

    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(inputFile, &input_comps, 1, &num_comps)) != 0){
        return err_code;
    }

    int width = input_comps[0].width, height = input_comps[0].height;

    // Allocate storage for the filtered output

    int outHeight = height;
    int outWidth = width;
    my_image_comp *output_comps = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        output_comps[n].init(outHeight, outWidth, 0); // Don't need a border for output
    }

    // Process the image, all in floating point (easy)
    for (int n=0; n < num_comps; n++){
       input_comps[n].perform_boundary_extension_symmetric();
    }

    for(int n=0; n < num_comps; n++){
        resampleImage(input_comps+n, output_comps+n, transform);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void resampleImage(my_image_comp *in, my_image_comp *out, Mat<float> transform){
    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->stride; c++){
            Mat<float> output(2, 1);
            output(0, 0) = r;
            output(1, 0) = c;
            Mat<float> input = transform*output;

            float inR = input(0, 0);
            float inC = input(1, 0);

            if(inC > -1 && inC < out->width && inR > -1 && inR < out->height){
                out->buf[r*out->stride+c] = bilinear_interpolation_2D(inC-((int)inC), inR-((int)inR), in->buf[((int)inR)*in->stride+((int)inC)], in->buf[((int)inR)*in->stride+((int)inC+1)], in->buf[((int)inR+1)*in->stride+((int)inC)], in->buf[((int)inR+1)*in->stride+((int)inC+1)]);
            }
        }
    }
}

float getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<float> offset){
    return 1;
}

/*
perform 2D bilinear interpolation, 
(x,y) are the normalised coordinates, 
val_xy are the values at the corners

0________1_->x
 |
 |
 |
1|
 |
 \/
 y
*/

float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11){
    float y0 = bilinear_interpolation_1D(x, val_00, val_01);
    float y1 = bilinear_interpolation_1D(x, val_10, val_11);
    return bilinear_interpolation_1D(y, y0, y1);
}

float bilinear_interpolation_1D(float x, float x0_val, float x1_val){
    if(x == 0.0){
        return x0_val;
    } else {
        return (1-x)*x0_val + x*x1_val;
    }
}