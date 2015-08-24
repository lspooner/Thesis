#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include "image_comps.h"
#include "io_bmp.h"
#include "affine.h"
#include "wavelet.h"
#include "resize.h"

#define PI 3.141592
#define MATCH_LEVELS 1

using namespace std;
using namespace arma;

Mat<float> matchImages(my_image_comp *LR, my_image_comp *HR, float *LRpoints, float *HRpoints);
Mat<float> computeAffineTransform(float *LRpoints, float *HRpoints);
Mat<float> scaleTransform(Mat<float> transform, int *scaleDiff);
void resampleImage(my_image_comp *in, my_image_comp *out, Mat<float> transform, Mat<float> offset);
double getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff);
double getMetric_Wavelet(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff);
double getMetric_Wavelet2(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff);
double getMetric_Intensity(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff);
float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11);
float bilinear_interpolation_1D(float x, float x0_val, float x1_val);

int matchImages(char *LRfilename, char *HRfilename, float* LRpoints, float *HRpoints){

    // Read the input image
    int err_code;
    my_image_comp *LR_input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(LRfilename, &LR_input_comps, 0, &num_comps)) != 0){
        return err_code;
    }

    my_image_comp *HR_input_comps = NULL;
    if ((err_code = readBMP(HRfilename, &HR_input_comps, 2, &num_comps)) != 0){
        return err_code;
    }

    for(int n=0; n < num_comps; n++){
        HR_input_comps[n].perform_boundary_extension_symmetric();
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
    int scaleDiff;
    transform = scaleTransform(transform, &scaleDiff);

    Mat<float> offset_resample(2, 1);
    offset_resample(0,0) = offset(0,0)-(int)offset(0,0);
    offset_resample(1,0) = offset(1,0)-(int)offset(1,0);
    offset(0,0) = (float)((int)offset(0,0));
    offset(1,0) = (float)((int)offset(1,0));

    my_image_comp *HRTransformed = new my_image_comp;
    //int spacing = (int)pow(2, scaleDiff+MATCH_LEVELS-1);
    //int boundary = spacing*2+1;
    int boundary = 10*scaleDiff;
    HRTransformed->init(HR->height, HR->width, boundary);
    resampleImage(HR, HRTransformed, transform, offset_resample);

    double bestMatchQuality = getMatchMetric(LR, HRTransformed, offset, scaleDiff);
    printf("Initial mse = %f\n", bestMatchQuality);

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 3; j++){
            for(float diff = -0.1; diff <= 0.1; diff+=0.2){
                printf("i = %d, j = %d, diff = %f\n", i, j, diff);
                cout << transform_offset << endl;
                transform_offset(i, j) += diff;
                cout << transform_offset << endl;
                transform = transform_offset.submat(0, 0, 1, 1);
                offset = transform_offset.submat(0, 2, 1, 2);

                //printf("Transform before scaling:\n");
                //cout << transform << endl;
                transform = scaleTransform(transform, &scaleDiff);
                //printf("Transform after scaling:\n");
                //cout << transform << endl;

                offset_resample(0,0) = offset(0,0)-(int)offset(0,0);
                offset_resample(1,0) = offset(1,0)-(int)offset(1,0);
                offset(0,0) = (float)((int)offset(0,0));
                offset(1,0) = (float)((int)offset(1,0));

                delete HRTransformed;
                HRTransformed = new my_image_comp;
                //int spacing = (int)pow(2, scaleDiff+MATCH_LEVELS-1);
                //boundary = spacing*2+1;
                boundary = 10*scaleDiff;
                HRTransformed->init(HR->height, HR->width, boundary);


                for(int r = 0; r < HR->height; r++){
                    for(int c = 0; c < HR->width; c++){
                        assert(HR->buf[r*HR->stride+c] >= -1);
                        assert(HR->buf[r*HR->stride+c] < 1);
                    }
                }
                cout << transform << endl;
                cout << offset_resample << endl;
                resampleImage(HR, HRTransformed, transform, offset_resample);
                double matchQuality = getMatchMetric(LR, HRTransformed, offset, scaleDiff);

                printf("new mse = %f, best mse = %f\n", matchQuality, bestMatchQuality);

                if(matchQuality < bestMatchQuality){
                    //restart looking into coordinates
                    i = 0;
                    j = -1;
                    diff = 1;
                    bestMatchQuality = matchQuality;
                } else {
                    transform_offset(i, j) -= diff;
                }
            }
        }
    }
    delete HRTransformed;
    return transform_offset;
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

Mat<float> scaleTransform(Mat<float> transform, int *scaleDiff){
    assert(transform.n_rows == 2);
    assert(transform.n_cols == 2);

    Mat<float> tinv = inv(transform);

    Mat<float> *corners = new Mat<float>[4];
    for(int i = 0; i < 4; i++){
        corners[i].resize(2, 1);
    }
    corners[0](0,0) =  PI;
    corners[0](1,0) =  PI;
    corners[1](0,0) =  PI;
    corners[1](1,0) = -PI;
    corners[2](0,0) = -PI;
    corners[2](1,0) =  PI;
    corners[3](0,0) = -PI;
    corners[3](1,0) = -PI;

    float maxCorner = PI;

    for(int i = 0; i < 4; i++){
        Mat<float> outCorner = tinv*corners[i];
        maxCorner = max(maxCorner, abs(outCorner(0,0)));
        maxCorner = max(maxCorner, abs(outCorner(1,0)));
    }

    maxCorner /= PI;

    //printf("MaxCorner = %f\n", maxCorner);

    *scaleDiff = pow(2, (int)(log2(maxCorner)+0.5));

    //printf("end of scaling function:, scaleDiff = %d, maxCorner = %f\n", *scaleDiff, maxCorner);
    //cout << transform << endl; 
    //cout << (1.0/(float)(*scaleDiff))*transform << endl; 
    delete[] corners;
    //return ((maxCorner*maxCorner)/(float)(*scaleDiff))*transform;
    return (1.0/(float)*scaleDiff)*tinv;
}

int resampleImage(char* inputFile, char* outputFile, arma::Mat<float> transform){

    assert(transform.n_rows == 2);
    assert(transform.n_cols == 3);

    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(inputFile, &input_comps, 2, &num_comps)) != 0){
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

    Mat<float> transform_only = transform.submat(0, 0, 1, 1);
    Mat<float> offset = transform.submat(0, 2, 1, 2);

    for(int n=0; n < num_comps; n++){
        resampleImage(input_comps+n, output_comps+n, transform_only, offset);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void resampleImage(my_image_comp *in, my_image_comp *out, Mat<float> transform, Mat<float> offset){
    assert(in->border >= 2);

    //printf("Resampling image with transform:\n");
    //cout << transform << endl;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->stride; c++){
            Mat<float> output(2, 1);
            output(0, 0) = r;
            output(1, 0) = c;
            Mat<float> input = transform*output+offset;

            float inR = input(0, 0);
            float inC = input(1, 0);

            if(inC > -2 && inC < out->width+1 && inR > -2 && inR < out->height+1){
                int r0 = (int)inR;
                int r1 = ((int)inR)+1;
                int c0 = (int)inC;
                int c1 = ((int)inC)+1;
                float v_00 = in->buf[r0*in->stride+c0];
                float v_01 = in->buf[r0*in->stride+c1];
                float v_10 = in->buf[r1*in->stride+c0];
                float v_11 = in->buf[r1*in->stride+c1];
                float x = inC-((int)inC);
                if(x < 0){
                    x = 1+x;
                }
                float y = inR-((int)inR);
                if(y < 0){
                    y = 1+y;
                }
                out->buf[r*out->stride+c] = bilinear_interpolation_2D(x, y, v_00, v_01, v_10, v_11);
            } else {
                //printf("NO OUTPUT, r = %d, inR = %f, c = %d, inC = %f\n", r, inR, c, inC);
                out->buf[r*out->stride+c] = INVALID;
            }
        }
    }
}

double getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff){
    //return getMetric_Intensity(LR, HR, offset, scaleDiff);
    return getMetric_Wavelet(LR, HR, offset, scaleDiff);
    //return getMetric_Wavelet2(LR, HR, offset, scaleDiff);
}

double getMetric_Intensity(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff){
    my_image_comp *HR_reduced = new my_image_comp;
    HR_reduced->init(HR->height/scaleDiff, HR->width/scaleDiff, 0);
    HR->perform_boundary_extension_symmetric();

    reductionFilter(HR, HR_reduced, scaleDiff*10, scaleDiff);

    double mse = 0.0;
    int numPixels = 0;

    for(int r = 0; r < HR_reduced->height; r++){
        for(int c = 0; c < HR_reduced->width; c++){
            if(HR_reduced->buf[r*HR_reduced->stride+c] != INVALID && LR->buf[(r+(int)offset(1,0))*LR->stride+c+(int)offset(0,0)] != INVALID){
                double err = HR_reduced->buf[r*HR_reduced->stride+c] - LR->buf[(r+(int)offset(1,0))*LR->stride+c+(int)offset(0,0)];
                mse += err*err;
                numPixels++;
                //printf("r = %d, c = %d, err = %f, mse = %f\n", r, c, err, mse);
            }
        }
    }
    delete HR_reduced;
    return mse/numPixels;
}

double getMetric_Wavelet(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff){
    my_image_comp *LR_offset = new my_image_comp;
    LR_offset->init((HR->height)/scaleDiff, (HR->width)/scaleDiff, MATCH_LEVELS*2+1);
    for(int r = 0; r < LR_offset->height; r++){
        for(int c = 0; c < LR_offset->width; c++){
            LR_offset->buf[r*LR_offset->stride+c] = LR->buf[(r+(int)offset(1,0))*LR->stride+(c+(int)offset(0,0))];
        }
    }

    Mat<float> analysis_offset(2,1);
    analysis_offset(0,0) = 0;
    analysis_offset(1,0) = 0;

    my_image_comp *LR_offset_wavelet = new my_image_comp;
    LR_offset_wavelet->init(LR_offset->height, LR_offset->width, 0);

    analysis_5_3(LR_offset, LR_offset_wavelet, MATCH_LEVELS, analysis_offset);
    delete LR_offset;

    my_image_comp *HR_wavelet = new my_image_comp;
    HR_wavelet->init(HR->height, HR->width, 0);

    int levels = log2(scaleDiff)+MATCH_LEVELS;
    analysis_5_3(HR, HR_wavelet, levels, analysis_offset);

    double mse = 0.0;
    int numPixels = 0;

    int LR_spacing = (int)pow(2, MATCH_LEVELS);
    for(int r = 0; r < LR_offset_wavelet->height; r++){
        for(int c = 0; c < LR_offset_wavelet->width; c++){
            if((r%LR_spacing != 0 || c%LR_spacing != 0) && HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff] != INVALID && LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] != INVALID){
                double err = LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] - HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff];
                mse += err*err;
                numPixels++;
                //printf("r = %d, c = %d, err = %f, mse = %f\n", r, c, err, mse);
            }
        }
    }
    delete LR_offset_wavelet;
    delete HR_wavelet;
    return mse/numPixels;
}

double getMetric_Wavelet2(my_image_comp *LR, my_image_comp *HR, Mat<float> offset, int scaleDiff){
    my_image_comp *LR_offset = new my_image_comp;
    LR_offset->init((HR->height)/scaleDiff, (HR->width)/scaleDiff, 2*2+1);
    for(int r = 0; r < LR_offset->height; r++){
        for(int c = 0; c < LR_offset->width; c++){
            LR_offset->buf[r*LR_offset->stride+c] = LR->buf[(r+(int)offset(1,0))*LR->stride+(c+(int)offset(0,0))];
        }
    }

    Mat<float> analysis_offset(2,1);
    analysis_offset(0,0) = 0;
    analysis_offset(1,0) = 0;

    my_image_comp *LR_offset_wavelet = new my_image_comp;
    LR_offset_wavelet->init(LR_offset->height, LR_offset->width, 0);

    analysis_5_3(LR_offset, LR_offset_wavelet, 2, analysis_offset);
    delete LR_offset;

    my_image_comp *HR_wavelet = new my_image_comp;
    HR_wavelet->init(HR->height, HR->width, 0);

    analysis_5_3(HR, HR_wavelet, log2(scaleDiff)+2, analysis_offset);

    double mse = 0.0;
    int numPixels = 0;

    int LR_spacing = 4;
    for(int r = 0; r < LR_offset_wavelet->height; r++){
        for(int c = 0; c < LR_offset_wavelet->width; c++){
            if((HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff] != INVALID && LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] != INVALID) && 
                ((r%LR_spacing == 1 && (c%LR_spacing == 0 || c%LR_spacing == 1)) || (c%LR_spacing == 1 && (r%LR_spacing == 0 || r%LR_spacing == 1)))){
                double err = LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] - HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff];
                mse += err*err;
                numPixels++;
                //printf("r = %d, c = %d, err = %f, mse = %f\n", r, c, err, mse);
            }
        }
    }
    delete LR_offset_wavelet;
    delete HR_wavelet;
    return mse/numPixels;
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