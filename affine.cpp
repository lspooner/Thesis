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
//Mat<float> scaleTransform(Mat<float> transform, int *scaleDiff);
void scaleTransform(Mat<float> transform_input, int waveletLevels, Mat<float> *sTransform, Mat<float> *HRoffset, Mat<int> *LRoffset, int *scaleDiff);
void resampleImage(my_image_comp *in, my_image_comp *out, Mat<float> transform, Mat<float> offset);
double getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff);
double getMetric_Wavelet(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff);
double getMetric_Wavelet2(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff);
double getMetric_MultipleWavelet(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff, int levels);
double getMetric_Intensity(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff);
float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11);
float bilinear_interpolation_1D(float x, float x0_val, float x1_val);

//TODO: make matching work if images don't entirely overlap

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

    int waveletLevels = 4;

    Mat<float> transform(2,2);
    Mat<float> offset_resample(2,1);
    Mat<int> offset(2,1);
    int scaleDiff;

    scaleTransform(transform_offset, waveletLevels, &transform, &offset_resample, &offset, &scaleDiff);

    my_image_comp *HRTransformed = new my_image_comp;
    int boundary = 2*256+1; //TODO: see below
    HRTransformed->init(HR->height, HR->width, boundary);
    resampleImage(HR, HRTransformed, transform, offset_resample);

    double bestMatchQuality = getMatchMetric(LR, HRTransformed, offset, scaleDiff);
    //printf("Initial mse = %f\n", bestMatchQuality);

    float tempHRpoints[6];
    for (int i = 0; i < 6; ++i){
        tempHRpoints[i] = HRpoints[i];
    }
    float tempBestMatchQuality = bestMatchQuality;
    bool changed = true;
    while(changed){
        changed = false;
        for(int i = 0; i < 6; i+=2){
            for(float diff_x = -0.5; diff_x <= 0.5; diff_x+=0.5){
                for(float diff_y = -0.5; diff_y <= 0.5; diff_y+=0.5){
                    if(!(diff_x == 0 && diff_y == 0)){
                        HRpoints[i] += diff_x;
                        HRpoints[i+1] += diff_y;
                        transform_offset = computeAffineTransform(LRpoints, HRpoints);

                        scaleTransform(transform_offset, waveletLevels, &transform, &offset_resample, &offset, &scaleDiff);

                        delete HRTransformed;
                        HRTransformed = new my_image_comp;
                        boundary = 2*256+1; //TODO: work out boundary dynamically
                        HRTransformed->init(HR->height, HR->width, boundary);

                        //cout << transform << endl;
                        //cout << offset_resample << endl;
                        resampleImage(HR, HRTransformed, transform, offset_resample);
                        double matchQuality = getMatchMetric(LR, HRTransformed, offset, scaleDiff);

                        //printf("new mse = %f, best mse = %f\n", matchQuality, bestMatchQuality);

                        if(matchQuality < tempBestMatchQuality){
                            changed = true;
                            tempBestMatchQuality = matchQuality;
                            for (int k = 0; k < 6; ++k){
                                tempHRpoints[k] = HRpoints[k];
                            }
                        }
                        HRpoints[i] -= diff_x;
                        HRpoints[i+1] -= diff_y;
                    }
                }
            }
            bestMatchQuality = tempBestMatchQuality;
            for (int k = 0; k < 6; ++k){
                HRpoints[k] = tempHRpoints[k];
            }
        }
    }
    /*printf("==========Matched Points==============\n");
    printf("                 LR                   \n");
    for (int i = 0; i < 6; ++i){
        printf("%f ", LRpoints[i]);
    }
    printf("\n");
    printf("                 HR                   \n");
    for (int i = 0; i < 6; ++i){
        printf("%f ", HRpoints[i]);
    }
    printf("\n");
    printf("======================================\n");*/
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

void scaleTransform(Mat<float> transform_input, int waveletLevels, Mat<float> *sTransform, Mat<float> *HRoffset, Mat<int> *LRoffset, int *scaleDiff){
    assert(transform_input.n_rows == 2);
    assert(transform_input.n_cols == 3);

    (*sTransform) = transform_input.submat(0, 0, 1, 1);

    Mat<float> tinv = inv((*sTransform));

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

    *scaleDiff = pow(2, (int)(log2(maxCorner)+0.5));
    int levelDiff = log2(*scaleDiff);

    float sx = 1.0/sqrt((*sTransform)(0,0)*(*sTransform)(0,0) + (*sTransform)(1,0)*(*sTransform)(1,0));
    float sy = 1.0/sqrt((*sTransform)(0,1)*(*sTransform)(0,1) + (*sTransform)(1,1)*(*sTransform)(1,1));

    //printf("sd = %d, sx = %f, sy = %f\n", (*scaleDiff), sx, sy);
    Mat<float> scaleT(2,2);
    scaleT(0,0) = sx*sx/(*scaleDiff);
    scaleT(0,1) = 0;
    scaleT(1,0) = 0;
    scaleT(1,1) = sy*sy/(*scaleDiff);

    (*HRoffset) = transform_input.submat(0, 2, 1, 2);

    //offset between LR image and HR image
    (*LRoffset)(0) = (int)(*HRoffset)(0);
    (*LRoffset)(1) = (int)(*HRoffset)(1);

    //offset that HR image will be warped by
    (*HRoffset)(0) = (*HRoffset)(0)-(*LRoffset)(0);
    (*HRoffset)(1) = (*HRoffset)(1)-(*LRoffset)(1);

    (*HRoffset)(0) += (*LRoffset)(0)%(int)pow(2, waveletLevels-levelDiff);
    (*HRoffset)(1) += (*LRoffset)(1)%(int)pow(2, waveletLevels-levelDiff);

    (*HRoffset)(0) *= -1;
    (*HRoffset)(1) *= -1;

    (*LRoffset)(0) -= (*LRoffset)(0)%(int)pow(2, waveletLevels-levelDiff);
    (*LRoffset)(1) -= (*LRoffset)(1)%(int)pow(2, waveletLevels-levelDiff);

    (*sTransform) = (*sTransform)*scaleT;
    (*HRoffset) = (*scaleDiff)*(*HRoffset);
}

/*Mat<float> scaleTransform(Mat<float> transform, int *scaleDiff){
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

    float sx = 1.0/sqrt(transform(0,0)*transform(0,0) + transform(1,0)*transform(1,0));
    float sy = 1.0/sqrt(transform(0,1)*transform(0,1) + transform(1,1)*transform(1,1));

    printf("sd = %d, sx = %f, sy = %f\n", (*scaleDiff), sx, sy);
    Mat<float> scaleT(2,2);
    scaleT(0,0) = sx*sx/(*scaleDiff);
    scaleT(0,1) = 0;
    scaleT(1,0) = 0;
    scaleT(1,1) = sy*sy/(*scaleDiff);

    //printf("end of scaling function:, scaleDiff = %d, maxCorner = %f\n", *scaleDiff, maxCorner);
    //cout << transform << endl; 
    //cout << (1.0/(float)(*scaleDiff))*transform << endl; 
    delete[] corners;
    //return ((maxCorner*maxCorner)/(float)(*scaleDiff))*transform;
    //return (1.0/(float)*scaleDiff)*tinv;
    return transform*scaleT;
}*/

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

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->stride; c++){
            Mat<float> output(2, 1);
            output(1, 0) = r;
            output(0, 0) = c;
            Mat<float> input = transform*output+offset;

            float inR = input(1, 0);
            float inC = input(0, 0);

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
                out->buf[r*out->stride+c] = INVALID;
            }
        }
    }
}

//TODO: add variable to change match metric
double getMatchMetric(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff){
    //return getMetric_Intensity(LR, HR, offset, scaleDiff);
    //return getMetric_Wavelet(LR, HR, offset, scaleDiff);
    //return getMetric_Wavelet2(LR, HR, offset, scaleDiff);
    return getMetric_MultipleWavelet(LR, HR, offset, scaleDiff, 6);
}

double getMetric_Intensity(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff){
    my_image_comp *HR_reduced = new my_image_comp;
    HR_reduced->init(HR->height/scaleDiff, HR->width/scaleDiff, 0);
    HR->perform_boundary_extension_symmetric();

    reductionFilter(HR, HR_reduced, scaleDiff*10, scaleDiff);

    double mse = 0.0;
    int numPixels = 0;

    for(int r = 0; r < HR_reduced->height; r++){
        for(int c = 0; c < HR_reduced->width; c++){
            if(HR_reduced->buf[r*HR_reduced->stride+c] != INVALID && LR->buf[(r+offset(1))*LR->stride+c+offset(0)] != INVALID){
                double err = HR_reduced->buf[r*HR_reduced->stride+c] - LR->buf[(r+offset(1))*LR->stride+c+offset(0)];
                mse += err*err;
                numPixels++;
            }
        }
    }
    delete HR_reduced;
    return mse/numPixels;
}

double getMetric_Wavelet(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff){
    my_image_comp *LR_offset = new my_image_comp;
    LR_offset->init((HR->height)/scaleDiff, (HR->width)/scaleDiff, MATCH_LEVELS*2+1);
    for(int r = 0; r < LR_offset->height; r++){
        for(int c = 0; c < LR_offset->width; c++){
            LR_offset->buf[r*LR_offset->stride+c] = LR->buf[(r+offset(1))*LR->stride+(c+offset(0))];
        }
    }

    my_image_comp *LR_offset_wavelet = new my_image_comp;
    LR_offset_wavelet->init(LR_offset->height, LR_offset->width, 0);

    analysis_5_3(LR_offset, LR_offset_wavelet, MATCH_LEVELS);
    delete LR_offset;

    my_image_comp *HR_wavelet = new my_image_comp;
    HR_wavelet->init(HR->height, HR->width, 0);

    int levels = log2(scaleDiff)+MATCH_LEVELS;
    analysis_5_3(HR, HR_wavelet, levels);

    double mse = 0.0;
    int numPixels = 0;

    int LR_spacing = (int)pow(2, MATCH_LEVELS);
    for(int r = 0; r < LR_offset_wavelet->height; r++){
        for(int c = 0; c < LR_offset_wavelet->width; c++){
            if((r%LR_spacing != 0 || c%LR_spacing != 0) && HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff] != INVALID && LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] != INVALID){
                double err = LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] - HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff];
                mse += err*err;
                numPixels++;
            }
        }
    }
    delete LR_offset_wavelet;
    delete HR_wavelet;
    return mse/numPixels;
}

double getMetric_Wavelet2(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff){
    my_image_comp *LR_offset = new my_image_comp;
    LR_offset->init((HR->height)/scaleDiff, (HR->width)/scaleDiff, 2*2+1);
    for(int r = 0; r < LR_offset->height; r++){
        for(int c = 0; c < LR_offset->width; c++){
            LR_offset->buf[r*LR_offset->stride+c] = LR->buf[(r+offset(1))*LR->stride+(c+offset(0))];
        }
    }

    my_image_comp *LR_offset_wavelet = new my_image_comp;
    LR_offset_wavelet->init(LR_offset->height, LR_offset->width, 0);

    analysis_5_3(LR_offset, LR_offset_wavelet, 2);
    delete LR_offset;

    my_image_comp *HR_wavelet = new my_image_comp;
    HR_wavelet->init(HR->height, HR->width, 0);

    analysis_5_3(HR, HR_wavelet, log2(scaleDiff)+2);

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
            }
        }
    }
    delete LR_offset_wavelet;
    delete HR_wavelet;
    return mse/numPixels;
}

double getMetric_MultipleWavelet(my_image_comp *LR, my_image_comp *HR, Mat<int> offset, int scaleDiff, int levels){
    my_image_comp *LR_offset = new my_image_comp;
    LR_offset->init((HR->height)/scaleDiff, (HR->width)/scaleDiff, 2*pow(2,levels)+1);
    for(int r = 0; r < LR_offset->height; r++){
        for(int c = 0; c < LR_offset->width; c++){
            LR_offset->buf[r*LR_offset->stride+c] = LR->buf[(r+offset(1))*LR->stride+(c+offset(0))];
        }
    }

    my_image_comp *LR_offset_wavelet = new my_image_comp;
    LR_offset_wavelet->init(LR_offset->height, LR_offset->width, 0);

    analysis_5_3(LR_offset, LR_offset_wavelet, levels);
    delete LR_offset;

    my_image_comp *HR_wavelet = new my_image_comp;
    HR_wavelet->init(HR->height, HR->width, 0);

    analysis_5_3(HR, HR_wavelet, log2(scaleDiff)+levels);

    double mse[levels-2] = {0.0};
    int numPixels[levels-2] = {0};

    int LR_spacing = (int)pow(2, levels);
    for(int r = 0; r < LR_offset_wavelet->height; r++){
        for(int c = 0; c < LR_offset_wavelet->width; c++){
            for(int l = 0; l < levels-2; l++){
                if((HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff] != INVALID && LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] != INVALID)){
                    //0 level: thing % lrspacing, values that are r == 1 && c <= 1, 
                    //1 level: thing % lrspacing, values that are r == 2, 3 && c <= 3
                    //n level: thing % lrspacing, values that are r >= 2^n, r < 2^n+1 && c < 2^n+1
                    if((r%LR_spacing >= pow(2, l) && r%LR_spacing < pow(2, l+1) && c%LR_spacing < pow(2, l+1)) || 
                       (c%LR_spacing >= pow(2, l) && c%LR_spacing < pow(2, l+1) && r%LR_spacing < pow(2, l+1))){
                        double err = LR_offset_wavelet->buf[r*LR_offset_wavelet->stride+c] - HR_wavelet->buf[r*scaleDiff*HR_wavelet->stride+c*scaleDiff];
                        mse[l] += err*err;
                        numPixels[l]++;
                    }
                }
            }
        }
    }
    double metric = 0.0;
    delete LR_offset_wavelet;
    delete HR_wavelet;
    for(int l = 0; l < levels-2; l++){
        metric += pow(5, l)*(mse[l]/numPixels[l]);
    }
    return metric;
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