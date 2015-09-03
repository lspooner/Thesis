#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "wavelet.h"
#include "combine.h"
#include "affine.h"

#include <iostream>

void combineWavelets(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int method, Mat<int>offset);
void maxValueCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset);
void maxValueCombinationRetainLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset);
void maxValueCombinationAverageLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset);
void maxValueCombinationBlendLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset);
void allBlendCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset);
void getLuminance(my_image_comp *in, my_image_comp *out);
void getHRinfo(my_image_comp *image, Mat<int> offset, int scaleDiff, int height, int width, float* HRinfo);
void getDiscontinuities(my_image_comp *image, Mat<int> offset, int height, int width, float* discontinuities);

int combineImages(char* LRinputFile, char* HRinputFile, char* outputFile, int waveletType, int combineType){
    // Read the input images
    int err_code;
    my_image_comp *input_comps_LR = NULL;
    my_image_comp *input_comps_HR = NULL;
    int num_comps;
    int waveletLevels = 4;
    //printf("about to read in images\n");
    if ((err_code = readBMP(LRinputFile, &input_comps_LR, HL_A97_LP*pow(2, waveletLevels)+HL_A97_HP, &num_comps)) != 0){
        return err_code;
    }
    //printf("Read in LR image\n");
    if ((err_code = readBMP(HRinputFile, &input_comps_HR, 2, &num_comps)) != 0){ //border need only be 1 when transform is implemented
        return err_code;
    }
    //printf("Read in HR image\n");

    int LRwidth = input_comps_LR[0].width, LRheight = input_comps_LR[0].height;
    int HRwidth = input_comps_HR[0].width, HRheight = input_comps_HR[0].height;

    /*printf("Enter LR coordinates: ");
    float LRpoints[6];
    assert(scanf("%f %f %f %f %f %f", (&LRpoints[0]), (&LRpoints[1]), (&LRpoints[2]), (&LRpoints[3]), (&LRpoints[4]), (&LRpoints[5])) == 6);

    printf("Enter HR coordinates: ");
    float HRpoints[6];
    assert(scanf("%f %f %f %f %f %f", (&HRpoints[0]), (&HRpoints[1]), (&HRpoints[2]), (&HRpoints[3]), (&HRpoints[4]), (&HRpoints[5])) == 6);*/

    //noteboard
    float LRpoints[6] = {1345, 959, 2055, 943, 1927, 1488};
    //float HRpoints[6] = {764, 413, 3070, 353, 2670, 2136};
    float HRpoints[6] = {763.5, 413, 3074, 353, 2669.5, 2136};


    //bookshelf
    //float LRpoints[6] = {1442, 1117, 1861, 1094, 1786, 1369};
    //float HRpoints[6] = {112, 628, 2983, 502, 2457, 2381};


    Mat<float> transform_offset;

    if(num_comps == 3){
        my_image_comp *LR_luminance = new my_image_comp;
        LR_luminance->init(LRheight, LRwidth, 0);
        getLuminance(input_comps_LR, LR_luminance);
        //assert(outputBMP((char*)"images/LR_luminance.bmp", LR_luminance, 1) == 0);

        my_image_comp *HR_luminance = new my_image_comp;
        HR_luminance->init(HRheight, HRwidth, 2);
        getLuminance(input_comps_HR, HR_luminance);
        //assert(outputBMP((char*)"images/HR_luminance.bmp", HR_luminance, 1) == 0);

        transform_offset = matchImages(LR_luminance, HR_luminance, LRpoints, HRpoints);
    } else {
        transform_offset = matchImages(input_comps_LR, input_comps_HR, LRpoints, HRpoints);
    }

    transform_offset = computeAffineTransform(LRpoints, HRpoints);

    Mat<float> transform(2,2);
    Mat<float> offset_resample(2,1);
    Mat<int> offset(2,1);
    int scaleDiff;

    scaleTransform(transform_offset, waveletLevels, &transform, &offset_resample, &offset, &scaleDiff);

    int levelDiff = log2(scaleDiff);

    /*printf("scaleDiff = %d, levelDiff = %d\n", scaleDiff, levelDiff);
    cout << transform << endl;
    cout << offset_resample << endl;
    cout << offset << endl;*/

    int Outwidth = LRwidth*scaleDiff, Outheight = LRheight*scaleDiff;

    my_image_comp *output_comps = new my_image_comp[num_comps];

    float totalHRinfo = 0.0;
    float totalDiscontinuities = 0.0;

    for(int n = 0; n < num_comps; n++){
        //printf("n = %d, about to transform HR\n", n);
        my_image_comp *inHR_transformed = new my_image_comp;
        inHR_transformed->init(HRheight, HRwidth, HL_A53_LP*pow(2, waveletLevels)+HL_A53_HP);
        input_comps_HR[n].perform_boundary_extension_symmetric();
        resampleImage(input_comps_HR+n, inHR_transformed, transform, offset_resample);
        delete[] input_comps_HR[n].handle;

        //printf("n = %d, about to analyse HR\n", n);
        my_image_comp *inHR_wavelet = new my_image_comp;
        inHR_wavelet->init(HRheight, HRwidth, 0);
        if(waveletType == WAVELET_53){
            analysis_5_3(inHR_transformed, inHR_wavelet, waveletLevels);
        } else {
            analysis_9_7(inHR_transformed, inHR_wavelet, waveletLevels);
        }
        delete inHR_transformed;

        //printf("n = %d, about to analyse LR\n", n);
        my_image_comp *inLR_wavelet_small = new my_image_comp;
        inLR_wavelet_small->init(LRheight, LRwidth, 0);
        if(waveletType == WAVELET_53){
            analysis_5_3(input_comps_LR+n, inLR_wavelet_small, waveletLevels-levelDiff);
        } else {
            analysis_9_7(input_comps_LR+n, inLR_wavelet_small, waveletLevels-levelDiff);
        }
        delete[] input_comps_LR[n].handle;

        //printf("n = %d, about to increase LR\n", n);
        my_image_comp *inLR_wavelet = new my_image_comp;
        inLR_wavelet->init(Outheight, Outwidth, 0);
        increaseWaveletLevel(inLR_wavelet_small, inLR_wavelet, scaleDiff);
        delete inLR_wavelet_small;

        //printf("n = %d, about to combine wavelets\n", n);
        my_image_comp *out_wavelet = new my_image_comp;
        out_wavelet->init(Outheight, Outwidth, HL_S97_HP*pow(2, waveletLevels)+HL_S97_LP);
        combineWavelets(inLR_wavelet, inHR_wavelet, out_wavelet, pow(2, waveletLevels), combineType, offset*scaleDiff);
        delete inLR_wavelet;
        delete inHR_wavelet;

        float HRinfo;
        getHRinfo(out_wavelet, offset*scaleDiff, scaleDiff, HRheight, HRwidth, &HRinfo);

        //printf("n = %d, about to synthesise\n", n);
        output_comps[n].init(Outheight, Outwidth, 0); // Don't need a border for output
        if(waveletType == WAVELET_53){
            synthesis_5_3(out_wavelet, output_comps+n, waveletLevels);
        } else {
            synthesis_9_7(out_wavelet, output_comps+n, waveletLevels);
        }
        delete out_wavelet;

        float discontinuities;
        getDiscontinuities(output_comps+n, offset*scaleDiff, HRheight, HRwidth, &discontinuities);

        totalHRinfo += HRinfo;
        totalDiscontinuities += discontinuities;
    }

    printf("HRinfo = %f, discontinuities = %f\n", totalHRinfo, totalDiscontinuities);

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] output_comps;
    return 0;
}

int combineWavelets(char* LRinputFile, char* HRinputFile, char* outputFile, int LRspacing, int method, int Roffset, int Coffset){
    // Read the input image
    int err_code;
    my_image_comp *inLR = NULL;
    my_image_comp *inHR = NULL;
    int num_comps;
    if ((err_code = readBMP(LRinputFile, &inLR, 0, &num_comps)) != 0){
        return err_code;
    }
    if ((err_code = readBMP(HRinputFile, &inHR, 0, &num_comps)) != 0){
        return err_code;
    }

    int width = inLR[0].width, height = inLR[0].height;

    // Allocate storage for the filtered output
    int outHeight = height*LRspacing;
    int outWidth = width*LRspacing;
    my_image_comp *output_comps = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        output_comps[n].init(outHeight, outWidth, 0); // Don't need a border for output
    }

    Mat<int> offset(2,1);
    offset(0) = Coffset;
    offset(1) = Roffset;

    for (int n=0; n < num_comps; n++){
        combineWavelets(inLR+n, inHR+n, output_comps+n, LRspacing, method, offset);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] inLR;
    delete[] output_comps;
    return 0;
}

void combineWavelets(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int method, Mat<int>offset){
    //TODO: fix this: note: currently assumes HR image is completely enclosed in LR images

    if(method == MAX_VALUE){
        maxValueCombination(inLR, inHR, out, LRspacing, offset);
    } else if(method == MAX_VALUE_RETAIN_LR){
        maxValueCombinationRetainLR(inLR, inHR, out, LRspacing, offset);
    } else if(method == MAX_VALUE_AVERAGE_LR){
        maxValueCombinationAverageLR(inLR, inHR, out, LRspacing, offset);
    } else if(method == MAX_VALUE_BLEND_LR){
        maxValueCombinationBlendLR(inLR, inHR, out, LRspacing, offset);
    } else if(method == BLEND_ALL){
        allBlendCombination(inLR, inHR, out, LRspacing, offset);
    }
}

void maxValueCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset){
    Mat<int> offset_high(2,1);
    offset_high(1) = offset(1)+inHR->height;
    offset_high(0) = offset(0)+inHR->width;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            if(r < offset(1) || r >= offset_high(1) || c < offset(0) || c >= offset_high(0) ){
                if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                if(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))] == INVALID){
                     out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    float maxVal = max(abs(inLR->buf[r*inLR->stride+c]), abs(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]));
                    if(maxVal != inLR->buf[r*inLR->stride+c] && maxVal != inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]){
                        maxVal *= -1;
                    }
                    out->buf[r*out->stride+c] = maxVal;
                }
            }
        }
    }
}

void maxValueCombinationRetainLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset){
    Mat<int> offset_high(2,1);
    offset_high(1) = offset(1)+inHR->height;
    offset_high(0) = offset(0)+inHR->width;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            if(r < offset(1) || r >= offset_high(1) || c < offset(0) || c >= offset_high(0) ){
                if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                if(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))] == INVALID || (r%LRspacing == 0 && c%LRspacing == 0)){
                     out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    float maxVal = max(abs(inLR->buf[r*inLR->stride+c]), abs(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]));
                    if(maxVal != inLR->buf[r*inLR->stride+c] && maxVal != inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]){
                        maxVal *= -1;
                    }
                    out->buf[r*out->stride+c] = maxVal;
                }
            }
        }
    }
}

void maxValueCombinationBlendLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset){
    Mat<int> offset_high(2,1);
    offset_high(1) = offset(1)+inHR->height;
    offset_high(0) = offset(0)+inHR->width;

    int blendBoundary = 20;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            if(r < offset(1) || r >= offset_high(1) || c < offset(0) || c >= offset_high(0) ){
                if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                if(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))] == INVALID){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else if(r%LRspacing == 0 && c%LRspacing == 0){
                    if((r + blendBoundary*LRspacing < offset_high(1) && r - blendBoundary*LRspacing > offset(1)) && (c + blendBoundary*LRspacing < offset_high(0) && c - blendBoundary*LRspacing > offset(0))){
                        out->buf[r*out->stride+c] = inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))];
                    } else {
                        float scaleR = min(((float)(r-offset(1))/(float)LRspacing)/(float)blendBoundary, ((float)(offset_high(1)-r)/(float)LRspacing)/(float)blendBoundary);
                        float scaleC = min(((float)(c-offset(0))/(float)LRspacing)/(float)blendBoundary, ((float)(offset_high(0)-c)/(float)LRspacing)/(float)blendBoundary);
                        float scale = min(scaleR, scaleC);
                        out->buf[r*out->stride+c] = (1-scale)*inLR->buf[r*inLR->stride+c] + scale*inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))];
                    }
                } else {
                    float maxVal = max(abs(inLR->buf[r*inLR->stride+c]), abs(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]));
                    if(maxVal != inLR->buf[r*inLR->stride+c] && maxVal != inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]){
                        maxVal *= -1;
                    }
                    out->buf[r*out->stride+c] = maxVal;
                }
            }
        }
    }
}


void maxValueCombinationAverageLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset){
    Mat<int> offset_high(2,1);
    offset_high(1) = offset(1)+inHR->height;
    offset_high(0) = offset(0)+inHR->width;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            if(r < offset(1) || r >= offset_high(1) || c < offset(0) || c >= offset_high(0) ){
                if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                if(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))] == INVALID){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = (inLR->buf[r*inLR->stride+c]+inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))])/2.0;
                } else {
                    float maxVal = max(abs(inLR->buf[r*inLR->stride+c]), abs(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]));
                    if(maxVal != inLR->buf[r*inLR->stride+c] && maxVal != inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))]){
                        maxVal *= -1;
                    }
                    out->buf[r*out->stride+c] = maxVal;
                }
            }
        }
    }
}

void allBlendCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, Mat<int>offset){
    Mat<int> offset_high(2,1);
    offset_high(1) = offset(1)+inHR->height;
    offset_high(0) = offset(0)+inHR->width;

    int blendBoundary = 20;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            if(r < offset(1) || r >= offset_high(1) || c < offset(0) || c >= offset_high(0) ){
                if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                if(inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))] == INVALID){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else if((r + blendBoundary*LRspacing < offset_high(1) && r - blendBoundary*LRspacing > offset(1)) && (c + blendBoundary*LRspacing < offset_high(0) && c - blendBoundary*LRspacing > offset(0))){
                    out->buf[r*out->stride+c] = inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))];
                } else {
                    float scaleR = min(((float)(r-offset(1))/(float)LRspacing)/(float)blendBoundary, ((float)(offset_high(1)-r)/(float)LRspacing)/(float)blendBoundary);
                    float scaleC = min(((float)(c-offset(0))/(float)LRspacing)/(float)blendBoundary, ((float)(offset_high(0)-c)/(float)LRspacing)/(float)blendBoundary);
                    float scale = min(scaleR, scaleC);
                    out->buf[r*out->stride+c] = (1-scale)*inLR->buf[r*inLR->stride+c] + scale*inHR->buf[(r-offset(1))*inHR->stride+(c-offset(0))];
                }
            }
        }
    }
}

void getLuminance(my_image_comp *in, my_image_comp *out){
    Mat<float> colourMatrix(3,3);
    colourMatrix(0,0) = 0.299;
    colourMatrix(0,1) = 0.587;
    colourMatrix(0,2) = 0.114;
    colourMatrix(1,0) = -0.169;
    colourMatrix(1,1) = -0.331;
    colourMatrix(1,2) = 0.500;
    colourMatrix(2,0) = 0.500;
    colourMatrix(2,1) = -0.419;
    colourMatrix(2,2) = -0.081;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            Mat<float> RGB(3,1);
            RGB(0,0) = in[0].buf[r*in[0].stride+c];
            RGB(1,0) = in[1].buf[r*in[1].stride+c];
            RGB(2,0) = in[2].buf[r*in[2].stride+c];

            Mat<float> YCbCr = colourMatrix*RGB;
            out->buf[r*out->stride+c] = YCbCr(0,0);
        }
    }
}

void getHRinfo(my_image_comp *image, Mat<int> offset, int scaleDiff, int height, int width, float* HRinfo){
    //calculate high res info by adding all the wavelet coefficients in the highest bands
    *HRinfo = 0.0;
    for(int r = offset(1); r < offset(1)+height; r++){
        for(int c = offset(0); c < offset(0) + width; c++){
            if(c%scaleDiff != 0 || r%scaleDiff != 0){
                *HRinfo += abs(image->buf[r*image->stride+c]);
            }
        }
    }
}

void getDiscontinuities(my_image_comp *image, Mat<int> offset, int height, int width, float* discontinuities){
    //calculate the discontinuities by adding up the step difference around the edge of the image
    *discontinuities = 0.0;
    for(int r = offset(1); r <= offset(1)+height; r+=height){
        for(int c = offset(0); c < offset(0)+width; c++){
            *discontinuities += abs(image->buf[r*image->stride+c] - image->buf[(r-1)*image->stride+c]);
        }
    }

    for(int c = offset(0); c <= offset(0)+width; c+=width){
        for(int r = offset(1); r < offset(1)+height; r++){
            *discontinuities += abs(image->buf[r*image->stride+c] - image->buf[r*image->stride+c-1]);
        }
    }
}