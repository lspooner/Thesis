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

void combineWavelets(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int method, int Roffset, int Coffset);
void maxValueCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset);

int combineImages(char* LRinputFile, char* HRinputFile, char* outputFile, int waveletType, int combineType){
    // Read the input images
    int err_code;
    my_image_comp *input_comps_LR = NULL;
    my_image_comp *input_comps_HR = NULL;
    int num_comps;
    int waveletLevels = 4;
    printf("about to read in images\n");
    if ((err_code = readBMP(LRinputFile, &input_comps_LR, HL_A97_LP*pow(2, waveletLevels)+HL_A97_HP, &num_comps)) != 0){
        return err_code;
    }
    printf("Read in LR image\n");
    if ((err_code = readBMP(HRinputFile, &input_comps_HR, 2, &num_comps)) != 0){ //border need only be 1 when transform is implemented
        return err_code;
    }
    printf("Read in HR image\n");

    /*printf("Enter LR coordinates: ");
    float LRpoints[6];
    assert(scanf("%f %f %f %f %f %f", (&LRpoints[0]), (&LRpoints[1]), (&LRpoints[2]), (&LRpoints[3]), (&LRpoints[4]), (&LRpoints[5])) == 6);

    printf("Enter HR coordinates: ");
    float HRpoints[6];
    assert(scanf("%f %f %f %f %f %f", (&HRpoints[0]), (&HRpoints[1]), (&HRpoints[2]), (&HRpoints[3]), (&HRpoints[4]), (&HRpoints[5])) == 6);*/

    float LRpoints[6] = {1337, 955, 2047, 939, 1919, 1484};
    float HRpoints[6] = {764, 413, 3077, 353, 2666, 2132};

    //TODO: work out borders!!!

    //find transform between images
    //hence find scale between images
    //scale hardcoded and transform assumed to be I for now

    Mat<float> transform_offset = matchImages(input_comps_LR, input_comps_HR, LRpoints, HRpoints);
    //Mat<float> transform_offset = computeAffineTransform(LRpoints, HRpoints);


    /*monument match
       1.6102e-01   2.5884e-03   1.9040e+03
  -3.2612e-03   1.6129e-01   1.4617e+03


    Mat<float> transform_offset(2, 3);
    transform_offset(0,0) = 0.16102;
    transform_offset(0,1) = 0.0025884;
    transform_offset(0,2) = 1904.0;
    transform_offset(1,0) = -0.0032612;
    transform_offset(1,1) = 0.16129;
    transform_offset(1,2) = 1461.7;*/

    Mat<float> transform = transform_offset.submat(0, 0, 1, 1);
    Mat<float> offset = transform_offset.submat(0, 2, 1, 2);
    int scaleDiff;
    transform = scaleTransform(transform, &scaleDiff);
    //transform = transform_offset.submat(0, 0, 1, 1)*3;

    int levelDiff = log2(scaleDiff);

    Mat<float> offset_resample(2, 1);
    offset_resample(0,0) = offset(0,0)-(int)offset(0,0);
    offset_resample(1,0) = offset(1,0)-(int)offset(1,0);
    offset(0,0) = (float)((int)offset(0,0));
    offset(1,0) = (float)((int)offset(1,0));

    offset_resample(0,0) += (int)offset(0,0)%(int)pow(2, waveletLevels-levelDiff);
    offset_resample(1,0) += (int)offset(1,0)%(int)pow(2, waveletLevels-levelDiff);

    offset_resample *= -1;

    offset(0,0) -= (int)offset(0,0)%(int)pow(2, waveletLevels-levelDiff);
    offset(1,0) -= (int)offset(1,0)%(int)pow(2, waveletLevels-levelDiff);

    printf("scaleDiff = %d, levelDiff = %d\n", scaleDiff, levelDiff);
    cout << transform << endl;
    cout << offset_resample << endl;
    cout << offset << endl;

    int LRwidth = input_comps_LR[0].width, LRheight = input_comps_LR[0].height;
    int HRwidth = input_comps_HR[0].width, HRheight = input_comps_HR[0].height;
    int Outwidth = LRwidth*scaleDiff, Outheight = LRheight*scaleDiff;

    my_image_comp *output_comps = new my_image_comp[num_comps];

    for(int n = 0; n < num_comps; n++){
        printf("n = %d, about to transform HR\n", n);
        my_image_comp *inHR_transformed = new my_image_comp;
        inHR_transformed->init(HRheight, HRwidth, HL_A53_LP*pow(2, waveletLevels)+HL_A53_HP);
        input_comps_HR[n].perform_boundary_extension_symmetric();
        resampleImage(input_comps_HR+n, inHR_transformed, transform, offset_resample);
        delete[] input_comps_HR[n].handle;

        Mat<float> offset_wavelet(2, 1);
        offset_wavelet(0,0) = 0;
        offset_wavelet(1,0) = 0;

        printf("n = %d, about to analyse HR\n", n);
        my_image_comp *inHR_wavelet = new my_image_comp;
        inHR_wavelet->init(HRheight, HRwidth, 0);
        if(waveletType == WAVELET_53){
            analysis_5_3(inHR_transformed, inHR_wavelet, waveletLevels, offset_wavelet);
        } else {
            analysis_9_7(inHR_transformed, inHR_wavelet, waveletLevels, offset_wavelet);
        }
        delete inHR_transformed;

        printf("n = %d, about to analyse LR\n", n);
        my_image_comp *inLR_wavelet_small = new my_image_comp;
        inLR_wavelet_small->init(LRheight, LRwidth, 0);
        if(waveletType == WAVELET_53){
            analysis_5_3(input_comps_LR+n, inLR_wavelet_small, waveletLevels-levelDiff, offset_wavelet);
        } else {
            analysis_9_7(input_comps_LR+n, inLR_wavelet_small, waveletLevels-levelDiff, offset_wavelet);
        }
        delete[] input_comps_LR[n].handle;

        printf("n = %d, about to increase LR\n", n);
        my_image_comp *inLR_wavelet = new my_image_comp;
        inLR_wavelet->init(Outheight, Outwidth, 0);
        increaseWaveletLevel(inLR_wavelet_small, inLR_wavelet, scaleDiff);
        delete inLR_wavelet_small;

        printf("n = %d, about to combine wavelets\n", n);
        my_image_comp *out_wavelet = new my_image_comp;
        out_wavelet->init(Outheight, Outwidth, HL_S97_HP*pow(2, waveletLevels)+HL_S97_LP);
        combineWavelets(inLR_wavelet, inHR_wavelet, out_wavelet, pow(2, waveletLevels), combineType, (int)offset(1,0)*scaleDiff, (int)offset(0,0)*scaleDiff);
        delete inLR_wavelet;
        delete inHR_wavelet;

        printf("n = %d, about to synthesise\n", n);
        output_comps[n].init(Outheight, Outwidth, 0); // Don't need a border for output
        if(waveletType == WAVELET_53){
            synthesis_5_3(out_wavelet, output_comps+n, waveletLevels, offset_wavelet);
        } else {
            synthesis_9_7(out_wavelet, output_comps+n, waveletLevels, offset_wavelet);
        }
        delete out_wavelet;
    }

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

    for (int n=0; n < num_comps; n++){
        combineWavelets(inLR+n, inHR+n, output_comps+n, LRspacing, method, Roffset, Coffset);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] inLR;
    delete[] output_comps;
    return 0;
}

void combineWavelets(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int method, int Roffset, int Coffset){
    //note: currently assumes HR image is completely enclosed in LR images
    //assert(out->height == inLR->height*LRspacing);
    //assert(out->width == inLR->width*LRspacing);

    printf("Combining wavelets now!\n");

    if(method == MAX_VALUE){
        maxValueCombination(inLR, inHR, out, LRspacing, Roffset, Coffset);
    }
}

void maxValueCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset){
    int Roffset_high = Roffset+inHR->height;
    int Coffset_high = Coffset+inHR->width;

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            if(r < Roffset || r >= Roffset_high || c < Coffset || c >= Coffset_high ){
                if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                if(inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)] == INVALID){
                     out->buf[r*out->stride+c] = inLR->buf[r*inLR->stride+c];
                } else {
                    float maxVal = max(abs(inLR->buf[r*inLR->stride+c]), abs(inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)]));
                    if(maxVal != inLR->buf[r*inLR->stride+c] && maxVal != inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)]){
                        maxVal *= -1;
                    }
                    out->buf[r*out->stride+c] = maxVal;//inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)];
                }
            }
        }
    }
}