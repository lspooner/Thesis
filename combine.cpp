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
void maxValueCombinationRetainLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset);
void maxValueCombinationAverageLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset);
void maxValueCombinationBlendLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset);
void allBlendCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset);
void getLuminance(my_image_comp *in, my_image_comp *out);
void getHRinfo(my_image_comp *image, Mat<float> offset, int scaleDiff, int height, int width, float* HRinfo);
void getDiscontinuities(my_image_comp *image, Mat<float> offset, int height, int width, float* discontinuities);

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

    int LRwidth = input_comps_LR[0].width, LRheight = input_comps_LR[0].height;
    int HRwidth = input_comps_HR[0].width, HRheight = input_comps_HR[0].height;

    /*printf("Enter LR coordinates: ");
    float LRpoints[6];
    assert(scanf("%f %f %f %f %f %f", (&LRpoints[0]), (&LRpoints[1]), (&LRpoints[2]), (&LRpoints[3]), (&LRpoints[4]), (&LRpoints[5])) == 6);

    printf("Enter HR coordinates: ");
    float HRpoints[6];
    assert(scanf("%f %f %f %f %f %f", (&HRpoints[0]), (&HRpoints[1]), (&HRpoints[2]), (&HRpoints[3]), (&HRpoints[4]), (&HRpoints[5])) == 6);*/

    //noteboard
    //float LRpoints[6] = {1337, 955, 2047, 939, 1919, 1484};
    //float HRpoints[6] = {764, 413, 3077, 353, 2666, 2132};
    float LRpoints[6] = {1345, 959, 2055, 943, 1927, 1488};
    //float HRpoints[6] = {764, 413, 3070, 353, 2670, 2136};
    float HRpoints[6] = {763.5, 413, 3074, 353, 2669.5, 2136};


    //1345.000000 959.000000 2055.000000 943.000000 1927.000000 1488.000000 
    //763.500000 413.000000 3074.000000 353.000000 2669.500000 2136.000000 


    //float LRpoints[6] = {1571, 1339, 1815, 874, 1239, 1356};
    //float HRpoints[6] = {1532, 1658, 2320, 142, 450, 1720};

    //float LRpoints[6] = {1571.19, 1339.3, 1618.92, 1093.47, 1815.09, 874.113};
    //float HRpoints[6] = {1532.19, 1657.67, 1683.03, 862.294, 2320.36, 141.963};

/*
-- Good Match [0] Keypoint 1: [1532.19, 1657.67] Keypoint 2: [1571.19, 1339.3]
-- Good Match [1] Keypoint 1: [1683.03, 862.294] Keypoint 2: [1618.92, 1093.47]
-- Good Match [2] Keypoint 1: [2320.36, 141.963] Keypoint 2: [1815.09, 874.113]

---- Good Match [0] Keypoint 1: [1532.19, 1657.67] Keypoint 2: [1571.19, 1339.3]
-- Good Match [1] Keypoint 1: [1683.03, 862.294] Keypoint 2: [1618.92, 1093.47]
---- Good Match [2] Keypoint 1: [2320.36, 141.963] Keypoint 2: [1815.09, 874.113]
---- Good Match [3] Keypoint 1: [449.823, 1719.87] Keypoint 2: [1238.89, 1356.32]
-- Good Match [4] Keypoint 1: [1694.43, 257.949] Keypoint 2: [1622.91, 908.82]
-- Good Match [5] Keypoint 1: [1963.79, 969.275] Keypoint 2: [1707.59, 1067.16]
-- Good Match [6] Keypoint 1: [1736.52, 565.579] Keypoint 2: [1635.58, 1002.47]
-- Good Match [7] Keypoint 1: [1991.23, 197.346] Keypoint 2: [1635.58, 1002.47]
*/


    //bookshelf
    //float LRpoints[6] = {1442, 1117, 1861, 1094, 1786, 1369};
    //float HRpoints[6] = {112, 628, 2983, 502, 2457, 2381};

    //TODO: work out borders!!!

    //find transform between images
    //hence find scale between images
    //scale hardcoded and transform assumed to be I for now

    Mat<float> transform_offset;

    if(num_comps == 3){
        my_image_comp *LR_luminance = new my_image_comp;
        LR_luminance->init(LRheight, LRwidth, 0);
        getLuminance(input_comps_LR, LR_luminance);
        assert(outputBMP((char*)"images/LR_luminance.bmp", LR_luminance, 1) == 0);

        my_image_comp *HR_luminance = new my_image_comp;
        HR_luminance->init(HRheight, HRwidth, 2);
        getLuminance(input_comps_HR, HR_luminance);
        assert(outputBMP((char*)"images/HR_luminance.bmp", HR_luminance, 1) == 0);

        transform_offset = matchImages(LR_luminance, HR_luminance, LRpoints, HRpoints);
    } else {
        transform_offset = matchImages(input_comps_LR, input_comps_HR, LRpoints, HRpoints);
    }

    transform_offset = computeAffineTransform(LRpoints, HRpoints);


    /*monument match
       1.6102e-01   2.5884e-03   1.9040e+03
  -3.2612e-03   1.6129e-01   1.4617e+03*/

    /*noteboard match
    [0.3155647850225765, -0.0002451978272692557, 1088.04861934722;
    -0.178128803497289, 0.283625685752317, 1144.337763792789]

    [0.3054285682198384, -0.002093241794122846, 1106.683261008859;
    0.006549088440108208, 0.3103152772793613, 814.8630937837338]
    */


    /*Mat<float> transform_offset(2, 3);
    transform_offset(0,0) = 0.3054285682198384;
    transform_offset(0,1) = -0.002093241794122846;
    transform_offset(0,2) = 1106.683261008859;
    transform_offset(1,0) = 0.006549088440108208;
    transform_offset(1,1) = 0.3103152772793613;
    transform_offset(1,2) = 814.8630937837338;*/

    /*Mat<float> transform = transform_offset.submat(0, 0, 1, 1);
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

    offset_resample(0,0) += scaleDiff*((int)offset(0,0)%(int)pow(2, waveletLevels-levelDiff));
    offset_resample(1,0) += scaleDiff*((int)offset(1,0)%(int)pow(2, waveletLevels-levelDiff));

    //offset_resample *= -1;

    offset(0,0) -= (int)offset(0,0)%(int)pow(2, waveletLevels-levelDiff);
    offset(1,0) -= (int)offset(1,0)%(int)pow(2, waveletLevels-levelDiff);*/

    Mat<float> transform(2,2);
    Mat<float> offset_resample(2,1);
    Mat<float> offset(2,1);
    int scaleDiff;

    scaleTransform(transform_offset, waveletLevels, &transform, &offset_resample, &offset, &scaleDiff);

    int levelDiff = log2(scaleDiff);

    //offset(0,0) = 1108.0;
    //offset(1,0) = 828.0;


    printf("scaleDiff = %d, levelDiff = %d\n", scaleDiff, levelDiff);
    cout << transform << endl;
    cout << offset_resample << endl;
    cout << offset << endl;

    int Outwidth = LRwidth*scaleDiff, Outheight = LRheight*scaleDiff;

    my_image_comp *output_comps = new my_image_comp[num_comps];

    float totalHRinfo = 0.0;
    float totalDiscontinuities = 0.0;

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

        float HRinfo;
        getHRinfo(out_wavelet, offset*scaleDiff, scaleDiff, HRheight, HRwidth, &HRinfo);

        printf("n = %d, about to synthesise\n", n);
        output_comps[n].init(Outheight, Outwidth, 0); // Don't need a border for output
        if(waveletType == WAVELET_53){
            synthesis_5_3(out_wavelet, output_comps+n, waveletLevels, offset_wavelet);
        } else {
            synthesis_9_7(out_wavelet, output_comps+n, waveletLevels, offset_wavelet);
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
    } else if(method == MAX_VALUE_RETAIN_LR){
        maxValueCombinationRetainLR(inLR, inHR, out, LRspacing, Roffset, Coffset);
    } else if(method == MAX_VALUE_AVERAGE_LR){
        maxValueCombinationAverageLR(inLR, inHR, out, LRspacing, Roffset, Coffset);
    } else if(method == MAX_VALUE_BLEND_LR){
        maxValueCombinationBlendLR(inLR, inHR, out, LRspacing, Roffset, Coffset);
    } else if(method == BLEND_ALL){
        allBlendCombination(inLR, inHR, out, LRspacing, Roffset, Coffset);
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

void maxValueCombinationRetainLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset){
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
                if(inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)] == INVALID || (r%LRspacing == 0 && c%LRspacing == 0)){
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

void maxValueCombinationBlendLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset){
    int Roffset_high = Roffset+inHR->height;
    int Coffset_high = Coffset+inHR->width;

    int blendBoundary = 20;

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
                } else if(r%LRspacing == 0 && c%LRspacing == 0){
                    if((r + blendBoundary*LRspacing < Roffset_high && r - blendBoundary*LRspacing > Roffset) && (c + blendBoundary*LRspacing < Coffset_high && c - blendBoundary*LRspacing > Coffset)){
                        out->buf[r*out->stride+c] = inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)];
                    } else {
                        float scaleR = min(((float)(r-Roffset)/(float)LRspacing)/(float)blendBoundary, ((float)(Roffset_high-r)/(float)LRspacing)/(float)blendBoundary);
                        float scaleC = min(((float)(c-Coffset)/(float)LRspacing)/(float)blendBoundary, ((float)(Coffset_high-c)/(float)LRspacing)/(float)blendBoundary);
                        float scale = min(scaleR, scaleC);
                        out->buf[r*out->stride+c] = (1-scale)*inLR->buf[r*inLR->stride+c] + scale*inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)];
                    }
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


void maxValueCombinationAverageLR(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset){
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
                } else if(r%LRspacing == 0 && c%LRspacing == 0){
                    out->buf[r*out->stride+c] = (inLR->buf[r*inLR->stride+c]+inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)])/2.0;
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

void allBlendCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset){
    int Roffset_high = Roffset+inHR->height;
    int Coffset_high = Coffset+inHR->width;

    int blendBoundary = 20;

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
                } else if((r + blendBoundary*LRspacing < Roffset_high && r - blendBoundary*LRspacing > Roffset) && (c + blendBoundary*LRspacing < Coffset_high && c - blendBoundary*LRspacing > Coffset)){
                    out->buf[r*out->stride+c] = inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)];
                } else {
                    float scaleR = min(((float)(r-Roffset)/(float)LRspacing)/(float)blendBoundary, ((float)(Roffset_high-r)/(float)LRspacing)/(float)blendBoundary);
                    float scaleC = min(((float)(c-Coffset)/(float)LRspacing)/(float)blendBoundary, ((float)(Coffset_high-c)/(float)LRspacing)/(float)blendBoundary);
                    float scale = min(scaleR, scaleC);
                    out->buf[r*out->stride+c] = (1-scale)*inLR->buf[r*inLR->stride+c] + scale*inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)];
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

void getHRinfo(my_image_comp *image, Mat<float> offset, int scaleDiff, int height, int width, float* HRinfo){
    //calculate high res info by adding all the wavelet coefficients in the highest bands
    *HRinfo = 0.0;
    for(int r = (int)offset(1); r < (int)offset(1)+height; r++){
        for(int c = (int)offset(0); c < (int)offset(0) + width; c++){
            if(c%scaleDiff != 0 || r%scaleDiff != 0){
                *HRinfo += abs(image->buf[r*image->stride+c]);
            }
        }
    }
}

void getDiscontinuities(my_image_comp *image, Mat<float> offset, int height, int width, float* discontinuities){
    //calculate the discontinuities by adding up the step difference around the edge of the image
    *discontinuities = 0.0;
    for(int r = (int)offset(1); r <= (int)offset(1)+height; r+=height){
        for(int c = (int)offset(0); c < (int)offset(0)+width; c++){
            *discontinuities += abs(image->buf[r*image->stride+c] - image->buf[(r-1)*image->stride+c]);
        }
    }

    for(int c = (int)offset(0); c <= (int)offset(0)+width; c+=width){
        for(int r = (int)offset(1); r < (int)offset(1)+height; r++){
            *discontinuities += abs(image->buf[r*image->stride+c] - image->buf[r*image->stride+c-1]);
        }
    }
}