#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "wavelet.h"
#include "combine.h"

void combineWavelets(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int method, int Roffset, int Coffset);
void maxValueCombination(my_image_comp *inLR, my_image_comp *inHR, my_image_comp *out, int LRspacing, int Roffset, int Coffset);

int combineImages(char* LRinputFile, char* HRinputFile, char* outputFile, int waveletType, int combineType){
    // Read the input images
    int err_code;
    my_image_comp *input_comps_LR = NULL;
    my_image_comp *input_comps_HR = NULL;
    int num_comps;
    int waveletLevels = 2;
    if ((err_code = readBMP(LRinputFile, &input_comps_LR, HL_A97_LP*pow(2, waveletLevels)+HL_A97_HP, &num_comps)) != 0){
        return err_code;
    }
    if ((err_code = readBMP(HRinputFile, &input_comps_HR, HL_A97_LP*pow(2, waveletLevels)+HL_A97_HP, &num_comps)) != 0){ //border need only be 1 when transform is implemented
        return err_code;
    }

    //TODO: work out borders!!!

    //find transform between images
    //hence find scale between images
    //scale hardcoded and transform assumed to be I for now

    int scaleDiff = 4;
    int levelDiff = log2(scaleDiff);

    int LRwidth = input_comps_LR[0].width, LRheight = input_comps_LR[0].height;
    int HRwidth = input_comps_HR[0].width, HRheight = input_comps_HR[0].height;
    int Outwidth = LRwidth*scaleDiff, Outheight = LRheight*scaleDiff;
    /*my_image_comp *inHR_transformed = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        inHR_transformed[n].init(HRheight, HRwidth, HL_A97_LP*pow(2, waveletLevels)+HL_A97_HP);
        input_comps_HR[n].perform_boundary_extension_symmetric();
        resampleImageComponents(input_comps_HR+n, inHR_transformed+n, transform);
    }
    delete[] input_comps_HR;*/

    my_image_comp *inLR_wavelet_small = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        inLR_wavelet_small[n].init(LRheight, LRwidth, 0);
        input_comps_LR[n].perform_boundary_extension_symmetric();
        if(waveletType == WAVELET_53){
            analysis_5_3(input_comps_LR+n, inLR_wavelet_small+n, waveletLevels-levelDiff);
        } else {
            analysis_9_7(input_comps_LR+n, inLR_wavelet_small+n, waveletLevels-levelDiff);
        }
    }
    delete[] input_comps_LR;

    my_image_comp *inLR_wavelet = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        inLR_wavelet[n].init(Outheight, Outwidth, 0);
        increaseWaveletLevel(inLR_wavelet_small+n, inLR_wavelet+n, scaleDiff);
        inLR_wavelet[n].perform_boundary_extension_symmetric();
    }
    delete[] inLR_wavelet_small;

    my_image_comp *inHR_wavelet = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        //should be inHR_transformed but will be input_comps_HR until this is implemented
        inHR_wavelet[n].init(HRheight, HRwidth, 0);
        input_comps_HR[n].perform_boundary_extension_symmetric();
        if(waveletType == WAVELET_53){
            analysis_5_3(input_comps_HR+n, inHR_wavelet+n, waveletLevels);
        } else {
            analysis_9_7(input_comps_HR+n, inHR_wavelet+n, waveletLevels);
        }
    }
    delete[] input_comps_HR;

    my_image_comp *out_wavelet = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        out_wavelet[n].init(Outheight, Outwidth, HL_S97_HP*pow(2, waveletLevels)+HL_S97_LP);
        combineWavelets(inLR_wavelet+n, inHR_wavelet+n, out_wavelet+n, pow(2, waveletLevels), combineType, 99, 99);
    }
    delete[] inLR_wavelet;
    delete[] inHR_wavelet;

    my_image_comp *output_comps = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        output_comps[n].init(Outheight, Outwidth, 0); // Don't need a border for output
        out_wavelet[n].perform_boundary_extension_symmetric();
        if(waveletType == WAVELET_53){
            synthesis_5_3(out_wavelet+n, output_comps+n, waveletLevels);
        } else {
            synthesis_9_7(out_wavelet+n, output_comps+n, waveletLevels);
        }
    }
    delete[] out_wavelet;

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
    assert(out->height == inLR->height*LRspacing);
    assert(out->width == inLR->width*LRspacing);

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
                    out->buf[r*out->stride+c] = inLR->buf[(r/LRspacing)*inLR->stride+(c/LRspacing)];
                } else {
                    out->buf[r*out->stride+c] = 0;
                }
            } else {
                out->buf[r*out->stride+c] = inHR->buf[(r-Roffset)*inHR->stride+(c-Coffset)];//std::max(inLR->buf[r*inLR->stride+c], inHR->buf[r*inHR->stride+c]);
            }
        }
    }
}