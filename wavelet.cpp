#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "wavelet.h"

#define HL_A53_LP 2
#define HL_A53_HP 1
#define HL_A97_LP 4
#define HL_A97_HP 3

#define HL_S53_LP 1
#define HL_S53_HP 2
#define HL_S97_LP 3
#define HL_S97_HP 4

void analysis(my_image_comp *in, my_image_comp *out, int spacing, int offset, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL);

int analysis_5_3(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_A53_LP*scale+HL_A53_HP, &num_comps)) != 0){
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

    float LPfilterhandle[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
    float HPfilterhandle[3] = {-0.25, 0.5, -0.25};

    float *LPfilter = LPfilterhandle + HL_A53_LP;
    float *HPfilter = HPfilterhandle + HL_A53_HP;

    for(int n=0; n < num_comps; n++){
        for (int spacing=1; spacing < scale; spacing*=2){
            //analysis_5_3_single(input_comps+n, output_comps+n, spacing, 0);
            analysis(input_comps+n, output_comps+n, spacing, 0, LPfilter, HPfilter, HL_A53_LP, HL_A53_HP);
        }
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

int analysis_9_7(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_A97_LP*scale+HL_A97_HP, &num_comps)) != 0){
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

    float LPfilterhandle[9] = {0.026749, -0.016864, -0.078223, 0.266864, 0.602949, 0.266864, -0.078223, -0.016864, 0.026749};
    float HPfilterhandle[7] = {0.045636, -0.028772, -0.295636, 0.557543, -0.295636, -0.028772, 0.045636};

    float *LPfilter = LPfilterhandle + HL_A97_LP;
    float *HPfilter = HPfilterhandle + HL_A97_HP;

    for(int n=0; n < num_comps; n++){
        for (int spacing=1; spacing < scale; spacing*=2){
            //analysis_5_3_single(input_comps+n, output_comps+n, spacing, 0);
            analysis(input_comps+n, output_comps+n, spacing, 0, LPfilter, HPfilter, HL_A97_LP, HL_A97_HP);
        }
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

int synthesis_5_3(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_S53_LP*scale+HL_S53_HP, &num_comps)) != 0){
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

    float LPfilterhandle[3] = {2*0.25, 2*0.5, 2*0.25};
    float HPfilterhandle[5] = {-2*0.125, -2*0.25, 2*0.75, -2*0.25, -2*0.125};

    float *LPfilter = LPfilterhandle + HL_S53_LP;
    float *HPfilter = HPfilterhandle + HL_S53_HP;

    for(int n=0; n < num_comps; n++){
        for (int spacing=scale/2; spacing >= 1; spacing/=2){
            //analysis_5_3_single(input_comps+n, output_comps+n, spacing, 0);
            analysis(input_comps+n, output_comps+n, spacing, 0, LPfilter, HPfilter, HL_S53_LP, HL_S53_HP);
        }
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

int synthesis_9_7(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_S97_LP*scale+HL_S97_HP, &num_comps)) != 0){
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

    float LPfilterhandle[7] = {-2*0.045636, -2*0.028772, 2*0.295636, 2*0.557543, 2*0.295636, -2*0.028772, -2*0.045636};
    float HPfilterhandle[9] = {2*0.026749, 2*0.016864, -2*0.078223, -2*0.266864, 2*0.602949, -2*0.266864, -2*0.078223, 2*0.016864, 2*0.026749};

    float *LPfilter = LPfilterhandle + HL_S97_LP;
    float *HPfilter = HPfilterhandle + HL_S97_HP;

    for(int n=0; n < num_comps; n++){
        for (int spacing=scale/2; spacing >= 1; spacing/=2){
            //analysis_5_3_single(input_comps+n, output_comps+n, spacing, 0);
            analysis(input_comps+n, output_comps+n, spacing, 0, LPfilter, HPfilter, HL_S97_LP, HL_S97_HP);
        }
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void analysis(my_image_comp *in, my_image_comp *out, int spacing, int offset, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL){
    // Check for consistent dimensions
    assert(in->border >= LP_HL*spacing+HP_HL);
    assert((out->height == in->height) && (out->width == in->width));

    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels];
    float *tempPicBuf = tempPic + out->width*in->border;

    //apply lowpass filter to the rows
    for(int r = -LP_HL; r < in->height + LP_HL; r++){
        for(int c = offset; c < out->width; c+=2*spacing){
            float sum = 0.0F;
            for(int n = -LP_HL; n <= LP_HL; n++){
                sum += LPfilter[n] * in->buf[r*in->stride+c+spacing*n];
            }
            tempPicBuf[r*out->width + c] = sum;
        }
    }

    //apply highpass filter to the rows
    for(int r = -LP_HL; r < in->height + LP_HL; r++){
        for(int c = offset+spacing; c < out->width; c+=2*spacing){
            float sum = 0.0F;
            for(int n = -HP_HL; n <= HP_HL; n++){
                sum += HPfilter[n] * in->buf[r*in->stride+c+spacing*n];
            }
            tempPicBuf[r*out->width + c] = sum;
        }
    }

    //apply lowpass filter to the columns
    for(int r = offset; r < out->height; r+=2*spacing){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            for(int n = -LP_HL; n <= LP_HL; n++){
                sum += LPfilter[n] * tempPicBuf[(r+spacing*n)*out->stride+c];
            }
            out->buf[r*out->stride+c] = sum;
        }
    }

    //apply highpass filter to the columns
    for(int r = offset+spacing; r < out->height; r+=2*spacing){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            for(int n = -HP_HL; n <= HP_HL; n++){
                sum += HPfilter[n] * tempPicBuf[(r+spacing*n)*out->stride+c];
            }
            out->buf[r*out->stride+c] = sum;
        }
    }

    delete[] tempPic;
}

void synthesis(my_image_comp *in, my_image_comp *out, int spacing, int offset, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL){
    // Check for consistent dimensions
    assert(in->border >= LP_HL*spacing+HP_HL);
    assert((out->height == in->height) && (out->width == in->width));

    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels]();
    float *tempPicBuf = tempPic + out->width*in->border;

    //apply lowpass filter to the rows
    for(int r = -LP_HL; r < in->height + LP_HL; r++){
        for(int n = -LP_HL; n <= LP_HL; n++){
            for(int c = offset+(spacing*n); c < out->width; c+=2*spacing){
            //TODO: work out how the fuck this indexing works
                if(c+n >=offset && c+n < out->width){
                    tempPicBuf[r*out->width + c + n] += LPfilter[n] * in->buf[r*in->stride+c+spacing*n];
                }
            }
        }
    }

    //apply highpass filter to the rows
    for(int r = -LP_HL; r < in->height + LP_HL; r++){
        for(int n = -HP_HL; n <= HP_HL; n++){
            for(int c = offset+spacing+(spacing*n); c < out->width; c+=2*spacing){
            //TODO: work out how the fuck this indexing works
                if(c+n >=offset && c+n < out->width){
                    tempPicBuf[r*out->width + c + n] += HPfilter[n] * in->buf[r*in->stride+c+spacing*n];
                }
            }
        }
    }

    //apply lowpass filter to the columns
    for(int c = 0; c < out->width; c++){
        for(int n = -LP_HL; n <= LP_HL; n++){
            for(int r = offset+(spacing*n); r < out->height; r+=2*spacing){
            //TODO: work out how the fuck this indexing works
                if(r+n >=offset && r+n < out->height){
                    out->buf[(r+n)*out->stride + c] += LPfilter[n] * in->buf[(r+spacing*n)*in->stride+c];
                }
            }
        }
    }

    //apply highpass filter to the columns
    for(int c = 0; c < out->width; c++){
        for(int n = -HP_HL; n <= HP_HL; n++){
            for(int r = offset+(spacing*n); r < out->height; r+=2*spacing){
            //TODO: work out how the fuck this indexing works
                if(r+n >=offset && r+n < out->height){
                    out->buf[(r+n)*out->stride + c] += HPfilter[n] * in->buf[(r+spacing*n)*in->stride+c];
                }
            }
        }
    }

    delete[] tempPic;
}