#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "wavelet.h"

#define HL_LP_5_3 2
#define HL_HP_5_3 1

void analysis_5_3_single(my_image_comp *in, my_image_comp *out, int spacing, int offset);

int analysis_5_3(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_LP_5_3*scale+HL_HP_5_3, &num_comps)) != 0){
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
        for (int spacing=1; spacing < scale; spacing*=2){
            analysis_5_3_single(input_comps+n, output_comps+n, spacing, 0);
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

void analysis_5_3_single(my_image_comp *in, my_image_comp *out, int spacing, int offset){
    // Check for consistent dimensions
    assert(in->border >= HL_LP_5_3);
    assert((out->height == in->height) && (out->width == in->width));

    float LPfilterhandle[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
    float HPfilterhandle[3] = {-0.25, 0.5, -0.25};

    float *LPfilter = LPfilterhandle + HL_LP_5_3;
    float *HPfilter = HPfilterhandle + HL_HP_5_3;

    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels];
    float *tempPicBuf = tempPic + out->width*in->border;

    /*printf("\n\nInput Image:\n");
    for(int r = 0; r < in->height; r++){
        for(int c = 0; c < in->width; c++){
            printf("%3f ", in->buf[r*in->stride+c]);
        }
        printf("\n");
    }*/

    //apply lowpass filter to the rows
    for(int r = -HL_LP_5_3; r < in->height + HL_LP_5_3; r++){
        for(int c = offset; c < out->width; c+=2*spacing){
            float sum = 0.0F;
            for(int n = -HL_LP_5_3; n <= HL_LP_5_3; n++){
                sum += LPfilter[n] * in->buf[r*in->stride+c+spacing*n];
            }
            tempPicBuf[r*out->width + c] = sum;
        }
    }

    //apply highpass filter to the rows
    for(int r = -HL_LP_5_3; r < in->height + HL_LP_5_3; r++){
        for(int c = offset+spacing; c < out->width; c+=2*spacing){
            float sum = 0.0F;
            for(int n = -HL_HP_5_3; n <= HL_HP_5_3; n++){
                sum += HPfilter[n] * in->buf[r*in->stride+c+spacing*n];
            }
            tempPicBuf[r*out->width + c] = sum;
        }
    }

    /*printf("\nIntermediate Image:\n");
    for(int r = 0; r < in->height; r++){
        for(int c = 0; c < out->width; c++){
            printf("%3f ", tempPicBuf[r*out->width+c]);
        }
        printf("\n");
    }*/

    //apply lowpass filter to the columns
    for(int r = offset; r < out->height; r+=2*spacing){
        for(int c = 0; c < out->width; c++){
            //printf("\nLow passing row %d column %d\n", r, c);
            float sum = 0.0F;
            for(int n = -HL_LP_5_3; n <= HL_LP_5_3; n++){
                sum += LPfilter[n] * tempPicBuf[(r+spacing*n)*out->stride+c];
                //printf("%f * %f = %f\n", LPfilter[n], tempPicBuf[(r+spacing*n)*out->stride+c], LPfilter[n] * tempPicBuf[(r+spacing*n)*out->stride+c]);
            }
            out->buf[r*out->stride+c] = sum;
        }
    }

    //apply highpass filter to the columns
    for(int r = offset+spacing; r < out->height; r+=2*spacing){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            for(int n = -HL_HP_5_3; n <= HL_HP_5_3; n++){
                sum += HPfilter[n] * tempPicBuf[(r+spacing*n)*out->stride+c];
            }
            out->buf[r*out->stride+c] = sum;
        }
    }

    /*printf("\nOutput Image:\n");
    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            printf("%3f ", out->buf[r*out->stride+c]);
        }
        printf("\n");
    }*/

    delete[] tempPic;
}