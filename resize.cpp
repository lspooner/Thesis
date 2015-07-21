#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "resize.h"

#define PI 3.1415
#define HL 3

void enlargementFilter(my_image_comp *in, my_image_comp *out, int halfLength, int scale);
void reductionFilter(my_image_comp *in, my_image_comp *out, int halfLength, int scale);
float* makeSincFilter(int halfLength, float centre, float scaling);
float* makeWindow(int halfLength, float centre);
float* makeSinc(int halfLenght, float centre, float scaling);
float sinc(float n);
float windowFunction(float n, int tau);

int enlargeImage(char* inputFile, char* outputFile, int H, int scale){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(inputFile, &input_comps, H, &num_comps)) != 0){
        return err_code;
    }

    int width = input_comps[0].width, height = input_comps[0].height;

    // Allocate storage for the filtered output
    int outHeight = height*scale;
    int outWidth = width*scale;
    my_image_comp *output_comps = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        output_comps[n].init(outHeight, outWidth, 0); // Don't need a border for output
    }

    // Process the image, all in floating point (easy)
    for (int n=0; n < num_comps; n++){
       input_comps[n].perform_boundary_extension_symmetric();
    }

    for (int n=0; n < num_comps; n++){
        enlargementFilter(input_comps+n, output_comps+n, H, scale);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

int reduceImage(char* inputFile, char* outputFile, int H, int scale){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(inputFile, &input_comps, H, &num_comps)) != 0){
        return err_code;
    }
    //printf("read in bmp ok\n");

    int width = input_comps[0].width, height = input_comps[0].height;

    // Allocate storage for the filtered output
    int outHeight = height/scale;
    int outWidth = width/scale;
    my_image_comp *output_comps = new my_image_comp[num_comps];
    int n;
    for (n=0; n < num_comps; n++){
        output_comps[n].init(outHeight, outWidth, 0); // Don't need a border for output
    }

    // Process the image, all in floating point (easy)
    for (n=0; n < num_comps; n++){
       input_comps[n].perform_boundary_extension_symmetric();
    }
    for (n=0; n < num_comps; n++){
        reductionFilter(input_comps+n, output_comps+n, H, scale);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void printFArray(float *Array, int size){
    printf("[ ");
    for(int i = 0; i < size; i++){
        printf("%f ", Array[i]);
    }
    printf("]\n");
}

void enlargementFilter(my_image_comp *in, my_image_comp *out, int halfLength, int scale){
    // Check for consistent dimensions
    assert(in->border >= halfLength);
    assert((out->height >= in->height) && (out->width >= in->width));

    //make filters
    //int filterLength = 2*halfLength + 1;
    float **filterhandles = new float*[scale];
    float **filter = new float*[scale];

    for(int i = 0; i < scale; i++){
        filterhandles[i] = makeSincFilter(halfLength, ((float)i)/scale, 1);
        filter[i] = filterhandles[i]+halfLength;
        
        printf("i = %d, scale = %d, centre = %f ", i, scale, ((float)-i)/scale);
        printFArray(filterhandles[i], halfLength*2+1);
    }

    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels];
    float *tempPicBuf = tempPic + out->width*in->border;

    //first apply filter to the rows
    for(int r = -halfLength; r < in->height + halfLength; r++){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            int filternum = c%scale;
            int inputC = c/scale;
            if(filternum != 0){
                inputC++;
            }
            //printf("inputC = %d, filternum = %d\n", inputC, filternum);
            //printFArray(filter[filternum], halfLength*2+1);
            for(int n = -halfLength; n <= halfLength; n++){
                sum += filter[filternum][n] * in->buf[r*in->stride+inputC+n];
            }
            tempPicBuf[r*out->width + c] = sum;
        }
    }

    //then apply filter to the columns
    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            int inputR = r/scale;
            int filternum = r%scale;
            if(filternum != 0){
                inputR++;
            }
            //printf("inputR = %d, filternum = %d\n", inputR, filternum);
            //printFArray(filter[filternum], halfLength*2+1);
            for(int n = -halfLength; n <= halfLength; n++){
                sum += filter[filternum][n] * tempPicBuf[(inputR+n)*out->width+c];
            }
            out->buf[r*out->stride+c] = sum;
        }
    }

    delete[] tempPic;
    delete[] filter;
    for(int i = 0; i < scale; i++){
        delete[] filterhandles[i];
    }
    delete[] filterhandles;
}

//reduces the image by integer scale
void reductionFilter(my_image_comp *in, my_image_comp *out, int halfLength, int scale){
    // Check for consistent dimensions
    assert(in->border >= halfLength);
    assert((out->height <= in->height) && (out->width <= in->width));
   
    //make filters
    float *filterhandle = makeSincFilter(halfLength, 0.0, (float)scale);
    float *filter = filterhandle+halfLength;
    
    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels];
    float *tempPicBuf = tempPic + out->width*in->border;

    //first apply filter to the rows
    for(int r = -halfLength; r < in->height + halfLength; r++){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            int inputC = c*scale;
            for(int n = -halfLength; n <= halfLength; n++){
                sum += filter[n] * in->buf[r*in->stride+inputC+n];
            }
            tempPicBuf[r*out->width + c] = sum;
        }
    }

    //then apply filter to the columns
    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            float sum = 0.0F;
            int inputR = r*scale;
            for(int n = -halfLength; n <= halfLength; n++){
                sum += filter[n] * tempPicBuf[(inputR+n)*out->width+c];
            }
            out->buf[r*out->stride+c] = sum;
        }
    }

    delete[] tempPic;
    delete[] filterhandle;
}

float* makeSincFilter(int halfLength, float centre, float scaling){
    float *filter;
    if(halfLength == 0){
        filter = new float[1];
        filter[0] = 1;
    } else {
        int size = 2*halfLength + 1;
        filter = new float[size];
        float *window = makeWindow(halfLength, centre);
        float *sincInterp = makeSinc(halfLength, centre, scaling);
        float sum = 0;

        for(int i = 0; i < size; i++){
            filter[i] = window[i] * sincInterp[i];
            sum += filter[i];
        }

        if(sum != 1){
            for(int i = 0; i < size; i++){
                filter[i] /= sum;
            }
        }

        delete[] window;
        delete[] sincInterp;
    }

    return filter;
}

float* makeWindow(int halfLength, float centre){
    int size = 2*halfLength + 1;
    float *window = new float[size];
    int tau = size/2+1;
    int i;
    for(i = 0; i < size; i++){
        window[i] = windowFunction(i-tau+1+centre, tau);
    }
    return window;
}

float* makeSinc(int halfLength, float centre, float scaling){
    int size = 2*halfLength + 1;
    float *sincInterp = new float[size];
    int tau = size/2;
    int i;
    for(i = 0; i < size; i++){
        sincInterp[i] = sinc((i-tau + centre)/scaling);
    }
    return sincInterp;
}

float sinc(float n){
    if(n == 0){
        return 1;
    } else {
        return (sin(PI*n)/(PI*n));
    }
}

float windowFunction(float n, int tau){
    return 0.5*(1+cos((PI*n)/tau));
}
