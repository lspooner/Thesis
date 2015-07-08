/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "io_bmp.h"
#include "image_comps.h"

void printUsage(char* fileName);
int differenceImage(char* inputFile1, char* inputFile2, char* outputFile);

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[]){

    if(argc == 5){
        
        int err_code=0;
        try {
            if((err_code = differenceImage(argv[2], argv[3], argv[4]))){
                throw err_code;
            }
        } catch (int exc) {
            if (exc == IO_ERR_NO_FILE){
                fprintf(stderr,"Cannot open supplied input or output file.\n");
            } else if (exc == IO_ERR_FILE_HEADER){
                fprintf(stderr,"Error encountered while parsing BMP file header.\n");
            } else if (exc == IO_ERR_UNSUPPORTED){
                fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                    "simple example supports only 8-bit and 24-bit data.\n");
            } else if (exc == IO_ERR_FILE_TRUNC){
                fprintf(stderr,"Input or output file truncated unexpectedly.\n");
            } else if (exc == IO_ERR_FILE_NOT_OPEN){
                fprintf(stderr,"Trying to access a file which is not open!(?)\n");
            } else {
                fprintf(stderr, "Unknown exception occured\n");
            }
            return -1;
        }
        return 0;
    } else {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    } 
}

void printUsage(char* fileName){
    fprintf(stderr, "Usage: %s <in bmp file 1> <in bmp file 2> <out bmp file> (finds the difference between 2 input files)\n", fileName);
}

int differenceImage(char* inputFile1, char* inputFile2, char* outputFile){
    // Read the input image1
    int err_code;
    my_image_comp *input_comps1 = NULL;
    int num_comps;
    if ((err_code = readBMP(inputFile1, &input_comps1, 0, &num_comps)) != 0){
        return err_code;
    }

    // Read the input image2
    my_image_comp *input_comps2 = NULL;
    if ((err_code = readBMP(inputFile2, &input_comps2, 0, &num_comps)) != 0){
        return err_code;
    }

    int width = input_comps1[0].width, height = input_comps1[0].height;

    // Allocate storage for the filtered output
    my_image_comp *output_comps = new my_image_comp[num_comps];
    for (int n=0; n < num_comps; n++){
        output_comps[n].init(height, width, 0); // Don't need a border for output
    }

    float ME = 0;
    float MSE = 0;
    //all images read in, time to find the difference
    for(int n = 0; n < num_comps; n++){
        my_image_comp *in1 = &input_comps1[n];
        my_image_comp *in2 = &input_comps2[n];
        my_image_comp *out = &output_comps[n];
        for(int r = 0; r < output_comps[n].height; r++){
            for(int c = 0; c < output_comps[n].width; c++){
                float diff = in1->buf[r*in1->stride+c] - in2->buf[r*in2->stride+c];
                ME += diff;
                MSE += (diff*diff);
                out->buf[r*out->stride+c] = 0.5*diff+DIFF_OFFSET;
            }
        }
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    //print stuff
    ME *= 1.0/(width*height);
    MSE *= 1.0/(width*height);
    float PSNR = 10*log10(255.0*255.0/MSE);

    printf("Mean Error = %f\n"
           "Mean Squared Error = %f\n"
           "Peak Signal to Noise Ratio = %f\n", ME, MSE, PSNR);

    delete[] input_comps1;
    delete[] input_comps2;
    delete[] output_comps;

    return 0;
}