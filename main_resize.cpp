/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <string.h>
#include "fixedPoint.h"
#include "io_bmp.h"
#include "image_comps.h"
#include "resize.h"

io_byte floatToByte(float num);
void printFloatArray(float *Array, int size);
void printUsage(char* fileName);

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[]){

    if(argc == 6 && (strcmp(argv[1], "-r") == 0 || strcmp(argv[1], "-e")==0)){
        
        int err_code=0;
        try {
            int scale;
            int H;

            sscanf(argv[4], "%d", &H);
            sscanf(argv[5], "%d", &scale);
            
            if(strcmp(argv[1], "-r") == 0){
                if((err_code = reduceImage(argv[2], argv[3], H, scale))){
                    throw err_code;
                }
            } else if(strcmp(argv[1], "-e") == 0){
                if((err_code = enlargeImage(argv[2], argv[3], H, scale))){
                    throw err_code;
                }
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

void printFloatArray(float *Array, int size){
    printf("[ ");
    for(int i = 0; i < size; i++){
        printf("%f ", Array[i]);
    }
    printf("]\n");
}

void printUsage(char* fileName){
    fprintf(stderr, "Usage: %s -r <in bmp file> <out bmp file> <H> <scale> (reduces image by integer scale using sinc interpolation)\n"
        "       %s -e <in bmp file> <out bmp file> <H> <scale> (enlarges image by integer scale using sinc interpolation)\n", fileName, fileName);
}
