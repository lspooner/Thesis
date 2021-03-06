/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <string.h>
#include <math.h>
#include "io_bmp.h"
#include "image_comps.h"
#include "combine.h"

void printUsage(char* fileName);

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[]){

    if((argc == 9 && strcmp(argv[1], "-ci") == 0) || (argc == 9 && strcmp(argv[1], "-cw")==0)){
        
        int err_code=0;
        try {
            if(strcmp(argv[1], "-ci") == 0){
                int waveletType;
                int combineType;
                int waveletLevels;
                float lamda;

                sscanf(argv[5], "%d", &waveletType);
                sscanf(argv[6], "%d", &combineType);
                sscanf(argv[7], "%d", &waveletLevels);
                sscanf(argv[8], "%f", &lamda);

                if((err_code = combineImages(argv[2], argv[3], argv[4], waveletType, combineType, waveletLevels, lamda))){
                    throw err_code;
                }
            } else {
                int LRspacing;
                int method;
                int Roffset;
                int Coffset;

                sscanf(argv[5], "%d", &LRspacing);
                sscanf(argv[6], "%d", &method);
                sscanf(argv[7], "%d", &Roffset);
                sscanf(argv[8], "%d", &Coffset);

                if((err_code = combineWavelets(argv[2], argv[3], argv[4], LRspacing, method, Roffset, Coffset))){
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

void printUsage(char* fileName){
    fprintf(stderr, "Usage: %s -ci <in LR bmp file> <in HR bmp file> <out bmp file> waveletType combineType waveletLevels\n"
                    "       %s -cw <in LR bmp file> <in HR bmp file> <out bmp file> LRspacing combineType Roffset Coffset"
                    "\n"
                    "waveletType: 0 = CDF5/3, 1 = CDF9/7\n"
                    "combineType: 0 = max coefficient\n"
                    "             1 = max coefficients, except LR\n"
                    "             2 = max coefficients, average LR\n"
                    "             3 = max coefficients, blend LR\n"
                    "             4 = blend all coefficients\n", fileName, fileName);
}
