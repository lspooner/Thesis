/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <string.h>
#include "io_bmp.h"
#include "image_comps.h"
#include "wavelet.h"

void printUsage(char* fileName);

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int main(int argc, char *argv[]){

    if(argc == 5 && (strcmp(argv[1], "-a53") == 0 || strcmp(argv[1], "-a97")==0 || strcmp(argv[1], "-s53")==0 || strcmp(argv[1], "-s97")==0)){
        
        int err_code=0;
        try {
            int levels;

            sscanf(argv[4], "%d", &levels);

            if(strcmp(argv[1], "-a53") == 0){
                if((err_code = analysis_5_3(argv[2], argv[3], levels))){
                    throw err_code;
                }
            } else if(strcmp(argv[1], "-a97") == 0){
                printf("Not implemented yet\n");
            } else if(strcmp(argv[1], "-s53") == 0){
                printf("Not implemented yet\n");
            } else if(strcmp(argv[1], "-s97") == 0){
                printf("Not implemented yet\n");
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
    fprintf(stderr, "Usage: %s -a53 <in bmp file> <out bmp file> <levels>\n"
        "       %s -a97 <in bmp file> <out bmp file> <levels>\n"
        "       %s -s53 <in bmp file> <out bmp file> <levels>\n"
        "       %s -s97 <in bmp file> <out bmp file> <levels>\n", fileName, fileName, fileName, fileName);
}
