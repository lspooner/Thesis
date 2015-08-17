/*****************************************************************************/
// File: image_comps.h
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include <armadillo>

#define INVALID 42

using namespace std;
using namespace arma;

/*****************************************************************************/
/* STRUCT                        my_image_comp                               */
/*****************************************************************************/

struct my_image_comp {
    // Data members: (these occupy space in the structure's block of memory)
    int width;
    int height;
    int stride; // width with border
    int border; // Extra rows/cols to leave around the boundary
    float *handle; // Points to start of allocated memory buffer
    float *buf; // Points to the first real image sample
    
    // Function members: (these do not occupy any space in memory)
    my_image_comp(){ 
        width = height = stride = border = 0;
        handle = buf = NULL;
    }

    ~my_image_comp(){ 
        if (handle != NULL){
            delete[] handle;
        }
    }

    void init(int height, int width, int border);

    void perform_boundary_extension();
    void perform_boundary_extension_zero_padding();
    void perform_boundary_extension_symmetric();
    void perform_boundary_extension_wavelet(int spacing, Mat<float> offset);
};

int readBMP(char* fileName, my_image_comp** in, int border, int* num_comps);
int outputBMP(char* fileName, my_image_comp* out_comps, int num_comps);