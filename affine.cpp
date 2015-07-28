#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "affine.h"

using namespace std;
using namespace arma;

void resampleImageComponents(my_image_comp *in, my_image_comp *out, Mat<float> transform);
float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11);
float bilinear_interpolation_1D(float x, float x0_val, float x1_val);

int resampleImage(char* inputFile, char* outputFile, arma::Mat<float> transform){

    assert(transform.n_rows == 2);
    assert(transform.n_cols == 2);

    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    if ((err_code = readBMP(inputFile, &input_comps, 1, &num_comps)) != 0){
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
        resampleImageComponents(input_comps+n, output_comps+n, transform);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void resampleImageComponents(my_image_comp *in, my_image_comp *out, Mat<float> transform){
    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->stride; c++){
            Mat<float> output(2, 1);
            output(0, 0) = r;
            output(1, 0) = c;
            Mat<float> input = transform*output;

            float inR = input(0, 0);
            float inC = input(1, 0);

            if(inC > -1 && inC < out->width && inR > -1 && inR < out->height){
                out->buf[r*out->stride+c] = bilinear_interpolation_2D(inC-((int)inC), inR-((int)inR), in->buf[((int)inR)*in->stride+((int)inC)], in->buf[((int)inR)*in->stride+((int)inC+1)], in->buf[((int)inR+1)*in->stride+((int)inC)], in->buf[((int)inR+1)*in->stride+((int)inC+1)]);
            }
        }
    }
}

/*
perform 2D bilinear interpolation, 
(x,y) are the normalised coordinates, 
val_xy are the values at the corners

0________1_->x
 |
 |
 |
1|
 |
 \/
 y
*/

float bilinear_interpolation_2D(float x, float y, float val_00, float val_01, float val_10, float val_11){
    float y0 = bilinear_interpolation_1D(x, val_00, val_01);
    float y1 = bilinear_interpolation_1D(x, val_10, val_11);
    return bilinear_interpolation_1D(y, y0, y1);
}

float bilinear_interpolation_1D(float x, float x0_val, float x1_val){
    if(x == 0.0){
        return x0_val;
    } else {
        return (1-x)*x0_val + x*x1_val;
    }
}