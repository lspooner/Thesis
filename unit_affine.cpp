#include <string.h>
#include <assert.h>
#include <iostream>
#include "io_bmp.h"
#include "image_comps.h"
#include "affine_unit.h"

#define EPSILON 0.00001

bool CompareFloats (float A, float B);

int main(int argc, char *argv[]){

    printf("Testing Affine transform functionality:\n");

    printf("Testing transform generator\n");
    float LRpoints[6] = {0, 0, 0, 1, 1, 1};
    float HRpoints[6] = {0, 0, 0, 1, 1, 1};

    Mat<float> transform = computeAffineTransform(LRpoints, HRpoints);

    assert(transform(0,0) == 1);
    assert(transform(0,1) == 0);
    assert(transform(0,2) == 0);
    assert(transform(1,0) == 0);
    assert(transform(1,1) == 1);
    assert(transform(1,2) == 0);

    float LRpoints2[6] = {0, 0, 0, 1, 1, 0};
    float HRpoints2[6] = {0.5, 0, 0, 1, 1.5, 0};

    transform = computeAffineTransform(LRpoints2, HRpoints2);

    assert(transform(0,0) == 1);
    assert(transform(0,1) == 0.5);
    assert(transform(0,2) == -0.5);
    assert(transform(1,0) == 0);
    assert(transform(1,1) == 1);
    assert(transform(1,2) == 0);

    printf("Test Passed!\n");

    printf("Testing Bilinear transform\n");
    assert(CompareFloats(bilinear_interpolation_2D(0.5, 0.2, 91, 210, 162, 95), 146.1));
    assert(CompareFloats(bilinear_interpolation_2D(0, 0.4, 5, 6, 10, 9), 7));
    assert(CompareFloats(bilinear_interpolation_2D(0.4, 0, 5, 10, 10, 5), 7));
    printf("Test Passed!\n");

    printf("Testing scaling transform\n");

    Mat<float> sTransform(2, 2);
    sTransform(0,0) = 0.5;
    sTransform(0,1) = 0;
    sTransform(1,0) = 0;
    sTransform(1,1) = 0.5;

    sTransform = scaleTransform(sTransform);

    assert(sTransform(0,0) == 1);
    assert(sTransform(0,1) == 0);
    assert(sTransform(1,0) == 0);
    assert(sTransform(1,1) == 1);

    printf("Test Passed!\n");

    printf("Generating 8x8 image, testing it's best match with itself is the identity matrix\n");

    int boundary = 0;
    my_image_comp *image = new my_image_comp;
    image->init(8, 8, boundary);

    float v[8] = {(0-128)/128.0, (36-128)/128.0, (73-128)/128.0, (109-128)/128.0, (146-128)/128.0, (182-128)/128.0, (219-128)/128.0, (255-128)/128.0};
    
    float buf[64] = {v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
                     v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[6],
                     v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5],
                     v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4],
                     v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3],
                     v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2],
                     v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[1],
                     v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0]};

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            image->buf[r*image->stride + c] = buf[r*8+c];
        }
    }

    transform = matchImages(image, image, LRpoints, HRpoints);
    assert(transform(0,0) == 1);
    assert(transform(0,1) == 0);
    assert(transform(0,2) == 0);
    assert(transform(1,0) == 0);
    assert(transform(1,1) == 1);
    assert(transform(1,2) == 0);

    printf("Test Passed!\n");

    printf("ALL TESTS PASSED! YOU ARE AWESOME!\n");
}

bool CompareFloats (float A, float B){
   float diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}