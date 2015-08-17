#include <string.h>
#include <assert.h>
#include <iostream>
#include "io_bmp.h"
#include "image_comps.h"
#include "affine_unit.h"
#include "resize.h"

#define EPSILON 0.00001

bool CompareFloats (float A, float B);

int main(int argc, char *argv[]){

    printf("Testing Affine transform functionality:\n");

    printf("Testing transform generator\n");
    float LRpoints[6] = {0, 0, 0, 2, 2, 2};
    float HRpoints[6] = {0, 0, 0, 2, 2, 2};

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
    sTransform(0,0) = 2;
    sTransform(0,1) = 0;
    sTransform(1,0) = 0;
    sTransform(1,1) = 2;

    int scaleDiff;

    sTransform = scaleTransform(sTransform, &scaleDiff);

    assert(sTransform(0,0) == 1);
    assert(sTransform(0,1) == 0);
    assert(sTransform(1,0) == 0);
    assert(sTransform(1,1) == 1);

    assert(scaleDiff == 2);

    printf("Test Passed!\n");

    printf("Generating 8x8 image, testing it's best match with itself is the identity matrix\n");

    int boundary = 2;
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

    /*printf("Generating 16x16 image and 14x14 image offset by 1 and checking affine match\n");

    float buf16[16*16] = {-1.0000000, -0.8671875, -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,
                          -0.8671875, -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,
                          -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,
                          -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,
                          -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,
                          -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,
                          -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,
                          -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000,
                           0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125,
                           0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250,
                           0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375,
                           0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500,
                           0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500, -0.6015625,
                           0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500, -0.6015625, -0.7343750,
                           0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500, -0.6015625, -0.7343750, -0.8671875,
                           0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500, -0.6015625, -0.7343750, -0.8671875, -1.0000000};

    float buf14[14*14] = {-0.8671875, -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,
                          -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,
                          -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,
                          -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,
                          -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,
                          -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,
                          -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,
                           0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,
                           0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000,
                           0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125,
                           0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250,
                           0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375,
                           0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500,
                           0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500, -0.6015625};

    my_image_comp *image16 = new my_image_comp;
    image16->init(16, 16, boundary);

    for(int r = 0; r < 16; r++){
        for(int c = 0; c < 16; c++){
            image16->buf[r*image16->stride + c] = buf16[r*16+c];
        }
    }

    my_image_comp *image14 = new my_image_comp;
    image14->init(14, 14, boundary);

    for(int r = 0; r < 14; r++){
        for(int c = 0; c < 14; c++){
            image14->buf[r*image14->stride + c] = buf14[r*14+c];
        }
    }

    //float LRpoints2[6] = {0, 0, 0, 1, 1, 0};
    //float HRpoints2[6] = {0, 0.5, 0, 1, 1.5, 0};

    transform = matchImages(image16, image14, LRpoints, HRpoints);
    cout << transform << endl;
    assert(transform(0,0) == 1);
    assert(transform(0,1) == 0);
    assert(transform(0,2) == 0);
    assert(transform(1,0) == 0);
    assert(transform(1,1) == 1);
    assert(transform(1,2) == 1);

    printf("Test Passed!\n");*/


    printf("Testing image matching on images of lenna\n");

    my_image_comp *LR_lenna = NULL;
    int num_comps;
    assert(readBMP((char*)"images/lenna_r4_h40.bmp", &LR_lenna, 2, &num_comps)==0);

    my_image_comp *HR_lenna = NULL;
    assert(readBMP((char*)"images/lenna_mono_eye.bmp", &HR_lenna, 40, &num_comps)==0);

    float LRpoints_lenna[6] = {50, 50, 50, 54, 54, 50};
    float HRpoints_lenna[6] = {0, 0, 0, 1, 1, 0};

    my_image_comp *HR_lenna_reduced = new my_image_comp;
    HR_lenna_reduced->init(HR_lenna->height/4, HR_lenna->width/4, 0);
    HR_lenna_reduced->perform_boundary_extension_symmetric();
    reductionFilter(HR_lenna, HR_lenna_reduced, 40, 4);
    // Write the image back out again
    assert(outputBMP((char*)"images/lenna_eye_r4_h40.bmp", HR_lenna_reduced, 1) == 0);

    my_image_comp *LR_lenna_subset = new my_image_comp;
    LR_lenna_subset->init(HR_lenna->height/4, HR_lenna->width/4, 0);
    for(int r = 0; r < LR_lenna_subset->height; r++){
        for(int c = 0; c < LR_lenna_subset->width; c++){
            LR_lenna_subset->buf[r*LR_lenna_subset->stride+c] = LR_lenna->buf[(r+50)*LR_lenna->stride+(c+50)];
        }
    }
    // Write the image back out again
    assert(outputBMP((char*)"images/lenna_r4_h40_subset.bmp", LR_lenna_subset, 1) == 0);


    transform = matchImages(LR_lenna, HR_lenna, LRpoints_lenna, HRpoints_lenna);
    cout << transform << endl;
    assert(transform(0,0) == 4);
    assert(transform(0,1) == 0);
    assert(transform(0,2) == 50);
    assert(transform(1,0) == 0);
    assert(transform(1,1) == 4);
    assert(transform(1,2) == 50);

    printf("Test Passed!\n");

    printf("ALL TESTS PASSED! YOU ARE AWESOME!\n");
}

bool CompareFloats (float A, float B){
   float diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}