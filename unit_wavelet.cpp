#include <string.h>
#include <assert.h>
#include "io_bmp.h"
#include "image_comps.h"
#include "wavelet_unit.h"

#define EPSILON 0.00001

bool CompareFloats (float A, float B);

int main(int argc, char *argv[]){

    printf("Testing Input Output functionality:\n");
    printf("Initialising 8x8 image\n");

    int boundary = 5;
    my_image_comp *image = new my_image_comp;
    image->init(8, 8, boundary);

    assert(image->width == 8);
    assert(image->height == 8);
    assert(image->border == boundary);

    printf("Test Passed!\n");

    printf("Generating Buffer\n");
    float v[8] = {(0-128)/128.0, (36-128)/128.0, (73-128)/128.0, (109-128)/128.0, (146-128)/128.0, (182-128)/128.0, (219-128)/128.0, (255-128)/128.0};
    
    float buf[64] = {v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
                     v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[6],
                     v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5],
                     v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4],
                     v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3],
                     v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2],
                     v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[1],
                     v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0]};

    printf("Writing buffer to image\n");
    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            image->buf[r*image->stride + c] = buf[r*8+c];
        }
    }

    char* filename = (char*)"images/square_unit.bmp";
    printf("Writing image to %s\n", filename);
    int err_code = outputBMP(filename, image, 1);
    assert(err_code == 0);

    printf("Test Passed!\n");

    printf("Reading Image from %s\n", filename);
    delete image;
    int num_comps;
    err_code = readBMP(filename, &image, boundary, &num_comps);
    assert(err_code == 0);

    assert(num_comps == 1);
    assert(image->width == 8);
    assert(image->height == 8);
    assert(image->border == boundary);

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            assert(image->buf[r*image->stride + c] == buf[r*8+c]);
        }
    }

    printf("Test Passed!\n");

    printf("Performing boundary symmetric extension\n");
    float buf_extended[(8+2*boundary)*(8+2*boundary)] = {
        v[4], v[5], v[6], v[7], v[6], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[3], v[4], v[5], v[6], v[7],
        v[5], v[6], v[7], v[6], v[5], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[4], v[5], v[6], v[7], v[6],
        v[6], v[7], v[6], v[5], v[4], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[5], v[6], v[7], v[6], v[5],
        v[7], v[6], v[5], v[4], v[3], v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[6], v[7], v[6], v[5], v[4],
        v[6], v[5], v[4], v[3], v[2], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[7], v[6], v[5], v[4], v[3],
        v[5], v[4], v[3], v[2], v[1], v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2],
        v[6], v[5], v[4], v[3], v[2], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[7], v[6], v[5], v[4], v[3],
        v[7], v[6], v[5], v[4], v[3], v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[6], v[7], v[6], v[5], v[4],
        v[6], v[7], v[6], v[5], v[4], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[5], v[6], v[7], v[6], v[5],
        v[5], v[6], v[7], v[6], v[5], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[4], v[5], v[6], v[7], v[6],
        v[4], v[5], v[6], v[7], v[6], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[3], v[4], v[5], v[6], v[7],
        v[3], v[4], v[5], v[6], v[7], v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[2], v[3], v[4], v[5], v[6],
        v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0], v[1], v[2], v[3], v[4], v[5],
        v[3], v[4], v[5], v[6], v[7], v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[2], v[3], v[4], v[5], v[6],
        v[4], v[5], v[6], v[7], v[6], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[2], v[3], v[4], v[5], v[6], v[7],
        v[5], v[6], v[7], v[6], v[5], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[3], v[4], v[5], v[6], v[7], v[6],
        v[6], v[7], v[6], v[5], v[4], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[4], v[5], v[6], v[7], v[6], v[5],
        v[7], v[6], v[5], v[4], v[3], v[2], v[3], v[4], v[5], v[6], v[7], v[6], v[5], v[6], v[7], v[6], v[5], v[4]};

    Mat<float> offset(2, 1);
    offset(0,0) = 0;
    offset(1,0) = 0;

    image->perform_boundary_extension_wavelet(1, offset);

    for(int r = 0; r < 8+2*boundary; r++){
        for(int c = 0; c < 8+2*boundary; c++){
            assert(image->handle[r*image->stride + c] == buf_extended[r*(8+2*boundary)+c]);
        }
    }
    printf("Test Passed!\n");

    printf("Performing single wavelet 5_3 analysis with general analysis function\n");

    float LPfilterhandle_a53[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
    float HPfilterhandle_a53[3] = {-0.25, 0.5, -0.25};
    float *LPfilter_a53 = LPfilterhandle_a53 + HL_A53_LP;
    float *HPfilter_a53 = HPfilterhandle_a53 + HL_A53_HP;

    float buf_a53_l1[8*8] = {-1.00195312500,  0.00000000000, -0.43164062500,  0.00000000000,  0.12084960940, -0.03564453125,  0.86938476560,  0.0712890625,
                              0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.08911132813, -0.0693359375,
                             -0.43164062500,  0.00000000000,  0.12976074220, -0.01782226563,  0.82482910160,  0.10693359380,  0.74462890630, -0.1782226563,
                              0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.0019531250,
                              0.12084960940, -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.05847167969, -0.1425781250,
                             -0.03564453125, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.0019531250,
                              0.86938476560,  0.08911132813,  0.74462890630, -0.01782226563,  0.05847167969,  0.00000000000, -0.57421875000, -0.1425781250,
                              0.07128906250, -0.06933593750, -0.17822265630,  0.00195312500, -0.14257812500,  0.00195312500, -0.14257812500,  0.0019531250};

    my_image_comp *image_a53_l1 = new my_image_comp;
    image_a53_l1->init(8, 8, boundary);

    analysis(image, image_a53_l1, 1, offset, LPfilter_a53, HPfilter_a53, HL_A53_LP, HL_A53_HP);

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            assert(CompareFloats(image_a53_l1->buf[r*image_a53_l1->stride+c], buf_a53_l1[r*8+c]));
        }
    }
    printf("Test Passed!\n");

    printf("Performing second level wavelet 5_3 analysis with general analysis function\n");
    float buf_a53_l2[8*8] = {-1.05486297600,  0.00000000000, -0.05736541748,  0.00000000000,  0.44638442990, -0.03564453125,  0.35644531250,  0.0712890625,
                              0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.08911132813, -0.0693359375,
                             -0.05736541748,  0.00000000000, -0.06182098389, -0.01782226563,  0.13617324830,  0.10693359380, -0.01782226563, -0.1782226563,
                              0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.0019531250,
                              0.44638442990, -0.01782226563,  0.13617324830,  0.10693359380,  0.60163307190, -0.01782226563, -0.37538146970, -0.1425781250,
                             -0.03564453125, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.0019531250,
                              0.35644531250,  0.08911132813, -0.01782226563, -0.01782226563, -0.37538146970,  0.00000000000,  0.03341674805, -0.1425781250,
                              0.07128906250, -0.06933593750, -0.17822265630,  0.00195312500, -0.14257812500,  0.00195312500, -0.14257812500,  0.0019531250};

    image_a53_l1->perform_boundary_extension_wavelet(2, offset);

    float buf_a53_l1_extended[(8+2*boundary)*(8+2*boundary)] = {
         0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.03564453125, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.0019531250, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.03564453125,
        -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.12084960940, -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.05847167969, -0.1425781250,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.12084960940,
        -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.0019531250,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,
         0.10693359380,  0.82482910160, -0.01782226563,  0.12976074220,  0.00000000000, -0.43164062500,  0.00000000000,  0.12976074220, -0.01782226563,  0.82482910160,  0.10693359380,  0.74462890630, -0.1782226563,  0.82482910160, -0.01782226563,  0.12976074220,  0.00000000000, -0.43164062500,
        -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.08911132813, -0.0693359375, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,
        -0.03564453125,  0.12084960940,  0.00000000000, -0.43164062500,  0.00000000000, -1.00195312500,  0.00000000000, -0.43164062500,  0.00000000000,  0.12084960940, -0.03564453125,  0.86938476560,  0.0712890625,  0.12084960940,  0.00000000000, -0.43164062500,  0.00000000000, -1.00195312500,
        -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.08911132813, -0.0693359375, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,
         0.10693359380,  0.82482910160, -0.01782226563,  0.12976074220,  0.00000000000, -0.43164062500,  0.00000000000,  0.12976074220, -0.01782226563,  0.82482910160,  0.10693359380,  0.74462890630, -0.1782226563,  0.82482910160, -0.01782226563,  0.12976074220,  0.00000000000, -0.43164062500,
        -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.0019531250,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,
        -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.12084960940, -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.05847167969, -0.1425781250,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.12084960940,
         0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.03564453125, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.0019531250, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.03564453125,
         0.00000000000,  0.05847167969, -0.01782226563,  0.74462890630,  0.08911132813,  0.86938476560,  0.08911132813,  0.74462890630, -0.01782226563,  0.05847167969,  0.00000000000, -0.57421875000, -0.1425781250,  0.05847167969, -0.01782226563,  0.74462890630,  0.08911132813,  0.86938476560,
         0.00195312500, -0.14257812500,  0.00195312500, -0.17822265630, -0.06933593750,  0.07128906250, -0.06933593750, -0.17822265630,  0.00195312500, -0.14257812500,  0.00195312500, -0.14257812500,  0.0019531250, -0.14257812500,  0.00195312500, -0.17822265630, -0.06933593750,  0.07128906250,
        -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.12084960940, -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.05847167969, -0.1425781250,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.12084960940,
        -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.0019531250,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,
         0.10693359380,  0.82482910160, -0.01782226563,  0.12976074220,  0.00000000000, -0.43164062500,  0.00000000000,  0.12976074220, -0.01782226563,  0.82482910160,  0.10693359380,  0.74462890630, -0.1782226563,  0.82482910160, -0.01782226563,  0.12976074220,  0.00000000000, -0.43164062500,
        -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.08911132813, -0.0693359375, -0.01782226563,  0.00195312500,  0.00000000000,  0.00195312500,  0.00000000000,
        -0.03564453125,  0.12084960940,  0.00000000000, -0.43164062500,  0.00000000000, -1.00195312500,  0.00000000000, -0.43164062500,  0.00000000000,  0.12084960940, -0.03564453125,  0.86938476560,  0.0712890625,  0.12084960940,  0.00000000000, -0.43164062500,  0.00000000000, -1.00195312500};

    for(int r = 0; r < 8+2*boundary; r++){
        for(int c = 0; c < 8+2*boundary; c++){
            assert(CompareFloats(image_a53_l1->handle[r*image_a53_l1->stride+c], buf_a53_l1_extended[r*(8+2*boundary)+c]));
        }
    }

    my_image_comp *image_a53_l2 = new my_image_comp;
    image_a53_l2->init(8, 8, boundary);

    analysis(image_a53_l1, image_a53_l2, 2, offset, LPfilter_a53, HPfilter_a53, HL_A53_LP, HL_A53_HP);

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            assert(CompareFloats(image_a53_l2->buf[r*image_a53_l2->stride+c], buf_a53_l2[r*8+c]));
        }
    }
    printf("Test Passed!\n");

    printf("Performing second level wavelet 5_3 synthesis\n");

    float LPfilterhandle_s53[3] = {2*0.25, 2*0.5, 2*0.25};
    float HPfilterhandle_s53[5] = {-2*0.125, -2*0.25, 2*0.75, -2*0.25, -2*0.125};

    float *LPfilter_s53 = LPfilterhandle_s53 + HL_S53_LP;
    float *HPfilter_s53 = HPfilterhandle_s53 + HL_S53_HP;

    image_a53_l2->perform_boundary_extension_wavelet(2, offset);

    float buf_a53_s53_l2[64] = {-1.00195312500,  0.00000000000, -0.43164062500,  0.00000000000,  0.12084960940, -0.03564453125,  0.86938476560,  0.0712890625,
                                 0.00000000000,  0.00195312500,  0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.08911132813, -0.0693359375,
                                -0.43164062500,  0.00000000000,  0.12976074220, -0.01782226563,  0.82482910160,  0.10693359380,  0.74462890630, -0.1782226563,
                                 0.00000000000,  0.00195312500, -0.01782226563, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.0019531250,
                                 0.12084960940, -0.01782226563,  0.82482910160,  0.10693359380,  0.82482910160, -0.01782226563,  0.05847167969, -0.1425781250,
                                -0.03564453125, -0.03369140625,  0.10693359380, -0.03369140625, -0.01782226563,  0.00195312500,  0.00000000000,  0.0019531250,
                                 0.86938476560,  0.08911132813,  0.74462890630, -0.01782226563,  0.05847167969,  0.00000000000, -0.57421875000, -0.1425781250,
                                 0.07128906250, -0.06933593750, -0.17822265630,  0.00195312500, -0.14257812500,  0.00195312500, -0.14257812500,  0.0019531250};

    my_image_comp *image_a53_s53_l2 = new my_image_comp;
    image_a53_s53_l2->init(8, 8, boundary);

    synthesis(image_a53_l2, image_a53_s53_l2, 2, offset, LPfilter_s53, HPfilter_s53, HL_S53_LP, HL_S53_HP);

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            assert(CompareFloats(image_a53_s53_l2->buf[r*image_a53_s53_l2->stride+c], buf_a53_s53_l2[r*8+c]));
        }
    }
    printf("Test Passed!\n");

    printf("Performing single wavelet 5_3 synthesis\n");

    image_a53_s53_l2->perform_boundary_extension_wavelet(1, offset);

    float buf_a53_s53_l1[64] = {-1.0000000, -0.7187500, -0.4296875, -0.1484375,  0.1406250,  0.4218750,  0.7109375,  0.9921875,
                                -0.7187500, -0.4296875, -0.1484375,  0.1406250,  0.4218750,  0.7109375,  0.9921875,  0.7109375,
                                -0.4296875, -0.1484375,  0.1406250,  0.4218750,  0.7109375,  0.9921875,  0.7109375,  0.4218750,
                                -0.1484375,  0.1406250,  0.4218750,  0.7109375,  0.9921875,  0.7109375,  0.4218750,  0.1406250,
                                 0.1406250,  0.4218750,  0.7109375,  0.9921875,  0.7109375,  0.4218750,  0.1406250, -0.1484375,
                                 0.4218750,  0.7109375,  0.9921875,  0.7109375,  0.4218750,  0.1406250, -0.1484375, -0.4296875,
                                 0.7109375,  0.9921875,  0.7109375,  0.4218750,  0.1406250, -0.1484375, -0.4296875, -0.7187500,
                                 0.9921875,  0.7109375,  0.4218750,  0.1406250, -0.1484375, -0.4296875, -0.7187500, -1.0000000};

    my_image_comp *image_a53_s53_l1 = new my_image_comp;
    image_a53_s53_l1->init(8, 8, boundary);

    synthesis(image_a53_s53_l2, image_a53_s53_l1, 1, offset, LPfilter_s53, HPfilter_s53, HL_S53_LP, HL_S53_HP);

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            assert(CompareFloats(image_a53_s53_l1->buf[r*image_a53_s53_l1->stride+c], buf_a53_s53_l1[r*8+c]));
        }
    }
    printf("Test Passed!\n");

    /*printf("Testing two level wavelet analysis with offset 7, 5\n");

    my_image_comp *image16 = new my_image_comp;
    image16->init(16, 16, boundary);

    float buf16 [16*16] = { 0.0000000, -0.8671875, -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,
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
                            0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.3281250,  0.1953125,  0.0625000, -0.0703125, -0.2031250, -0.3359375, -0.4687500, -0.6015625, -0.7343750, -0.8671875,  0.0000000};

    for(int r = 0; r < 16; r++){
        for(int c = 0; c < 16; c++){
            image16->buf[r*image16->stride + c] = buf16[r*16+c];
        }
    }

    offset(0,0) = 7;
    offset(1,0) = 5;

    float buf16_a53_l1[16*16] = { 0.0000000, -0.8671875, -0.7343750, -0.6015625, -0.4687500, -0.33593750000, -0.20312500000, -0.07031250000,  0.06250000000,  0.19531250000,  0.32812500000,  0.4609375000,  0.5937500,  0.7265625,  0.8593750,  0.9921875,
                                 -0.8671875, -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.20312500000, -0.07031250000,  0.06250000000,  0.19531250000,  0.32812500000,  0.46093750000,  0.5937500000,  0.7265625,  0.8593750,  0.9921875,  0.8593750,
                                 -0.7343750, -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.07031250000,  0.06250000000,  0.19531250000,  0.32812500000,  0.46093750000,  0.59375000000,  0.7265625000,  0.8593750,  0.9921875,  0.8593750,  0.7265625,
                                 -0.6015625, -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.06250000000,  0.19531250000,  0.32812500000,  0.46093750000,  0.59375000000,  0.72656250000,  0.8593750000,  0.9921875,  0.8593750,  0.7265625,  0.5937500,
                                 -0.4687500, -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.19531250000,  0.32812500000,  0.46093750000,  0.59375000000,  0.72656250000,  0.85937500000,  0.9921875000,  0.8593750,  0.7265625,  0.5937500,  0.4609375,
                                 -0.3359375, -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.32812500000,  0.46093750000,  0.59375000000,  0.72656250000,  0.85937500000,  0.99218750000,  0.8593750000,  0.7265625,  0.5937500,  0.4609375,  0.3281250,
                                 -0.2031250, -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.46093750000,  0.59375000000,  0.72656250000,  0.85937500000,  0.99218750000,  0.85937500000,  0.7265625000,  0.5937500,  0.4609375,  0.3281250,  0.1953125,
                                 -0.0703125,  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.58959960940, -0.00830078125,  0.91333007810,  0.04980468750,  0.91333007810, -0.00830078125,  0.5895996094,  0.0000000,  0.3281250,  0.0000000,  0.0625000,
                                  0.0625000,  0.1953125,  0.3281250,  0.4609375,  0.5937500, -0.00830078125, -0.01660156250,  0.04980468750, -0.01660156250, -0.00830078125,  0.00000000000,  0.0000000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                                  0.1953125,  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.91333007810,  0.04980468750,  0.91333007810, -0.00830078125,  0.58959960940,  0.00000000000,  0.3281250000,  0.0000000,  0.0625000,  0.0000000, -0.2031250,
                                  0.3281250,  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.04980468750, -0.01660156250, -0.00830078125,  0.00000000000,  0.00000000000,  0.00000000000,  0.0000000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                                  0.4609375,  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.91333007810, -0.00830078125,  0.58959960940,  0.00000000000,  0.32812500000,  0.00000000000,  0.0625000000,  0.0000000, -0.2031250,  0.0000000, -0.4687500,
                                  0.5937500,  0.7265625,  0.8593750,  0.9921875,  0.8593750, -0.00830078125,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.0000000000,  0.0000000,  0.0000000,  0.0000000,  0.0000000,
                                  0.7265625,  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.58959960940,  0.00000000000,  0.32812500000,  0.00000000000,  0.06250000000,  0.00000000000, -0.2031250000,  0.0000000, -0.4531250,  0.0312500, -0.8281250,
                                  0.8593750,  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.00000000000,  0.0000000000,  0.0000000,  0.0312500,  0.0625000, -0.1875000,
                                  0.9921875,  0.8593750,  0.7265625,  0.5937500,  0.4609375,  0.32812500000,  0.00000000000,  0.06250000000,  0.00000000000, -0.20312500000,  0.00000000000, -0.4687500000,  0.0000000, -0.8281250, -0.1875000, -0.4375000};

    my_image_comp *image_a53_l1 = new my_image_comp;
    image_a53_l1->init(8, 8, boundary);

    analysis(image, image_a53_l1, 1, offset, LPfilter_a53, HPfilter_a53, HL_A53_LP, HL_A53_HP);

    for(int r = 0; r < 8; r++){
        for(int c = 0; c < 8; c++){
            assert(CompareFloats(image_a53_l1->buf[r*image_a53_l1->stride+c], buf_a53_l1[r*8+c]));
        }
    }

    printf("Test Passed!\n");*/

    printf("ALL TESTS PASSED! YOU ARE AWESOME!\n");
}

bool CompareFloats (float A, float B){
   float diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}
