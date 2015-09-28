#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "io_bmp.h"
#include "wavelet.h"

using namespace std;
using namespace arma;

void analysis(my_image_comp *in, my_image_comp *out, int spacing, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL);
void synthesis(my_image_comp *in, my_image_comp *out, int spacing, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL);

void copyBuffer(my_image_comp *copyFrom, my_image_comp *copyTo);

int analysis_5_3(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_A53_LP*scale+HL_A53_HP, &num_comps)) != 0){
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

    for(int n=0; n < num_comps; n++){
        analysis_5_3(input_comps+n, output_comps+n, levels);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void analysis_5_3(my_image_comp *in, my_image_comp *out, int levels){
    if(levels != 0){
        float LPfilterhandle[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
        float HPfilterhandle[3] = {-0.25, 0.5, -0.25};

        float *LPfilter = LPfilterhandle + HL_A53_LP;
        float *HPfilter = HPfilterhandle + HL_A53_HP;

        int scale = pow(2.0, (float) levels);
        for (int spacing=1; spacing < scale; spacing*=2){
            in->perform_boundary_extension_wavelet(spacing);
            analysis(in, out, spacing, LPfilter, HPfilter, HL_A53_LP, HL_A53_HP);
            if(spacing < scale/2){
                copyBuffer(out, in);
            }
        }
    } else {
        copyBuffer(in, out);
    }
}

int analysis_9_7(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_A97_LP*scale+HL_A97_HP, &num_comps)) != 0){
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

    for(int n=0; n < num_comps; n++){
        analysis_9_7(input_comps+n, output_comps+n, levels);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void analysis_9_7(my_image_comp *in, my_image_comp *out, int levels){
    float LPfilterhandle[9] = {0.026749, -0.016864, -0.078223, 0.266864, 0.602949, 0.266864, -0.078223, -0.016864, 0.026749};
    float HPfilterhandle[7] = {0.045636, -0.028772, -0.295636, 0.557543, -0.295636, -0.028772, 0.045636};

    float *LPfilter = LPfilterhandle + HL_A97_LP;
    float *HPfilter = HPfilterhandle + HL_A97_HP;

    int scale = pow(2.0, (float) levels);
    for (int spacing=1; spacing < scale; spacing*=2){
        in->perform_boundary_extension_wavelet(spacing);
        analysis(in, out, spacing, LPfilter, HPfilter, HL_A97_LP, HL_A97_HP);
        if(spacing < scale/2){
            copyBuffer(out, in);
        }
    }
}

int synthesis_5_3(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_S53_HP*scale+HL_S53_LP, &num_comps)) != 0){
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

    for(int n=0; n < num_comps; n++){
        synthesis_5_3(input_comps+n, output_comps+n, levels);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void synthesis_5_3(my_image_comp *in, my_image_comp *out, int levels){
    float LPfilterhandle[3] = {2*0.25, 2*0.5, 2*0.25};
    float HPfilterhandle[5] = {-2*0.125, -2*0.25, 2*0.75, -2*0.25, -2*0.125};

    float *LPfilter = LPfilterhandle + HL_S53_LP;
    float *HPfilter = HPfilterhandle + HL_S53_HP;

    int scale = pow(2.0, (float) levels);
    for (int spacing=scale/2; spacing >= 1; spacing/=2){
        in->perform_boundary_extension_wavelet(spacing);
        synthesis(in, out, spacing, LPfilter, HPfilter, HL_S53_LP, HL_S53_HP);
        if(spacing > 1){
            copyBuffer(out, in);
        }
    }
}

int synthesis_9_7(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, HL_S97_HP*scale+HL_S97_LP, &num_comps)) != 0){
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

    for(int n=0; n < num_comps; n++){
        synthesis_9_7(input_comps+n, output_comps+n, levels);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void synthesis_9_7(my_image_comp *in, my_image_comp *out, int levels){
    float LPfilterhandle[7] = {-2*0.045636, -2*0.028772, 2*0.295636, 2*0.557543, 2*0.295636, -2*0.028772, -2*0.045636};
    float HPfilterhandle[9] = {2*0.026749, 2*0.016864, -2*0.078223, -2*0.266864, 2*0.602949, -2*0.266864, -2*0.078223, 2*0.016864, 2*0.026749};

    float *LPfilter = LPfilterhandle + HL_S97_LP;
    float *HPfilter = HPfilterhandle + HL_S97_HP;

    int scale = pow(2.0, (float) levels);
    for (int spacing=scale/2; spacing >= 1; spacing/=2){
        in->perform_boundary_extension_wavelet(spacing);
        synthesis(in, out, spacing, LPfilter, HPfilter, HL_S97_LP, HL_S97_HP);
        if(spacing > 1){
            copyBuffer(out, in);
        }
    }
}

void analysis(my_image_comp *in, my_image_comp *out, int spacing, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL){
    // Check for consistent dimensions
    assert(in->border >= LP_HL*spacing+HP_HL);
    assert((out->height == in->height) && (out->width == in->width));

    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels];
    float *tempPicBuf = tempPic + out->width*in->border;

    //apply lowpass filter to the rows
    bool invalid = false;
    for(int r = -LP_HL*spacing; r < in->height + LP_HL*spacing; r+=spacing){
        for(int c = 0; c < out->width; c+=2*spacing){
            float sum = 0.0F;
            invalid = false;
            for(int n = -LP_HL; n <= LP_HL && !invalid; n++){
                if(in->buf[r*in->stride+c+spacing*n] == INVALID){
                    invalid = true;
                }
                sum += LPfilter[n] * in->buf[r*in->stride+c+spacing*n];
            }
            if(invalid){
                tempPicBuf[r*out->width + c] = INVALID;
            } else {
                tempPicBuf[r*out->width + c] = sum;
            }
        }
    }

    //apply highpass filter to the rows
    for(int r = -LP_HL*spacing; r < in->height + LP_HL*spacing; r+=spacing){
        for(int c = spacing; c < out->width; c+=2*spacing){
            float sum = 0.0F;
            invalid = false;
            for(int n = -HP_HL; n <= HP_HL && !invalid; n++){
                if(in->buf[r*in->stride+c+spacing*n] == INVALID){
                    invalid = true;
                }
                sum += HPfilter[n] * in->buf[r*in->stride+c+spacing*n];
            }
            if(invalid){
                tempPicBuf[r*out->width + c] = INVALID;
            } else {
                tempPicBuf[r*out->width + c] = sum;
            }
        }
    }

    //apply lowpass filter to the columns
    for(int r = 0; r < out->height; r+=2*spacing){
        for(int c = 0; c < out->width; c+=spacing){
            float sum = 0.0F;
            invalid = false;
            for(int n = -LP_HL; n <= LP_HL && !invalid; n++){
                if(tempPicBuf[(r+spacing*n)*out->width+c] == INVALID){
                    invalid = true;
                }
                sum += LPfilter[n] * tempPicBuf[(r+spacing*n)*out->width+c];
            }
            if(invalid){
                out->buf[r*out->stride+c] = INVALID;
            } else {
                out->buf[r*out->stride+c] = sum;
            }
        }
    }

    //apply highpass filter to the columns
    for(int r = spacing; r < out->height; r+=2*spacing){
        for(int c = 0; c < out->width; c+=spacing){
            float sum = 0.0F;
            invalid = false;
            for(int n = -HP_HL; n <= HP_HL && !invalid; n++){
                if(tempPicBuf[(r+spacing*n)*out->width+c] == INVALID){
                    invalid = true;
                }
                sum += HPfilter[n] * tempPicBuf[(r+spacing*n)*out->width+c];
            }
            if(invalid){
                out->buf[r*out->stride+c] = INVALID;
            } else {
                out->buf[r*out->stride+c] = sum;
            }
        }
    }

    if(spacing > 1){
        //fill in other values
        for(int r = 0; r < out->height; r++){
            for(int c = 0; c < out->width; c++){
                if(r%spacing != 0 || c%spacing != 0){
                    out->buf[r*out->stride+c] = in->buf[r*in->stride+c];
                }
            }
        }
    }

    delete[] tempPic;
}

void synthesis(my_image_comp *in, my_image_comp *out, int spacing, float* LPfilter, float* HPfilter, int LP_HL, int HP_HL){
    // Check for consistent dimensions
    assert(in->border >= HP_HL*spacing+LP_HL);
    assert((out->height == in->height) && (out->width == in->width));

    int numPixels = (in->height+2*in->border)*out->width;
    float *tempPic = new float[numPixels]();
    float *tempPicBuf = tempPic + out->width*in->border;

    for(int i = 0; i < out->stride*(out->height+2*out->border); i++){
        out->handle[i] = 0;
    }

    //apply lowpass filter to the rows
    for(int r = -HP_HL*spacing; r < in->height + HP_HL*spacing; r+=spacing){
        for(int c = -HP_HL*spacing; c < in->stride + HP_HL*spacing; c+=2*spacing){
            for(int n = -LP_HL; n <= LP_HL; n++){
                if(c+spacing*n >= 0 && c+spacing*n < out->width){
                    if(in->buf[r*in->stride+c] == INVALID || tempPicBuf[r*out->width+c+spacing*n] == INVALID){
                        tempPicBuf[r*out->width+c+spacing*n] = INVALID;
                    } else {
                        tempPicBuf[r*out->width+c+spacing*n] += LPfilter[n] * in->buf[r*in->stride+c];
                    }
                }
            }
        }
    }

    //apply highpass filter to the rows
    for(int r = -HP_HL*spacing; r < in->height + HP_HL*spacing; r+=spacing){
        for(int c = -HP_HL*spacing+spacing; c < in->stride + HP_HL*spacing; c+=2*spacing){
            for(int n = -HP_HL; n <= HP_HL; n++){
                if(c+spacing*n >= 0 && c+spacing*n < out->width){
                    if(in->buf[r*in->stride+c] == INVALID || tempPicBuf[r*out->width+c+spacing*n] == INVALID){
                        tempPicBuf[r*out->width+c+spacing*n] = INVALID;
                    } else {
                        tempPicBuf[r*out->width+c+spacing*n] += HPfilter[n] * in->buf[r*in->stride+c];
                    }
                }
            }
        }
    }

    //apply lowpass filter to the columns
    for(int r = -HP_HL*spacing; r < out->height + HP_HL*spacing; r+=2*spacing){
        for(int c = 0; c < out->width; c+=spacing){
            for(int n = -LP_HL; n <= LP_HL; n++){
                if(r+spacing*n >= 0 && r+spacing*n < out->height){
                    if(tempPicBuf[r*out->width+c] == INVALID || out->buf[(r+spacing*n)*out->stride+c] == INVALID){
                        out->buf[(r+spacing*n)*out->stride+c] = INVALID;
                    } else {
                        out->buf[(r+spacing*n)*out->stride+c] += LPfilter[n] * tempPicBuf[r*out->width+c];
                    }
                }
            }
        }
    }

    //apply highpass filter to the columns
    for(int r = -HP_HL*spacing+spacing; r < out->height + HP_HL*spacing; r+=2*spacing){
        for(int c = 0; c < out->width; c+=spacing){
            for(int n = -HP_HL; n <= HP_HL; n++){
                if(r+spacing*n >= 0 && r+spacing*n < out->height){
                    if(tempPicBuf[r*out->width+c] == INVALID || out->buf[(r+spacing*n)*out->stride+c] == INVALID){
                        out->buf[(r+spacing*n)*out->stride+c] = INVALID;
                    } else {
                        out->buf[(r+spacing*n)*out->stride+c] += HPfilter[n] * tempPicBuf[r*out->width+c];
                    }
                }
            }
        }
    }

    if(spacing > 1){
        //fill in other values
        for(int r = 0; r < out->height; r++){
            for(int c = 0; c < out->width; c++){
                if(r%spacing != 0 || c%spacing != 0){
                    out->buf[r*out->stride+c] = in->buf[r*in->stride+c];
                }
            }
        }
    }

    delete[] tempPic;
}

int increaseWaveletLevel(char* inputFile, char* outputFile, int levels){
    // Read the input image
    int err_code;
    my_image_comp *input_comps = NULL;
    int num_comps;
    int scale = pow(2.0, (float) levels);
    if ((err_code = readBMP(inputFile, &input_comps, 0, &num_comps)) != 0){
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

    for(int n=0; n < num_comps; n++){
        increaseWaveletLevel(input_comps+n, output_comps+n, scale);
    }

    // Write the image back out again
    if((err_code = outputBMP(outputFile, output_comps, num_comps)) != 0){
        return err_code;
    }

    delete[] input_comps;
    delete[] output_comps;
    return 0;
}

void increaseWaveletLevel(my_image_comp *in, my_image_comp *out, int scale){
    assert(out->height == in->height*scale);
    assert(out->width == in->width*scale);

    for(int r = 0; r < out->height; r++){
        for(int c = 0; c < out->width; c++){
            out->buf[r*out->stride+c] = 0;
        }
    }

    for(int r = 0; r < out->height; r+=scale){
        for(int c = 0; c < out->width; c+=scale){
            out->buf[r*out->stride+c] = in->buf[(r/scale)*in->stride+(c/scale)];
        }
    }
}

void copyBuffer(my_image_comp *copyFrom, my_image_comp *copyTo){
    assert(copyFrom->height == copyTo->height);
    assert(copyFrom->width == copyTo->width);

    for(int r = 0; r < copyTo->height; r++){
        for(int c = 0; c < copyTo->width; c++){
            copyTo->buf[r*copyTo->stride + c] = copyFrom->buf[r*copyFrom->stride+c];
        }
    }
}
