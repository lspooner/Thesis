#include "io_bmp.h"
#include "image_comps.h"

io_byte floatToByte(float num);

void my_image_comp::init(int height, int width, int border){
    this->width = width;  
    this->height = height;  
    this->border = border;
    stride = width + 2*border;
    if (handle != NULL){
        delete[] handle; // Delete mem allocated by any previous `init' call
    }
    printf("Getting %lu bytes\n", sizeof(float)*stride*height+2*border);
    handle = new float[stride*(height+2*border)];
    buf = handle + (border*stride) + border;
}

void my_image_comp::perform_boundary_extension(){
    int r, c;

    // First extend upwards
    float *first_line = buf;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            first_line[-r*stride+c] = first_line[c];
        }
    }

    // Now extend downwards
    float *last_line = buf+(height-1)*stride;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            last_line[r*stride+c] = last_line[c];
        }
    }

    // Now extend all rows to the left and to the right
    float *left_edge = buf-border*stride;
    float *right_edge = left_edge + width - 1;
    for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride){
        for (c=1; c <= border; c++){
            left_edge[-c] = left_edge[0];
            right_edge[c] = right_edge[0];
        }
    }
}

void my_image_comp::perform_boundary_extension_zero_padding(){
    int r, c;

    // First extend upwards
    float *first_line = buf;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            first_line[-r*stride+c] = 0;
        }
    }

    // Now extend downwards
    float *last_line = buf+(height-1)*stride;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            last_line[r*stride+c] = 0;
        }
    }

    // Now extend all rows to the left and to the right
    float *left_edge = buf-border*stride;
    float *right_edge = left_edge + width - 1;
    for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride){
        for (c=1; c <= border; c++){
            left_edge[-c] = 0;
            right_edge[c] = 0;
        }
    }
}

void my_image_comp::perform_boundary_extension_symmetric(){
    //not technically the right way to do, just assuming our image is bigger than the extension
    int r, c;

    // First extend upwards
    float *first_line = buf;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            first_line[-r*stride+c] = first_line[r*stride+c];
        }
    }

    // Now extend downwards
    float *last_line = buf+(height-1)*stride;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            last_line[r*stride+c] = last_line[-r*stride+c];
        }
    }

    // Now extend all rows to the left and to the right
    float *left_edge = buf-border*stride;
    float *right_edge = left_edge + width - 1;
    for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride){
        for (c=1; c <= border; c++){
            left_edge[-c] = left_edge[c];
            right_edge[c] = right_edge[-c];
        }
    }
}

void my_image_comp::perform_boundary_extension_wavelet(int spacing, Mat<float> offset){
    //not technically the right way to do, just assuming our image is bigger than the extension
    int r, c;

    int extendColumn = width+(int)offset(0,0)%spacing - 2*spacing;
    int extendRow = height+(int)offset(1,0)%spacing - 2*spacing;

    // First extend upwards
    float *first_line = buf;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            first_line[-r*stride+c] = first_line[r*stride+c];
        }
    }

    // Now extend downwards
    float *last_line = buf+(height-1)*stride;
    int r_extended = extendRow;
    for (r=1; r <= border; r++){
        for (c=0; c < width; c++){
            last_line[r*stride+c] = buf[r_extended*stride+c];
        }
        r_extended--;
    }

    // Now extend all rows to the left and to the right
    float *left_edge = buf-border*stride;
    float *right_edge = left_edge + width - 1;
    for (r=-border; r < height+border; r++, left_edge+=stride, right_edge+=stride){
        int c_extended = extendColumn;
        for (c=1; c <= border; c++){
            left_edge[-c] = left_edge[c];
            right_edge[c] = buf[r*stride + c_extended];
            c_extended--;
        }
    }
}

int readBMP(char* fileName, my_image_comp** input_comps, int border, int* num_comps){
    // Read the input image
    bmp_in in;
    int err_code;
    if ((err_code = bmp_in__open(&in, fileName)) != 0){
        return err_code;
    }

    int width = in.cols, height = in.rows;
    int n;
    *num_comps = in.num_components;
    *input_comps = new my_image_comp[*num_comps];
    for (n=0; n < *num_comps; n++){
        (*input_comps)[n].init(height,width,border); // Leave a border
    }

    /*Mat<float> convertColour(3,3);
    convertColour(0,0) = 0.299;
    convertColour(0,1) = 0.587;
    convertColour(0,2) = 0.114;
    convertColour(1,0) = -0.169;
    convertColour(1,1) = -0.331;
    convertColour(1,2) = 0.500;
    convertColour(2,0) = 0.500;
    convertColour(2,1) = -0.419;
    convertColour(2,2) = -0.081;

    Mat<float> colourOffset(3,1);
    colourOffset(0,0) = 0;
    colourOffset(1,0) = 128;
    colourOffset(2,0) = 128;*/
    
    //printf("\nReading in BMP\n");
    int r; // Declare row index
    io_byte *line = new io_byte[width*(*num_comps)];
    for (r=height-1; r >= 0; r--){
        // "r" holds the true row index we are reading, since the image is
        // stored upside down in the BMP file.
        if ((err_code = bmp_in__get_line(&in,line)) != 0){
            return err_code;
        }
        //if(*num_comps == 1){
            for (n=0; n < *num_comps; n++){
                io_byte *src = line+n; // Points to first sample of component n
                float *dst = (*input_comps)[n].buf + r * (*input_comps)[n].stride;
                for (int c=0; c < width; c++, src+=(*num_comps)){
                    dst[c] = (((float) *src)-128)/128.0; 
                    // The cast to type "float" is not
                    // strictly required here, since bytes can always be
                    // converted to floats without any loss of information.
                }
            }
        /*} else if(*num_comps == 3){

            io_byte *src = line; // Points to first sample of component n
            //float *dst = (*input_comps)[n].buf + r * (*input_comps)[n].stride;
            for (int c=0; c < width; c++, src+=(*num_comps)){
                Mat<float> RGB(3,1);
                RGB(0,0) = (float) src[0];
                RGB(1,0) = (float) src[1];
                RGB(2,0) = (float) src[2];

                Mat<float> YCbCr = convertColour*RGB + colourOffset;

                (*input_comps)[0].buf[r*(*input_comps)[0].stride+c] = ((YCbCr(0,0) - 128)/128.0);
                (*input_comps)[1].buf[r*(*input_comps)[1].stride+c] = ((YCbCr(1,0) - 128)/128.0);
                (*input_comps)[2].buf[r*(*input_comps)[2].stride+c] = ((YCbCr(2,0) - 128)/128.0);

                //dst[c] = (((float) *src)-128)/128.0; 
                // The cast to type "float" is not
                // strictly required here, since bytes can always be
                // converted to floats without any loss of information.
            }
        } else {
            printf("I can't handle this many components!\n");
            assert(0);
        }*/
    }
    bmp_in__close(&in);
    delete[] line;
    return 0;
}

int outputBMP(char* fileName, my_image_comp* out_comps, int num_comps){
        // Write the image back out again
    bmp_out out;
    int err_code;
    int width = out_comps[0].width;
    int height = out_comps[0].height;
    if ((err_code = bmp_out__open(&out, fileName, width, height, num_comps)) != 0){
        return err_code;
    }

    /*Mat<float> convertColour(3,3);
    convertColour(0,0) = 1.000;
    convertColour(0,1) = 0.000;
    convertColour(0,2) = 1.400;
    convertColour(1,0) = 1.000;
    convertColour(1,1) = -0.343;
    convertColour(1,2) = -0.711;
    convertColour(2,0) = 1.000;
    convertColour(2,1) = 1.765;
    convertColour(2,2) = 0.000;

    Mat<float> colourOffset(3,1);
    colourOffset(0,0) = 0;
    colourOffset(1,0) = -128;
    colourOffset(2,0) = -128;*/

    io_byte *line = new io_byte[width*num_comps];
    for (int r=height-1; r >= 0; r--){
        // "r" holds the true row index we are writing, since the image is
        // written upside down in BMP files.
        //if(num_comps == 1){
            for (int n=0; n < num_comps; n++){
                io_byte *dst = line+n; // Points to first sample of component n
                float *src = out_comps[n].buf + r * out_comps[n].stride;
                for (int c=0; c < width; c++, dst+=num_comps){
                    *dst = floatToByte(src[c]); 
                    // The cast to type "io_byte" is
                    // required here, since floats cannot generally be
                    // converted to bytes without loss of information.  The
                    // compiler will warn you of this if you remove the cast.
                    // There is in fact not the best way to do the
                    // conversion.  You should fix it up in the lab.
                }
            }
        /*} else if (num_comps == 3){
            io_byte *dst = line; // Points to first sample of component n
            //float *src = out_comps[n].buf + r * out_comps[n].stride;
            for (int c=0; c < width; c++, dst+=num_comps){
                Mat<float> YCbCr(3,1);
                YCbCr(0,0) = floatToByte(out_comps[0].buf[r*out_comps[0].stride+c]);
                YCbCr(1,0) = floatToByte(out_comps[1].buf[r*out_comps[1].stride+c]);
                YCbCr(2,0) = floatToByte(out_comps[2].buf[r*out_comps[2].stride+c]);

                Mat<float> RGB = convertColour*(YCbCr + colourOffset);

                dst[0] = RGB(0,0);
                dst[1] = RGB(1,0);
                dst[2] = RGB(2,0);
            }

        } else {
            printf("I can't handle this many components!\n");
            assert(0);
        }*/
        bmp_out__put_line(&out,line);
    }
    bmp_out__close(&out);
    delete[] line;
    return 0;
}

io_byte floatToByte(float num){
    //num += 0.5; //round to nearest int
    if(num >= 1){
        return 255;
    } else if(num <= -1){
        return 0;
    } else {
        return (io_byte) (num*128+128.5);
    }
}