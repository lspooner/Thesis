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
    
    //printf("\nReading in BMP\n");
    int r; // Declare row index
    io_byte *line = new io_byte[width*(*num_comps)];
    for (r=height-1; r >= 0; r--){
        // "r" holds the true row index we are reading, since the image is
        // stored upside down in the BMP file.
        if ((err_code = bmp_in__get_line(&in,line)) != 0){
            return err_code;
        }
        for (n=0; n < *num_comps; n++){
            io_byte *src = line+n; // Points to first sample of component n
            float *dst = (*input_comps)[n].buf + r * (*input_comps)[n].stride;
            for (int c=0; c < width; c++, src+=(*num_comps)){
                //printf("%d ", *src);
                dst[c] = (((float) *src)-128)/128.0; 
                // The cast to type "float" is not
                // strictly required here, since bytes can always be
                // converted to floats without any loss of information.
            }
            //printf("\n");
        }
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

    io_byte *line = new io_byte[width*num_comps];
    for (int r=height-1; r >= 0; r--){
        // "r" holds the true row index we are writing, since the image is
        // written upside down in BMP files.
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