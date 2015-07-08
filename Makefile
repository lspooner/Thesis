CC = g++
CFLAGS = -Wall -Werror -O2 -o 
SOURCES = image_comps.cpp io_bmp.cpp resize.cpp

all: resample diff affine

resample: clean $(SOURCES) main_resample.cpp
	$(CC) $(CFLAGS) resample $(SOURCES) main_resample.cpp
    
diff: clean $(SOURCES) main_diff.cpp
	$(CC) $(CFLAGS) diff $(SOURCES) main_diff.cpp
	
affine: clean affine.cpp main_affine.cpp
	$(CC) $(CFLAGS) affine affine.cpp main_affine.cpp
	
affine_armadillo: clean affine_armadillo.cpp main_affine_armadillo.cpp
	$(CC) $(CFLAGS) affine_armadillo -I /home/lauren/armadillo-5.100.1/include -lblas -llapack affine_armadillo.cpp main_affine_armadillo.cpp

clean:
	@rm -f resample diff affine affine_armadillo *.o core
