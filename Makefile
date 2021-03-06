CC = g++
CFLAGS = -Wall -Werror -O2 -o 
UFLAGS = -Wall -Werror -g -o 
SOURCES = image_comps.cpp io_bmp.cpp

all: resize diff wavelet affine combine

unit_all: unit_wavelet unit_affine

resize: $(SOURCES) resize.cpp main_resize.cpp
	$(CC) $(CFLAGS) resize $(SOURCES) resize.cpp main_resize.cpp -larmadillo

diff: $(SOURCES) main_diff.cpp
	$(CC) $(CFLAGS) diff $(SOURCES) main_diff.cpp -larmadillo

wavelet: $(SOURCES) wavelet.cpp main_wavelet.cpp
	$(CC) $(CFLAGS) wavelet $(SOURCES) wavelet.cpp main_wavelet.cpp -larmadillo

affine: $(SOURCES) affine.cpp main_affine.cpp wavelet.cpp resize.cpp
	$(CC) $(CFLAGS) affine $(SOURCES) affine.cpp main_affine.cpp wavelet.cpp resize.cpp -larmadillo

combine: $(SOURCES) combine.cpp main_combine.cpp wavelet.cpp affine.cpp resize.cpp
	$(CC) $(CFLAGS) combine $(SOURCES) combine.cpp main_combine.cpp wavelet.cpp affine.cpp resize.cpp -larmadillo

unit_wavelet: $(SOURCES) wavelet.cpp unit_wavelet.cpp
	$(CC) $(UFLAGS) unit_wavelet $(SOURCES) wavelet.cpp unit_wavelet.cpp

unit_affine: $(SOURCES) affine.cpp unit_affine.cpp wavelet.cpp resize.cpp
	$(CC) $(UFLAGS) unit_affine $(SOURCES) affine.cpp unit_affine.cpp wavelet.cpp resize.cpp -larmadillo

clean:
	@rm -f resize diff wavelet affine combine unit_wavelet unit_affine *.o core
