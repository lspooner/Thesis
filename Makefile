CC = g++
CFLAGS = -Wall -Werror -O2 -o 
UFLAGS = -Wall -Werror -g -o 
SOURCES = image_comps.cpp io_bmp.cpp

all: resize diff wavelet affine combine

unit_all: unit_wavelet

resize: clean $(SOURCES) resize.cpp main_resize.cpp
	$(CC) $(CFLAGS) resize $(SOURCES) resize.cpp main_resize.cpp

diff: clean $(SOURCES) main_diff.cpp
	$(CC) $(CFLAGS) diff $(SOURCES) main_diff.cpp

wavelet: clean $(SOURCES) wavelet.cpp main_wavelet.cpp
	$(CC) $(CFLAGS) wavelet $(SOURCES) wavelet.cpp main_wavelet.cpp

affine: clean $(SOURCES) affine.cpp main_affine.cpp
	$(CC) $(CFLAGS) affine $(SOURCES) affine.cpp main_affine.cpp -larmadillo

combine: clean $(SOURCES) combine.cpp main_combine.cpp wavelet.cpp
	$(CC) $(CFLAGS) combine $(SOURCES) combine.cpp main_combine.cpp wavelet.cpp

unit_wavelet: clean $(SOURCES) wavelet.cpp unit_wavelet.cpp
	$(CC) $(UFLAGS) unit_wavelet $(SOURCES) wavelet.cpp unit_wavelet.cpp

clean:
	@rm -f resize diff wavelet affine combine unit_wavelet *.o core
