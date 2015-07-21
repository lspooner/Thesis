CC = g++
CFLAGS = -Wall -Werror -O2 -o 
SOURCES = image_comps.cpp io_bmp.cpp

all: resize diff wavelet

resize: clean $(SOURCES) main_resize.cpp
	$(CC) $(CFLAGS) resize $(SOURCES) resize.cpp main_resize.cpp

diff: clean $(SOURCES) main_diff.cpp
	$(CC) $(CFLAGS) diff $(SOURCES) main_diff.cpp

wavelet: clean $(SOURCES) main_wavelet.cpp
	$(CC) $(CFLAGS) wavelet $(SOURCES) wavelet.cpp main_wavelet.cpp

clean:
	@rm -f resize diff wavelet *.o core
