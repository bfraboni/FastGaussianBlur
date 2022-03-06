# make
fastblur: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -O3 -fopenmp

all: fastblur

clean:
	rm fastblur
