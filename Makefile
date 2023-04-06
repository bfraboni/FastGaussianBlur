# make
fastblur: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -O3 -fopenmp -std=c++17

debug: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -Og -g -std=c++17

single: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -O3 -std=c++17

all: fastblur

clean:
	rm fastblur
