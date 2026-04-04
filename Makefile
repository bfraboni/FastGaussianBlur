# make
# make
fastblur: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -O3 -fopenmp -std=c++17 -mavx2

debug: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -Og -g -std=c++17 -mavx2

single: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur -O3 -std=c++17 -mavx2

all: fastblur

clean:
	rm fastblur