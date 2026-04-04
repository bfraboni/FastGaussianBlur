# make
# make
fastblur: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur_opt -O3 -fopenmp -std=c++17 -mavx2

debug: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur_opt -Og -g -std=c++17 -mavx2

single: main.cpp fast_gaussian_blur_template.h
	g++ main.cpp -o fastblur_opt -O3 -std=c++17 -mavx2

all: fastblur

clean:
	rm fastblur