# make
blur_float: blur_float.cpp
	g++ blur_float.cpp -o blur_float -O3

blur_float_rgb: blur_float_rgb.cpp
	g++ blur_float_rgb.cpp -o blur_float_rgb -O3

blur_int: blur_int.cpp
	g++ blur_int.cpp -o blur_int -O3

blur_int_rgb: blur_int_rgb.cpp
	g++ blur_int_rgb.cpp -o blur_int_rgb -O3

all: blur_float blur_float_rgb blur_int blur_int_rgb

clean:
	rm blur_float blur_float_rgb blur_int blur_int_rgb
