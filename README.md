# Fast Gaussian Blur

C++ implementation of a fast gaussian blur approximation in linear time. It is based on a blog post by Ivan Kutskir: [blog](http://blog.ivank.net/fastest-gaussian-blur.html). Which refers to a presentation by Wojciech Jarosz: [slides](http://elynxsdk.free.fr/ext-docs/Blur/Fast_box_blur.pdf). Which itself describes an algorithm from the paper **Fast Almost-Gaussian Filtering** by Peter Kovesi: [site](https://www.peterkovesi.com/matlabfns/#integral), [paper](https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf). The code uses STB_IMAGE and STB_IMAGE_WRITE by stb for image manipulation: [stb github](https://github.com/nothings/stb). 

## Details

There are several implementations, from slowest to fastest:

- integer buffer per channel blur in `blur_int.cpp`
- floating point buffer per channel blur in `blur_float.cpp`
- integer buffer all channels blur in `blur_int_rgb.cpp`
- floating point buffer all channels blur in `blur_float_rgb.cpp`
- unsigned char buffer all channels blur in `blur_uchar_rgb.cpp`
- `fast_gaussian_blur.h` is a WIP header that regroup the implementations, ands the last cache coherent version that has not its own example file for now.

Integer versions are slower due to additional float -> int and int -> float casts. The fastest version blurs 2000k pixels in ~15ms on all cores of a Ryzen 7 2700X CPU with OpenMP. Hence it may be used for real-time applications with reasonable image resolutions. A SIMD vectorized or a GPU version of this algorithm could be significantly faster (but may be painful for the developper for arbitrary channels number / data sizes).

## Compilation

In a Unix term or WSL term you can use the provided makefile; use `make all` to build all targets without dependencies.

## Usage

Run the program with the following command:

`./program <input_image_filename> <sigma> <output_image_filename>`

- program can be any of `blur_float_rgb`, `blur_float`, `blur_int_rgb`, `blur_int`, `blur_uchar_rgb`.
- input_image_filename should be any of [.jpg, .png, .bmp, .tga, .psd, .gif, .hdr, .pic, .pnm].
- sigma is the desired Gaussian blur standard deviation (should be positive).
- output_image_filename should be any of [.png, .jpg, .bmp]  (unknown extensions will be saved as .png by default).

## Results

The fast Gaussian blur approx is linear in time regarding the size of the input image, but independent of sigma hence sigma = 5 is equally fast as sigma = 50.
|Original|sigma = 2|sigma = 5|sigma = 10|sigma = 30|sigma = 50|
|:---:|:---:|:---:|:---:|:---:|:---:|
![](data/demo.png)|![](data/blur2.png)|![](data/blur5.png)|![](data/blur10.png)|![](data/blur30.png)|![](data/blur50.png)|

## Performance

![](data/time.png)

The above graph shows the average exectution time of blur algorithm w.r.t pixel number on Ryzen 7 2700X. The dashed blue line highlights the fact that column major traversal of large image buffer may result in cache incohenrency. Hence we can perform image buffer transpositions and only cache coherent row major traversals to mitigate the problem. However the transposition step is not a cache friendly operation thus on large image buffer we observe the slope w/ transpose slowly increasing w.r.t image size. Performing the image transposition with fixed squared blocks per thread helps preserving the cache coherency and results in the fastest version of the algortihm (flip bloc).   

## Licence

You may use, distribute and modify this code under the terms of the MIT license. For further details please refer to : https://mit-license.org/

## TODO

- [ ] to read [recursive gaussian filters](https://software.intel.com/content/dam/develop/external/us/en/documents/cwp546-181134.pdf)
- [ ] to read [Fast O(1) bilateral filtering using trigonometric range kernels](http://bigwww.epfl.ch/chaudhury/Fast%20bilateral%20filtering.pdf)
