# Fast Gaussian Blur

C++ implementation of a fast gaussian blur approximation in linear time. It is based on a blog post by Ivan Kutskir: [blog](http://blog.ivank.net/fastest-gaussian-blur.html). Which refers to a presentation by Wojciech Jarosz: [slides](http://elynxsdk.free.fr/ext-docs/Blur/Fast_box_blur.pdf). Which itself describes an algorithm from the paper **Fast Almost-Gaussian Filtering** by Peter Kovesi: [site](https://www.peterkovesi.com/matlabfns/#integral), [paper](https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf). The code uses STB_IMAGE and STB_IMAGE_WRITE by stb for image manipulation: [stb github](https://github.com/nothings/stb). 

**Note**: the fast gaussian blur algorithm is not accurate on image boundaries. It performs a diffusion of the signal with several independant passes, each pass depending of the preceding one. Some of the diffused signal is lost near borders and results in a slight loss of accuracy for next pass. This problem can be solved by increasing the image support of half the box kernel extent at each pass of the algorithm. The added padding would in this case capture the diffusion and make the next pass accurate. On contrary true Gaussian blur does not suffer this problem since it is performed in one pass only.

## Details

There are several implementations, from slowest to fastest:

- integer buffer per channel blur in `blur_int.cpp`
- floating point buffer per channel blur in `blur_float.cpp`
- integer buffer all channels blur in `blur_int_rgb.cpp`
- floating point buffer all channels blur in `blur_float_rgb.cpp`
- unsigned char buffer all channels blur in `blur_uchar_rgb.cpp`
- `fast_gaussian_blur.h` is a WIP header that regroup the implementations, and a cache coherent version that has not its own example file for now.
- `fast_gaussian_blur_template.h` is a WIP header that contains the latest templated cache coherent version that has not its own example file for now.

Integer versions are slower due to additional float -> int and int -> float casts. The fastest version blurs 2000k pixels in ~7ms on all cores of a Ryzen 7 2700X CPU with OpenMP. Hence it may be used for real-time applications with reasonable image resolutions. A SIMD vectorized or a GPU version of this algorithm could be significantly faster (but may be painful for the developper for arbitrary channels number / data sizes).

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

## Acknowledgments

Special thanks to Jean-Claude Iehl (@jciehl) for our insightful discussions and his passion for making code simple and fast. 

## Licence

You may use, distribute and modify this code under the terms of the MIT license. For further details please refer to : https://mit-license.org/

## References

- [Recursive gaussian filters](https://software.intel.com/content/dam/develop/external/us/en/documents/cwp546-181134.pdf)
- [Fast O(1) bilateral filtering using trigonometric range kernels](http://bigwww.epfl.ch/chaudhury/Fast%20bilateral%20filtering.pdf)
- [A Survey of Gaussian Convolution Algorithms](http://www.ipol.im/pub/art/2013/87/)
- [Filtering by repeated integration](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.72.4795)
- [Fast Filter Spreading and its Applications](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2009/EECS-2009-54.pdf)
- [Fast Almost-Gaussian Filtering](https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf)
- [Fast image convolutions](http://elynxsdk.free.fr/ext-docs/Blur/Fast_box_blur.pdf)
