#include <iostream>
#include <chrono>

// image io
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// fast blur
#include "fast_gaussian_blur_template.h"
#include "avx2_blur.h" // <--- Our new file

typedef unsigned char uchar;

int main(int argc, char * argv[])
{   
    if( argc < 4 )
    {
        printf("%s [input] [output] [sigma] [order - optional] [border - optional]\n", argv[0]);
        exit(1);
    }

    // load image
    int width, height, channels;
    uchar * image_data = stbi_load(argv[1], &width, &height, &channels, 0);
    printf("Source image: %s %dx%d (%d channels)\n", argv[1], width, height, channels);

    // read parameters
    const float sigma = std::atof(argv[3]);
    const int passes = argc > 4 ? std::atoi(argv[4]) : 3;
    const std::string policy = argc > 5 ? std::string(argv[5]) : "mirror";
    Border border;
    if (policy == "mirror")         border = Border::kMirror;
    else if (policy == "extend")    border = Border::kExtend;
    else if (policy == "crop")      border = Border::kKernelCrop;
    else if (policy == "wrap")      border = Border::kWrap;
    else                            border = Border::kMirror;
    
    std::size_t size = width * height * channels;
    uchar * new_image = new uchar[size];
    uchar * old_image = new uchar[size];
    
    for(std::size_t i = 0; i < size; ++i)
    {
        old_image[i] = image_data[i];
    }
    
    // stats
    auto start = std::chrono::system_clock::now(); 
    
    // ---------------------------------------------------------
    // THE ROUTER
    // ---------------------------------------------------------
    if (channels == 3) {
        printf("Routing to Custom AVX2 RGB Blur...\n");
        
        fast_gaussian_blur_avx2(old_image, new_image, width, height, sigma);
    } else {
        printf("Routing to standard template blur...\n");
        fast_gaussian_blur(old_image, new_image, width, height, channels, sigma, passes, border);
    }
    
    auto end = std::chrono::system_clock::now();
    float elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    printf("Time %.4fms\n", elapsed);

    // convert result
    for(std::size_t i = 0; i < size; ++i)
    {
        image_data[i] = (uchar)(new_image[i]);
    }

    // save image
    std::string file(argv[2]), ext = file.substr(file.size()-3);
    if( ext == "bmp" )
        stbi_write_bmp(argv[2], width, height, channels, image_data);
    else if( ext == "jpg" )
        stbi_write_jpg(argv[2], width, height, channels, image_data, 90);
    else
    {
        if( ext != "png" ) file = file.substr(0, file.size()-4) + std::string(".png");
        stbi_write_png(file.c_str(), width, height, channels, image_data, channels*width);
    }
    
    stbi_image_free(image_data);
    delete[] new_image;
    delete[] old_image;

    return 0;
}