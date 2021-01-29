// Copyright (C) 2017 Basile Fraboni
// Copyright (C) 2014 Ivan Kutskir
// All Rights Reserved
// You may use, distribute and modify this code under the
// terms of the MIT license. For further details please refer 
// to : https://mit-license.org/
//

//!
//! \file blur.cpp
//! \author Basile Fraboni
//! \date 2017
//!
//! \brief The software is a C++ implementation of a fast  
//! Gaussian blur algorithm by Ivan Kutskir. For further details 
//! please refer to : 
//! http://blog.ivank.net/fastest-gaussian-blur.html
//!
//! Floating point single precision version
//!

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>

//!
//! \fn void std_to_box(int boxes[], float sigma, int n)  
//!
//! \brief this function converts the standard deviation of 
//! Gaussian blur into dimensions of boxes for box blur. For 
//! further details please refer to :
//! https://www.peterkovesi.com/matlabfns/#integral
//! https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf
//!
//! \param[out] boxes   boxes dimensions
//! \param[in] sigma    Gaussian standard deviation
//! \param[in] n        number of boxes
//!
void std_to_box(int boxes[], float sigma, int n)  
{
    // ideal filter width
    float wi = std::sqrt((12*sigma*sigma/n)+1); 
    int wl = std::floor(wi);  
    if(wl%2==0) wl--;
    int wu = wl+2;
                
    float mi = (12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
    int m = std::round(mi);
                
    for(int i=0; i<n; i++) 
        boxes[i] = ((i < m ? wl : wu) - 1) / 2;
}

//!
//! \fn void horizontal_blur_rgb(float * in, float * out, int w, int h, int c, int r)   
//!
//! \brief this function performs the horizontal blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void horizontal_blur_rgb(float * in, float * out, int w, int h, int c, int r) 
{
    float iarr = 1.f / (r+r+1);
    for(int i=0; i<h; i++) 
    {
        int ti = i*w; 
        int li = ti;  
        int ri = ti+r;

        float fv[3] = { in[ti*c+0], in[ti*c+1], in[ti*c+2] };                  
        float lv[3] = { in[(ti+w-1)*c+0], in[(ti+w-1)*c+1], in[(ti+w-1)*c+2] };
        float val[3] = { (r+1)*fv[0], (r+1)*fv[1], (r+1)*fv[2] };              

        for(int j=0; j<r; j++) 
        { 
            val[0] += in[(ti+j)*c+0]; 
            val[1] += in[(ti+j)*c+1]; 
            val[2] += in[(ti+j)*c+2]; 
        }

        for(int j=0; j<=r; j++) 
        { 
            val[0] += in[ri*c+0] - fv[0]; 
            val[1] += in[ri*c+1] - fv[1]; 
            val[2] += in[ri*c+2] - fv[2]; 
            out[ti*c+0] = val[0]*iarr; 
            out[ti*c+1] = val[1]*iarr; 
            out[ti*c+2] = val[2]*iarr; 
            ++ri;
            ++ti;
        }

        for(int j=r+1; j<w-r; j++) 
        { 
            val[0] += in[ri*c+0] - in[li*c+0]; 
            val[1] += in[ri*c+1] - in[li*c+1]; 
            val[2] += in[ri*c+2] - in[li*c+2]; 
            out[ti*c+0] = val[0]*iarr; 
            out[ti*c+1] = val[1]*iarr; 
            out[ti*c+2] = val[2]*iarr; 
            ++ri;
            ++li; 
            ++ti;
        }

        for(int j=w-r; j<w; j++) 
        { 
            val[0] += lv[0] - in[li*c+0]; 
            val[1] += lv[1] - in[li*c+1]; 
            val[2] += lv[2] - in[li*c+2]; 
            out[ti*c+0] = val[0]*iarr; 
            out[ti*c+1] = val[1]*iarr; 
            out[ti*c+2] = val[2]*iarr; 
            ++li; 
            ++ti;
        }
    }
}

//!
//! \fn void total_blur_rgb(float * in, float * out, int w, int h, int c, int r)   
//!
//! \brief this function performs the total blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void total_blur_rgb(float * in, float * out, int w, int h, int c, int r) 
{
    // radius range on either side of a pixel + the pixel itself
    float iarr = 1.f / (r+r+1);
    for(int i=0; i<w; i++) 
    {
        int ti = i;
        int li = ti;
        int ri = ti+r*w;

        float fv[3] = {in[ti*c+0], in[ti*c+1], in[ti*c+2] };
        float lv[3] = {in[(ti+w*(h-1))*c+0], in[(ti+w*(h-1))*c+1], in[(ti+w*(h-1))*c+2] };
        float val[3] = {(r+1)*fv[0], (r+1)*fv[1], (r+1)*fv[2] };

        for(int j=0; j<r; j++) 
        { 
            val[0] += in[(ti+j*w)*c+0]; 
            val[1] += in[(ti+j*w)*c+1]; 
            val[2] += in[(ti+j*w)*c+2]; 
        }

        for(int j=0; j<=r; j++) 
        { 
            val[0] += in[ri*c+0] - fv[0]; 
            val[1] += in[ri*c+1] - fv[1]; 
            val[2] += in[ri*c+2] - fv[2]; 
            out[ti*c+0] = val[0]*iarr; 
            out[ti*c+1] = val[1]*iarr; 
            out[ti*c+2] = val[2]*iarr; 
            ri+=w;
            ti+=w;
        }

        for(int j=r+1; j<h-r; j++) 
        { 
            val[0] += in[ri*c+0] - in[li*c+0]; 
            val[1] += in[ri*c+1] - in[li*c+1]; 
            val[2] += in[ri*c+2] - in[li*c+2]; 
            out[ti*c+0] = val[0]*iarr; 
            out[ti*c+1] = val[1]*iarr; 
            out[ti*c+2] = val[2]*iarr; 
            li+=w; 
            ri+=w;
            ti+=w;
        }
        
        for(int j=h-r; j<h; j++) 
        { 
            val[0] += lv[0] - in[li*c+0]; 
            val[1] += lv[1] - in[li*c+1]; 
            val[2] += lv[2] - in[li*c+2]; 
            out[ti*c+0] = val[0]*iarr; 
            out[ti*c+1] = val[1]*iarr; 
            out[ti*c+2] = val[2]*iarr; 
            li+=w; 
            ti+=w;
        }
    }
}

//!
//! \fn void box_blur_rgb(float * in, float * out, int w, int h, int c, int r)   
//!
//! \brief this function performs a box blur pass. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void box_blur_rgb(float *& in, float *& out, int w, int h, int c, int r) 
{
    std::swap(in, out);
    horizontal_blur_rgb(out, in, w, h, c, r);
    total_blur_rgb(in, out, w, h, c, r);
    // Note to myself : 
    // here we could go anisotropic with different radiis rx,ry in HBlur and TBlur
}

//!
//! \fn void fast_gaussian_blur_rgb(float * in, float * out, int w, int h, int c, float sigma)   
//!
//! \brief this function performs a fast Gaussian blur. Applying several
//! times box blur tends towards a true Gaussian blur. Three passes are sufficient
//! for good results. For further details please refer to :  
//! http://blog.ivank.net/fastest-gaussian-blur.html
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] sigma        gaussian std dev
//!
void fast_gaussian_blur_rgb(float *& in, float *& out, int w, int h, int c, float sigma) 
{
    // sigma conversion to box dimensions
    int boxes[3];
    std_to_box(boxes, sigma, 3);
    box_blur_rgb(in, out, w, h, c, boxes[0]);
    box_blur_rgb(out, in, w, h, c, boxes[1]);
    box_blur_rgb(in, out, w, h, c, boxes[2]);
}

//! \code{.cpp}
int main(int argc, char * argv[])
{   
    if( argc < 2 ) exit(1);
    const char * image_file = argv[1];
    const float sigma = argc > 2 ? std::atof(argv[2]) : 3.;
    const char * output_file = argc > 3 ? argv[3] : "blur.png";
    
    // image loading
    int width, height, channels;
    unsigned char * image_data = stbi_load(argv[1], &width, &height, &channels, 0);
    std::cout << "Source image: " << width<<"x" << height << " ("<<channels<<")" << std::endl;
    if(channels < 3)
    {
        std::cout<< "Input images must be RGB images."<<std::endl;
        exit(1);
    }

    // copy data
    int size = width * height * channels;
    
    // output channels r,g,b
    float * new_image = new float[size];
    float * old_image = new float[size];

    // channels copy r,g,b
    for(int i = 0; i < size; ++i)
        old_image[i] = image_data[i] / 255.f;

    // per channel filter
    auto start = std::chrono::system_clock::now();
    fast_gaussian_blur_rgb(old_image, new_image, width, height, channels, sigma);
    auto end = std::chrono::system_clock::now();

    // stats
    float elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "time " << elapsed << "ms" << std::endl;

    // channels copy r,g,b
    for(int i = 0; i < size; ++i)
        image_data[i] = (unsigned char) std::min(255.f, std::max(0.f, 255.f * new_image[i]));

    // save
    std::string file(output_file);
    std::string ext = file.substr(file.size()-3);
    if( ext == "bmp" )
        stbi_write_bmp(output_file, width, height, channels, image_data);
    else if( ext == "jpg" )
        stbi_write_jpg(output_file, width, height, channels, image_data, 90);
    else
    {
        if( ext != "png" )
        {
            std::cout << "format '" << ext << "' not supported writing default .png" << std::endl; 
            file = file.substr(0, file.size()-4) + std::string(".png");
        }
        stbi_write_png(file.c_str(), width, height, channels, image_data, channels*width);
    }
    stbi_image_free(image_data);

    // clean memory
    delete[] new_image;
    delete[] old_image;
    return 0;
}
//! \endcode