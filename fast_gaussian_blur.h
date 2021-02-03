// Copyright (C) 2017-2021 Basile Fraboni
// Copyright (C) 2014 Ivan Kutskir
// All Rights Reserved
// You may use, distribute and modify this code under the
// terms of the MIT license. For further details please refer 
// to : https://mit-license.org/
//

#include <omp.h>

//!
//! \file fast_gaussian_blur.h
//! \author Basile Fraboni
//! \date 2021
//!
//! \brief The software is a C++ implementation of a fast  
//! Gaussian blur algorithm by Ivan Kutskir. For further details 
//! please refer to : 
//! http://blog.ivank.net/fastest-gaussian-blur.html
//!

typedef unsigned char uchar;

void sigma_to_box_radius(int boxes[], float sigma, int n);  
int sigma_to_kernel_radius(float sigma);
float gaussian( float x, float mu, float sigma );
void make_gaussian_kernel(float sigma, float *& kernel, int& k);
void transpose(uchar * in, uchar * out, int w, int h, int c);

void horizontal_blur(float * in, float * out, int w, int h, int c, int r); 
void horizontal_blur(uchar * in, uchar * out, int w, int h, int c, int r); 

void total_blur(float * in, float * out, int w, int h, int c, int r); 
void total_blur(uchar * in, uchar * out, int w, int h, int c, int r); 

void fast_gaussian_blur(float *& in, float *& out, int w, int h, int c, float sigma); 
void fast_gaussian_blur(uchar *& in, uchar *& out, int w, int h, int c, float sigma); 
void fast_gaussian_blur_transpose(uchar *& in, uchar *& out, int w, int h, int c, float sigma); 

void gaussian_blur(float * in, float * out, int w, int h, int c, float * kernel, int k);
void gaussian_blur(uchar * in, uchar * out, int w, int h, int c, float * kernel, int k);

//!
//! \fn void sigma_to_box_radius(float boxes[], float sigma, int n)  
//!
//! \brief this function converts the standard deviation of 
//! Gaussian blur into box radius for each box blur pass. For 
//! further details please refer to :
//! https://www.peterkovesi.com/matlabfns/#integral
//! https://www.peterkovesi.com/papers/FastGaussianSmoothing.pdf
//!
//! \param[out] boxes   box radiis
//! \param[in] sigma    Gaussian standard deviation
//! \param[in] n        number of box blur pass
//!
void sigma_to_box_radius(int boxes[], float sigma, int n)  
{
    // ideal filter width
    float wi = std::sqrt((12*sigma*sigma/n)+1); 
    int wl = wi; // no need std::floor  
    if(wl%2==0) wl--;
    int wu = wl+2;
                
    float mi = (12*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
    int m = mi+0.5f; // replaced std::round by adding 0.5f + cast ton integer
                
    for(int i=0; i<n; i++) 
        boxes[i] = ((i < m ? wl : wu) - 1) / 2;
}

//! \fn int sigma_to_kernel_radius(float sigma)
int sigma_to_kernel_radius(float sigma)
{
    // capture almost all the gaussian extent
    return std::ceil(2.57f*sigma);
}

//! \fn float gaussian( float x, float mu, float sigma ) 
float gaussian( float x, float mu, float sigma ) 
{
    const float a = ( x - mu ) / sigma;
    return std::exp( -0.5f * a * a );
}

//! \fn void make_gaussian_kernel(float sigma, float *& kernel, int& k) 
void make_gaussian_kernel(float sigma, float *& kernel, int& k) 
{
    // kernel radius from sigma
    k = sigma_to_kernel_radius(sigma);
    // kernel alloc
    kernel = new float[(2*k+1)*(2*k+1)];
    // unnormalized kernel weights 
    for (int row = 0, index = 0; row < 2*k+1; row++)
        for (int col = 0; col < 2*k+1; col++, index++) 
            kernel[index] = gaussian(row, k, sigma) * gaussian(col, k, sigma);
    // we do not need to normalize since it will be done during the convolution
}

//!
//! \fn void horizontal_blur(float * in, float * out, int w, int h, int c, int r)   
//! \fn void horizontal_blur(uchar * in, uchar * out, int w, int h, int c, int r)   
//! \fn void horizontal_blur_no_round(uchar * in, uchar * out, int w, int h, int c, int r)   
//!
//! \brief These functions perform the horizontal blur pass for box blur. 
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] r            box dimension
//!
void horizontal_blur(float * in, float * out, int w, int h, int c, int r) 
{
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<h; i++) 
    {
        int ti = i*w; 
        int li = ti;  
        int ri = ti+r;

        float fv[4], lv[4], val[4]; // fixed max to 4 channels 
        for(int ch = 0; ch < c; ++ch)
        {
            fv[ch] = in[ti*c+ch];
            lv[ch] = in[(ti+w-1)*c+ch];
            val[ch] = (r+1)*fv[ch];
        }

        for(int j=0; j<r; j++) 
        for(int ch = 0; ch < c; ++ch)
        {
            val[ch] += in[(ti+j)*c+ch]; 
        }

        for(int j=0; j<=r; j++, ri++, ti++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - fv[ch]; 
            out[ti*c+ch] = val[ch]*iarr; 
        }

        for(int j=r+1; j<w-r; j++, ri++, ti++, li++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr; 
        }

        for(int j=w-r; j<w; j++, ti++, li++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += lv[ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr; 
        }
    }
}

void horizontal_blur(uchar * in, uchar * out, int w, int h, int c, int r) 
{
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<h; i++) 
    {
        int ti = i*w; 
        int li = ti;  
        int ri = ti+r;

        int fv[4], lv[4], val[4]; // fixed max to 4 channels 
        for(int ch = 0; ch < c; ++ch)
        {
            fv[ch] = in[ti*c+ch];
            lv[ch] = in[(ti+w-1)*c+ch];
            val[ch] = (r+1)*fv[ch];
        }

        for(int j=0; j<r; j++) 
        for(int ch = 0; ch < c; ++ch)
        {
            val[ch] += in[(ti+j)*c+ch]; 
        }

        for(int j=0; j<=r; j++, ri++, ti++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - fv[ch]; 
            out[ti*c+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to uchar
        }

        for(int j=r+1; j<w-r; j++, ri++, ti++, li++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to uchar
        }

        for(int j=w-r; j<w; j++, ti++, li++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += lv[ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to uchar
        }
    }
}

//!
//! \fn void total_blur(float * in, float * out, int w, int h, int c, int r)   
//! \fn void total_blur(uchar * in, uchar * out, int w, int h, int c, int r)   
//! \fn void total_blur_no_round(uchar * in, uchar * out, int w, int h, int c, int r)   
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
void total_blur(float * in, float * out, int w, int h, int c, int r) 
{
    // radius range on either side of a pixel + the pixel itself
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<w; i++) 
    {
        int ti = i;
        int li = ti;
        int ri = ti+r*w;

        float fv[4], lv[4], val[4]; // fixed max to 4 channels 
        for(int ch = 0; ch < c; ++ch)
        {
            fv[ch] = in[ti*c+ch];
            lv[ch] = in[(ti+w*(h-1))*c+ch];
            val[ch] = (r+1)*fv[ch];
        }

        for(int j=0; j<r; j++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[(ti+j*w)*c+ch]; 
        }

        for(int j=0; j<=r; j++, ri+=w, ti+=w) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - fv[ch]; 
            out[ti*c+ch] = val[ch]*iarr; 
        }

        for(int j=r+1; j<h-r; j++, ri+=w, ti+=w, li+=w) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr; 
        }
        
        for(int j=h-r; j<h; j++, ti+=w, li+=w) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += lv[ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr; 
        }
    }
}

void total_blur(uchar * in, uchar * out, int w, int h, int c, int r) 
{
    // radius range on either side of a pixel + the pixel itself
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<w; i++) 
    {
        int ti = i;
        int li = ti;
        int ri = ti+r*w;

        int fv[4], lv[4], val[4]; // fixed max to 4 channels 
        for(int ch = 0; ch < c; ++ch)
        {
            fv[ch] = in[ti*c+ch];
            lv[ch] = in[(ti+w*(h-1))*c+ch];
            val[ch] = (r+1)*fv[ch];
        }

        for(int j=0; j<r; j++) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[(ti+j*w)*c+ch]; 
        }

        for(int j=0; j<=r; j++, ri+=w, ti+=w) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - fv[ch]; 
            out[ti*c+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to uchar
        }

        for(int j=r+1; j<h-r; j++, ri+=w, ti+=w, li+=w) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += in[ri*c+ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to uchar
        }
        
        for(int j=h-r; j<h; j++, ti+=w, li+=w) 
        for(int ch = 0; ch < c; ++ch)
        { 
            val[ch] += lv[ch] - in[li*c+ch]; 
            out[ti*c+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to uchar
        }
    }
}

//!
//! \fn void fast_gaussian_blur(float *& in, float *& out, int w, int h, int c, float sigma)   
//! \fn void fast_gaussian_blur(uchar *& in, uchar *& out, int w, int h, int c, float sigma)   
//! \fn void fast_gaussian_blur_no_round(uchar * in, uchar * out, int w, int h, int c, float sigma)   
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
void fast_gaussian_blur(float *& in, float *& out, int w, int h, int c, float sigma) 
{
    // sigma conversion to box dimensions
    int boxes[3];
    sigma_to_box_radius(boxes, sigma, 3);

    std::swap(in, out);
    // Note: Hblur or Tblur passes order does not matter 
    // Note: Anisotropy is possible with different boxes size for TBlur and HBlur 
    // Note; Maybe it is possible to perform all HBlur passes at once, same for TBlur  
    horizontal_blur(out, in, w, h, c, boxes[0]);
    horizontal_blur(in, out, w, h, c, boxes[1]);
    horizontal_blur(out, in, w, h, c, boxes[2]);

    total_blur(in, out, w, h, c, boxes[0]);
    total_blur(out, in, w, h, c, boxes[1]);
    total_blur(in, out, w, h, c, boxes[2]);
}

void fast_gaussian_blur(uchar *& in, uchar *& out, int w, int h, int c, float sigma) 
{
    // sigma conversion to box dimensions
    int boxes[3];
    sigma_to_box_radius(boxes, sigma, 3);
    
    std::swap(in, out);
    // Note: Hblur or Tblur passes order does not matter 
    // Note: Anisotropy is possible with different boxes size for TBlur and HBlur 
    // Note; Maybe it is possible to perform all HBlur passes at once, same for TBlur  
    horizontal_blur(out, in, w, h, c, boxes[0]);
    horizontal_blur(in, out, w, h, c, boxes[1]);
    horizontal_blur(out, in, w, h, c, boxes[2]);

    total_blur(in, out, w, h, c, boxes[0]);
    total_blur(out, in, w, h, c, boxes[1]);
    total_blur(in, out, w, h, c, boxes[2]);
}

void transpose(uchar * in, uchar * out, int w, int h, int c)
{
    const std::size_t b = w * h - 1;
    #pragma omp parallel for
    for(std::size_t i = 0; i < b+1; i++)
    {
        std::size_t o = i == b ? b : h * i % b; // transpose index
        for(int ch = 0; ch < c; ch++)
            out[o*c+ch] = in[i*c+ch];
    }
}

void fast_gaussian_blur_transpose(uchar *& in, uchar *& out, int w, int h, int c, float sigma) 
{
    // sigma conversion to box dimensions for 3 passes
    int n = 3;
    int boxes[3];
    sigma_to_box_radius(boxes, sigma, n);

    horizontal_blur(in, out, w, h, c, boxes[0]);
    horizontal_blur(out, in, w, h, c, boxes[1]);
    horizontal_blur(in, out, w, h, c, boxes[2]);

    transpose(out, in, w, h, c);

    horizontal_blur(in, out, h, w, c, boxes[0]);
    horizontal_blur(out, in, h, w, c, boxes[1]);
    horizontal_blur(in, out, h, w, c, boxes[2]);

    transpose(out, in, h, w, c);
    // swap pointers    
    std::swap(in, out);
}

//!
//! \fn void gaussian_blur(float * in, float * out, int w, int h, int c, float * kernel, int k)
//! \fn void gaussian_blur(uchar * in, uchar * out, int w, int h, int c, float * kernel, int k)
//!
//! \brief this function performs a Gaussian blur.
//!
//! \param[in,out] in       source channel
//! \param[in,out] out      target channel
//! \param[in] w            image width
//! \param[in] h            image height
//! \param[in] c            image channels
//! \param[in] kernel       kernel weights array (square)
//! \param[in] k            kernel radius (kernel size is 2*k+1)
template<typename T>
void gaussian_blur(T * in, T * out, int w, int h, int c, float * kernel, int k)
{
    #pragma omp parallel for
    for (int row = 0; row < h; row++)
    for (int col = 0; col < w; col++) 
    {
        float val[4] = {0, 0, 0, 0}; // fixed max to 4 channels 
        float sum = 0;
        for(int irow = row-k, kindex = 0; irow <= row+k; ++irow)
        for(int icol = col-k; icol <= col+k; ++icol, ++kindex)
        {
            if(irow < 0 || icol < 0 || irow > h-1 || icol > w-1) 
                continue;
            
            int iindex = (irow*w+icol)*c;
            float weight = kernel[kindex];
            sum += weight;
            for(int ch = 0; ch < c; ++ch)
                val[ch] += in[iindex+ch] * weight;
        }
        for(int ch = 0; ch < c; ++ch)
            out[(row*w+col)*c+ch] = val[ch]/sum;
    }
}