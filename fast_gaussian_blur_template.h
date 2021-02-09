#pragma once

#include <type_traits>
#include <algorithm>
#include <cmath>

template<typename T, int C>
void horizontal_blur_template(T * in, T * out, const int w, const int h, const int r)
{
    float iarr = 1.f / (r+r+1);
    #pragma omp parallel for
    for(int i=0; i<h; i++) 
    {
        int ti = i*w; 
        int li = ti;  
        int ri = ti+r;

        std::conditional_t<std::is_integral<T>::value, int, float> fv[C], lv[C], val[C];

        for(int ch = 0; ch < C; ++ch)
        {
            fv[ch] = in[ti*C+ch];
            lv[ch] = in[(ti+w-1)*C+ch];
            val[ch] = (r+1)*fv[ch];
        }

        for(int j=0; j<r; j++) 
        for(int ch = 0; ch < C; ++ch)
        {
            val[ch] += in[(ti+j)*C+ch]; 
        }

        for(int j=0; j<=r; j++, ri++, ti++) 
        for(int ch = 0; ch < C; ++ch)
        { 
            val[ch] += in[ri*C+ch] - fv[ch]; 
            if( std::is_integral<T>::value )
                out[ti*C+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to integer type
            else
                out[ti*C+ch] = val[ch]*iarr;
        }

        for(int j=r+1; j<w-r; j++, ri++, ti++, li++) 
        for(int ch = 0; ch < C; ++ch)
        { 
            val[ch] += in[ri*C+ch] - in[li*C+ch]; 
            if( std::is_integral<T>::value )
                out[ti*C+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to integer type
            else
                out[ti*C+ch] = val[ch]*iarr;

        }

        for(int j=w-r; j<w; j++, ti++, li++) 
        for(int ch = 0; ch < C; ++ch)
        { 
            val[ch] += lv[ch] - in[li*C+ch]; 
            if( std::is_integral<T>::value )
                out[ti*C+ch] = val[ch]*iarr+0.5f; // avoid std::round by adding 0.5f and cast to integer type
            else
                out[ti*C+ch] = val[ch]*iarr;
        }
    }
}

template<typename T>
void horizontal_blur_template_dispatch(T * in, T * out, const int w, const int h, const int c, const int r)
{
    switch(c)
    {
        case 1: horizontal_blur_template<T,1>(in, out, w, h, r); break;
        case 2: horizontal_blur_template<T,2>(in, out, w, h, r); break;
        case 3: horizontal_blur_template<T,3>(in, out, w, h, r); break;
        case 4: horizontal_blur_template<T,4>(in, out, w, h, r); break;
        default: printf("%d channels is not supported yet. Add a specific case if possible or fall back to the generic version.", c); break;
    }
}

template<typename T, int C>
void flip_block_template(T * in, T * out, const int w, const int h)
{
    constexpr int block = 256;
    #pragma omp parallel for collapse(2)
    for(int x= 0; x < w; x+= block)
    {
        for(int y= 0; y < h; y+= block)
        {
            T * p = in + y*w*C + x*C;
            T * q = out + y*C + x*h*C;
            
            const int blockx= std::min(w, x+block) - x;
            const int blocky= std::min(h, y+block) - y;
            for(int xx= 0; xx < blockx; xx++)
            {
                for(int yy= 0; yy < blocky; yy++)
                {
                    for(int k= 0; k < C; k++)
                        q[k]= p[k];
                    //~ std::memcpy(q, p, C);
                    p+= w*C;
                    q+= C;                    
                }
                // repositionne les pointeurs sur le prochain pixel
                p+= -blocky*w*C + C;
                q+= -blocky*C + h*C;
            }
        }
    }
}

template<typename T>
void flip_block_template_dispatch(T * in, T * out, const int w, const int h, const int c)
{
    switch(c)
    {
        case 1: flip_block_template<T,1>(in, out, w, h); break;
        case 2: flip_block_template<T,2>(in, out, w, h); break;
        case 3: flip_block_template<T,3>(in, out, w, h); break;
        case 4: flip_block_template<T,4>(in, out, w, h); break;
        default: printf("%d channels is not supported yet. Add a specific case if possible or fall back to the generic version.", c); break;
    }
}

template<typename T>
void fast_gaussian_blur_template(T *& in, T *& out, int w, int h, int c, float sigma) 
{
    int n = 3;
    int boxes[3];
    sigma_to_box_radius(boxes, sigma, n);

    horizontal_blur_template_dispatch(in, out, w, h, c, boxes[0]);
    horizontal_blur_template_dispatch(out, in, w, h, c, boxes[1]);
    horizontal_blur_template_dispatch(in, out, w, h, c, boxes[2]);
    flip_block_template_dispatch(out, in, w, h, c);
    
    horizontal_blur_template_dispatch(in, out, h, w, c, boxes[0]);
    horizontal_blur_template_dispatch(out, in, h, w, c, boxes[1]);
    horizontal_blur_template_dispatch(in, out, h, w, c, boxes[2]);
    flip_block_template_dispatch(out, in, h, w, c);
    
    std::swap(in, out);    
}