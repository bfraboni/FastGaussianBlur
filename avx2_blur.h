#pragma once
#include <immintrin.h>
#include <algorithm>
#include <cstring>
#include "fast_gaussian_blur_template.h"

typedef unsigned char uchar;

// ---------------------------------------------------------
inline void blur_row_avx2(const uchar* in_row, short* out_row, int width, int radius) {
    int window_size = 2 * radius + 1;
    __m256i acc = _mm256_setzero_si256();

    for (int i = -radius; i <= radius; i++) {
        int safe_x = std::max(0, std::min(i, width - 1));
        int offset = safe_x * 3;
        // Safe 3-byte read to prevent edge-case segfaults
        int val = in_row[offset] | (in_row[offset+1] << 8) | (in_row[offset+2] << 16);
        __m128i vec = _mm_setr_epi32(val, 0, 0, 0); 
        __m256i px16 = _mm256_cvtepu8_epi16(vec);
        acc = _mm256_add_epi16(acc, px16);
    }

    for (int x = 0; x < width; x++) {
        short r = _mm256_extract_epi16(acc, 0);
        short g = _mm256_extract_epi16(acc, 1);
        short b = _mm256_extract_epi16(acc, 2);

        out_row[x * 3 + 0] = r;
        out_row[x * 3 + 1] = g;
        out_row[x * 3 + 2] = b;

        int left_x = std::max(0, x - radius);
        int right_x = std::min(width - 1, x + radius + 1);

        int out_val = in_row[left_x * 3] | (in_row[left_x * 3 + 1] << 8) | (in_row[left_x * 3 + 2] << 16);
        int in_val  = in_row[right_x * 3] | (in_row[right_x * 3 + 1] << 8) | (in_row[right_x * 3 + 2] << 16);

        __m128i out_vec = _mm_setr_epi32(out_val, 0, 0, 0);
        __m128i in_vec = _mm_setr_epi32(in_val, 0, 0, 0);

        acc = _mm256_sub_epi16(acc, _mm256_cvtepu8_epi16(out_vec));
        acc = _mm256_add_epi16(acc, _mm256_cvtepu8_epi16(in_vec));
    }
}

// ---------------------------------------------------------
void avx2_fused_pass(const uchar* in, uchar* out, int width, int height, int radius) {
    int window_size = 2 * radius + 1;
    float total_inverse = 1.0f / (window_size * window_size);
    __m256 f_mult = _mm256_set1_ps(total_inverse);

    // Pad memory allocations to the nearest multiple of 8 to prevent SIMD out-of-bounds reads
    int padded_channels = ((width * 3 + 7) / 8) * 8;
    
    int ring_size = window_size; 
    short** ring_buffer = new short*[ring_size];
    for (int i = 0; i < ring_size; i++) {
        ring_buffer[i] = new short[padded_channels]();
    }
    int* col_acc = new int[padded_channels](); 

    // --- PHASE 1: PRE-FILL THE ACCUMULATOR ---
    // Safely clamp and accumulate the top rows before the master loop starts
    for (int i = -radius; i <= radius; i++) {
        int safe_y = std::max(0, std::min(i, height - 1));
        int ring_idx = i + radius; 
        
        blur_row_avx2(in + safe_y * width * 3, ring_buffer[ring_idx], width, radius);
        
        for (int x = 0; x < width * 3; x += 8) {
            __m256i current_acc = _mm256_loadu_si256((__m256i*)(col_acc + x));
            __m128i in_16 = _mm_loadu_si128((__m128i*)(ring_buffer[ring_idx] + x));
            __m256i in_32 = _mm256_cvtepi16_epi32(in_16);
            current_acc = _mm256_add_epi32(current_acc, in_32);
            _mm256_storeu_si256((__m256i*)(col_acc + x), current_acc);
        }
    }

    // --- PHASE 2: THE MASTER LOOP ---
    for (int y = 0; y < height; y++) {
        int ring_idx = y % ring_size;

        // Step 2A: Write current accumulator to image AND subtract the outgoing row
        for (int x = 0; x < width * 3; x += 8) {
            __m256i current_acc = _mm256_loadu_si256((__m256i*)(col_acc + x));
            
            // 1. Math formatting
            __m256 f_acc = _mm256_cvtepi32_ps(current_acc);
            __m256 f_blurred = _mm256_mul_ps(f_acc, f_mult);
            __m256i i_blurred = _mm256_cvtps_epi32(f_blurred);
            
            // 2. Linear SSE Packing
            __m128i lo = _mm256_castsi256_si128(i_blurred);
            __m128i hi = _mm256_extracti128_si256(i_blurred, 1);
            __m128i pack16 = _mm_packus_epi32(lo, hi);
            __m128i pack8 = _mm_packus_epi16(pack16, pack16);
            
            // 3. Safe SIMD Tail Write (Prevents segfaults on unaligned image widths)
            if (x + 8 <= width * 3) {
                _mm_storel_epi64((__m128i*)(out + y * width * 3 + x), pack8);
            } else {
                uchar temp[8];
                _mm_storel_epi64((__m128i*)temp, pack8);
                int remaining = width * 3 - x;
                for (int j = 0; j < remaining; j++) {
                    out[y * width * 3 + x + j] = temp[j];
                }
            }

            // 4. Subtract the oldest row (which is leaving the window)
            __m128i out_16 = _mm_loadu_si128((__m128i*)(ring_buffer[ring_idx] + x));
            __m256i out_32 = _mm256_cvtepi16_epi32(out_16);
            current_acc = _mm256_sub_epi32(current_acc, out_32);
            
            _mm256_storeu_si256((__m256i*)(col_acc + x), current_acc);
        }

        // Step 2B: Horizontally blur the new incoming row, overwriting the old one
        int incoming_y = std::min(y + radius + 1, height - 1);
        blur_row_avx2(in + incoming_y * width * 3, ring_buffer[ring_idx], width, radius);

        // Step 2C: Add the newly blurred row to the accumulator
        for (int x = 0; x < width * 3; x += 8) {
            __m256i current_acc = _mm256_loadu_si256((__m256i*)(col_acc + x));
            __m128i in_16 = _mm_loadu_si128((__m128i*)(ring_buffer[ring_idx] + x));
            __m256i in_32 = _mm256_cvtepi16_epi32(in_16);
            current_acc = _mm256_add_epi32(current_acc, in_32);
            _mm256_storeu_si256((__m256i*)(col_acc + x), current_acc);
        }
    }

    // Cleanup
    for (int i = 0; i < ring_size; i++) delete[] ring_buffer[i];
    delete[] ring_buffer;
    delete[] col_acc;
}

// ---------------------------------------------------------
void fast_gaussian_blur_avx2(uchar*& in, uchar*& out, int width, int height, float sigma) {
    int boxes[3];
    sigma_to_box_radius(boxes, sigma, 3);
    for (int i = 0; i < 3; i++) {
        avx2_fused_pass(in, out, width, height, boxes[i]);
        std::swap(in, out);
    }
}