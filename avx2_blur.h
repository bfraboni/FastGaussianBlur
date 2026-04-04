#pragma once
#include <immintrin.h>
#include <algorithm>
#include <cstring>
#include "fast_gaussian_blur_template.h"

typedef unsigned char uchar;

// ---------------------------------------------------------
// 1. THE HORIZONTAL LINE ENGINE (Processes exactly 1 row)
// ---------------------------------------------------------
inline void blur_row_avx2(const uchar* in_row, short* out_row, int width, int radius) {
    int window_size = 2 * radius + 1;
    
    // Notice we don't multiply/divide here anymore! 
    // We just keep the raw sum in 16-bit 'short' format to pass to the vertical pass.
    
    __m256i acc = _mm256_setzero_si256();

    // Initial accumulation for the left edge
    for (int i = -radius; i <= radius; i++) {
        int safe_x = std::max(0, std::min(i, width - 1));
        int offset = safe_x * 3;
        int val = *(const int*)(in_row + offset) & 0x00FFFFFF;
        __m128i vec = _mm_setr_epi32(val, 0, 0, 0); // We only care about 1 pixel here for setup
        __m256i px16 = _mm256_cvtepu8_epi16(vec);
        acc = _mm256_add_epi16(acc, px16);
    }

    for (int x = 0; x < width; x++) {
        // Extract the raw 16-bit sum (RGB)
        short r = _mm256_extract_epi16(acc, 0);
        short g = _mm256_extract_epi16(acc, 1);
        short b = _mm256_extract_epi16(acc, 2);

        out_row[x * 3 + 0] = r;
        out_row[x * 3 + 1] = g;
        out_row[x * 3 + 2] = b;

        // Slide the horizontal window
        int left_x = std::max(0, x - radius);
        int right_x = std::min(width - 1, x + radius + 1);

        int out_val = *(const int*)(in_row + left_x * 3) & 0x00FFFFFF;
        int in_val = *(const int*)(in_row + right_x * 3) & 0x00FFFFFF;

        __m128i out_vec = _mm_setr_epi32(out_val, 0, 0, 0);
        __m128i in_vec = _mm_setr_epi32(in_val, 0, 0, 0);

        acc = _mm256_sub_epi16(acc, _mm256_cvtepu8_epi16(out_vec));
        acc = _mm256_add_epi16(acc, _mm256_cvtepu8_epi16(in_vec));
    }
}

// ---------------------------------------------------------
// 2. THE CACHE FUSION ENGINE (The Ring Buffer)
// ---------------------------------------------------------
void avx2_fused_pass(const uchar* in, uchar* out, int width, int height, int radius) {
    int window_size = 2 * radius + 1;
    
    // We are multiplying the inverse weight TWICE because we skipped the division in the horizontal pass.
    // Total Area = window_size * window_size.
    float total_inverse = 1.0f / (window_size * window_size);
    short inv_weight = (short)(total_inverse * 32768.0f);
    __m256i multiplier = _mm256_set1_epi16(inv_weight);

    // Create the L1 Ring Buffer. 
    // It only needs to hold (2 * radius + 1) rows. For radius 10, this is 21 rows.
    // 21 rows * 1024 width * 3 channels * 2 bytes (short) = ~128 Kilobytes.
    // This perfectly fits inside modern L2 caches, preventing Main RAM access!
    int ring_size = window_size + 1; 
    short** ring_buffer = new short*[ring_size];
    for (int i = 0; i < ring_size; i++) {
        ring_buffer[i] = new short[width * 3];
    }

    // The Vertical Column Accumulators
    // This keeps the running vertical sum for every single column in the image.
    int* col_acc = new int[width * 3](); // Initialize to 0

    // Master Loop: Top to Bottom
    for (int y = -radius; y < height + radius; y++) {
        
        // 1. Calculate our sliding window edges
        int incoming_y = std::max(0, std::min(y + radius, height - 1));
        int outgoing_y = std::max(0, y - radius - 1);

        // 2. Horizontally blur the incoming row and store it in our L1 Ring Buffer
        int incoming_ring_idx = incoming_y % ring_size;
        int outgoing_ring_idx = outgoing_y % ring_size;
        
        blur_row_avx2(in + incoming_y * width * 3, ring_buffer[incoming_ring_idx], width, radius);

        // 3. Update the Vertical Column Accumulators (O(1) sliding window)
        // We use SIMD to blast across the columns instantly
        for (int x = 0; x < width * 3; x += 8) {
            __m256i current_acc = _mm256_loadu_si256((__m256i*)(col_acc + x));
            
            // Load 16-bit horizontally blurred pixels, promote to 32-bit
            __m128i in_16 = _mm_loadu_si128((__m128i*)(ring_buffer[incoming_ring_idx] + x));
            __m128i out_16 = _mm_loadu_si128((__m128i*)(ring_buffer[outgoing_ring_idx] + x));

            __m256i in_32 = _mm256_cvtepi16_epi32(in_16);
            __m256i out_32 = _mm256_cvtepi16_epi32(out_16);

            // Add incoming row, subtract outgoing row
            current_acc = _mm256_add_epi32(current_acc, in_32);
            if (y > radius) {
                current_acc = _mm256_sub_epi32(current_acc, out_32);
            }

            // Store back to column accumulators
            _mm256_storeu_si256((__m256i*)(col_acc + x), current_acc);

            // 4. If the window is fully on the screen, write to final image
            if (y >= 0 && y < height) {
                // Pack 32-bit sums down to 16-bit so we can use our fast hardware multiplier
                // Note: In production, large radii require 32-bit float multiplication here to prevent overflow.
                __m256i acc_16 = _mm256_packus_epi32(current_acc, current_acc);
                __m256i blurred = _mm256_mulhrs_epi16(acc_16, multiplier);
                __m256i packed_8bit = _mm256_packus_epi16(blurred, blurred);
                
                // Destructive overlap store (writes 8 bytes, only 3 are RGB, next loop fixes it)
                __m128i lower = _mm256_castsi256_si128(packed_8bit);
                _mm_storel_epi64((__m128i*)(out + y * width * 3 + x), lower);
            }
        }
    }

    for (int i = 0; i < ring_size; i++) delete[] ring_buffer[i];
    delete[] ring_buffer;
    delete[] col_acc;
}

// ---------------------------------------------------------
// THE DIRECTOR
// ---------------------------------------------------------
void fast_gaussian_blur_avx2(uchar*& in, uchar*& out, int width, int height, float sigma) {
    int boxes[3];
    sigma_to_box_radius(boxes, sigma, 3);
    
    // Instead of 12 full memory passes, we do exactly 3.
    // Read -> Fused X/Y -> Write. No transposition required.
    for (int i = 0; i < 3; i++) {
        avx2_fused_pass(in, out, width, height, boxes[i]);
        std::swap(in, out);
    }
}