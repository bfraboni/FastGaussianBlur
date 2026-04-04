#include <iostream>
#include <immintrin.h>

int main() {
    const int width = 10;
    
    // Simulate 8 rows of an image, each 10 pixels wide. 
    // We initialize them with 1.0 just for testing.
    float rows[8][width];
    float output[8][width];
    
    for(int i = 0; i < 8; i++) {
        for(int j = 0; j < width; j++) {
            rows[i][j] = 1.0f;
        }
    }

    // This register holds the running sum for all 8 rows. Initialize it to 0.
    __m256 accumulators = _mm256_setzero_ps();

    // Slide horizontally across the width
    for(int x = 0; x < width; x++) {
        
        // 1. Manually gather 1 pixel from each of the 8 rows at the current column (x)
        __m256 incoming = _mm256_set_ps(
            rows[7][x], rows[6][x], rows[5][x], rows[4][x],
            rows[3][x], rows[2][x], rows[1][x], rows[0][x]
        );

        // 2. Add the incoming pixels to our accumulators
        accumulators = _mm256_add_ps(accumulators, incoming);

        // 3. Extract the accumulated values and store them in the output
        // (Normally we would multiply by an inverse radius here for the blur)
        float temp_out[8];
        _mm256_storeu_ps(temp_out, accumulators);
        
        for(int i = 0; i < 8; i++) {
            output[i][x] = temp_out[i];
        }
    }

    // Print Row 0 to see the sliding sum working
    std::cout << "Row 0 Accumulated Output:\n";
    for(int x = 0; x < width; x++) {
        std::cout << output[0][x] << " ";
    }
    std::cout << std::endl;

    return 0;
}