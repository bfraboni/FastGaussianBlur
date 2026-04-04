#!/bin/bash

# Ensure output directory exists
mkdir -p data/edge_cases

echo "====================================="
echo " Generating Test Images..."
echo "====================================="

# Standard Edge Cases
convert -size 1001x1001 xc:red -type TrueColor data/edge_cases/unaligned.png
convert -size 25x25 xc:yellow -type TrueColor data/edge_cases/tiny.png
convert -size 1920x10 xc:green -type TrueColor data/edge_cases/wide.png
convert -size 10x1920 xc:blue -type TrueColor data/edge_cases/tall.png
convert -size 3840x2160 xc:magenta -type TrueColor data/edge_cases/4k.png

# Nightmare Edge Cases
convert -size 1024x1024 xc:cyan -type TrueColor data/edge_cases/overflow.png
convert -size 2x2 xc:black -type TrueColor data/edge_cases/micro.png
convert -size 7680x4320 xc:orange -type TrueColor data/edge_cases/8k.png

# Safety check
if [ ! -f "data/edge_cases/unaligned.png" ]; then
    echo "ERROR: Images were not generated!"
    exit 1
fi

echo ""
echo "====================================="
echo " FINAL PRODUCTION BENCHMARKS"
echo "====================================="

run_test() {
    TITLE=$1
    FILE=$2
    SIGMA=$3

    echo "----------------------------------------------------------------"
    echo "$TITLE"
    echo "----------------------------------------------------------------"
    
    echo ">>> ORIGINAL (fastblur):"
    ./fastblur "$FILE" "out_orig_${FILE##*/}" "$SIGMA" | grep -E "Time" || echo "CRASHED OR FAILED"
    
    echo ">>> OPTIMIZED (fastblur_opt):"
    ./fastblur_opt "$FILE" "out_opt_${FILE##*/}" "$SIGMA" | grep -E "Time|Routing" || echo "CRASHED OR FAILED"
    echo ""
}

# Standard tests
run_test "1. UNALIGNED WIDTH (1001x1001, Sigma 10)" "data/edge_cases/unaligned.png" 10
run_test "2. EXTREME SIGMA (25x25, Sigma 50)" "data/edge_cases/tiny.png" 50
run_test "3. EXTREME WIDE (1920x10, Sigma 10)" "data/edge_cases/wide.png" 10
run_test "4. EXTREME TALL (10x1920, Sigma 10)" "data/edge_cases/tall.png" 10
run_test "5. 4K CACHE BUSTER (3840x2160, Sigma 15)" "data/edge_cases/4k.png" 15

# Nightmare tests
run_test "6. THE 16-BIT OVERFLOW TRAP (1024x1024, Sigma 100)" "data/edge_cases/overflow.png" 100
run_test "7. THE MICRO-MEMORY SEGFAULT TRAP (2x2, Sigma 5)" "data/edge_cases/micro.png" 5
run_test "8. THE ZERO-SIGMA PARADOX (1001x1001, Sigma 0.1)" "data/edge_cases/unaligned.png" 0.1
run_test "9. 8K L3 CACHE SPILLER (7680x4320, Sigma 10)" "data/edge_cases/8k.png" 10

echo "Done."