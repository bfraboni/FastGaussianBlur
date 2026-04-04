#!/bin/bash

mkdir -p data/real_images

echo "Downloading reliable test images from Picsum..."
# 1. Nature/Landscape
wget -q "https://picsum.photos/id/1018/1024/768.jpg" -O data/real_images/landscape.jpg
# 2. Animal (Pug in a blanket)
wget -q "https://picsum.photos/id/1025/1024/768.jpg" -O data/real_images/macaw.jpg
# 3. Architecture (Castle)
wget -q "https://picsum.photos/id/1040/1024/768.jpg" -O data/real_images/architecture.jpg

# Safety check: Make sure they actually downloaded as JPEG image data, not HTML errors
if ! file data/real_images/landscape.jpg | grep -q "JPEG"; then
    echo "ERROR: Download failed. The files are not valid JPEGs."
    exit 1
fi

echo "Images downloaded successfully! Running Blur Tests (Sigma = 10)..."

for img in landscape macaw architecture; do
    echo "Processing $img..."
    
    # Run Original
    ./fastblur data/real_images/${img}.jpg out_orig_${img}.png 0 > /dev/null
    
    # Run Optimized
    ./fastblur_opt data/real_images/${img}.jpg out_opt_${img}.png 0 > /dev/null
done

echo "Done! Check your folder and visually compare out_orig_*.png with out_opt_*.png"