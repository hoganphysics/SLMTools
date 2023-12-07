import numpy as np
from PIL import Image

def save_as_bmp(data, output_file):
    # Convert the data to a NumPy array and ensure it's in 8-bit format
    img_data = np.array(data).astype(np.uint8)

    # Create an image object from the array
    img = Image.fromarray(img_data, 'L')  # 'L' mode for grayscale

    # Save as 8-bit BMP
    img.save(output_file, format='BMP')
