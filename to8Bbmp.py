import numpy as np

from PIL import Image
# if this doesn't work becaue of missing dependencies, you need to install the dependencies in the environment julia/pycall is using to run python
# you can find the environment by running 
# using PyCall
# PyCall.python
# on my ms windows machine it is C:\Users\micha\.julia\conda\3\x86_64
# Now i have to activate that environment and install pillow. 
# I use conda for creating packages, and vscode for editing, and the julia extension for vscode to run julia code, 
# so i just click the python version in the bottom left right of vscode and select the environment i want to use.
# then i open a terminal and type pip install pillow
# now it should work

def save_as_bmp(data, output_file):
    # Convert the data to a NumPy array and ensure it's in 8-bit format
    img_data = np.array(data).astype(np.uint8)

    # Create an image object from the array
    img = Image.fromarray(img_data, 'L')  # 'L' mode for grayscale

    # Save as 8-bit BMP
    img.save(output_file, format='BMP')


# data = [[0, 0, 0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0, 0, 0]]

# save_as_bmp(data, 'test.bmp')