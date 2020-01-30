from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# image = fits.open("dat/b62_g_5s_end_001.fits")    # makes a HDU list object
#
# image.info()
#
# for key in image[0].header:
#     print(key, image[0].header[key])

from reduction import write_out_fits

test_array = np.zeros((500,500))
test_array[0:250,0:250] = 1.0
plt.imshow(test_array)
plt.show()

write_out_fits(test_array, "test_array.fits")
