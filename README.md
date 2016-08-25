# aia
Place for sharing python code written for doing things with AIA data.

## Sun Intensity

The SunIntensity package's primary usage is for normalizing the brightness of AIA images so that they are at the standard brightness from the beginning of the mission. The sun_intensity.py modules contains all the code necessary for collecting and processing the original data, however this is rather computationally expensive and will take some time. That code is here mostly to show how the values are obtained. The best way to get the normalization values is to use the get_intensity function in the sun_intensity_lw module, which only contains the code needed for accessing the saved lightcurve in the [wavelength]regression.pkl file. I recommend saving the module somewhere and using imp to import it into your program.

```
import imp
sun_intensity = imp.load_source('sun_intensity', '/path/to/sun_intensity_lw.py')
```
