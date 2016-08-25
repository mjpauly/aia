# aia
Place for sharing python code written for doing things with AIA data.

## Sun Intensity

The SunIntensity package's primary usage is for normalizing the brightness of AIA images so that they are at the standard brightness from the beginning of the mission. The sun_intensity.py modules contains all the code necessary for collecting and processing the original data, however this is rather computationally expensive and will take some time. That code is here mostly to show how the values are obtained. The best way to get the normalization values is to use the get_intensity function in the sun_intensity_lw module, which only contains the code needed for accessing the saved lightcurve in the regression file. I recommend saving the module somewhere and using `imp` to import it into your program.

```
import imp
sun_intensity = imp.load_source('sun_intensity', '/path/to/sun_intensity_lw.py')
```

Use the get_intensity functions to get the brightness value for a given day. This value is the median brightness of the disk of the Sun, normaized for exposure.

```
from datetime import datetime
wavelength = '304'
mission_start_brightness = sun_intensity.get_intensity(datetime(2010,5,1), wavelength)
today_brightness = sun_intensity.get_intensity(datetime.utcnow(), wavelength)
dim_factor = mission_start_brightness / today brightness
```

This dim factor can then be used to set the maximun and minimum limits for pixel values of an AIA image. I will upload the python version of the AIA image producing routine soon.
