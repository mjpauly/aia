"""image producing routine for AIA
"""

from datetime import datetime
import matplotlib.colors as colors
import numpy as np
from scipy import misc

import sunpy.map
import sun_intensity

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


STANDARD_INT = {'131': 6.99685, '171': 4.99803, '193': 2.9995,
		'211': 4.99801, '304': 4.99941, '335': 6.99734, '94': 4.99803}
SQRT_NORM = {'131': False, '171': True, '193': False, '211': False,
				'304': False, '335': False, '94': True}
MINMAX = {'131': (7, 1200),'171': (10, 6000),'193': (120, 6000),
		'211': (30, 13000), '304': (50, 2000),'335': (3.5, 1000),
		'94': (1.5, 50)}

def process_img(fits_file, fname=None):
	"""Produces a png image of the Sun from a fits file
	Optional kwarg fname determining whether to save the image directly
	or to return a pil image
	"""
	hdr = sun_intensity.getFitsHdr(fits_file)
	wavelength = str(hdr['wavelnth'])
	exptime = hdr['EXPTIME']
	cmap = sunpy.cm.get_cmap('sdoaia' + wavelength)
	cmap.set_bad()
	imin, imax = MINMAX[wavelength]

	themap = sunpy.map.Map(fits_file)
	data = themap.data / exptime #  normalize for exposure
	norm_scale = STANDARD_INT[wavelength]
        dim_factor = sun_intensity.get_dim_factor(themap.date, wavelength)
	data = data * norm_scale
	data = np.flipud(data)

	imin = imin / dim_factor #  adjust min limit
	imax = imax / dim_factor
	data[0,0] = imin
	data[0,1] = imax

	if SQRT_NORM[wavelength]:
		norm = colors.PowerNorm(1)
		data = np.sqrt(np.clip(data, imin, imax))
	else:
		norm = colors.LogNorm(vmin=imin, vmax=imax, clip=True)
	pil_img = misc.toimage(cmap(norm(data)))

	draw = ImageDraw.Draw(pil_img)
	font_height = int(len(data) / 64)
	font = ImageFont.truetype('/Library/Fonts/Arial.ttf', font_height)
	draw.text((font_height, len(data) - (2 * font_height)),
			'SDO/AIA- ' + wavelength + ' ' +
			themap.date.strftime('%Y-%m-%d %H:%M:%S'), font=font)
	if fname:
		pil_img.save(fname)
	else:
		return pil_img
