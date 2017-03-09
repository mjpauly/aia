"""image producing routine for AIA
"""

from datetime import datetime
import matplotlib.colors as colors
import numpy as np
from scipy import misc
from skimage.transform import downscale_local_mean
import os
from subprocess import Popen
from shutil import rmtree
from matplotlib.cm import get_cmap

from sunpy.map import Map
from sunpy.instr.aia import aiaprep

import sun_intensity

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


STANDARD_INT = {'131': 6.99685, '171': 4.99803, '193': 2.9995,
                '211': 4.99801, '304': 4.99941, '335': 6.99734,
                '94': 4.99803}
SQRT_NORM = {'131': False, '171': True, '193': False, '211': False,
             '304': False, '335': False, '94': True}
MINMAX = {'131': (7, 1200),'171': (10, 6000),'193': (120, 6000),
          '211': (30, 13000), '304': (50, 2000),'335': (3.5, 1000),
          '94': (1.5, 50)}

def process_img(fits_file, fname=None, downscale=None,
                rescale_brightness=True, side_by_side=False,
                timestamp=True):
        """Produces a png image of the Sun from a fits file

        fits_file: name of fits file to process

        fname: determines whether to save the image directly
        or to return a pil image

        downscale: if filled, it downscales the data by (x, y) factor.
        E.g. (8,8)

        rescale_brightness: determines if brightness correction is done.

        side_by_side: make a side-by-side comparison of the scaled
        and not brightness scaled images

        timestamp: show timestamp or not
        """
        hdr = sun_intensity.getFitsHdr(fits_file)
        wavelength = str(hdr['wavelnth'])
        exptime = hdr['EXPTIME']
        cmap = get_cmap('sdoaia' + wavelength)
        cmap.set_bad()
        imin, imax = MINMAX[wavelength]

        themap = Map(fits_file)
        if hdr['lvl_num'] != 1.5:
                # perform aiaprep if data not at level 1.5
                themap = aiaprep(themap)
        data = themap.data / exptime #  normalize for exposure
        norm_scale = STANDARD_INT[wavelength]
        dim_factor = sun_intensity.get_dim_factor(themap.date, wavelength)
        data = data * norm_scale
        data = np.flipud(data)
        if downscale:
                data = downscale_local_mean(data, downscale)

        if rescale_brightness or side_by_side:
                imin = imin / dim_factor # brightness correction
                imax = imax / dim_factor
        data[0,0] = imin # first pixel set to min
        data[0,1] = imax # second pixel sit to max

        if SQRT_NORM[wavelength]:
                norm = colors.PowerNorm(1)
                data = np.sqrt(np.clip(data, imin, imax))
        else:
                norm = colors.LogNorm(vmin=imin, vmax=imax, clip=True)
        pil_img = misc.toimage(cmap(norm(data)))

        width, height = pil_img.size
        if side_by_side:
                new_img = Image.new('RGB', (width * 2, height))
                new_img.paste(pil_img, (0, 0))
                second_image = process_img(fits_file, downscale=downscale,
                                           rescale_brightness=False,
                                           timestamp=False)
                new_img.paste(second_image, (width, 0))
                pil_img = new_img

        if timestamp:
                draw = ImageDraw.Draw(pil_img)
                font_height = int(height / 64)
                font = ImageFont.truetype('/Library/Fonts/Arial.ttf',
                                          font_height)
                draw.text((font_height, height - (2 * font_height)),
                          'SDO/AIA- ' + wavelength + ' ' +
                          themap.date.strftime('%Y-%m-%d %H:%M:%S'),
                          font=font)

        if fname:
                pil_img.save(fname)
        else:
                return pil_img


def process_hmi(fits_file, rsun_obs, cdelt, fname=None,
                downscale=None, timestamp=True):
        """Produces a png image of the Sun from a fits file

        fits_file: name of fits file to process

        fname: determines whether to save the image directly
        or to return a pil image

        downscale: if filled, it downscales the data by (x, y) factor.
        E.g. (8,8)

        timestamp: show timestamp or not
        """
        hdr = sun_intensity.getFitsHdr(fits_file)
        # cmap = get_cmap('hmimag')
        # cmap.set_bad()
        themap = Map(fits_file)
        data = themap.data
        data = np.flipud(data)
        r_pix = (rsun_obs / cdelt) - 10 # -10 for testing
        import ipdb
        ipdb.set_trace()
        mask = sun_intensity.get_disk_mask(data.shape, r_pix)
        data[mask] = 0
        if downscale:
                data = downscale_local_mean(data, downscale)

        norm = colors.LogNorm(1)
        cmap = mpl_cmap('Greys_r')
        pil_img = misc.toimage(cmap(norm(data)))
        #pil_img = misc.toimage(cmap(data))
        

        width, height = pil_img.size
        if timestamp:
                draw = ImageDraw.Draw(pil_img)
                font_height = int(height / 64)
                font = ImageFont.truetype('/Library/Fonts/Arial.ttf',
                                          font_height)
                draw.text((font_height, height - (2 * font_height)),
                          'SDO/HMI- Magnetogram' +
                          themap.date.strftime('%Y-%m-%d %H:%M:%S'),
                          font=font)

        if fname:
                pil_img.save(fname)
        else:
                return pil_img


def make_movie(fits_list, movname='outfile.mov', framerate=60, **kwargs):
        """Produces a movie from the list of fits files provided
        **kwargs pases args for each frame to mov_img
        """
        if not os.path.exists('/tmp/'):
                os.mkdir('/tmp/')
        if not os.path.exists('/tmp/aia_movie/'):
                os.mkdir('/tmp/aia_movie/')
        for i, fits_file in enumerate(fits_list):
                process_img(fits_file,
                            fname='/tmp/aia_movie/{:05d}.jpg'.format(i),
                            **kwargs)
        pop = Popen(['ffmpeg', '-y', '-framerate', str(framerate), '-i',
                     '/tmp/aia_movie/%05d.jpg', movname])
        pop.wait()
        rmtree('/tmp/aia_movie/')
