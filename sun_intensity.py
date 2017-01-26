"""Module for finding the brightness values for AIA images using the date
"""

import mov_img
from datetime import datetime, timedelta
from fetcher import fetch
from sunpy.time import parse_time
from astropy.io import fits
from sunpy.map import Map
from sunpy.instr.aia import aiaprep
import numpy as np
import multiprocessing as mp
from scipy.signal import savgol_filter
import os
import imp
pb0r = imp.load_source('pb0r', os.path.join(os.path.dirname(__file__), 'pb0r.py'))
import pandas as pd


csv_path = os.path.join(os.path.dirname(__file__), 'aia_rescaling_data.csv')
wavelengths = ['131','171','193','211','304','335','94']


def create_date_series(tstart):
        if isinstance(tstart, tuple):
                dt = datetime(*tstart)
        else:
                dt = tstart
        end = datetime.utcnow()
        step = timedelta(days=1)
        result = []
        while dt < end:
                result.append(dt)
                dt += step
        return result


# create list of datetimes and date strings
# starts at minute 1 to prevent hitting a leap second
datetime_list = create_date_series((2010,5,1, 0, 1))
date_list = [str(date.date()) for date in datetime_list]


def getFitsHdr(fle):
        """get the header for a fits file
        """
        f = fits.open(fle)
        hdr = f[-1].header
        f.close()
        return hdr


def get_disk_mask(data, r_pix):
        """returns the array mask for only the disk of the Sun
        """
        x, y = np.meshgrid(*map(np.arange, data.shape), indexing='ij')
        return (np.sqrt((x - 2047.5)**2 + (y - 2047.5)**2) > r_pix)


def process_med_int(fle):
        """Processes 1 image and extracts the median intensity on the disk
        normalized for exposure time.
        """
        amap = Map(fle)
        amap = aiaprep(amap)
        data = amap.data

        date = amap.date
        radius = pb0r.pb0r(date, arcsec=True)['sd'].value
        hdr = getFitsHdr(fle)
        exp_time = hdr['exptime']

        r_pix = radius / 0.6 # radius of sun in pixels
        disk_mask = get_disk_mask(data, r_pix)
        disk_data = np.ma.array(data, mask=disk_mask)
        med_int = np.ma.median(disk_data) # np.median doesn't support masking
        
        return med_int / exp_time


def process_wave(wave):
        """Gets the median intensities for a wavelength, the filtered regression,
        and all the paths.

        If no *good* data is found in first 6 hours of day at 15 minutes steps,
        then the value is replaced with NaN in the series.
        Good images are those that have a "quality" rating of 0

        At the end, all NaNs are filled with the last known value
        Unkown values in the beginning are filled from the next known value
        """
        paths = pd.Series(index=date_list)
        raw = pd.Series(index=date_list)
        for date in datetime_list:
                fles = fetch(date, date + timedelta(minutes=1), wave)
                missing_data = False
                while (len(fles) == 0) or (getFitsHdr(fles[0])['quality'] != 0):
                        date += timedelta(minutes=15)
                        fles = fetch(date, date + timedelta(minutes=1), wave)
                        if date.hour >= 6:
                                missing_data = True
                                break
                print(date)
                if not missing_data:
                        index = [str(date.date())]
                        fle = fles[0]
                        med_int = process_med_int(fle)
                        paths.loc[index] = fle
                        raw.loc[index] = med_int
        paths = paths.ffill() # propagate missing values forwards
        paths = paths.bfill() # backwards. (if initial dates lack data)
        raw = raw.ffill()
        raw = raw.bfill()
        return [wave, paths, raw]
        """
        except:
                with open('/Users/pauly/Desktop/{}.pkl'.format(wave), 'wb') as f:
                        pickle.dump([wave, paths, raw])
        """


def main(compute_regression=True):
        """Gets all the sun intensities for all wavelengths.
        Uses multiprocessing.
        """
        csv_dict = {}
        with mp.Pool(processes=12) as pool:
                for r in pool.imap_unordered(process_wave, wavelengths):
                        wave = r[0]
                        csv_dict[wave + '_paths'] = r[1]
                        csv_dict[wave + '_raw'] = r[2]
                        if compute_regressoin:
                                csv_dict[wave + '_filtered'] = pd.Series(
                                        savgol_filter(r[2], 301, 2), index=date_list)
        df = pd.DataFrame(csv_dict)
        return df


def update_csv():
        """Updates the csv
        """
        if os.path.exists(csv_path):
                df = open_csv()
                latest_date = parse_time(df.index[-1])
                global datetime_list
                global date_list
                datetime_list = create_date_series(latest_date
                                                          + timedelta(days=1)
                                                          + timedelta(minutes=1))
                date_list = [str(date.date()) for date in datetime_list]
                df2 = main(compute_regression=False)
                df = df.append(df2)
                for wave in wavelengths:
                        df[wave + '_filtered'] = pd.Series(
                                savgol_filter(df[wave + '_raw'], 301, 2),
                                index=df.index)
        else:
                df = main()
        return df


def get_dim_factor(date, wave):
        """Gets the intensity scale factor for a day
        Unlike before, this outputs the scale factor an image's data
        should be multiplied by. Not the actual intensity.

        date: datetime object or date string
        Datetime objects are truncated to the current day.
        wave: string denoting the wavelength (eg '171')
        """
        df = open_csv()
        if not isinstance(date, str):
                date = str(date.date())
        if date > df.index[-1]:
                # Use latest known time if date ahead of csv
                date = df.index[-1]
        if date < df.index[0]:
                # Use first known time if date before start of csv data
                date = df.index[0]
        scale_factor = (df.loc['2010-05-01', wave + '_filtered']
                        / df.loc[date, wave + '_filtered'])
        return scale_factor


def get_today_factors():
        """Outputs a dictionary containing today's dim factors
        """
        factors = {}
        date = datetime.utcnow()
        for wave in wavelengths:
                factors[wave] = get_dim_factor(date, wave)
        return factors


def make_movie_frames(downscale=(8,8), rescale_brightness=True,
                      side_by_side=True, sampling_rate=7):
        """Produces the frames for a time-series movie of AIA data

        downscale: if filed, downscales the data by (x, y) factor

        rescale: determines if brightness correction takes place
        """
        df = open_csv()
        out_path = os.path.join(os.path.dirname(__file__), 'movies/')
        if not os.path.exists(out_path):
                os.mkdir(out_path)
        for wave in wavelengths:
                out_path_wave = os.path.join(out_path, wave + '/')
                if not os.path.exists(out_path_wave):
                        os.mkdir(out_path_wave)
                for i, img_path in enumerate(df[wave + '_paths'][::sampling_rate]):
                        mov_img.process_img(img_path,
                                            '{}{:03d}.jpg'.format(out_path_wave, i),
                                            downscale=downscale,
                                            rescale_brightness=rescale_brightness,
                                            side_by_side=side_by_side)


def make_animation(fname, frame_dir):
        """Makes a gif from a set of frames

        fname: name of output file

        frame_dir: dir containing frames formatted as dir/*.png

        gif or movie format determined from fname
        """
        from subprocess import Popen
        if fname[-4:] == '.gif':
                Popen(['convert', '-loop', '1', '-delay', '2', '-dispose', 'previous',
                       frame_dir, fname])
        else:
                Popen(['~/ffmpeg', '-y', '-i', frame_dir, '-filter:v',
                       '"setpts=1.0*PTS"', '-c:v', 'libx264', '-crf', '15',
                       '-pix_fmt', 'yuv420p', fname])


def make_all_animations(gif=False):
        """Makes all the gifs

        gif: output a gif or movie if False
        """
        if gif:
                for wave in wavelengths:
                        make_animation('movies/{}.gif'.format(wave),
                                       'movies/{}/*'.format(wave))
        else:
                for wave in wavelengths:
                        make_animation('movies/{}.mov'.format(wave),
                                       '\'movies/{}/{}.jpg\''.format(wave, '%03d'))


def open_csv():
        if not os.path.exists(csv_path):
                raise FileNotFoundError('Can not find AIA data CSV file.')
        return pd.read_csv(csv_path, index_col=0)


if __name__ == '__main__':
        df = update_csv()
        df.to_csv(csv_path)
        # remove unnecessary columns for json file creation
        drop_columns = [x for x in df.columns if 'raw' in x]
        drop_columns.extend([x for x in df.columns if 'paths' in x])
        df.drop(drop_columns, inplace=True, axis=1)
        df.to_json('{}json'.format(csv_path[:-3]), orient='index')
