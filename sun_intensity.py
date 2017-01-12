"""module for finding the brightness values for AIA images using the date

THIS IS OLD CODE. CURRENT DEVELOPMENT MOVED TO THE REPOS DIRECTORY
"""

from datetime import datetime, timedelta
from fetcher import fetch
from sunpy.time import parse_time
from astropy.io import fits
import pickle
from sunpy.map import Map
from sunpy.instr.aia import aiaprep
import numpy as np
import imp
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date
from scipy.signal import savgol_filter
import os
import pytz
pb0r = imp.load_source('pb0r',
                       os.path.join(os.path.dirname(__file__), 'pb0r.py'))
import pandas as pd
"""
import logging # next 3 lines make mp pass exceptions stderr
mpl = mp.log_to_stderr()
mpl.setLevel(logging.INFO)
"""


WAVELENGTHS = ['131','171','193','211','304','335','94']


def storePaths(wavelength):
        """stores the file location in a txt file
        """
        now = datetime.utcnow()
        start = datetime(2010,5,1,0,0,0)
        halfyear = timedelta(days=182)
        lst = []
        while (start + halfyear < now):
                lst.extend(fetch(start, start + halfyear, wavelength, td=86400))
                start += halfyear
        lst.extend(fetch(start, now, wavelength, td=86400))

        qualCheck(lst, wavelength)

        f = open(wavelength + 'file_locations.txt','a')
        for fle in lst:
                f.write(fle + '\n')
        f.close()


def getFitsHdr(fle):
        """get the header for a fits file
        """
        f = fits.open(fle)
        hdr = f[-1].header
        f.close()
        return hdr


def qualCheck(lst, wavelength):
        """checks the quality of each file and calls qualReplace if != 0
        """
        for idx, fle in enumerate(lst):
                if os.path.isfile(fle):
                        hdr = getFitsHdr(fle)
                        if (hdr['quality'] != 0):
                                date = parse_time(hdr['date-obs'])
                                qualReplace(fle, date, lst, wavelength)
                else:
                        lst.remove(fle)
                        print(fle + ' does not exist. File removed.')


def qualReplace(fle, date, lst, wavelength):
        """replaces images with quality != 0 with one taken 1 minute later
        Checks that image, too, and keeps going recursively until
        finding a good image
        """
        if date.hour < 1:
                fivemin = timedelta(minutes=5)
                poss_replace = fetch(date, date + fivemin, wavelength)[-1]
                hdr = getFitsHdr(poss_replace)
                if (hdr['quality'] == 0):
                        idx = lst.index(fle)
                        lst.remove(fle)
                        lst.insert(idx, poss_replace)
                else:
                        qualReplace(fle, date + fivemin, lst, wavelength)
        else:
                lst.remove(fle)


def process(fle):
        """processes 1 images and extracts the following features
        data format: [DATETIME, SUN_RADIUS, HEADER, QUALITY, EXP_TIME,
        MEAN_INT, MED_INT, 10_INT, 90_INT, N_HI_LAT_INT, S_HI_LAT_INT,
        LO_LAT_INT, FILENAME]
        """
        amap = Map(fle)
        amap = aiaprep(amap)
        data = amap.data

        date = amap.date
        radius = pb0r.pb0r(date, arcsec=True)['sd'].value
        hdr = getFitsHdr(fle)
        qual = hdr['quality']
        exp_time = hdr['exptime']

        r_pix = radius / 0.6 # radius of sun in pixels
        disk_mask = get_disk_mask(data, r_pix)
        disk_data = np.ma.array(data, mask=disk_mask)
        mean_int = np.mean(disk_data)
        med_int = np.ma.median(disk_data) # np.median doesn't support masking
        per_10_int = np.percentile(disk_data, 10)
        per_90_int = np.percentile(disk_data, 90)
        
        n_hi_lat_int, s_hi_lat_int = hi_lat_int(data, r_pix, disk_mask)
        lo_lat_int = lo_lat_int_func(data, r_pix, disk_mask)

        return [date, radius, hdr, qual, exp_time, mean_int,
                        med_int, per_10_int, per_90_int, n_hi_lat_int,
                        s_hi_lat_int, lo_lat_int, fle]


def get_disk_mask(data, r_pix):
        """returns the array mask for only the disk of the Sun
        """
        x, y = np.meshgrid(*map(np.arange, data.shape), indexing='ij')
        return (np.sqrt((x - 2047.5)**2 + (y - 2047.5)**2) > r_pix)


def hi_lat_int(data, r_pix, disk_mask):
        """returns the intensity of the Sun at the North and South pole at
        latitudes greater than 70 degrees
        """
        x, y = np.meshgrid(*map(np.arange, data.shape), indexing='ij')
        # mask for the southern pole first
        mask = np.logical_or((x < 2047.5 + r_pix * np.sin(70 * np.pi / 180)),
                                                 disk_mask)
        masked_data = np.ma.array(data, mask=mask)
        north_int = np.mean(masked_data)
 
        # rotate mask 180 degrees to apply it on the south pole
        mask = np.rot90(mask, k=2)
        masked_data = np.ma.array(data, mask=mask)
        south_int = np.mean(masked_data)

        return north_int, south_int


def lo_lat_int_func(data, r_pix, disk_mask):
        x, y = np.meshgrid(*map(np.arange, data.shape), indexing='ij')
        # number of pixels between the center of Sun and 10 degrees latitude 
        lat_pix_10 = r_pix * np.sin(10 * np.pi / 180)
        # mask pixels that are too low, too high, and too far from the center of the sun
        mask = np.logical_or((x < 2047.5 - lat_pix_10), (x > 2047.5 + lat_pix_10), disk_mask)
        masked_data = np.ma.array(data, mask=mask)
        return np.mean(masked_data)


def process_wrapper(fle):
        """processes one file and prints out the extracted features
        data format: [DATETIME, SUN_RADIUS, HEADER, QUALITY, EXP_TIME, MEAN_INT, MED_INT, 10_INT, 90_INT, N_HI_LAT_INT, S_HI_LAT_INT, LO_LAT_INT, FILENAME]
        """
        lst = process(fle)
        print('date: ' + lst[0].strftime('%Y-%m-%dT%H:%M:%S'))
        print('radius: ' + str(lst[1]) + ' arcseconds')
        print('quality: ' + str(lst[3]))
        print('exposure time: ' + str(lst[4]))
        print('mean intensity: ' + str(lst[5]))
        print('median intensity: ' + str(lst[6]))
        print('10th percentile intensity: ' + str(lst[7]))
        print('90th percentile intensity: ' + str(lst[8]))
        print('north high lat intensity: ' + str(lst[9]))
        print('south high lat intensity: ' + str(lst[10]))
        print('low lat intensity: ' + str(lst[11]))
        print('filename: ' + str(lst[12]))


def analysis(fle_name):
        data = open_norm(fle_name)
        wavelength = fle_name[:-8]
        
        dates = date2num(data[:,0])
        radius = data[:,1]
        exp_time = data[:,4]
        mean_int = data[:,5]
        med_int = data[:,6]
        per_10_int = data[:,7]
        per_90_int= data[:,8]
        n_hi_lat_int = data[:,9]
        s_hi_lat_int = data[:,10]
        lo_lat_int = data[:,11]
        plt.subplot(2,1,1)
        plt.plot_date(dates, mean_int, 'k-', xdate=True, ydate=False, lw=2)
        plt.plot_date(dates, n_hi_lat_int, 'g-', xdate=True, ydate=False, lw=2)
        plt.plot_date(dates, s_hi_lat_int, 'r-', xdate=True, ydate=False, lw=2)
        plt.plot_date(dates, lo_lat_int, 'y-', xdate=True, ydate=False, lw=2)
        plt.legend(['mean intensity', 'north pole intensity',
                                'south pole intensity', 'low latitude intensity'])
        plt.title('Intensity in particular regions of the Sun')

        plt.subplot(2,1,2)
        # plt.plot_date(dates, radius, 'k-', xdate=True, ydate=False, lw=2)
        # plt.plot_date(dates, exp_time, 'k-', xdate=True, ydate=False, lw=2)
        
        plt.plot_date(dates, med_int, 'k-', xdate=True, ydate=False)
        sav = savgol_filter(med_int, 301, 2)
        plt.plot_date(dates, sav, 'r-', xdate=True, ydate=False, lw=2)
        
        dumpreg(wavelength, [dates, sav])

        
def patchProcessGaps(fle_name, wavelength):
        """patches gaps the the data left by a crash with the multiprocessing
        not always necessary
        """
        have_lst = []
        with open(os.path.join(os.path.dirname(__file__), fle_name),'rb') as f:
                data_lst = pickle.load(f)
                for pt in data_lst:
                        have_lst.append(pt[-1])
        to_get_lst = [] 
        with open(wavelength + 'file_locations.txt','r') as f:
                for fle in f.readlines():
                        fle = fle[0:-1] # remove newline character
                        if fle not in have_lst:
                                to_get_lst.append(fle)
        print('Number of data points missing:',len(to_get_lst))

        # save along the way
        with mp.Pool(processes=12) as pool:
                for idx, r in enumerate(pool.imap_unordered(process, to_get_lst)):
                        data_lst.append(r)
                        if (idx % 80 == 0):
                                with open(os.path.join(os.path.dirname(__file__),
                                                                           str(idx) + '.pkl'), 'wb') as pkl_f:
                                        pickle.dump(data_lst, pkl_f)

        with open(os.path.join(os.path.dirname(__file__),
                                                   wavelength + 'data.pkl'),'wb') as f:
                pickle.dump(data_lst, f)


def replace_pts(wavelength):
        """takes in a wavelength and replaces the points in the data file with
        bad attributes with ones with better exposure times or qualities
        """
        with open(os.path.join(os.path.dirname(__file__),
                                                   wavelength + 'data.pkl'), 'rb') as f:
                data = pickle.load(f)
        onemin = timedelta(minutes=1)

        for pt in data:
                hdr = pt[2]
                if (hdr['quality'] != 0):
                        datatime = pt[0]
                        poss_replace = fetch(datatime, datatime + onemin,
                                                                 wavelength)[-1]
                        hdr = getFitsHdr(poss_replace)
                        while (hdr['quality'] != 0):
                                datatime += onemin
                                poss_replace = fetch(datatime, datatime + onemin,
                                                                         wavelength)[-1]
                                hdr = getFitsHdr(poss_replace)
                        data.remove(pt)
                        data.append(process(poss_replace))
        with open(os.path.join(os.path.dirname(__file__),
                                                   wavelength + 'data.pkl'),'wb') as f:
                pickle.dump(data, f)


def open_norm(fle_name):
        """opens and normalizes for the exposure time given the data in a
        given pickle file
        used when dealing with the raw, extracted features, not the regression
        """
        with open(os.path.join(os.path.dirname(__file__), fle_name), 'rb') as f:
                data = pickle.load(f)
        data = np.array(sorted(data, key= lambda x: x[0]))
        # normalize for exposure time
        data[:,5:-1] = data[:,5:-1] / np.reshape(data[:,4],(len(data),1))
        return data


def get_intensity(date, wavelength):
        """takes in a datetime instance and gets the intensity of the Sun
        based on its apparent brightness on that day
        """
        date = date2num(date)
        reg = openreg(wavelength)
        idx = find_nearest(reg[0],date)
        return reg[1][idx]


def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx


def process_multi(wavelength):
        datalst = []
        flelst2 = []
        with open(wavelength + 'file_locations.txt','r') as f:
                for fle in f.readlines():
                        flelst2.append(fle[0:-1]) # remove newline character

        with mp.Pool(processes=12) as pool:
                for idx, r in enumerate(pool.imap_unordered(process, flelst2)):
                        datalst.append(r)
                        if (idx % 79 == 0):
                                with open(os.path.join(os.path.dirname(__file__),
                                                       str(idx) + '.pkl'),'wb') as pkl_f:
                                        pickle.dump(datalst, pkl_f)

        with open(os.path.join(os.path.dirname(__file__),
                               wavelength + 'data.pkl'),'wb') as pkl_f:
                pickle.dump(datalst,pkl_f)

        replace_pts(wavelength)


def regressionComparison():
        lines = ['b-','g-','r-','c-','m-','y-','k-']
        for idx, wave in enumerate(WAVELENGTHS):
                reg = openreg(wave)
                plt.plot_date(reg[0], reg[1], lines[idx], xdate=True, ydate=False,
                                          lw=2)
        plt.legend(['131','171','193','211','304','335','94'])
        plt.title('Different Wavelength Regression Comparisons')


def openreg(wave):
        with open(os.path.join(os.path.dirname(__file__),
                                                   wave + 'regression.pkl'), 'rb') as f:
                return pickle.load(f)


def dumpreg(wave, data):
        with open(os.path.join(os.path.dirname(__file__),
                                                   wave + 'regression.pkl'), 'wb') as f:
                pickle.dump(data, f)


def main(wavelength):
        storePaths(wavelength)
        process_multi(wavelength)


def update():
        """Looks at existing regression to see what the latest image dates are,
        Fetches new files locations for more recent images,
        And then processes them and updates the regression.
        Regression intensity normalized for exposure time (using division)
        """
        utc = pytz.utc
        now = utc.localize(datetime.utcnow())
        for wave in WAVELENGTHS:
                dates = openreg(wave)[0]
                latest_date = num2date(dates[-1])
                new = fetch(latest_date + timedelta(days=1), now, wave, td=86400)
                newdata = []
                with mp.Pool(processes=12) as pool:
                        for idx, r in enumerate(pool.imap_unordered(process, new)):
                                newdata.append(r)
                newdata = np.array(newdata)
                data = open_norm(wave + 'data.pkl')
                dates = data[:,0]
                newdates = newdata[:,0]
                dates = date2num(np.append(dates, newdates))

                med_int = data[:,6]
                newmed_int = newdata[:,6] / newdata[:,4]
                med_int = np.append(med_int, newmed_int)
                sav = savgol_filter(med_int, 301, 2)
                dumpreg(wave, [dates, sav])
        print('Regression updated. {}'.format(now))


def get_dim_factors(date):
        dim_factors = []
        for wave in WAVELENGTHS:
                mission_start_brightness = get_intensity(datetime(2010,5,1), wave)
                today_brightness = get_intensity(date, wave)
                dim_factor = mission_start_brightness / today_brightness
                dim_factors.append(dim_factor)
        return dict(zip(WAVELENGHTS, dim_factors))


def roundTime(dt=None, dateDelta=timedelta(minutes=1)):
        """Round a datetime object to a multiple of a timedelta
        dt : datetime.datetime object, default now.
        dateDelta : timedelta object, we round to a multiple of this, default 1 minute.
        Author: Thierry Husson 2012 - Use it as you want but don't blame me.
        Stijn Nevens 2014 - Changed to use only datetime objects as variables
        """
        roundTo = dateDelta.total_seconds()

        if dt == None : dt = datetime.now()
        seconds = (dt - dt.min).seconds
        # // is a floor division, not a comment on following line:
        rounding = (seconds+roundTo/2) // roundTo * roundTo
        return dt + timedelta(0,rounding-seconds,-dt.microsecond)


def make_csv(fname):
        d = {}
        for wave in WAVELENGTHS:
                data = open_norm(wave + 'data.pkl')
                dates = data[:,0]
                indices = [str(roundTime(dt=date, dateDelta=timedelta(days=1)).date()) for date in dates]
                med_int = np.array([num.item() for num in data[:,6]])
                sav = savgol_filter(med_int, 301, 2)
                d[wave + 'raw'] = pd.Series(med_int, index=indices)
                d[wave + 'filtered'] = pd.Series(sav, index=indices)
        df = pd.DataFrame(d)
        df.to_csv(path_or_buf=fname)

# -----------start of new csv-based scaling tracking--------------
CSV_PATH = os.path.join(os.path.dirname(__file__), 'aia_rescaling_data.csv')


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
DATETIME_LIST = create_date_series((2010,5,1, 0, 1))
DATE_LIST = [str(date.date()) for date in DATETIME_LIST]


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
        paths = pd.Series(index=DATE_LIST)
        raw = pd.Series(index=DATE_LIST)
        for date in DATETIME_LIST:
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
        sav = pd.Series(savgol_filter(raw, 301, 2), index=DATE_LIST)
        return [wave, paths, raw, sav]
        """
        except:
                with open('/Users/pauly/Desktop/{}.pkl'.format(wave), 'wb') as f:
                        pickle.dump([wave, paths, raw])
        """


def updated_main():
        """Gets all the sun intensities for all wavelengths.
        Uses multiprocessing.
        """
        csv_dict = {}
        with mp.Pool(processes=12) as pool:
                for r in pool.imap_unordered(process_wave, WAVELENGTHS):
                        wave = r[0]
                        csv_dict[wave + '_paths'] = r[1]
                        csv_dict[wave + '_raw'] = r[2]
                        csv_dict[wave + '_filtered'] = r[3]
        df = pd.DataFrame(csv_dict)
        return df


def update_csv():
        """Updates the csv
        """
        if os.path.exists(CSV_PATH):
                df = pd.read_csv(CSV_PATH)
                latest_date = parse_time(df.index[-1])
                DATETIME_LIST = create_date_series(latest_date + timedelta(days=1) + timedelta(minutes=1))
                DATE_LIST = [str(date.date()) for date in DATETIME_LIST]
                df2 = updated_main(fname=None)
                df = df.append(df2)
        else:
                df = updated_main()
        return df


def get_csv_intensity(date, wave):
        """Gets the intensity scale factor for a day
        Unlike before, this outputs the scale factor an image's data
        should be multiplied by. Not the actual intensity.

        date: datetime object or date string
        Datetime objects are truncated to the current day.
        wave: string denoting the wavelength (eg '171')
        """
        df = pd.read_csv(CSV_PATH)
        if not isinstance(date, str)
                date = str(date.date())
        scale_factor = df.loc['2010-05-01', wave] / df.loc[date, wave]
        return scale_factor


if __name__ == '__main__':
        df = update_csv()
        df.to_csv(CSV_PATH)
