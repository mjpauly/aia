"""Program for getting the adjusted image scaling factors for the current date
"""

import pickle
import sun_intensity
from datetime import datetime

if __name__ == '__main__':
	waves = ['94', '131', '171', '193', '211', '304', '335']
	start = datetime(2010,5,1)
	end = datetime.utcnow()
	for wave in waves:
		print(wave + ':')
		print('median intensity: {}'.format(sun_intensity.get_intensity(start, wave) / 
											sun_intensity.get_intensity(end, wave)))
		print()
