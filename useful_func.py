import numpy as np
from scipy import signal
import scipy.io.wavfile as wav



def bp(y, fs=40/100, lowcut=.05, highcut=.1, order=3) :
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return signal.lfilter(b, a, y)

def get_mov_rms(y, fs=40/100, window=150):
    return np.sqrt(signal.convolve(y**2, np.ones(int(fs*window))/int(fs*window), mode='same'))

def get_mov_mean(y, fs=40/100, window=150):
    return signal.convolve(y, np.ones(int(fs*window))/int(fs*window), mode='same')

def listen(y, sample_rate=50000, filename='temp.wav'):
    y = y / np.max(np.abs(y))
    y = np.int16(y * 32767)
    wav.write(filename, sample_rate, y)