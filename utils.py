import netCDF4
import numpy as np
from numpy.fft import fft, ifft
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.dates as mpldts
from scipy.stats import binned_statistic
import datetime# as dt
import urllib
import time
import calendar
import csv

# Find nearest value in numpy array
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# Convert to unix timestamp
def get_unix_timestamp(human_time,dateFormat):
    unix_timestamp = int(calendar.timegm(datetime.datetime.strptime(human_time, dateFormat).timetuple()))
    return unix_timestamp

# Convert to human readable timestamp
def get_human_timestamp(unix_timestamp, dateFormat):
    human_timestamp = datetime.datetime.utcfromtimestamp(int(unix_timestamp)).strftime(dateFormat)
    return human_timestamp


def compute_1d_spectrum(eta,dx):
    dk = 1/len(eta)
    etak = fft(eta)
    Nfft = len(eta)
    k = np.fft.fftfreq(Nfft)/(dx)
    # compute spectrum
    PSD=np.abs(etak)**2/Nfft**2/dk # power spectral density ( m^2 / cpm )
    # keep only positive frequencies
    Npos=int(Nfft/2)
    kp = k[0:Npos]
    PSD=2*PSD[0:Npos] # times two to account for the variance of the ambigous negative frequencies
    return PSD,kp  # smoothed spectrum PSDs at center wavenumber ksc

def spectrum_smooth(PSD,w,Nave):
    # Smooth the spectum by a fcator Nave (compute a bin average by a factor Nave)
    # w is the wavenumber or frquency array
    # PSD is the power spectrum
    PSDs=binned_statistic(w,PSD,bins=len(PSD)/Nave).statistic
    ws=binned_statistic(w,PSD,bins=len(PSD)/Nave).bin_edges
    wsc=(ws[0:-1:]+ws[1::])/2 # center bin
    return PSDs,wsc  # smoothed spectrum PSDs at center wavenumber wsc

def filt1d(data, window):
    """
    Smooths a 1D spectrum using a moving average filter.
    
    Parameters:
    data (array-like): The 1D spectral data to be smoothed.
    window_size (int): The size of the moving average window.
    
    Returns:
    numpy.ndarray: The smoothed 1D spectrum, with length reduced by (window_size - 1).
    
    Author: UCONN AirSea Lab, 2022
    """

    return np.convolve(data, window, mode='same')

def compute_auto_spectrum(data1,dx):
    dk = 1/len(data1)
    etak = fft(data1)
    Nfft = len(data1)
    k = np.fft.fftfreq(Nfft)/(dx)
    # compute spectrum
    PSD= (etak * np.conj(etak))/Nfft**2/dk # power spectral density ( m^2 / cpm )
    # keep only positive frequencies
    Npos=int(Nfft/2)
    kp = k[0:Npos]
    PSD=2*PSD[0:Npos] # times two to account for the variance of the ambigous negative frequencies
    return PSD,kp  # smoothed spectrum PSDs at center wavenumber ksc

def compute_cross_spectrum(data1,data2,dx):
    dk = 1/len(data1)
    data1k = fft(data1)
    data2k = fft(data2)
    Nfft = len(data1)
    k = np.fft.fftfreq(Nfft)/(dx)
    # compute spectrum
    PSD= (data1k * np.conj(data2k))/Nfft**2/dk # power spectral density ( m^2 / cpm )
    # keep only positive frequencies
    Npos=int(Nfft/2)
    kp = k[0:Npos]
    PSD=2*PSD[0:Npos] # times two to account for the variance of the ambigous negative frequencies

    return PSD,kp  # smoothed spectrum PSDs at center wavenumber ksc


def fourier_calculations(x,y,z,dt):
    Nfft = len(x)
    dt = 0.78125
    dk = 1/len(x)
    
    # Hann taper
    taper=.5*(1 - np.cos(2*np.pi*(np.arange(0,Nfft)/(Nfft-1))))
    
    ####################### apply taper x ###########################
    x_taper= x*taper
    np.where(x_taper < -4)
    ind =np.where(x_taper < -4)
    x_taper[ind] = 0.#None
    Fac_x=np.var(x)/np.var(x_taper)
    
    
    [PSDt_x,w]=compute_1d_spectrum(x_taper,dt)
    dw=np.mean(np.diff(w));
    PSDt_x=PSDt_x*Fac_x # correct variance lost due to taper
    
    Nave=10
    # smooth the spectra by coarsening the spectrum by a factor of 10 and averaging the data in each bin.
    [PSDts_x,wsc]=spectrum_smooth(PSDt_x,w,Nave)
    DOFs=2*Nave #degress of freedom
    
    
    # add addional smoothing with a running filter
    
    wfilt=np.ones(5);#/window_size;
    DOFf=sum(wfilt);
    # normalize filter
    wfilt=wfilt/sum(wfilt)
    
    PSD_x=filt1d(PSDts_x, wfilt);
    
    DOFs2=DOFs*DOFf
    
    ################################# apply taper y ###############################
    y_taper= y*taper
    np.where(y_taper < -4)
    ind =np.where(y_taper < -4)
    y_taper[ind] = 0.#None
    Fac_y=np.var(y)/np.var(y_taper)
    
    
    [PSDt_y,w]=compute_1d_spectrum(y_taper,dt)
    dw=np.mean(np.diff(w));
    PSDt_y=PSDt_y*Fac_y # correct variance lost due to taper
    
    Nave=10
    # smooth the spectra by coarsening the spectrum by a factor of 10 and averaging the data in each bin.
    [PSDts_y,wsc]=spectrum_smooth(PSDt_y,w,Nave)
    DOFs=2*Nave #degress of freedom
    
    
    # add addional smoothing with a running filter
    
    wfilt=np.ones(5);#/window_size;
    DOFf=sum(wfilt);
    # normalize filter
    wfilt=wfilt/sum(wfilt)
    
    PSD_y=filt1d(PSDts_y, wfilt);
    
    DOFs2=DOFs*DOFf
    DOFs2
    
    ####################### apply taper z #########################
    z_taper= z*taper
    np.where(z_taper < -4)
    ind =np.where(z_taper < -4)
    z_taper[ind] = 0.#None
    Fac_z=np.var(z)/np.var(z_taper)
    
    [PSDt_z,w]=compute_1d_spectrum(z_taper,dt)
    dw=np.mean(np.diff(w));
    PSDt_z=PSDt_z*Fac_z # correct variance lost due to taper
    
    Nave=10
    # smooth the spectra by coarsening the spectrum by a factor of 10 and averaging the data in each bin.
    [PSDts_z,wsc]=spectrum_smooth(PSDt_z,w,Nave)
    DOFs=2*Nave #degress of freedom
    
    
    # add addional smoothing with a running filter
    
    wfilt=np.ones(5);#/window_size;
    DOFf=sum(wfilt);
    # normalize filter
    wfilt=wfilt/sum(wfilt)
    
    PSD_z=filt1d(PSDts_z, wfilt);
    
    DOFs2=DOFs*DOFf
    DOFs2
    

    
    #finding Qyz
    [PSD_yz,kp] = compute_cross_spectrum(y,z,dt)
    [PSD_yz,wsc]=spectrum_smooth(PSD_yz,w,Nave)
    DOFs=2*Nave #degress of freedom
    
    # add addional smoothing with a running filter
    
    wfilt=np.ones(5);#/window_size;
    DOFf=sum(wfilt);
    # normalize filter
    wfilt=wfilt/sum(wfilt)
    
    PSD_yz=filt1d(PSD_yz, wfilt);
    Qyz = np.imag(PSD_yz)
    
    
    #finding Qxz
    [PSD_xz,kp] = compute_cross_spectrum(x,z,dt)
    np.shape(PSD_xz)
    
    [PSD_xz,wsc]=spectrum_smooth(PSD_xz,w,Nave)
    DOFs=2*Nave #degress of freedom

    # add addional smoothing with a running filter
    
    wfilt=np.ones(5);#/window_size;
    DOFf=sum(wfilt);
    # normalize filter
    wfilt=wfilt/sum(wfilt)
    
    PSD_xz=filt1d(PSD_xz, wfilt);
    Qxz = np.imag(PSD_xz)
    
    
    #finding Cxy
    [PSD_xy,kp] = compute_cross_spectrum(x,y,dt)
    [PSD_xy,wsc]=spectrum_smooth(PSD_xy,w,Nave)
    DOFs=2*Nave #degress of freedom
    
    # add addional smoothing with a running filter
    
    wfilt=np.ones(5);#/window_size;
    DOFf=sum(wfilt);
    # normalize filter
    wfilt=wfilt/sum(wfilt)
    
    PSD_xy=filt1d(PSD_xy, wfilt);
    Cxy = np.real(PSD_xy)
    
    
    #calculating fourier coefficients
    a1 = Qxz/((PSD_z*(PSD_x+PSD_y))**(1/2))
    b1 = Qyz/((PSD_z*(PSD_x+PSD_y))**(1/2))
    b2 = (2*Cxy)/(PSD_x+PSD_y)
    a2 = (PSD_x - PSD_y)/(PSD_x+PSD_y)

    return PSD_z,wsc,a1,b1,a2,b2


def plot_spec(PSD_z,wsc,a1,b1,a2,b2):
    theta = np.linspace(0,2*np.pi,115)
    
    D_f_theta = 1/(2*np.pi) + 1/(np.pi)*((a1[:,np.newaxis]*np.cos(1*theta) 
                                      + b1[:,np.newaxis]*np.sin(1*theta)) 
                                     + (a2[:,np.newaxis]*np.cos(2*theta) 
                                                + b2[:,np.newaxis]*np.sin(2*theta)))

    Psi = np.multiply(PSD_z,D_f_theta)


    fig = plt.figure(figsize=(6,6))
    
    #radii = Fq[0:35] # Only call frequency bands > 0.3 Hz
    #thetas = Dmean_rad[:]
    thetas = theta[:]
    
    ## Color-scale
    contour_levels = np.arange(-0.02,0.1,0.001) # Manually set colorbar min/max values
    #contour_levels = np.arange(np.min(Psi),np.max(Psi),0.0001) # Manually set colorbar min/max values
    #contour_levels = np.arange(np.max(Psi)/1000,np.max(Psi)/100,0.0001) # Manually set colorbar min/max values
    #contour_levels = np.arange(Edfloat/1000,Edfloat/2,0.0001) # Automatically set colorbar min/max values based on 'Ed'
    
    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    ylabels = ([20,10,6.7,5,4])
    ax.set_yticklabels(ylabels)
    
    #colorax = ax.contourf(thetas, radii, Ednew[0:35],contour_levels)
    #colorax = ax.contourf(thetas, wsc, Psi,contour_levels)
    colorax = ax.contourf(theta, wsc, Psi,contour_levels, extend="both")
    
    
    ## Set titles and colorbar
    #suptitle('Polar Spectrum at Stn. ' + stn, fontsize=22, y=0.95, x=0.45)
    #title(startdate, fontsize=16, y=1.11)
    
    cbar = fig.colorbar(colorax)
    cbar.set_label('Energy Density (m*m/Hz/deg)', rotation=270, fontsize=16)
    cbar.ax.get_yaxis().labelpad = 30
    
    degrange = range(0,360,30)
    #lines, labels = plt.thetagrids(degrange, labels=None, frac = 1.07)

    print(np.max(Psi), "is the max value of energy density")