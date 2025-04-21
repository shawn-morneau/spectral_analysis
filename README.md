# spectral_analysis
///////////////////////////////////
/////Spectral Analysis Project/////
///////////////////////////////////

This polar plot of the energy density spectrum is an example of the result of Spec_function_final

<div align="center">
    <img src="figs/energy_density_polar_plot.png" width="600">
</div>

This plot shows the difference between the 1d spectrum of our calculations vs the pre-calculated cdip data

<div align="center">
    <img src="figs/Ed_calc_vs_cdip.png" width="600">
</div>

//Data//
    The data folder is used for housing downloaded data.
    The data includes .nc files from the CDIP website, .csv files holding values from THREDDS, and .npz files holding values from THREDDS.
    Data in this folder is not yet added due to large file sizes.

//figs//
    The figs folder holds two example figures created from the Spec_function_final and 1d_comparison scripts.
    These two figures will be saved from the scripts into the figs folder.

//Example_scripts//
    The example scripts folder houses many secondary scripts that are works in progress or are not particularly useful.


    fetch_cdip_directional_example
        This script takes the x,y,z variables from THREDDS and saves them in a .npz file (or optionally a .csv).
        This file will be called for script xyz_spectrum.
        * Saving data is not useful if you are pulling for data straight from THREDDS and using it in the same script, skip the middle-man

    Realtime_test
        This script is used for testing of the CDIP realtime dataset


//Scripts//
    The scripts folder houses the main, completed scripts.

    xyz_spectrum
        This script takes data from saved values and makes 1d plots of x,y,z.
        * The script takes saved x,y,z values from a .npz file 
        * This script is used with fetch_cdip_directional_example, as this is how you will obtain your variables

    1d_comparison  
        This script goes into further detail in the comparison between our calculations and their precalculated cdip data.
        Comparisons are made to see how close our calculations are to their calculations.

    Spec_function_final
        Formulas used in this script are pulled from utils.py        
        Calculations are made for the Fourier coefficients and plotted.
        The energy density we find is made into a polar plot. and is compared against the cdip plot.
        * Note: Currently the plot of energy density does not line up with CDIP website, WIP
