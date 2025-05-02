# Ocean Waves Spectral Analysis


<div align="center">
    <img src="figs/energy_density_polar_plot.png" width="600">
</div>



The `data/` folder contains various datasets used for analysis:

- `.nc` files downloaded from the [CDIP website](https://cdip.ucsd.edu/)
- `.csv` files from THREDDS
- `.npz` files storing x, y, z values retrieved from THREDDS


## Scripts

The `scripts/` folder includes the main, completed scripts for wave spectral analysis.

### `xyz_spectrum.py`
- Loads x, y, z data from a `.npz` file and generates 1D plots.
- Requires output from `fetch_cdip_directional_example.py`.

### `1d_comparison.py`
- Compares custom spectral calculations against CDIP’s precomputed data.
- Evaluates the accuracy of custom methods.

### `spec_function_final.py`
- Uses functions from `utils.py` to compute Fourier coefficients and energy densities.
- Generates a polar plot of energy density for comparison with CDIP’s visual output.

The `example_scripts/` folder contains work-in-progress or secondary scripts.
### `fetch_cdip_directional_example.py`
- Retrieves x, y, z variables from THREDDS and saves them in `.npz` or optionally `.csv` format.
- Intended to be used with `xyz_spectrum.py`.

### `realtime_test.py`
- A simple script for testing access to the CDIP real-time dataset.

