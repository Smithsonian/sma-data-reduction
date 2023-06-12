
This is an in-progress CASA-based calibration pipeline that uses a config
and flagging file format to be able to use a consistent calibration script
across different tracks.

To use this pipeline:

1. Append the `casa_reduction_pipeline` to your `.casa/config.py` CASA file: `sys.path.append("/PATH/TO/casa_reduction_pipeline/")`
2. Use a CASA version that has astropy installed (the pipeline versions seem to; e.g., `casa-6.4.1-12-pipeline-2022.2.0.64`)


There are two flavors of pipeline scripts:

1. `casa_sma_reduction_solsyssetjy.py` will use the Butler-JPL-Horizons 2012 Solar System models for the flux calibration.
2. `casa_sma_reduction_manualsetjy.py` requires the flux to be manually set for one of the calibrators (see `flux_stokesI` in `example_cfg_files/221023_05:03:22_bin32_bllacflux.cfg`).


.. note:: Both of these scripts currently work with dual tuning only. See the SMA CASA tutorial in this repo for an example using split tuning.

Each of these scripts has an additional version that creates extra quality assurance
plots in an HTML weblog (similar to the ALMA pipeline). These are the scripts ending
with `_with_quicklook_sma.py` and require the `quicklook-sma <https://github.com/e-koch/quicklook-sma>`_
to be installed for:

* accessible to the CASA environment by adding `sys.path.append("/PATH/TO/quicklook-sma/")` to the `.casa/config.py` file;
* a normal python environment that creates the HTML weblog and plotly figures.


.. warning:: The weblog remains experimental. Please check results closely.

