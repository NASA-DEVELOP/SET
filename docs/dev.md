# Development tips

If you want to work on SET's source code we suggest `pip install -e skyglow/`, so the package tracks code changes.
Use version control as the project is on Github already.
Install a PEP8 linter like `pip install pep8ify`. This is a  development dependency and should not be included in `geospatial.yml`, the conda environment file.
Add unit tests to `skyglow/skyglow/test.py` and test with `python -m unittest skyglow.test` from inside the outer skyglow folder.
To do validation with NPS data, `skyglow/skyglow/validation.py` contains all the functions necessary to convert from SET units (microcandelas/m^2) to NPS units (magnitudes/arc second^2) and do a validation on the datasets. Make sure you've generated skyglow maps for the angles you want to compare. __!!!__ NPS data is in a strange format, so import first into ArcMap and Project Raster (to something like Authalic Sphere). Save that raster as a `.tif` and that's the input for validation.

---

Once the project goes through code release the documentation website written with Sphinx can be updated and tested (and should be).

### Outstanding issues
Though a definite correlation between NPS on-site data and SET result data exists, -180/180 azimuth angles return unusually high values, which causes these values to be outliers. It is most likely the result of some logic/arithmetic error in the code that the 2018 Summer US Urban Development could not debug. We hope you find a fix for this, though excluding -180/180 azimuths does work as a quick fix for better hemispheres.

There are issues on Github:
- The light emission model used in SET is a __weighted__ average of three models. Currently, the weights are hardcoded, which may cause the model to favor urban environments over rural ones (where the light emitted follows a different pattern).

There are also a few issues outlined in `darksky.py` at the top, which have either been resolved or are not as crucial.
