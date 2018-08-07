# Development tips

If you want to work on SET's source code we suggest `pip install -e skyglow/`, so the package tracks code changes.
Use version control as the project is on Github already.
Install a PEP8 linter like `pip install pep8ify`. This is a  development dependency and should not be included in `geospatial.yml`, the conda environment file.
Add unit tests to `skyglow/skyglow/test.py` and test with `python -m unittest skyglow.test` from inside the outer skyglow folder.

---

Once the project goes through code release the documentation website written with Sphinx can be updated and tested (and should be).

Though a functional relationship between NPS on-site data and SET result data exists, -180/180 azimuth angles return unusually high values, which causes these values to be outliers. It is most likely the result of some logic/arithmetic error in the code that the 2018 Summer US Urban Development team inherited from previous terms. We hope you find a fix for this, though excluding -180/180 azimuths does work as a quick fix for better hemispheres.
