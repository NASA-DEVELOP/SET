# Installation instructions

__Dependencies__:
- Anaconda/Miniconda  
- Python >= 2.7

---
1. `git clone https://github.com/NASA-DEVELOP/SET.git`
2. `cd set`
3. `conda env create -f geospatial.yml`
4. `source activate geospatial` or `activate geospatial`

You now have a `conda` environment called `geospatial`, into which all the packages (`numpy`, `scipy`, `gdal`, etc. described in `geospatial.yml`) are installed, so if you ever open a new Terminal prompt, make sure to activate the environment again before using SET.

You can now also use SET through its command-line interface `skyglow`.

`skyglow gui` opens up the GUI interface.
`skyglow --help` lists all the options and arguments for the CLI.

Ex: `skyglow sgmap_single ~/Downloads/Imagery/gi.tif -l 30.39 -kam 1.2 -z 80 -a 180 --debug` runs the sgmap_single action with a latitude of 30.39, atmospheric clarity parameter 1.2, zenith 80, azimuth 180, and sets the log level to debug.
