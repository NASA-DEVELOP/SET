# Code breakdown

The toolbox is packaged in the `skyglow/` folder. The main entry into the CLI is `skyglow/__main__.py`. It contains the CLI, built with `argsparse`. Depending on what action is requested _(gui, sgmap_single, kernel_lib, sgmap_multiple,_ or _hemisphere)_ it calls a different function:

`gui`: `skyglow.main()`. Calls the entry point to the GUI.

`sgmap_single`: `darksky.sgmapper`. Creates a single artificial brightness map for a given input VIIRS image, regional latitude, regional atmospheric clarity, and zenith and azimuth viewing angles.

`kernel_lib`: `darksky.krn_unpack` (asynchronous) or `darksky.krnlibber` (synchronous). Generates a series of light propagation kernels for the angle combinations listed in a given `.csv` file.

`sgmap_multiple`: `darksky.multisgmapper`. Creates an artificial brightness map for each kernel in a given folder (the bulk version of `sgmap_single`).

`hemisphere`: `darksky.generate_hem`. Creates a hemispherical visualization using values extracted from artificial brightness maps in a given folder.

---

`darksky.sgmapper` and `darsky.krnlibber` perform the same first step. They generate a light propagation kernel from the given parameters:
- Regional latitude: the approximate latitude of the study area (used to calculate kernel size and Earth radius at that latitude)
- Atmospheric clarity ratio: the ratio of aerosols to small particles in the air in the study area (affects light scattering; look up Rayleigh and Mie scattering for more information). >1 means air polluted with aerosols, <1 means clearer air.
- Zenith angle: vertical viewing angle from zenith (opposite of altitude)
- Azimuth angle: viewing angle from North (0 degrees)
- Input VIIRS image: VIIRS DNB data of study area with a __400km buffer__ for proper kernel calculation

The parameters are given to `darksky.fsum_2d`, an internal function that first calculates the distance of each pixel in the kernel to the central pixel of the kernel (stored as one matrix), the angle between each pixel and the central pixel (beta angle), etc. `fsum_2d` then repeatedly calls `darsky.fsum_single` to calculate total light propagation along the vector from each pixel to the central pixel (`fsum_single` is called once for every pixel in the kernel, so this step takes a long time). The result is saved as a faux `geotiff` with geographic projection and reference the same as the input VIIRS image.

The purpose of the generated light propagation kernel is to be a summary of light propagation for the given viewing angle (import the kernel to ArcMap and play with the Symbology tab in the layer settings to reveal the kernel). It can then be applied to __every__ pixel on the input VIIRS image by making that pixel the center pixel (observation site) and calculating total light pollution at that pixel for the kernel's viewing angle by taking the weighted average of the pixels under the kernel (light sources) (more important pixels for the angle get higher values).

![convolution illustration](https://i.imgur.com/DETberU.png)

`kernel_lib` can generate kernels in bulk for a list of angles given in a `.csv` file. The results are saved to the directory the command was executed from. Then, `sgmap_multiple` can be used to create artificial skyglow maps for all the kernels. The resulting maps can be used to generate a hemispherical visualization of skyglow at any pixel (geographic point) of the study area by specifying the point in latitude/longitude and giving the folder path where all of the skyglow maps are contained. The pixel values from each skyglow map at the geographic point given are interpolated between to create the hemisphere.

![hemisphere compilation](https://i.imgur.com/6Z2qrKt.png)
