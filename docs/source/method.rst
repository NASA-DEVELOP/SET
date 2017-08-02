===============
**Methodology**
===============

This page gives insight into the way the computer model applies the physic and geometry of light propagation. Although not necessary to run the program, understanding the methodology offers benefits for anyone looking to improve SET or use it for academic purposes.

**Kernel Creation**
-------------------

All of SET's methods for applying its light propagation model are in the file ``Itest.py``, which executes the ``main`` function. In this function, SET first checks if a kernel is already created and, if no kernel is supplied, creates a kernel. A kernel is essentially a matrix with numerical coefficients representing the weight of light scattering from the source to the line of sight of the observer. Each pixel in the kernel thus contains a unitless value that comapres the amount of light emitted from the ground to the amount of light scattered.

The kernel creation process is executed in ``main`` in lines 30-38::

    if os.path.isfile(kerneltiffpath) is False:
        # Estimate the 2d propagation function
        # bottom bottom_lat = 40.8797
        # top lat= 46.755666
        propkernel, totaltime = fsum_2d(regionlat_arg, k_am_arg, zen_arg, azimuth_arg)
        logger.debug('propagation array: %s', propkernel)
        kerneltiffpath = 'kernel_' + str(regionlat_arg) + '_' +  str(k_am_arg) + '_' + str(zen_arg) + '_' + str(azimuth_arg)
        array_to_geotiff(propkernel, kerneltiffpath, filein)
        logger.info("time for prop function ubreak 10: %s", totaltime)

The block calls ``fsum_2d``, which goes from line 109-230 and produces the array that serves as the light propagation kernel. ``fsum_2d`` begins by gathering several inputs: regional latitude, Earth radius, and the array of beta angles, which depends on the observer's azimuth angle.
