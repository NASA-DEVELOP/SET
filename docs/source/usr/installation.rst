================
**Installation**
================

.. image:: ../_static/adult_install.jpeg
   :scale: 7%
   :align: right

SET is still in its early stages of development and does not yet offer a double-click
installer. As a result, running the program requires you to setup the development environment,
which can be tricky for those not familiar with programming in the command line. We've tried to be as
descriptive as possible in our installation guide for our users. This step-by-step guide should be enough to setup the program correctly, but if you have any questions or trouble installing, feel free to ask for help from the :doc:`support page <trouble/contact>` and describe your issue to us.

**Windows**
-----------

1. Download `Anaconda - Python 2.7 version <https://www.continuum.io/downloads>`_

.. warning::

   When installing Anaconda, you must change Advanced Options to add Anaconda to your PATH. Disregard the "Not recommended" warning and check the "Add Anaconda to my PATH environment variable" option.

2. Download `Git <https://git-scm.com/downloads>`_

3. Open Anaconda Prompt and clone the SET repository by typing::
	
	``git clone <SET Repository Link> skyglow``

.. figure:: ../_static/install_git_clone.PNG
   :figwidth: 750

   *Example of cloning the repository from Bitbucket, a we-based hsoting service for Git similar to Github. SET's repository link can be found on its Github page, under the green "Clone or download" button. Your link will be unique to your account.*

4. Switch into skyglow and install the geospatial environment, which contains all the necessary packages::

	``cd skyglow``
	``conda env create -f geospatial.yml``

5. Activate the environment (must be done every time you open the prompt)::

	``activate geospatial``

6. Run the program::

	``python Skyglow.py``

**Required Packages**
----------------------

The following dependencies are needed for the program to run properly. If creating the environment from the ``geospatial.yml`` file does not work, you may need to isntall the following packages manually. For help with installation, check the FAQ or contact the developers on the :doc:`support page <trouble/contact>`. 

* `anaconda`_
* `numpy`_ via ``conda install numpy``
* `GDAL`_ via ``conda install gdal``
* `matplotlib`_ via ``conda install matplotlib``
* matplotlib_scalebar via ``pip install matplotlib_scalebar``
* `Pillow`_ via ``conda install pillow``
* playsound via ``pip install playsound``


.. _anaconda: http://continuum.io/downloads
.. _numpy: http://www.numpy.org/
.. _gdal: http://www.gdal.org/
.. _matplotlib: https://matplotlib.org/
.. _Pillow: https://pypi.python.org/pypi/Pillow/2.7.0
