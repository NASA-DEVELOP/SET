================
**Installation**
================

SET is still in its early stages of development, and does not yet offer a double-click
installer. As a result, running the program requires you to setup the development environment,
which can be tricky for those not well-versed in programming. We've tried to be as
descriptive as possible in our installation guide for our users. If you have any questions or trouble installing, feel free to ask for help from the :doc:`support page <trouble/contact>` and describe your issue to us.

**Required Packages**
----------------------

- bzip2=1.0.6=vc9_3
- curl=7.52.1=vc9_0
- cycler=0.10.0=py27_0
- decorator=4.0.11=py27_0
- freetype=2.5.5=vc9_2
- functools32=3.2.3.2=py27_0
- gdal=2.1.0=py27_0
- geos=3.5.0=vc9_0
- hdf4=4.2.12=vc9_0
- hdf5=1.8.15.1=vc9_4
- icu=57.1=vc9_0
- jpeg=8d=vc9_2
- kealib=1.4.6=vc9_0
- libgdal=2.1.0=vc9_0
- libnetcdf=4.3.3.1=vc9_4
- libpng=1.6.27=vc9_0
- libtiff=4.0.6=vc9_2
- matplotlib=2.0.2=np112py27_0
- mkl=2017.0.1=0
- networkx=1.11=py27_0
- numpy=1.12.1=py27_0
- olefile=0.44=py27_0
- openssl=1.0.2k=vc9_0
- pillow=3.4.2=py27_0
- pip=9.0.1=py27_1
- proj4=4.9.2=vc9_0
- pyparsing=2.1.4=py27_0
- pyqt=5.6.0=py27_2
- python=2.7.13=1
- python-dateutil=2.6.0=py27_0
- pytz=2017.2=py27_0
- pywavelets=0.5.2=np112py27_0
- qt=5.6.2=vc9_0
- scikit-image=0.13.0=np112py27_0
- scipy=0.19.0=np112py27_0
- setuptools=27.2.0=py27_1
- sip=4.18=py27_0
- six=1.10.0=py27_0
- subprocess32=3.2.7=py27_0
- tk=8.5.18=vc9_0
- vs2008_runtime=9.00.30729.5054=0
- wheel=0.29.0=py27_0
- xerces-c=3.1.4=vc9_0
- zlib=1.2.8=vc9_3
- pip:
  - archook==1.1.0
prefix: C:\Users\rbavery\AppData\Local\Continuum\Anaconda3\envs\geospatial

**Windows**
-----------

1. Download `Python Anaconda 2.7 <https://www.continuum.io/downloads>`_
2. Download `Git Bash <https://git-scm.com/downloads>`_
3. Open a terminal, navigate to the directory you would like to contain the files for the Toolbox
   ``git clone <insert link to SkyglowGithubhere>``

4. Run (need to double check that this is right)
    ``conda create env geospatial.yml``

5. Run
    ``python Itest.py``
