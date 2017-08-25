=====================
**Updating The Site**
=====================

To build the docs, you'll need a couple of packages:

   * sphinx (installed with anaconda)
   * `sphinx-bootstrap-theme`_

.. note::

   This site is generated with ``sphinx`` . **do not** attempt to modify the html files of this
   site, instead follow these instructions for writing docs in reStructured Markdown and then building the html files using Sphinx! Sphinx makes it easy to write and create documentation websites.

The easiest way to built the site would be keeping two repositories of *SET* active, and keep
one for developing, the other for docs. For this example the development is done on the *master*
branch, and the website will always be on the *gh-pages* branch.

.. rubric:: Cloning the Repositories

1. First, go ahead and clone two repositories into the workplace of your choice.::
   
   > mkdir SET
   > cd SET
   > git clone https://Github.com/NASA-DEVELOP/SET .
   > cd ..
   > mkdir docs
   > cd docs
   > git clone https://Github.com/NASA-DEVELOP/SET .

2. Now, we have two folders, your development folder ``SET`` and your documentation site
   ``docs``. You'll need to switch the docs repository to *gh-pages*.::
   
   > cd docs
   > git checkout gh-pages

3. So you can now go on and write code, develop and create additional docs in your
   ``SET`` workspace. When you're finally ready to rebuild the site , this will require you
   to clear your ``docs`` branch and populate it with your new generated docs.

   .. warning::

      It is **very** important you specify a directory when using the ``rm`` command. *trust me*

   .. code-block:: none

      > rm -rf docs/* 
      > cd SET

4. Simply build the sphinx docs now with your build directory as docs!::

   > sphinx-build -b html docs/source docs/build

5. Head inside your docs folder, force commit the changes and the site should be live. All done.

   > cd docs
   > git push --force

.. _sphinx-bootstrap-theme: https://ryan-roemer.github.io/sphinx-bootstrap-theme/