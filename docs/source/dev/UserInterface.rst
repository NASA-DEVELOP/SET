==================
**User Interface**
==================

SET's user interface is created nearly entirely using Tkinter, which is Python's de facto standard GUI. The menu is created by initializing a class `SkyglowEstimationToolbox`. Following the initialization are plenty of functions that make use of Tkinter's library of features, such as a menu button, canvas, etc. More information, guides, and tutorials on Tkinter can be found online.

The `generate_map` function runs `Itest.py` and makes use of the `threading` library to run the UI window and Itest program simultaneously (preventing one or the other from freezing). The `LogRedirector` and `StderrRedirector` simply redirect logging lines from the terminal to the Tkinter progress window.
