"""
Name: Skyglow Estimation Toolbox
Authors: Ryan Avery, Dr. Kenton Ross, Stanley Yu
Description: Python toolbox that generates artificial skyglow maps using data from NASA and NOAA's
             Visible Infrared Imaging Radiometer Suite (VIIRS) satellite sensor.
References:
(1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh, B.A. Portnov, N.A.
	  Rybnikov and R. Furgoni, 2016. The new workd atlas of artificial night sky brightness.
	  Sci. Adv. 2.
(2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial night sky
	  brightness mapped from DMSP satellite Operational Linescan System measurements. Mon.
	  Not. R. Astron. Soc. 318.
(3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. Astron. Soc.
	  Pac. 101.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# python2
try:
    from Tkinter import (Tk, Toplevel, PanedWindow, Frame, Label, Entry, Button,
                         Canvas, Scrollbar, Text, Menubutton, Menu, Checkbutton,
                         BOTH, VERTICAL, HORIZONTAL, CENTER, NE, E, W, Y, N, S,
                         StringVar, IntVar, Radiobutton, IntVar)
    import ttk
    import tkFileDialog as filedialog
    from tkFont import Font
# python3
except:
    from tkinter import (Tk, Toplevel, PanedWindow, Frame, Label, Entry, Button,
                         Canvas, Scrollbar, Text, Menubutton, Menu, Checkbutton,
                         BOTH, VERTICAL, HORIZONTAL, CENTER, NE, E, W, Y, N, S,
                         StringVar, IntVar, ttk, filedialog, Radiobutton, IntVar)
    from tkinter.font import Font
import matplotlib
matplotlib.use('TkAgg')
from PIL import Image, ImageTk, ImageOps
from numpy import loadtxt
import webbrowser
import threading
from multiprocessing import Pool
import sys, os
import time
import logging
logger = logging.getLogger()

import skyglow.constants as constants
import skyglow.darksky as darksky


# Main class of the software. Establishes GUI.
class SkyglowEstimationToolbox:

    ##################################################
    # Initialization
    def __init__(self, root):
        self.root = root

        # Radio action buttons
        self.action = None
        self.sgmap_single_btn, self.krn_lib_btn, self.multi_map_btn = None, None, None

        self.file_log_var = StringVar()
        self.csv_file_var = StringVar()
        self.krn_folder_var = StringVar()
        self.output_folder_var = StringVar()
        self.sgmap_folder_var = StringVar()

        self.krn_ent_var = StringVar()
        self.krn_var, self.hem_var = IntVar(), IntVar()
        self.img, self.cdiag = None, None
        self.lat_lbl, self.lat_entry = None, None
        self.k_lbl, self.k_entry = None, None
        self.zen_lbl, self.zen_entry = None, None
        self.azi_lbl, self.azi_entry = None, None
        self.krn_lvl, self.krn_entry, self.krn_btn = None, None, None
        self.txt_redir, self.prg_log = None, None
        self.map_btn, self.gen_krn_btn = None, None

        # Sets window title, size, and icon on screen.
        self.root.title("Skyglow Estimation Toolbox (SET)")
        self.root.geometry('%dx%d+%d+%d' % (constants.SW*0.75, constants.SH*0.75, 25, 25))
        self.root.iconbitmap(os.path.join(os.getcwd(), constants.ICO))
        self.root.resizable(False, False)
        self.root.update_idletasks()

        # Creates three paned windows for the main screen.
        base = PanedWindow()
        base.pack(fill=BOTH, expand=1)
        sub1 = PanedWindow(base, orient=VERTICAL, height=self.root.winfo_height()*3/4)
        base.add(sub1)
        sub2 = PanedWindow(sub1, orient=HORIZONTAL, height=self.root.winfo_height()/5)
        sub1.add(sub2)

        # Creates frame for holding inputs.
        self.input_frame = Frame(sub2)
        sub2.add(self.input_frame)

        # Creates frame for bottom half of main screen.
        self.img_frame = Frame(sub1, bd=2, bg='white', relief="sunken")
        sub1.add(self.img_frame)

        # Creates canvas for displaying images.
        self.img_canvas = Canvas(self.img_frame, bd=2, relief="groove",
                                 width=constants.SW*0.6, height=self.root.winfo_height()*3/4*0.9)
        self.img_canvas.place(relx=.5, rely=.5, anchor=CENTER)

        # Creates help button for link to documentation, instructions, and about.
        self.help_btn = Menubutton(self.input_frame, text="Help", relief="raised",
                                   bd=2, width=8, pady=1)
        #self.help_btn.place(relx=1, rely=0, anchor=NE)
        self.help_btn.grid(column=4,columnspan=1,row=0)
        self.help_btn_menu = Menu(self.help_btn, tearoff=0)
        doc = 'https://github.com/NASA-DEVELOP'
        self.help_btn_menu.add_command(label="Documentation", command=lambda: self.open_url(doc))
        self.help_btn_menu.add_command(label="Instructions", command=self.instructions)
        self.help_btn_menu.add_separator()
        self.help_btn_menu.add_command(label="About", command=self.about)
        self.help_btn["menu"] = self.help_btn_menu

    # Sets up input GUI and image display screen.
    def main_screen(self):
        self.action = IntVar()

        btn_width = int(constants.SW/60)
        file_width = int(constants.SW/18)
        lbl_width = int(constants.SW/60)
        gen_width = int(constants.SW/42)
        radio_font = Font(family='TkDefaultFont', size=12)
        self.sgmap_single_btn = Radiobutton(self.input_frame, text="Generate Artificial Skyglow Map", font=radio_font,
                                            width=btn_width, variable=self.action, value='sng',
                                            command=self.sng_popup)
        self.krn_lib_btn = Radiobutton(self.input_frame, text="Generate Kernel Library", font=radio_font,
                                       width=btn_width, variable=self.action, value='krn',
                                       command=self.krn_popup)
        self.multi_map_btn = Radiobutton(self.input_frame, text="Generate Maps from Multiple Kernels", font=radio_font,
                                         width=btn_width, variable=self.action, value='mul',
                                         command=self.mul_popup)
        self.hem_map_btn = Radiobutton(self.input_frame, text="Generate Hemispherical Visualization", font=radio_font,
                                         width=btn_width, variable=self.action, value='hem',
                                         command=self.hem_popup)
        #Place widget
        self.sgmap_single_btn.grid(column=0,columnspan=1, row=0)
        self.krn_lib_btn.grid(column=1, columnspan=1, row=0)
        self.multi_map_btn.grid(column=2, columnspan=1, row=0)
        self.hem_map_btn.grid(column=3, columnspan=1, row=0)

        # VIIRS Image Reference File
        self.file_lbl = Label(self.input_frame, text="Image File:", width=lbl_width, anchor=E)
        self.file_log = Entry(self.input_frame, width=file_width, bd=2, relief="sunken",
                         textvariable=self.file_log_var)
        self.browse_btn = Button(self.input_frame, text="Browse", command=self.import_viirs)

        # Angles CSV File
        self.csv_file_lbl = Label(self.input_frame, text="Angles CSV File:", width=lbl_width, anchor=E)
        self.csv_file_log = Entry(self.input_frame, width=file_width, bd=2, relief="sunken",
                         textvariable=self.csv_file_var)
        self.csv_browse_btn = Button(self.input_frame, text="Browse", command=self.import_csv)

        # Multiple Maps form Kernel library
        self.mul_file_lbl = Label(self.input_frame, text="Kernel Folder:", width=lbl_width, anchor=E)
        self.mul_file_log = Entry(self.input_frame, width=file_width, bd=2, relief="sunken",
                         textvariable=self.krn_folder_var)
        self.mul_browse_btn = Button(self.input_frame, text="Browse", command=self.import_krn_folder)

        # MultiKrn Map Output Location
        self.output_lbl = Label(self.input_frame, text="Output Location:", width=lbl_width, anchor=E)
        self.output_log = Entry(self.input_frame, width=file_width, bd=2, relief="sunken",
                 textvariable=self.output_folder_var)
        self.output_btn = Button(self.input_frame, text="Browse", command=self.import_out_folder)

        # Hemisphere Output Location
        self.sgmap_folder_lbl = Label(self.input_frame, text="Skyglow Map Location:", width=lbl_width, anchor=E)
        self.sgmap_folder_log = Entry(self.input_frame, width=file_width, bd=2, relief="sunken",
                 textvariable=self.sgmap_folder_var)
        self.sgmap_folder_btn = Button(self.input_frame, text="Browse", command=self.import_sgmap_folder)

        # Import Kernel Checkbutton
        self.check_lbl = Label(self.input_frame, text="Import Kernel:", width=lbl_width, anchor=E)

        self.krn_chk = Checkbutton(self.input_frame, anchor=W, variable=self.krn_var,
                              command=self.checkbtn_val)

        self.hem_chk_lbl = Label(self.input_frame, text="Generate kernels for hemisphere:", width=lbl_width, anchor=E)

        self.hem_chk = Checkbutton(self.input_frame, anchor=W, variable=self.hem_var)

        # Region Latitude (deg), Grand Teton National park = 43.7904 degrees N
        self.lat_lbl = Label(self.input_frame, text="Latitude (deg):", width=lbl_width, anchor=E)
        self.lat_entry = Entry(self.input_frame, width=btn_width, bd=2, relief="sunken")
        self.lon_lbl = Label(self.input_frame, text="Longitude (deg):", width=lbl_width, anchor=E)
        self.lon_entry = Entry(self.input_frame, width=btn_width, bd=2, relief="sunken")

        # Atmospheric Clarity Parameter, REF 2, Eq. 12, p. 645
        self.k_lbl = Label(self.input_frame, text="Atmospheric Clarity Parameter:",
        				   width=btn_width, anchor=E)
        self.k_entry = Entry(self.input_frame, width=btn_width, bd=2, relief="sunken")

        # Zenith angle (deg), z, REF 2, Fig. 6, p.648
        self.zen_lbl = Label(self.input_frame, text="Zenith Angle (deg):", width=lbl_width, anchor=E)
        self.zen_entry = Entry(self.input_frame, width=btn_width, bd=2, relief="sunken")

        # Azimuth angle (deg)
        self.azi_lbl = Label(self.input_frame, text="Azimuth Angle (deg):", width=lbl_width, anchor=E)
        self.azi_entry = Entry(self.input_frame, width=btn_width, bd=2, relief="sunken")

        self.krn_lbl = Label(self.input_frame, text="Kernel File:", width=lbl_width, anchor=E)
        self.krn_ent = Entry(self.input_frame, width=file_width, bd=2, relief="sunken",
                             textvariable=self.krn_ent_var)
        self.krn_btn = Button(self.input_frame, text="Browse", command=self.import_krn)

        # Generate Artificial Skyglow Map Button
        self.map_btn = Button(self.input_frame, text="Generate Artificial Skyglow Map",
                                 width=gen_width, command=self.generate_map)
        # Generate Kernal library button for SET
        self.gen_krn_btn = Button(self.input_frame, text="Generate Kernel Library",
                                    width=gen_width, command=self.generate_krn)
        # Generate Map of Multiple Kernals(word better later on)
        self.mul_map_btn = Button(self.input_frame, text="Generate Maps from Multiple Kernels",
                                    width=gen_width, command=self.generate_mmap)
        # Generate Hemispherical Visualization Display of Skyglow
        self.hem_gen_btn = Button(self.input_frame, text="Generate Hemisphere",
                                    width=gen_width, command=self.generate_hem)


    # Imports a TIFF file for referencing sky brightness in the region of interest.
    def import_viirs(self):
        # Allows user to search through his directory for VIIRS Image file.
        file_types = [('TIFF Files', '*.tif'), ('All files', '*')]
        file_name = filedialog.askopenfilename(initialdir='/', title="Select file", filetypes=file_types)
        self.file_log_var.set(file_name)

        # Checks to see if file is empty. If not, displays image on canvas.
        if file_name != '':
            pilimg = Image.open(file_name)
            pilimg_width, pilimg_height = pilimg.size
            pilimg.tile = [t for t in pilimg.tile if t[1][2] < pilimg_width and t[1][3] < pilimg_height]
            canvas_size = (self.img_canvas.winfo_width(), self.img_canvas.winfo_height())
            pilimg_r = pilimg.resize(canvas_size, Image.ANTIALIAS)
            pilimg_col = ImageOps.colorize(ImageOps.grayscale(pilimg_r), (0,0,0), (255,255,255))
            pilimg_cont = ImageOps.autocontrast(pilimg_col, cutoff=.4, ignore=None)
            self.img = ImageTk.PhotoImage(pilimg_cont)
            self.img_canvas.create_image(canvas_size[0]/2, canvas_size[1]/2, image=self.img)
        else:
            print('File is empty.')

    def import_csv(self):
        # Allows user to search through his directory for VIIRS Image file.
        file_types = [('CSV Files', '*.csv'), ('All files', '*')]
        file_name = filedialog.askopenfilename(initialdir='/', title="Select file", filetypes=file_types)
        self.csv_file_var.set(file_name)

        if file_name is '':
            print('File is empty.')

    def import_krn_folder(self):
        # Allows user to search through his directory for VIIRS Image file.
        krn_dir = filedialog.askdirectory(initialdir='/', title="Select kernel folder")
        self.krn_folder_var.set(krn_dir)

        if krn_dir is '':
            print('Directory is empty.')

    def import_out_folder(self):
        # Allows user to search through his directory for VIIRS Image file.
        output_dir = filedialog.askdirectory(initialdir='/', title="Select output folder")
        self.output_folder_var.set(output_dir)

        if output_dir is '':
            print('Directory is empty.')

    # Imports a TIFF file containing the kernel data from previous input parameters.
    def import_krn(self):
        # Allows user to search through his directory for kernel file.
        file_types = [('TIFF Files', '*.tif'), ('All files', '*')]
        file_name = filedialog.askopenfilename(initialdir='/', title="Select file", filetypes=file_types)
        self.krn_ent_var.set(file_name)

    # Selects the skyglow map output folder for hemisphere
    def import_sgmap_folder(self):
        sgmap_dir = filedialog.askdirectory(initialdir='/', title="Select skyglow map folder")
        self.sgmap_folder_var.set(sgmap_dir)

        if sgmap_dir is '':
            print('Directory is empty.')

    def sng_popup(self):
        self.remove_all()

        self.check_lbl.grid(column=0, row=2)
        self.krn_chk.place(relx=.22, rely=.41, anchor=CENTER)

        self.file_lbl.grid(column=0, row=1)
        self.file_log.grid(column=1, columnspan=3, row=1)
        self.browse_btn.grid(column=4, row=1, sticky=W, padx=3)

        self.lat_lbl.grid(column=0, row=3)
        self.lat_entry.grid(column=1, row=3)

        self.k_lbl.grid(column=2, row=3)
        self.k_entry.grid(column=3, row=3)

        self.zen_lbl.grid(column=0, row=4)
        self.zen_entry.grid(column=1, row=4)

        self.azi_lbl.grid(column=2, row=4)
        self.azi_entry.grid(column=3, row=4)

        self.map_btn.grid(column=1, columnspan=3, row=5, sticky=N+S+E+W)

    def krn_popup(self):
        self.remove_all()

        # latitude
        self.lat_lbl.grid(column=0, row=3)
        self.lat_entry.grid(column=1, row=3)

        # atmospheric clarity
        self.k_lbl.grid(column=2, row=3)
        self.k_entry.grid(column=3, row=3)

        # angles file
        self.csv_file_lbl.grid(column=0, row=1)
        self.csv_file_log.grid(column=1, columnspan=3, row=1)
        self.csv_browse_btn.grid(column=4, row=1,sticky=W, padx=3)

        # input VIIRS image
        self.file_lbl.grid(column=0, row=2)
        self.file_log.grid(column=1, columnspan=3, row=2)
        self.browse_btn.grid(column=4, row=2, sticky=W, padx=3)

        self.hem_chk_lbl.grid(column=0, row=4)
        self.hem_chk.place(relx=.21, rely=.69)

        self.gen_krn_btn.grid(column=1, columnspan=3, row=5, sticky=N+S+E+W)

    def hem_popup(self):
        self.remove_all()

        # Skyglow Map Folder
        self.sgmap_folder_lbl.grid(column=0, row=1)
        self.sgmap_folder_log.grid(column=1, columnspan=3, row=1)
        self.sgmap_folder_btn.grid(column=4, row=1, sticky=W, padx=3)

        # Latitude entry
        self.lat_lbl.grid(column=0, row=3)
        self.lat_entry.grid(column=1, row=3)

        # Longitude entry
        self.lon_lbl.grid(column=2, row=3)
        self.lon_entry.grid(column=3, row=3)

        # Generate Hemispherical Visualization button
        self.hem_gen_btn.grid(column=1, columnspan=3, row=4, sticky=N+S+E+W)

    def mul_popup(self):
        self.remove_all()

        # Kernel folder location
        self.mul_file_lbl.grid(column=0, row=1)
        self.mul_file_log.grid(column=1, columnspan=3, row=1)
        self.mul_browse_btn.grid(column=4, row=1, sticky=W, padx=3)

        # input VIIRS image
        self.file_lbl.grid(column=0, row=2)
        self.file_log.grid(column=1, columnspan=3, row=2)
        self.browse_btn.grid(column=4, row=2, sticky=W, padx=3)

        # Choose output location
        self.output_lbl.grid(column=0, row=3)
        self.output_log.grid(column=1, columnspan=3, row=3)
        self.output_btn.grid(column=4, row=3, sticky=W, padx=3)

        # Generate map from kernel folder
        self.mul_map_btn.grid(column=1, columnspan=3, row=4, sticky=N+S+E+W)

    def remove_all(self):
        self.check_lbl.grid_remove()
        self.krn_chk.place_forget()
        self.hem_chk.place_forget()
        self.hem_chk_lbl.grid_remove()
        self.file_lbl.grid_remove()
        self.file_log.grid_remove()
        self.browse_btn.grid_remove()
        self.krn_lbl.grid_remove()
        self.krn_ent.grid_remove()
        self.krn_btn.grid_remove()
        self.lat_lbl.grid_remove()
        self.lat_entry.grid_remove()
        self.k_lbl.grid_remove()
        self.k_entry.grid_remove()
        self.zen_lbl.grid_remove()
        self.zen_entry.grid_remove()
        self.azi_lbl.grid_remove()
        self.azi_entry.grid_remove()
        self.map_btn.grid_remove()
        self.gen_krn_btn.grid_remove()
        self.mul_map_btn.grid_remove()
        self.csv_file_lbl.grid_remove()
        self.csv_file_log.grid_remove()
        self.csv_browse_btn.grid_remove()
        self.mul_file_lbl.grid_remove()
        self.mul_file_log.grid_remove()
        self.mul_browse_btn.grid_remove()
        self.output_lbl.grid_remove()
        self.output_log.grid_remove()
        self.output_btn.grid_remove()
        self.hem_gen_btn.grid_remove()
        self.lat_lbl.grid_remove()
        self.lat_entry.grid_remove()
        self.lon_lbl.grid_remove()
        self.lon_entry.grid_remove()
        self.sgmap_folder_lbl.grid_remove()
        self.sgmap_folder_log.grid_remove()
        self.sgmap_folder_btn.grid_remove()

    # Changes interface based on whether Kernel Checkbutton is selected.
    def checkbtn_val(self):
        # Import Kernel File widgets when Kernel Checkbutton is marked.
        if self.krn_var.get():
            self.lat_lbl.grid_remove()
            self.lat_entry.grid_remove()
            self.k_lbl.grid_remove()
            self.k_entry.grid_remove()
            self.zen_lbl.grid_remove()
            self.zen_entry.grid_remove()
            self.azi_lbl.grid_remove()
            self.azi_entry.grid_remove()
            self.krn_lbl.grid(column=0, row=2)
            self.krn_ent.grid(column=1, columnspan=3, row=2)
            self.krn_btn.grid(column=4, row=2, sticky=W, padx=3)
            self.krn_chk.place_forget()
            self.krn_chk.place(relx=0.19, rely=.5)
        # Input parameter widgets when Kernel Checkbuttton is unmarked
        else:
            self.krn_lbl.grid_remove()
            self.krn_ent.grid_remove()
            self.krn_btn.grid_remove()
            self.lat_lbl.grid(column=0, row=3)
            self.lat_entry.grid(column=1, row=3)
            self.k_lbl.grid(column=2, row=3)
            self.k_entry.grid(column=3, row=3)
            self.zen_lbl.grid(column=0, row=4)
            self.zen_entry.grid(column=1, row=4)
            self.azi_lbl.grid(column=2, row=4)
            self.azi_entry.grid(column=3, row=4)
            self.krn_chk.place_forget()
            self.krn_chk.place(relx=0.22, rely=.41, anchor=CENTER)
    # Simple function for opening URLs.
    @staticmethod
    def open_url(url):
        webbrowser.open_new(url)

    # Opens an instruction window that guides the user through the inputs.
    def instructions(self):
        # Instantiates separate Toplevel instruction window.
        instr_window = Toplevel(self.root)
        instr_window.geometry('550x575+25+25')
        instr_window.title('Instructions')
        instr_window.wm_iconbitmap(constants.ICO)
        instr_window.resizable(False, False)

        # Creatse Scrollbar and Frame for containing other widgets.
        instr_scroll = Scrollbar(instr_window)
        instr_scroll.pack(fill=Y, side="right")
        instr_frame = Frame(instr_window, bg='white')
        instr_frame.pack(fill=BOTH, side="left")

        # Adds instruction text from constants and adds image of Cinzano's diagram.
        instr = Text(instr_frame, width=65, height=40, padx=10, pady=5, bd=0, wrap="word")
        instr.insert("end", constants.INSTR)
        cdiagram_file = Image.open("./static/cinzano_diagram.PNG")
        cdiagram_file = cdiagram_file.resize((500, 450), Image.ANTIALIAS)
        self.cdiag = ImageTk.PhotoImage(cdiagram_file)
        instr.image_create("end", image=self.cdiag)
        instr.tag_add("top", "1.0", "4.10")
        instr.tag_config("top", font='Times 12 bold')
        instr.tag_add("body", "5.0", "19.20")
        instr.tag_config("body", font='Times 12')
        instr.insert("end", constants.CDIAG)
        instr.pack()
        instr_scroll.config(command=instr.yview)

    # Opens an about window that gives authors, SET version number, and icon credit.
    def about(self):
        # Instantiates a new Toplevel about window.
        about_window = Toplevel(self.root)
        about_window.geometry('350x335+25+25')
        about_window.title('About')
        about_window.wm_iconbitmap(constants.ICO)
        about_window.resizable(False, False)

        # Adds text to about window.
        about = Text(about_window, width=50, height=30, padx=10, pady=3)
        about.insert("end", constants.ABOUT)
        about.tag_add("abt", "1.0", "21.30")
        about.tag_config("abt", font='Times 10 bold', justify=CENTER)
        about.pack()

    # Constructs a progress window to monitor darkskypy.
    def progress(self):
        # Instantiates a new Toplevel window and frame for progress bar and loading log.
        prg_window = Toplevel(self.root)
        prg_window.geometry('650x325+250+250')
        prg_window.title('Generating Artificial Skyglow Map...')
        prg_window.iconbitmap(constants.ICO)
        prg_window.resizable(False, False)
        prg_frame = Frame(prg_window)
        prg_frame.pack(fill=BOTH)

        # Creates Scrollbar, Progressbar, and Label for checking progress..
        prg_scroll = Scrollbar(prg_frame)
        prg_scroll.pack(fill=Y, side="right")
        self.prg_bar = ttk.Progressbar(prg_frame, orient=HORIZONTAL, length=750,
                                  mode='indeterminate')
        self.prg_bar.pack()
        self.prg_bar.start()
        prg_lbl_txt = StringVar()
        prg_lbl = Label(prg_frame, textvariable=prg_lbl_txt)
        prg_lbl.pack()

        # Displays message log that prints from log file and starts darkskypy.
        self.prg_log = Text(prg_frame, width=90, padx=5, pady=5, relief="sunken")
        self.prg_log.pack()
        self.prg_log.insert("end", "*****Progress Log*****\n=======================\n")
        self.prg_log.tag_add("abt", "1.0", "3.0")
        self.prg_log.tag_config("abt", font='Courier 12 bold', justify=CENTER)
        self.txt_redir = LogRedirector(self.prg_log)
        logger.addHandler(self.txt_redir)
        sys.stderr = StderrRedirector(self.prg_log)
        prg_lbl_txt.set("Start time: " + str(time.asctime()))

    # Updates progress window to prevent it from freezing.
    def update_progress(self):
        self.prg_log.update()
        if len(threading.enumerate()) == 1:
            self.prg_bar.stop()
        else:
            self.prg_bar.start()
        self.root.after(1000, self.update_progress)

    # Generates artificial skyglow map based on VIIRS reference and local parameters.
    def generate_map(self):
        # Acquires input arguments.
        lat_in, k_in, zen_in, azi_in, file_in, krn_file_in = 0, 0, 0, 0, '', ''
        if self.krn_var.get():
            krn_file_in = self.krn_ent_var.get()
        else:
            lat_in = float(self.lat_entry.get())
            k_in = float(self.k_entry.get())
            zen_in = float(self.zen_entry.get())
            azi_in = float(self.azi_entry.get())
        file_in = self.file_log_var.get()

        self.progress()

        # Create new threads to run light propagation model simultaneously.
        p_thread = threading.Thread(target=self.update_progress())
        t_thread = threading.Thread(target=darksky.sgmapper,
                                    args=(lat_in, k_in, zen_in, azi_in, file_in, krn_file_in))
        t_thread.setDaemon(True)
        p_thread.start()
        t_thread.start()

    # Generate Kernels
    def generate_krn(self):
        # Acquires input arguments
        csv_in, file_in, lat_in, k_in, hem = '', '', 0, 0, False
        csv_in = self.csv_file_var.get()
        file_in = self.file_log_var.get()
        lat_in = float(self.lat_entry.get())
        k_in = float(self.k_entry.get())
        hem = self.hem_var.get()

        self.progress()

        # Create new threads to run light propagation model simultaneously.
        p_thread = threading.Thread(target=self.update_progress())
        with open(csv_in, "rb") as f:
            angle_list = loadtxt(f, delimiter=",", skiprows=1)
        p_thread.start()
        for angle_set in angle_list:
            t_thread = threading.Thread(target=darksky.generate_krn,
                                        args=(lat_in, k_in, angle_set[0],
                                              angle_set[1], file_in, hem))
            t_thread.setDaemon(True)

            t_thread.start()

    # Generate Multi Map from Kernel Library
    def generate_mmap(self):
        # Acquires input arguments
        krn_folder_in, file_in, output_in, = '', '', ''
        krn_folder_in = self.krn_folder_var.get()
        file_in = self.file_log_var.get()
        output_in = self.output_folder_var.get()

        self.progress()

        # Create new threads to run light propagation model simultaneously.
        p_thread = threading.Thread(target=self.update_progress())
        t_thread = threading.Thread(target=darksky.multisgmapper,
                                    args=(file_in, krn_folder_in, output_in))
        t_thread.setDaemon(True)
        p_thread.start()
        t_thread.start()

    #Generate Hemispherical Visualization of data from Skyglow map folder
    def generate_hem(self):
        sgmap_folder_in, lat_in, lon_in, = '', 0, 0
        sgmap_folder_in = self.sgmap_folder_var.get()
        lat_in = float(self.lat_entry.get())
        lon_in = float(self.lon_entry.get())
        darksky.generate_hem(lat_in, lon_in, sgmap_folder_in)

# Redirects formatted lines from the log file to a widget.
class LogRedirector(logging.Handler):
    def __init__(self, widget):
        logging.Handler.__init__(self)
        self.widget = widget

    def emit(self, record):
        self.widget.config(state='normal')
        self.widget.insert("end", self.format(record) + '\n')
        self.widget.see("end")
        self.widget.config(state='disabled')


class StderrRedirector(object):
    def __init__(self, widget):
        self.widget = widget

    def write(self, msg):
        self.widget.insert("end", msg)
        self.widget.see("end")


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    # Creates Tkinter root and initializes SET.
    logger.info('Starting SET program')
    root = Tk()
    tbox = SkyglowEstimationToolbox(root)

    # Sets up main screen of SET window.
    logger.info('Setting up main screen')
    tbox.main_screen()

    # Runs program.
    root.mainloop()
    logger.removeHandler(tbox.txt_redir)
    logger.info('Terminated SET Program')

if __name__ == '__main__':
    main()
