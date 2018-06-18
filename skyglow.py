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
try:
    from Tkinter import (Tk, Toplevel, PanedWindow, Frame, Label, Entry, Button,
                         Canvas, Scrollbar, Text, Menubutton, Menu, Checkbutton,
                         BOTH, VERTICAL, HORIZONTAL, CENTER, NE, E, W, Y,
                         StringVar, IntVar)
    import ttk
    import tkFileDialog as filedialog
except:
    from tkinter import (Tk, Toplevel, PanedWindow, Frame, Label, Entry, Button,
                         Canvas, Scrollbar, Text, Menubutton, Menu, Checkbutton,
                         BOTH, VERTICAL, HORIZONTAL, CENTER, NE, E, W, Y,
                         StringVar, IntVar, ttk, filedialog)
import matplotlib
matplotlib.use('TkAgg')
from PIL import Image, ImageTk, ImageOps
import webbrowser
import threading
import sys, os
import time
import constants
import darkskypy
import logging
logger = logging.getLogger()


# Main class of the software. Establishes GUI.
class SkyglowEstimationToolbox:

    ##################################################
    # Initialization
    def __init__(self, root):
        self.root = root
        self.file_log_var = StringVar()
        self.krn_ent_var = StringVar()
        self.krn_var = IntVar()
        self.img, self.cdiag = None, None
        self.lat_lbl, self.lat_entry = None, None
        self.k_lbl, self.k_entry = None, None
        self.zen_lbl, self.zen_entry = None, None
        self. azi_lbl, self.azi_entry = None, None
        self. krn_lvl, self.krn_entry, self.krn_btn = None, None, None
        self.txt_redir, self.prg_log = None, None
        self.map_btn = None

        # Sets window title, size, and icon on screen.
        self.root.title("Skyglow Estimation Toolbox (SET)")
        self.root.geometry('%dx%d+%d+%d' % (constants.SW*0.75, constants.SH*0.75, 25, 25))
        self.root.iconbitmap(os.path.join(os.getcwd(), constants.ICO))
        self.root.resizable(False, False)

        # Creates three paned windows for the main screen.
        base = PanedWindow()
        base.pack(fill=BOTH, expand=1)
        sub1 = PanedWindow(base, orient=VERTICAL)
        base.add(sub1)
        sub2 = PanedWindow(sub1, orient=HORIZONTAL)
        sub1.add(sub2)

        # Creates frame for holding inputs.
        self.input_frame = Frame(sub2)
        sub2.add(self.input_frame)

        # Creates frame for bottom half of main screen.
        self.img_frame = Frame(sub1, bd=2, bg='white', relief="sunken")
        sub1.add(self.img_frame)

        # Creates canvas for displaying images.
        self.img_canvas = Canvas(self.img_frame, bd=2, relief="groove",
                                 width=constants.SW*0.6, height=constants.SH*0.6)
        print(constants.SW, constants.SH)
        self.img_canvas.place(relx=.5, rely=.5, anchor=CENTER)

        # Creates help button for link to documentation, instructions, and about.
        self.help_btn = Menubutton(self.input_frame, text="Help", relief="raised",
                                   bd=2, width=8, pady=1)
        self.help_btn.place(relx=1, rely=0, anchor=NE)
        self.help_btn_menu = Menu(self.help_btn, tearoff=0)
        doc = 'https://github.com/NASA-DEVELOP'
        self.help_btn_menu.add_command(label="Documentation", command=lambda: self.open_url(doc))
        self.help_btn_menu.add_command(label="Instructions", command=self.instructions)
        self.help_btn_menu.add_separator()
        self.help_btn_menu.add_command(label="About", command=self.about)
        self.help_btn["menu"] = self.help_btn_menu

    # Sets up input GUI and image display screen.
    def main_screen(self):
        # VIIRS Image Reference File
        file_lbl = Label(self.input_frame, text="Image File:", width=15, anchor=E)
        file_log = Entry(self.input_frame, width=93, bd=2, relief="sunken",
                         textvariable=self.file_log_var)
        browse_btn = Button(self.input_frame, text="Browse", command=self.import_viirs)
        file_lbl.grid(column=0, row=0)
        file_log.grid(column=1, columnspan=3, row=0)
        browse_btn.grid(column=4, row=0, sticky=W, padx=3)

        # Import Kernel Checkbutton
        check_lbl = Label(self.input_frame, text="Import Kernel:", width=15, anchor=E)
        check_lbl.grid(column=0, row=1)
        krn_chk = Checkbutton(self.input_frame, anchor=W, variable=self.krn_var,
                              command=self.checkbtn_val)
        krn_chk.place(relx=.1, rely=.31, anchor=CENTER)

        # Region Latitude (deg), Grand Teton National park = 43.7904 degrees N
        self.lat_lbl = Label(self.input_frame, text="Latitude (deg):", width=15, anchor=E)
        self.lat_entry = Entry(self.input_frame, width=30, bd=2, relief="sunken")

        # Atmospheric Clarity Parameter, REF 2, Eq. 12, p. 645
        self.k_lbl = Label(self.input_frame, text="Atmospheric Clarity Parameter:",
        				   width=25, anchor=E)
        self.k_entry = Entry(self.input_frame, width=30, bd=2, relief="sunken")

        # Zenith angle (deg), z, REF 2, Fig. 6, p.648
        self.zen_lbl = Label(self.input_frame, text="Zenith Angle (deg):", width=15, anchor=E)
        self.zen_entry = Entry(self.input_frame, width=30, bd=2, relief="sunken")

        # Azimuth angle (deg)
        self.azi_lbl = Label(self.input_frame, text="Azimuth Angle (deg):", width=25, anchor=E)
        self.azi_entry = Entry(self.input_frame, width=30, bd=2, relief="sunken")

        self.krn_lbl = Label(self.input_frame, text="Kernel File:", width=15, anchor=E)
        self.krn_ent = Entry(self.input_frame, width=93, bd=2, relief="sunken",
                             textvariable=self.krn_ent_var)
        self.krn_btn = Button(self.input_frame, text="Browse", command=self.import_krn)

        # Generate Artificial Skyglow Map Button
        self.map_btn = Button(self.input_frame, text="Generate Artificial Skyglow Map",
                                 width=79, command=self.generate_map)

        # Places widgets.
        self.lat_lbl.grid(column=0, row=2)
        self.lat_entry.grid(column=1, row=2)
        self.k_lbl.grid(column=2, row=2)
        self.k_entry.grid(column=3, row=2)
        self.zen_lbl.grid(column=0, row=3)
        self.zen_entry.grid(column=1, row=3)
        self.azi_lbl.grid(column=2, row=3)
        self.azi_entry.grid(column=3, row=3)
        self.map_btn.grid(column=1, columnspan=3, row=4)
    ##################################################

    # Imports a TIFF file for referencing sky brightness in the region of interest.
    def import_viirs(self):
        # Allows user to search through his directory for VIIRS Image file.
        file_types = [('TIFF Files', '*.tif'), ('All files', '*')]
        file_name = filedialog.Open(initialdir='/', title="Select file", filetypes=file_types)
        file_name = file_name.show()
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

    # Imports a TIFF file containing the kernel data from previous input parameters.
    def import_krn(self):
        # Allows user to search through his directory for kernel file.
        file_types = [('TIFF Files', '*.tif'), ('All files', '*')]
        file_name = tkFileDialog.Open(initialdir='/', title="Select file", filetypes=file_types)
        file_name = file_name.show()
        self.krn_ent_var.set(file_name)

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
            self.map_btn.grid(column=1, columnspan=3, row=3)
        # Input parameter widgets when Kernel Checkbuttton is unmarked
        else:
            self.krn_lbl.grid_remove()
            self.krn_ent.grid_remove()
            self.krn_btn.grid_remove()
            self.lat_lbl.grid(column=0, row=2)
            self.lat_entry.grid(column=1, row=2)
            self.k_lbl.grid(column=2, row=2)
            self.k_entry.grid(column=3, row=2)
            self.zen_lbl.grid(column=0, row=3)
            self.zen_entry.grid(column=1, row=3)
            self.azi_lbl.grid(column=2, row=3)
            self.azi_entry.grid(column=3, row=3)
            self.map_btn.grid(column=1, columnspan=3, row=4)

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
        t_thread = threading.Thread(target=darkskypy.sgmapper,
                                    args=(lat_in, k_in, zen_in, azi_in, file_in, krn_file_in))
        t_thread.setDaemon(True)
        p_thread.start()
        t_thread.start()


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
    logging.basicConfig(level=logging.DEBUG)
    main()
