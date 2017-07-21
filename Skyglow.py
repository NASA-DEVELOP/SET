"""
Name: Skyglow Estimation Toolbox
Authors: Ryan Avery, Dr. Kenton Ross, Stanley Yu
Description: Python toolbox that generates artificial skyglow maps using data from NASA and NOAA's
             Visible Infrared Imaging Radiometer Suite (VIIRS) satellite sensor.
"""
from Tkinter import Tk, Toplevel, PanedWindow, Frame, Label, Entry, Button, Canvas, Scrollbar, \
    Text, Menubutton, Menu, BOTH, VERTICAL, HORIZONTAL, CENTER, LEFT, RIGHT, NE, E, W, Y, \
    WORD, END, SUNKEN, GROOVE, RAISED
from PIL import Image, ImageTk
import tkFileDialog
import webbrowser
import threading
import constants
import Itest
import logging
logger = logging.getLogger()


# Main class of the software. Establishes GUI.
class SkyglowEstimationToolbox:

    ##################################################
    # Initialization
    def __init__(self, root):
        self.root = root
        self.file_dialog = None
        self.img = None
        self.cdiag = None
        self.lat_entry = None
        self.ubr_entry = None
        self.zen_entry = None
        self.azi_entry = None

        # Sets window title, size, and icon on screen.
        self.root.title("Skyglow Estimation Toolbox (SET)")
        self.root.geometry('%dx%d+%d+%d' % (constants.SW*0.75, constants.SH*0.75, 25, 25))
        self.root.wm_iconbitmap(constants.ICO)
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
        self.img_frame = Frame(sub1, bd=2, bg='white', relief=SUNKEN)
        sub1.add(self.img_frame)

        # Creates canvas for displaying images.
        self.img_canvas = Canvas(self.img_frame, bd=2, relief=GROOVE,
                                 width=constants.SW*0.6, height=constants.SH*0.6)
        self.img_canvas.place(relx=.5, rely=.5, anchor=CENTER)

        # Creates help button for link to documentation, instructions, and about.
        self.help_btn = Menubutton(self.input_frame, text="Help", relief=RAISED,
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
        file_label = Label(self.input_frame, text="Image File:", width=15, anchor=E)
        self.file_dialog = Entry(self.input_frame, width=109, bd=2, relief=SUNKEN)
        browse_button = Button(self.input_frame, text="Browse", command=self.import_file)
        file_label.grid(column=0, row=0)
        self.file_dialog.grid(column=1, columnspan=3, row=0)
        browse_button.grid(column=4, row=0, sticky=W, padx=3)

        # Region Latitude (deg), Grand Teton National park = 43.7904 degrees N
        lat_label = Label(self.input_frame, text="Latitude (deg):", width=15, anchor=E)
        lat_label.grid(column=0, row=1)
        self.lat_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        self.lat_entry.grid(column=1, row=1)

        # Distance of Increasing Integration Increment (km), ubr, REF 2, Fig. 6, p.648
        ubr_label = Label(self.input_frame, text="Distance at which Integration Speed Increases"
                                                 " (km):", width=40, anchor=E)
        ubr_label.grid(column=2, row=1)
        self.ubr_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        self.ubr_entry.grid(column=3, row=1)

        # Zenith angle (deg), z, REF 2, Fig. 6, p.648
        zen_label = Label(self.input_frame, text="Zenith Angle (deg):", width=15, anchor=E)
        zen_label.grid(column=0, row=2)
        self.zen_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        self.zen_entry.grid(column=1, row=2)

        # Azimuth angle (deg)
        azi_label = Label(self.input_frame, text="Azimuth Angle (deg):", width=40, anchor=E)
        azi_label.grid(column=2, row=2)
        self.azi_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        self.azi_entry.grid(column=3, row=2)

        # Generate Artificial Skyglow Map
        map_button = Button(self.input_frame, text="Generate Artificial Skyglow Map", width=93,
                            command=self.generate_map)
        map_button.grid(column=1, columnspan=3, row=3)
    ##################################################

    # Imports a TIFF file for referencing sky brightness in the region of interest.
    def import_file(self):
        file_types = [('VIIRS Images', '*.tif'), ('All files', '*')]
        file_name = tkFileDialog.Open(initialdir='/', title="Select file", filetypes=file_types)
        file_name = file_name.show()
        self.file_dialog.insert(0, file_name)
        if file_name != '':
            pilimg = Image.open(file_name)
            canvas_size = (self.img_canvas.winfo_width(), self.img_canvas.winfo_height())
            pilimg_r = pilimg.resize(canvas_size, Image.ANTIALIAS)
            self.img = ImageTk.PhotoImage(pilimg_r)
            self.img_canvas.create_image(canvas_size[0]/2, canvas_size[1]/2, image=self.img)
        else:
            print('File is empty.')

    @staticmethod
    def open_url(url):
        webbrowser.open_new(url)

    def instructions(self):
        instr_window = Toplevel(self.root)
        instr_window.geometry('550x575+25+25')
        instr_window.title('Instructions')
        instr_window.wm_iconbitmap(constants.ICO)
        instr_window.resizable(False, False)

        instr_scroll = Scrollbar(instr_window)
        instr_scroll.pack(fill=Y, side=RIGHT)

        instr_frame = Frame(instr_window, bg='white')
        instr_frame.pack(fill=BOTH, side=LEFT)

        instr = Text(instr_frame, width=65, height=40, padx=10, pady=5, bd=0, wrap=WORD)
        instr.insert(END, constants.INSTR)
        cdiagram_file = Image.open("cinzano_diagram.PNG")
        cdiagram_file = cdiagram_file.resize((500, 450), Image.ANTIALIAS)
        self.cdiag = ImageTk.PhotoImage(cdiagram_file)
        instr.image_create(END, image=self.cdiag)
        instr.tag_add("top", "1.0", "4.10")
        instr.tag_config("top", font='Times 12 bold')
        instr.tag_add("body", "5.0", "19.20")
        instr.tag_config("body", font='Times 12')
        instr.insert(END, constants.CDIAG)
        instr.pack()
        instr_scroll.config(command=instr.yview)

    def about(self):
        about_window = Toplevel(self.root)
        about_window.geometry('350x335+25+25')
        about_window.title('About')
        about_window.wm_iconbitmap(constants.ICO)
        about_window.resizable(False, False)

        about = Text(about_window, width=50, height=30, padx=10, pady=3)
        about.insert(END, constants.ABOUT)
        about.tag_add("abt", "1.0", "21.30")
        about.tag_config("abt", font='Times 10 bold', justify=CENTER)
        about.pack()

    def update_root(self):
        self.root.update()
        self.root.after(1000, self.update_root)

    # Generates artificial skyglow map based on VIIRS reference and local parameters.
    def generate_map(self):
        # Acquires input arguments.
        lat_in = float(self.lat_entry.get())
        ubr_in = float(self.ubr_entry.get())
        zen_in = float(self.zen_entry.get())
        azi_in = float(self.azi_entry.get())
        file_in = self.file_dialog.get()

        # Create new threads
        r_thread = threading.Thread(target=self.update_root)
        t_thread = threading.Thread(target=Itest.main, 
                                    args=(lat_in, ubr_in, zen_in, azi_in, file_in))

        # Start new threads
        logger.info('Starting threads, running Itest')
        r_thread.start()
        t_thread.start()


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
    logger.info('Terminated SET Program')

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main()