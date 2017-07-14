"""
Name: Skyglow Estimation Toolbox
Authors: Ryan Avery, Dr. Kenton Ross, Stanley Yu
Description: Python toolbox that generates artificial skyglow maps using data from NASA and NOAA's
             Visible Infrared Imaging Radiometer Suite (VIIRS) satellite sensor.
"""
from Tkinter import Tk, PanedWindow, Frame, Label, Entry, Button, Canvas, \
    BOTH, VERTICAL, HORIZONTAL, CENTER, E, W, \
    SUNKEN, GROOVE
import tkFileDialog
from PIL import Image, ImageTk


# Main class of the software. Establishes GUI.
class SkyglowEstimationToolbox:

    ##################################################
    # Initialization
    def __init__(self, parent):
        self.parent = parent
        self.file_dialog = None
        self.img = None
        self.lat_entry = None
        self.ubr_entry = None
        self.zen_entry = 0
        self.azi_entry = 0

        # Sets window title and size on screen.
        self.parent.title("Skyglow Estimation Toolbox (SET)")
        sw = self.parent.winfo_screenwidth()
        sh = self.parent.winfo_screenheight()
        self.parent.geometry('%dx%d' % (sw*0.75, sh*0.75))

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

        # Canvas for displaying images.
        self.img_canvas = Canvas(self.img_frame, bd=2, relief=GROOVE, width=sw*0.6, height=sh*0.6)
        self.img_canvas.place(relx=.5, rely=.5, anchor=CENTER)

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
            # ### FIGURE OUT HOW TO MAKE TIFF DISPLAY MORE VISIBLE??? ###
            # pilimg_alt = pilimg_r.convert("RGB", palette='WEB')
            self.img = ImageTk.PhotoImage(pilimg_r)
            self.img_canvas.create_image(canvas_size[0]/2, canvas_size[1]/2, image=self.img)
        else:
            print('File is empty.')

    # Generates artificial skyglow map based on VIIRS reference and local parameters.
    def generate_map(self):
        # Acquires input arguments.
        lat_in = self.lat_entry.get()
        ubr_in = self.ubr_entry.get()
        zen_in = self.zen_entry.get()
        azi_in = self.azi_entry.get()
        # file_in = self.file_dialog.get()

        print(lat_in, ubr_in, zen_in, azi_in)


def main():
    # Creates Tkinter root and initializes SET.
    root = Tk()
    tbox = SkyglowEstimationToolbox(root)

    # Sets up main screen of SET window.
    tbox.main_screen()

    # Runs program.
    root.mainloop()

if __name__ == '__main__':
    main()