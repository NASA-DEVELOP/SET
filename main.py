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
        self.file = None            # *** Make sure to declare variables here that will be
        self.img = None             # used throughout the program ***

        # Creates window on screen.
        self.parent.title("Skyglow Estimation Toolbox (SET)")
        self.parent.geometry('1250x750')

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

        # Creates frame for displaying images.
        self.img_frame = Frame(sub1, bd=2, bg='white', relief=SUNKEN)
        sub1.add(self.img_frame)

        # Canvas for displaying images.
        # *** You should always declare higher level widgets like canvases and frames in the init
        # if possible. Also, no need to bind the canvas to an action, we'll call the action in
        # import_tiff ***
        self.img_canvas = Canvas(self.img_frame, bd=2, relief=GROOVE, width=750, height=500)
        self.img_canvas.place(relx=.5, rely=.5, anchor=CENTER)

    # Sets up input GUI and image display screen.
    def main_screen(self):
        # VIIRS Image Reference File
        file_label = Label(self.input_frame, text="Image File:", width=15, anchor=E)
        self.file_dialog = Entry(self.input_frame, width=110, bd=2, relief=SUNKEN)
        browse_button = Button(self.input_frame, text="Browse", command=self.import_tiff)
        file_label.grid(column=0, row=0)
        self.file_dialog.grid(column=1, columnspan=3, row=0)
        browse_button.grid(column=4, row=0, sticky=W, padx=3)

        # Region Latitude (deg), Grand Teton National park = 43.7904 degrees N
        lat_label = Label(self.input_frame, text="Latitude (deg):", width=15, anchor=E)
        lat_label.grid(column=0, row=1)
        lat_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        lat_entry.grid(column=1, row=1)

        # Distance of Increasing Integration Increment (km), ubr, REF 2, Fig. 6, p.648
        ubr_label = Label(self.input_frame, text="Distance at which Integration Speed Increases"
                                                 " (km):", width=40, anchor=E)
        ubr_label.grid(column=2, row=1)
        ubr_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        ubr_entry.grid(column=3, row=1)

        # Zenith angle (deg), z, REF 2, Fig. 6, p.648
        zen_label = Label(self.input_frame, text="Zenith Angle (deg):", width=15, anchor=E)
        zen_label.grid(column=0, row=2)
        zen_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        zen_entry.grid(column=1, row=2)

        # Azimuth angle (deg)
        azi_label = Label(self.input_frame, text="Azimuth Angle (deg):", width=40, anchor=E)
        azi_label.grid(column=2, row=2)
        azi_entry = Entry(self.input_frame, width=30, bd=2, relief=SUNKEN)
        azi_entry.grid(column=3, row=2)

        # Generate Artificial Skyglow Map
        map_button = Button(self.input_frame, text="Generate Artificial Skyglow Map", width=93)
        map_button.grid(column=1, columnspan=3, row=3)
    ##################################################

    # Imports a TIFF file for referencing sky brightness in the region of interest.
    # *** Its always a good idea to implement some sort of logging for sanity checking and
    # debugging. I just use print functions here ***
    def import_tiff(self):
        print('Importing VIIRS Image...')
        file_types = [('VIIRS Images', '*.tif'), ('All files', '*')]
        file_name = tkFileDialog.Open(initialdir='/', title="Select file", filetypes=file_types)
        file_name = file_name.show()
        print('Grabbed ' + str(file_name))
        if file_name != '':                         # *** Making sure an empty file wasn't chosen
            self.file = file_name                   # *** Set the permanent variable
            pilimg = Image.open(self.file)          # *** Open the image
            print('pilimg = ' + str(pilimg))
            self.img = ImageTk.PhotoImage(pilimg)   # *** Set the image as a class variable
            print('img = ' + str(self.img))
            self.img_canvas.create_image(375, 250, image=self.img)
        else:
            print('Fill is empty')
        # pilimg = Image.open(self.file_dialog.get())
        # img = ImageTk.PhotoImage(pilimg)
        # img_display = self.img_canvas.create_image(image=img)


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
