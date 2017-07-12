"""
Name: Skyglow Estimation Toolbox
Authors: Ryan Avery, Dr. Kenton Ross, Stanley Yu
Description: Python toolbox that generates artificial skyglow maps using data from NASA and NOAA's
             Visible Infrared Imaging Radiometer Suite (VIIRS) satellite sensor.
"""
import logging
from Tkinter import Tk, PanedWindow, Frame, Label, Entry, Button, \
    BOTH, VERTICAL, HORIZONTAL, E, W

# Main class of the software. Establishes GUI.
class SkyglowEstimationToolbox:

    ##################################################
    # Initialization

    def __init__(self,parent):
        self.parent = parent

        # Creates window on screen.
        self.parent.title("Skyglow Estimation Toolbox (SET)")
        self.parent.geometry('1250x750')

        # Creates three paned windows for the main screen.
        base = PanedWindow()
        base.pack(fill=BOTH,expand=1)
        sub1 = PanedWindow(base,orient=VERTICAL)
        base.add(sub1)
        sub2 = PanedWindow(sub1,orient=HORIZONTAL)
        sub1.add(sub2)

        # Creates frame for holding inputs.
        self.input_frame = Frame(sub2)
        sub2.add(self.input_frame)

        # Creates frame for displaying images.
        self.img_frame = Frame(sub1,bg='blue')
        sub1.add(self.img_frame)

    # Sets up input GUI and image dispay screen.
    def main_screen(self):
        # VIIRS Image Reference File
        file_label = Label(self.input_frame,text="Image File:",width=15,anchor=E)
        file_entry = Entry(self.input_frame,width=97)
        browse_button = Button(self.input_frame,text="Browse")
        file_label.grid(column=0,row=0)
        file_entry.grid(column=1,columnspan=3,row=0)
        browse_button.grid(column=4,row=0,sticky=W,padx=3)

        # Region Latitude (deg), Grand Teton National park = 43.7904 degrees N
        lat_label = Label(self.input_frame,text="Latitude (deg):",width=15,anchor=E)
        lat_label.grid(column=0,row=1)
        lat_entry = Entry(self.input_frame,width=40)
        lat_entry.grid(column=1,row=1)

        # Distance from site to point in sky, relaxing increment (km), u, REF 2, Fig. 6, p.648
        ubr_label = Label(self.input_frame,text="Distance from site to night sky (km):",width=30,anchor=E)
        ubr_label.grid(column=2,row=1)
        ubr_entry = Entry(self.input_frame)
        ubr_entry.grid(column=3,row=1)

        # Zenith angle (deg), z, REF 2, Fig. 6, p.648
        zen_label = Label(self.input_frame,text="Zenith Angle (deg):",width=15,anchor=E)
        zen_label.grid(column=0,row=2)
        zen_entry = Entry(self.input_frame,width=40)
        zen_entry.grid(column=1,row=2)

        # Azimuth angle (deg)
        azi_label = Label(self.input_frame,text="Azimuth Angle (deg):",width=30,anchor=E)
        azi_label.grid(column=2,row=2)
        azi_entry = Entry(self.input_frame)
        azi_entry.grid(column=3,row=2)

        # Temporary label for image.
        temp_label = Label(self.img_frame,text="Image goes here")
        temp_label.pack()

    ##################################################

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