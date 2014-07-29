__author__ = 'kgulliks'

import PL_Display
import numpy as np
import sys


utilities = PL_Display.CDisplay("", "")

#TODO: Fix this! I need to make a new badpixcor it looks like!
def FlatReduce(files, band):
    """
    Function to reduce the flat frames.
    """

    #Read in the file
    imgsA, hdrsA = utilities.readechellogram(files.Flats["ON"])
    imgsB, hdrsB = utilities.readechellogram(files.Flats["OFF"])

    #Bad pixel corrections
    imgsA, hdrsA = utilities.barebadpixcor(imgsA, band, headers=hdrsA)
    imgsB, hdrsB = utilities.barebadpixcor(imgsB, band, headers=hdrsB)

    #Combine
    on, hdr_on = utilities.barecombine(imgsA, method="median", hdr=hdrsA[0])
    off, hdr_off = utilities.barecombine(imgsB, method="median", hdr=hdrsB[0])

    #Subtract on from off
    flat, hdr = utilities.baresubtract(on, off, hdr=hdr_on)

    # Do cosmic ray correction
    flat = PL_Display.ip.cosmicrays(flat, threshold=2000, flux_ratio=0.2, wbox=3)

    # Normalize
    flat, hdr = utilities.barenormalize(flat, hdr, band=band)

    #Save the file
    PL_Display.savefits("%s/FLAT.fits" %files.WorkingDirectory, flat, hdr)




