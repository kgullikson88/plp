__author__ = 'kgulliks'

import PL_Display
import numpy as np
import sys
from scipy.signal import argrelmax, argrelmin


utilities = PL_Display.CDisplay("", "")


def FlatReduce(files, band):
    """
    Function to reduce the flat frames.
    """
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
    """

    flat, hdr = PL_Display.readfits("%s/FLAT.fits" %files.WorkingDirectory)

    #Find the apertures
    apertures, tops, bottoms = FindApertures(flat, hdr)

    #Trace the apertures




def FindApertures(image, header, threshold=0.1):
    """
      This function will find the apertures in a given image
    :param image: the image the find the apertures in
    :param header: the header for that image. The coefficients will be added to the header
    :param threshold: The fraction of the peak flux to include in the aperture
    :return: polynomial coefficients for the upper and lower bound
    """

    #Start by finding the apertures
    nx, ny = image.shape
    icol = np.median(image[:,(nx/2-10):(nx/2+10)], axis=1)
    pixels = np.arange(icol.size)
    peaks = argrelmax(icol, order=20)[0]
    peaks = peaks[icol[peaks] > 0.1]
    minima = argrelmin(icol, order=20)[0]

    #Find the edges of each aperture
    top = []
    bottom = []
    goodpeaks = []
    for p in peaks:
        m_less = minima[minima < p]
        m_greater = minima[minima > p]
        if m_less.size == 0 or m_greater.size == 0:
            continue
        m_less = m_less[-1]
        m_greater = m_greater[0]
        ap = np.where(np.logical_and(icol > threshold*icol[p],
                                     np.logical_and(pixels > m_less,
                                                    pixels < m_greater)))[0]
        goodpeaks.append(p)
        top.append(ap[-1])
        bottom.append(ap[0])

    goodpeaks = np.array(goodpeaks)
    top = np.array(top)
    bottom = np.array(bottom)
    #import pylab
    #pylab.plot(icol)
    #pylab.errorbar(goodpeaks, icol[goodpeaks], xerr=(goodpeaks-bottom, top-goodpeaks), fmt='k-')
    #pylab.show()
    return goodpeaks, top, bottom


def TraceApertures(image, header, apertures, tops, bottoms):
    pass







