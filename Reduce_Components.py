__author__ = 'kgulliks'

import PL_Display
import numpy as np
import sys
from scipy.signal import argrelmax, argrelmin
import matplotlib.pyplot as plt


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
    #apertures, tops, bottoms = FindApertures(flat, hdr)

    #Trace the apertures
    ap_coeffs1, ap_coeffs2 = TraceApertures(flat, hdr, band=band, target_path="%s/aperture_mapping/" %files.WorkingDirectory)





def FindApertures(icol, threshold=0.25, medsize=20):
    """
      This function will find the apertures in a given image
    :param icol: A 1d array containing the collapsed chip (see definition in TraceApertures)
    :param threshold: The fraction of the peak flux to include in the aperture
    :param medsize: The number of pixels to median-flatten to find the aperture
    :return: polynomial coefficients for the upper and lower bound
    """

    #Start by finding the apertures
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
        slope = (icol[m_greater] - icol[m_less])/(m_greater - m_less)
        icol2 = icol - slope*(pixels - m_less) + icol[m_less]
        ap = np.where(np.logical_and(icol2 > threshold*icol2[p],
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


def TraceApertures(image, header, band, startcol=1024, stepsize=40, medsize=20, ap_thresh=10, degree=3, target_path="./"):
    """
    Trace the apertures in the given image

    :param image:
    :param header:
    :param startcol: The column to start at (default is the middle)
    :param stepsize: The stepsize to take in finding the apertures
    :param medsize:  The number of pixels to median-flatten to find the aperture
    :param ap_thresh: Threshold for the allowed change in the aperture position in each step
    :param degree: The polynomial degree to fit each aperture to
    :return:
    """

    nx, ny = image.shape

    #Start at startcol
    icol = np.median(image[:,(startcol-medsize/2):(startcol+medsize/2)], axis=1)
    peaks, mtops, mbottoms = FindApertures(icol,medsize=medsize)
    n_ap = len(peaks)

    # make set of (x, y) positions for each aperture
    # add the (x, y) positions at default column
    apx1, apx2 = [], []
    apy1, apy2 = [], []
    for top, bottom in zip(mtops, mbottoms):
        apx1.append([startcol])
        apx2.append([startcol])
        apy1.append([top])
        apy2.append([bottom])

    # define the tracing direction (left, right)
    xsteps = np.arange(medsize,nx-medsize,stepsize)
    xlefts = xsteps[(xsteps < startcol)]
    xrights = xsteps[(xsteps > startcol)]

    # the left side columns
    ptops, pbottoms = np.copy(mtops), np.copy(mbottoms)  #Previous tops and bottoms
    for xpos in reversed(xlefts):
        icol = np.mean(image[:,(xpos-medsize/2):(xpos+medsize/2)],axis=1)
        peaks, tops, bottoms = np.array(FindApertures(icol, medsize=medsize))

        # check the aperture number based on original positions
        for i, ptop, pbottom in zip(range(n_ap), ptops, pbottoms):
            top = tops[(tops < ptop+ap_thresh) & ( tops > ptop-ap_thresh)]
            if len(top) == 1:
                # save this position to the (x, y) position list
                apx1[i].append(xpos)
                apy1[i].append(top[0])
                # save this position for the next tracing loop
                ptops[i] = top[0]
            elif len(top) == 0:
                print 'ap[%d](x=%d) : The matching aperture position was not found' \
                   % (i, xpos)
            else:
                print 'ap[%d](x=%d) : The matching aperture position was too many(%d)' \
                   % (i, xpos, len(top))

            bottom = bottoms[(bottoms < pbottom+ap_thresh) & ( bottoms > pbottom-ap_thresh)]
            if len(bottom) == 1:
                # save this position to the (x, y) position list
                apx2[i].append(xpos)
                apy2[i].append(bottom[0])
                # save this position for the next tracing loop
                pbottoms[i] = bottom[0]
            elif len(bottom) == 0:
                print 'ap[%d](x=%d) : The matching aperture position was not found' \
                   % (i, xpos)
            else:
                print 'ap[%d](x=%d) : The matching aperture position was too many(%d)' \
                   % (i, xpos, len(tcap2))


    # the right side columns
    ptops, pbottoms = np.copy(mtops), np.copy(mbottoms)  #Previous tops and bottoms
    for xpos in xrights:
        icol = np.mean(image[:,(xpos-medsize/2):(xpos+medsize/2)],axis=1)
        peaks, tops, bottoms = np.array(FindApertures(icol, medsize=medsize))

        # check the aperture number based on original positions
        for i, ptop, pbottom in zip(range(n_ap), ptops, pbottoms):
            top = tops[(tops < ptop+ap_thresh) & ( tops > ptop-ap_thresh)]
            if len(top) == 1:
                # save this position to the (x, y) position list
                apx1[i].append(xpos)
                apy1[i].append(top[0])
                # save this position for the next tracing loop
                ptops[i] = top[0]
            elif len(top) == 0:
                print 'ap[%d](x=%d) : The matching aperture position was not found' \
                   % (i, xpos)
            else:
                print 'ap[%d](x=%d) : The matching aperture position was too many(%d)' \
                   % (i, xpos, len(top))

            bottom = bottoms[(bottoms < pbottom+ap_thresh) & ( bottoms > pbottom-ap_thresh)]
            if len(bottom) == 1:
                # save this position to the (x, y) position list
                apx2[i].append(xpos)
                apy2[i].append(bottom[0])
                # save this position for the next tracing loop
                pbottoms[i] = bottom[0]
            elif len(bottom) == 0:
                print 'ap[%d](x=%d) : The matching aperture position was not found' \
                   % (i, xpos)
            else:
                print 'ap[%d](x=%d) : The matching aperture position was too many(%d)' \
                   % (i, xpos, len(tcap2))


    # sorting the (x, y) positions along x-direction for each aperture
    # polynomial fitting with degree
    ap_coeffs1, ap_coeffs2 = [], []
    z1, z2 = PL_Display.zscale(image)
    plt.figure(figsize=(10,10))
    plt.imshow(image, cmap='gray', vmin=z1, vmax=z2)
    for k in range(len(apx1)): #n_ap):
        tapx1 = np.array(apx1[k],dtype=np.float)
        tapx2 = np.array(apx2[k],dtype=np.float)
        tapy1 = np.array(apy1[k],dtype=np.float)
        tapy2 = np.array(apy2[k],dtype=np.float)

        tsort1 = np.argsort(tapx1)
        tsort2 = np.argsort(tapx2)
        coeff1 = np.polynomial.polynomial.polyfit(tapx1[tsort1],tapy1[tsort1],degree)
        coeff2 = np.polynomial.polynomial.polyfit(tapx2[tsort2],tapy2[tsort2],degree)

        # save the fitting coefficients
        ap_coeff = np.zeros([2,degree+1])
        ap_coeff[0,:] = coeff1
        ap_coeff[1,:] = coeff2
        np.savetxt(target_path+'apmap_%s_%02d.%03d.dat' % (band, degree, k), ap_coeff)

        ap_coeffs1.append(coeff1)
        ap_coeffs2.append(coeff2)

        yfit1 = np.polynomial.polynomial.polyval(tapx1[tsort1],coeff1)
        yfit2 = np.polynomial.polynomial.polyval(tapx2[tsort2],coeff2)

        plt.plot(tapx1[tsort1],tapy1[tsort1], 'go')
        plt.plot(tapx2[tsort2],tapy2[tsort2], 'ro')

        print 'ap[%d]: delta_y1 = %12.8f %12.8f ' \
          % (k, np.mean(tapy1[tsort1]-yfit1), np.std(tapy1[tsort1]-yfit1))
        print 'ap[%d]: delta_y2 = %12.8f %12.8f ' \
          % (k, np.mean(tapy2[tsort2]-yfit2), np.std(tapy2[tsort2]-yfit2))

        xx = np.arange(nx)
        yfit1 = np.polynomial.polynomial.polyval(xx,coeff1)
        yfit2 = np.polynomial.polynomial.polyval(xx,coeff2)

        plt.plot(xx,yfit1,'g-', linewidth=3, alpha=0.6)
        plt.plot(xx,yfit2,'r-', linewidth=3, alpha=0.6)

    plt.xlim(0,nx)
    plt.ylim(0,ny)
    plt.show()
    #plt.savefig(target_path+name+'_aptracing2.pdf')
    plt.close('all')

    return ap_coeffs1, ap_coeffs2










