__author__ = 'kgulliks'

import PL_Display
import numpy as np
import sys
from scipy.signal import argrelmax, argrelmin
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits


utilities = PL_Display.CDisplay("", "")


def FlatReduce(files, band, debug=False):
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
    print "Saving master flat to %s/FLAT.fits" %files.WorkingDirectory
    PL_Display.savefits("%s/FLAT.fits" %files.WorkingDirectory, flat, hdr)

    """
    flat, hdr = PL_Display.readfits("%s/FLAT.fits" %files.WorkingDirectory)

    #Find the apertures
    startcol = 1024
    medsize = 20
    icol = np.median(flat[:,(startcol-medsize/2):(startcol+medsize/2)], axis=1)
    apertures, tops, bottoms = FindApertures(icol)
    #sys.exit()

    #Trace the apertures
    ap_coeffs1, ap_coeffs2 = TraceApertures(flat,
                                            hdr,
                                            band=band,
                                            target_path="%s/aperture_mapping/" %files.WorkingDirectory,
                                            debug=debug)

    return ap_coeffs1, ap_coeffs2


def ArcReduce(files, band, aptops=None, apbottoms=None, debug=False):
    """
    Reduce the arc file: dark-correct, flat-correct, extract,
    and fit the wavelength solution

    :param files: The files structure
    :param band: Either H or K
    :param aptops: The coefficients for the tops of the apertures
    :param apbottoms: The coefficients for the bottoms of the apertures
    :param debug:
    :return:
    """

    # Read in the arc files
    imgsA, hdrsA = utilities.readechellogram(files.Arcs["ON"])
    imgsB, hdrsB = utilities.readechellogram(files.Arcs["OFF"])

    #Bad pixel corrections
    imgsA, hdrsA = utilities.barebadpixcor(imgsA, band, headers=hdrsA)
    imgsB, hdrsB = utilities.barebadpixcor(imgsB, band, headers=hdrsB)

    #Combine
    on, hdr_on = utilities.barecombine(imgsA, method="median", hdr=hdrsA[0])
    off, hdr_off = utilities.barecombine(imgsB, method="median", hdr=hdrsB[0])

    #Subtract on from off
    arc, hdr = utilities.baresubtract(on, off, hdr=hdr_on)

    #Divide by the flat
    flat, fhdr = PL_Display.readfits("%s/FLAT.fits" %files.WorkingDirectory)
    arc, hdr = utilities.flatcorrect([arc], flat, headers=[hdr])
    arc = arc[0]
    hdr = hdr[0]

    #Save the file
    print "Saving master arc to %s/ARC.fits" %files.WorkingDirectory
    PL_Display.savefits("%s/ARC.fits" %files.WorkingDirectory, arc, hdr)

    #Extract the apertures using the apertures defined by the flat
    if aptops is None or apbottoms is None:
        aptops, apbottoms = TraceApertures(flat,
                                           fhdr,
                                           band=band,
                                           target_path="%s/aperture_mapping/" %files.WorkingDirectory,
                                           debug=debug)
    spectra, ypos = ExtractFromApertures(arc, hdr, aptops, apbottoms)

    #Fit the Arc spectrum to determine the 2d wavelength calibration
    wave_coeffs = FitDispersion(spectra,
                                ypos,
                                band=band,
                                outpath="%s/aperture_mapping/" %files.WorkingDirectory)






    return None




def FindApertures(icol, threshold=0.25, medsize=20, debug=False):
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
    if debug:
      plt.plot(icol)
      plt.errorbar(goodpeaks, icol[goodpeaks], xerr=(goodpeaks-bottom, top-goodpeaks), fmt='ko')
      plt.show()
    return goodpeaks, top, bottom


def TraceApertures(image, header, band, startcol=1024, stepsize=40, medsize=20, ap_thresh=10, degree=3, target_path="./", debug=False):
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
                   % (i, xpos, len(bottom))


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
                   % (i, xpos, len(bottom))


    # sorting the (x, y) positions along x-direction for each aperture
    # polynomial fitting with degree
    ap_coeffs1, ap_coeffs2 = [], []
    if debug:
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

        if debug:
          plt.plot(tapx1[tsort1],tapy1[tsort1], 'go')
          plt.plot(tapx2[tsort2],tapy2[tsort2], 'ro')

        print 'ap[%d]: delta_y1 = %12.8f %12.8f ' \
          % (k, np.mean(tapy1[tsort1]-yfit1), np.std(tapy1[tsort1]-yfit1))
        print 'ap[%d]: delta_y2 = %12.8f %12.8f ' \
          % (k, np.mean(tapy2[tsort2]-yfit2), np.std(tapy2[tsort2]-yfit2))

        xx = np.arange(nx)
        yfit1 = np.polynomial.polynomial.polyval(xx,coeff1)
        yfit2 = np.polynomial.polynomial.polyval(xx,coeff2)

        if debug:
          plt.plot(xx,yfit1,'g-', linewidth=3, alpha=0.6)
          plt.plot(xx,yfit2,'r-', linewidth=3, alpha=0.6)

    if debug:
      plt.xlim(0,nx)
      plt.ylim(0,ny)
      plt.show()
      plt.savefig("%sFLAT_aptracing.pdf" %target_path)
      plt.close('all')

    return ap_coeffs1, ap_coeffs2




def ExtractFromApertures(img, hdr, aptops, apbottoms, startpixel=300, endpixel=100, debug=False, apsize=6):
    """
    Function to extract a spectrum directly from pre-defined apertures
    This function does not try to adjust the apertures at all!
    :param img: The image to extract
    :param hdr: The header of that image
    :param aptops: The coefficients for the tops of each aperture
    :param apbottoms: The coefficients for the bottoms of each aperture
    :param startpixel: The first pixel to use (because the blaze function cuts off)
    :param endpixel: The last pixel to use
    :return: a list of extracted spectra
    """

    nx, ny = img.shape
    pixels = np.arange(nx)

    spectra = []
    ypos = []

    for apnum, (topcoeffs, bottomcoeffs) in enumerate(zip(aptops, apbottoms)):
        top = np.polynomial.polynomial.polyval(pixels,topcoeffs)
        bottom = np.polynomial.polynomial.polyval(pixels,bottomcoeffs)

        ypos.append((top + bottom)[startpixel:-endpixel]/2.0)

        bottom = (top + bottom)/2.0 - apsize/2.0
        strip = np.zeros([apsize,nx], dtype=np.double)
        yy, xx = np.indices(strip.shape)
        ap_yy = np.array(yy, dtype=np.double) + bottom
        iyy = np.array(np.round(ap_yy), dtype=np.int)
        # calculate the fraction subtracted by integer pixel index
        fyy = ap_yy - iyy
        # find the valid points
        vv = np.where((iyy >= 0) & (iyy <= ny-2))
        # input the value into the strip
        #we use "round" 2013-09-03 meeting with SPak & Huynh Anh
        #and we don't scale or interpolate pixel intensities 2013-09-03 meeting with SPak & Huynh Anh
        #strip[yy[vv],xx[vv]] = \
        #   img[iyy[vv],xx[vv]] * (1.0 - fyy[yy[vv],xx[vv]]) + \
        #   img[(iyy[vv]+1),xx[vv]] * fyy[yy[vv],xx[vv]]
        strip[yy[vv],xx[vv]] = img[iyy[vv],xx[vv]]
        spectrum = np.mean(strip, axis=0)[startpixel:-endpixel]
        spectrum2 = np.median(strip, axis=0)[startpixel:-endpixel]

        """
        z1, z2 = PL_Display.zscale(strip)
        fig = plt.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.imshow(strip, vmin=z1, vmax=z2, aspect='auto')
        ax2.plot(np.arange(spectrum.size)+startpixel, spectrum, 'k-')
        ax2.plot(np.arange(spectrum2.size)+startpixel, spectrum2, 'b-')
        ax2.set_xlim((0, nx))
        ylim = ax2.get_ylim()
        ylim = [max(0, ylim[0]), ylim[1]]
        ax2.set_ylim(ylim)
        ax1.set_title("Aperture %i" %apnum)
        plt.show()
        """

        spectra.append(spectrum2)



    return spectra, ypos




def FitDispersion(spectra, ypos, band="H", outpath="./"):
    """
    This function will read in the spectra
    :param spectra: A list of arc lamp spectra for each order
    :param ypos: The y position on the chip, in the same format as spectra
    :return: fitting coefficients for a 2d polynomial
    """

    theoretical_wave = fits.open("mapping/IGRINS_%s_MAPPING_May.fits" %band)[0].data[0]

    th_wave, th_flux = np.loadtxt("ThArlines.dat", unpack=True)
    th_wave *= units.angstrom.to(units.micron)
    oh_wave, oh_flux = np.loadtxt("ohlines.dat", unpack=True)
    oh_wave *= units.angstrom.to(units.micron)

    def MakeSpec(x, wavelist, fluxlist):
        y = np.zeros(x.size)
        for w, f in zip(wavelist, fluxlist):
            i = np.argmin(np.abs(x - w))
            y[i] = f
        return y

    for i, (spec, y) in enumerate(zip(spectra, ypos)):
        np.savetxt("%sspec_ap%.3i.dat" %(outpath, i), spec)
        np.savetxt("%sypos_ap%.3i.dat" %(outpath, i), y)

        pixels = np.arange(spec.size, dtype=np.int)
        wave_est = theoretical_wave[np.round(y).astype(np.int),pixels+300]

        x = np.linspace(wave_est[0], wave_est[-1], spec.size)
        y1 = MakeSpec(x, th_wave, th_flux)
        y2 = MakeSpec(x, oh_wave, oh_flux)

        plt.plot(wave_est, spec/spec.max(), 'k-')
        plt.plot(x, y1/y1.max(), 'r--')
        plt.plot(x, y2/y2.max(), 'b--')
    plt.show()









