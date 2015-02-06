import numpy as np

from stsci_helper import stsci_median

class Destriper(object):
    def __init__(self):
        self.dy = dy = 128
        self.n_dy = 2048/dy
        self.dy_slices = [slice(iy*dy, (iy+1)*dy) for iy in range(2048/dy)]

    def _remove_vertical_smooth_bg(self, d, mask=None):
        ny, nx = d.shape
        iy = np.arange(ny)
        import scipy.ndimage as ni
        md = ni.median_filter(d, [1,7])
        mask2 = mask & np.isfinite(md)
        p_list = [np.polyfit(iy[~mask2[:,ix]],
                             md[:,ix][~mask2[:,ix]], 4) for ix in range(nx)]
        dd = np.vstack([np.polyval(p, iy) for p in p_list])
        return d - dd.T

    def get_stripe_pattern64(self, d, mask=None, concatenate=True):
        """
        if concatenate is True, return 2048x2048 array.
        Otherwise, 128x2048 array.
        """
        dy = 64
        n_dy = 2048/dy
        dy_slices = [slice(iy*dy, (iy+1)*dy) for iy in range(2048/dy)]
        from itertools import cycle
        if mask is not None:
            d = self._remove_vertical_smooth_bg(d, mask=mask)
            alt_sign = cycle([1, -1])
            dd = [d[sl][::next(alt_sign)] for sl in dy_slices]
            alt_sign = cycle([1, -1])
            msk = [mask[sl][::next(alt_sign)] for sl in dy_slices]
            ddm = stsci_median(dd, badmasks=msk)
            #dd1 = np.ma.array(dd, mask=msk)
            #ddm = np.ma.median(dd1, axis=0)
        else:
            alt_sign = cycle([1, -1])
            dd = [d[sl][::next(alt_sign)] for sl in dy_slices]
            ddm = np.median(dd, axis=0)


        if concatenate:
            return np.concatenate([ddm, ddm[::-1]] * (n_dy/2))
        else:
            return ddm

    def get_stripe_pattern128(self, d, mask=None, concatenate=True):
        """
        if concatenate is True, return 2048x2048 array.
        Otherwise, 128x2048 array.
        """
        dy_slices=self.dy_slices
        if mask is not None:
            dd = [d[sl] for sl in dy_slices]
            msk = [mask[sl] for sl in dy_slices]
            dd1 = np.ma.array(dd, mask=msk)
            ddm = np.ma.median(dd1, axis=0)
        else:
            dd = [d[sl] for sl in dy_slices]
            ddm = np.median(dd, axis=0)


        if concatenate:
            return np.concatenate([ddm] * self.n_dy)
        else:
            return ddm

    def get_stripe_pattern128_flat(self, d, mask=None, concatenate=True):
        """
        if concatenate is True, return 2048x2048 array.
        Otherwise, 128x2048 array.
        """
        dy_slices=self.dy_slices
        if mask is not None:
            dd = [d[sl] for sl in dy_slices]
            msk = [mask[sl] for sl in dy_slices]
            dd1 = np.ma.array(dd, mask=msk)
            dd2 = dd1.transpose([1, 0, 2]).reshape((self.dy, -1))
            ddm = np.ma.median(dd2, axis=1)
        else:
            dd = np.array([d[sl] for sl in dy_slices])
            dd2 = dd.transpose([1, 0, 2]).reshape((self.dy, -1))
            ddm = np.median(dd2, axis=1)


        if concatenate:
            return np.concatenate([ddm] * self.n_dy)
        else:
            return ddm


    def get_stripe_pattern2048(self, d, mask=None):
        """
        if concatenate is True, return 2048x2048 array.
        Otherwise, 128x2048 array.
        """
        if mask is not None:
            dd1 = np.ma.array(d, mask=mask)
            ddm = np.ma.median(dd1, axis=1)
        else:
            ddm = np.median(d, axis=1)

        return ddm

    def get_destriped(self, d, mask=None, hori=None, pattern=128):
        # if hori:
        #     s_hori = np.median(d, axis=0)
        #     d = d - s_hori
        if hori is None:
            if pattern == 2048:
                hori = True

        if pattern == 64:
            ddm = self.get_stripe_pattern64(d, mask=mask, concatenate=True)
            d_ddm = d - ddm
        elif pattern == 128:
            ddm = self.get_stripe_pattern128(d, mask=mask, concatenate=True)
            d_ddm = d - ddm
        elif pattern == 2048:
            ddm = self.get_stripe_pattern2048(d, mask=mask)
            #ddm = self.get_stripe_pattern128_flat(d, mask=mask)
            d_ddm = d - ddm[:,np.newaxis]

        if hori:
            d_ddm_masked = np.ma.array(d_ddm, mask=mask)
            s_hori = np.ma.median(d_ddm_masked, axis=0)
            d_ddm = d_ddm - s_hori


        return d_ddm

    def get_destriped_naive(self, d):
        """
        if concatenate is True, return 2048x2048 array.
        Otherwise, 128x2048 array.
        """
        s_vert = np.median(d, axis=0)
        d_vert = d - s_vert[ np.newaxis,:]

        s_hori = np.median(d_vert, axis=1)
        d_hori = d_vert - s_hori[:,np.newaxis]

        return d_hori

destriper = Destriper()

    # if 0:
    #     destriper = Destriper()

    #     hdu_list = igrins_log.get_cal_hdus(band, "flat_off")
    #     hdu = hdu_list[0]

    #     from mpl_toolkits.axes_grid1 import ImageGrid
    #     fig = plt.figure(figsize=(13, 5))
    #     grid = ImageGrid(fig, 111, (1, 3),
    #                      share_all=True,
    #                      label_mode="1",axes_pad=0.1)
    #     from matplotlib.colors import Normalize
    #     norm = Normalize(vmin=-15, vmax=30)
    #     im = grid[0].imshow(hdu.data, norm=norm,
    #                         origin="lower", cmap="gray_r",
    #                         interpolation="none")
    #     hh1 = destriper.get_destriped(hdu.data, pattern=2048)
    #     im2 = grid[1].imshow(hh1, norm=norm,
    #                          origin="lower", cmap="gray_r",
    #                          interpolation="none")

    #     hh2 = destriper.get_destriped(hdu.data, pattern=64)
    #     im3 = grid[2].imshow(hh2, norm=norm,
    #                          origin="lower", cmap="gray_r",
    #                          interpolation="none")

    #     grid[0].set_xlabel("x-pixel")
    #     grid[0].set_ylabel("y-pixel")
    #     fig.tight_layout()



    # if 0:
    #     destriper = Destriper()

    #     hdu_list = igrins_log.get_cal_hdus(band, "flat_off")
    #     hdu = hdu_list[0]


    #     hh1 = destriper.get_destriped(hdu.data, mask=bias_mask,
    #                                   pattern=2048)
    #     hh2 = destriper.get_destriped(hdu.data,
    #                                   pattern=64, mask=bias_mask)

    #     from matplotlib.colors import Normalize

    #     from mpl_toolkits.axes_grid1 import ImageGrid
    #     fig = plt.figure(figsize=(13, 5))
    #     grid = ImageGrid(fig, 111, (1, 3),
    #                      share_all=True,
    #                      label_mode="1",axes_pad=0.1,
    #                      cbar_mode="single")

    #     norm = Normalize(vmin=-15, vmax=30)
    #     for ax, dd in zip(grid, [hdu.data, hh1, hh2]):
    #         im = ax.imshow(dd, norm=norm,
    #                        origin="lower", cmap="gray_r",
    #                        interpolation="none")
    #     plt.colorbar(im, cax=grid[0].cax)

    #     grid[0].set_xlabel("x-pixel")
    #     grid[0].set_ylabel("y-pixel")
    #     fig.tight_layout()

    #     fig = plt.figure()
    #     grid = Grid(fig, 111, (3, 1),
    #                 share_x=True, share_y=True, share_all=True,
    #                 label_mode="1",axes_pad=0.1)
    #     for ax, dd in zip(grid, [hdu.data, hh1, hh2]):
    #         s = np.median(dd[:,1024-128:1024+128], axis=1)
    #         ax.plot(s)
    #         ax.axhline(0, color="0.5")
    #     grid[2].set_xlim(0, 2048)



def check_readout_pattern(fig, hdu_list, axis=1):

    from mpl_toolkits.axes_grid1 import Grid
    grid = Grid(fig, 111, (len(hdu_list), 1),
                share_x=True, share_y=True, share_all=True,
                label_mode="1")

    smin, smax = np.inf, -np.inf
    for hdu, ax in zip(hdu_list, grid):
        p_horizontal = np.median(hdu.data, axis=axis)
        p_horizontal0 = np.median(p_horizontal)
        p_horizontal -= p_horizontal0
        smin = min(smin, p_horizontal[100:-100].min())
        smax = max(smax, p_horizontal[100:-100].max())
        ax.plot(p_horizontal)
        ax.axhline(0, color="0.5")
        ax.locator_params(axis="y", nbins=5)

    ax = grid[-1]
    ax.set_xlim(0, len(p_horizontal))
    sminmax = max(-smin, smax)
    ax.set_ylim(-sminmax, sminmax)
    if axis == 1:
        ax.set_xlabel("y-pixel")
    elif axis == 0:
        ax.set_xlabel("x-pixel")

    ax.set_ylabel("ADUs")


def check_destriper(hdu, bias_mask, vmin, vmax):
    destriper = Destriper()



    hh1 = destriper.get_destriped(hdu.data, mask=bias_mask,
                                  pattern=2048)
    hh2 = destriper.get_destriped(hdu.data,
                                  pattern=64, mask=bias_mask)

    from matplotlib.colors import Normalize

    from mpl_toolkits.axes_grid1 import ImageGrid, Grid
    fig1 = plt.figure(figsize=(13, 5))
    grid = ImageGrid(fig1, 111, (1, 3),
                     share_all=True,
                     label_mode="1",axes_pad=0.1,
                     cbar_mode="single")

    norm = Normalize(vmin=vmin, vmax=vmax)
    for ax, dd in zip(grid, [hdu.data, hh1, hh2]):
        im = ax.imshow(dd, norm=norm,
                       origin="lower", cmap="gray_r",
                       interpolation="none")
    plt.colorbar(im, cax=grid[0].cax)

    grid[0].set_xlabel("x-pixel")
    grid[0].set_ylabel("y-pixel")
    fig1.tight_layout()

    fig2 = plt.figure()
    grid = Grid(fig2, 111, (3, 1),
                share_x=True, share_y=True, share_all=True,
                label_mode="1",axes_pad=0.1)
    for ax, dd in zip(grid, [hdu.data, hh1, hh2]):
        s = np.median(dd[:,1024-128:1024+128], axis=1)
        ax.plot(s)
        ax.axhline(0, color="0.5")
    grid[2].set_xlim(0, 2048)
    grid[2].set_xlabel("y-pixel")
    fig2.tight_layout()

    return fig1, fig2


if __name__ == "__main__":
    from igrins_log import IGRINSLog

    log_20140525 = dict(flat_off=range(64, 74),
                        flat_on=range(74, 84),
                        thar=range(3, 8))


    igrins_log = IGRINSLog("20140525", log_20140525)

    band = "K"


    if 1: # check ro pattern
        hdu_list = igrins_log.get_cal_hdus(band, "flat_off")

        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(9, 8))
        check_readout_pattern(fig, hdu_list, axis=1)
        fig.tight_layout()

        fig2 = plt.figure(figsize=(9, 8))
        check_readout_pattern(fig2, hdu_list, axis=0)
        fig2.tight_layout()

        fig.savefig("readout_horizontal_pattern_%s.pdf" % band)
        fig2.savefig("readout_vertical_pattern_%s.pdf" % band)

    if 1: # check destriper
        # with mask
        import pickle
        r = pickle.load(open("flat_info_20140316_%s.pickle" % band))

        from trace_flat import get_mask_bg_pattern
        bias_mask = get_mask_bg_pattern(r["flat_mask"],
                                        r["bottomup_solutions"])

        if 1:
            fig = plt.figure()
            ax=fig.add_subplot(111)
            ax.imshow(bias_mask, origin="lower", cmap="gray_r",
                       vmin=0, vmax=1)
            fig.savefig("destriper_mask_%s.pdf"%band, bbox_inches="tight")

        hdu_list = igrins_log.get_cal_hdus(band, "flat_off")
        hdu = hdu_list[0]

        fig11, fig12 = check_destriper(hdu, None, vmin=-25, vmax=30)
        fig13, fig14 = check_destriper(hdu, bias_mask, vmin=-25, vmax=30)

        fig11.savefig("destriper_test_image_flat_off_wo_mask_%s.pdf" % band)
        fig12.savefig("destriper_test_cut_flat_off_wo_mask_%s.pdf" % band)
        fig13.savefig("destriper_test_image_flat_off_w_mask_%s.pdf" % band)
        fig14.savefig("destriper_test_cut_flat_off_w_mask_%s.pdf" % band)


        hdu_list = igrins_log.get_cal_hdus(band, "thar")
        hdu = hdu_list[-1]

        fig15, fig16 = check_destriper(hdu, bias_mask, vmin=-30, vmax=80)

        fig15.savefig("destriper_test_image_thar_w_mask_%s.pdf" % band)
        fig16.savefig("destriper_test_cut_thar_w_mask_%s.pdf" % band)
