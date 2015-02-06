def figlist_to_pdfs(rootname, figlist, postfixes=None):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from itertools import izip, count

    if postfixes is None:
        postfixes = ("fig%02d" % i for i in count(1))

    for postfix, fig in izip(postfixes, figlist):
        FigureCanvasAgg(fig)
        fig.savefig("%s_%s.pdf" % (rootname, postfix))
        #fig2.savefig("align_zemax_%s_fig2_fit.pdf" % postfix)
        #fig3.savefig("align_zemax_%s_fig3_hist_dlambda.pdf" % postfix)

def figlist_to_json(rootname, figlist, postfixes=None):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from itertools import izip, count
    from mpld3 import save_json

    if postfixes is None:
        postfixes = ("fig%02d" % i for i in count(1))

    for postfix, fig in izip(postfixes, figlist):
        FigureCanvasAgg(fig)
        save_json(fig, "%s_%s.json" % (rootname, postfix))
