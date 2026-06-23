import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator, NullFormatter

 
plt.rcParams['text.usetex'] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}\usepackage{xcolor}"

 
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["DejaVu Serif"]
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["legend.fontsize"] = 9

def pad_yticks_by_digits(ax, n=1):
    """Add n digit-widths of left pad to the y-tick labels.
    Survives bbox_inches='tight' because pad moves real ink (unlike \\phantom).
    Call AFTER set_yticklabels, on the single-digit plot only."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fs = ax.get_yticklabels()[0].get_fontsize()
    probe = ax.text(0, 0, "$" + "0" * n + "$", fontsize=fs)   # measure n digits
    fig.canvas.draw()
    extra_pt = probe.get_window_extent(r).width * 72.0 / fig.dpi
    probe.remove()
    base = ax.yaxis.majorTicks[0].get_pad()
    ax.tick_params(axis="y", pad=base + extra_pt)
 
def main():
    ap = argparse.ArgumentParser(description="Plot condition number vs h.")
    ap.add_argument("files", nargs="+",
                    help="one or more .dat files; columns: h, #cell, L2_error, #iter, cond, assemble, solve")
    ap.add_argument("-o", "--output", help="save figure to this file instead of showing it")
    ap.add_argument("-n", "--last-n", type=int, default=None,
                    help="use only the last N refinements (smallest h) for the fit; "
                         "default: use all points")
    ap.add_argument("--no-fit", action="store_true",
                    help="don't draw the dashed power-law fit or show the slope in the legend")
    ap.add_argument("--linear", action="store_true",
                    help="use linear axes instead of log-log (default: log-log)")
    args = ap.parse_args()
 
    fig, ax = plt.subplots(figsize=(5, 4))

    colors = list(plt.cm.Set1.colors)
    colors[5] = (0.4, 0.4, 0.4)          # replace yellow with gray
    ax.set_prop_cycle(color=colors)

    cset = np.array(colors[1])  
    n = 7  
    factors = np.linspace(1.0, 0.01, n)   
    cset_shades = [tuple(cset * f) for f in factors]
    cset_shades =  ['#7fb9da', '#5da5d1', '#3f8fc5', '#2777b8', '#135fa7', '#08488e', '#08306b']
    #cset_shades = ['#86cc85', '#63bc6e', '#3fa85b', '#29914a', '#107a37', '#006227', '#00441b']

    blue = ['#9ecae1', '#6baed6', '#3182bd', '#08519c', '#08306b', '#011b3f']
    green = ['#a1d99b', '#74c476', '#31a354', '#006d2c', '#00441b', '#001f0f']
    purple = ['#bcbddc', '#9e9ac8', '#756bb1', '#54278f', '#3f007d', '#1f003f']
    orange = ['#fdbf7a', '#fdae6b', '#fd8d3c', '#d94801', '#8c2d04', '#3f1202']

    color = ['royalblue', 'navy']
    #color = ['palevioletred', 'crimson']
    color = ['royalblue', 'navy', colors[0], colors[3]]

    ax.set_prop_cycle(color=color)


    #'''
    labels = ["test case 2",
              "test case 2\_pc",
              "$\mathbf{T}$-$\Omega$",
              "$\mathbf{T}$-$\Omega$\_pc"]
    #'''

    '''
    labels = ["CG",
              "CG\_pc"]
    #'''

    '''
    labels = ["CG",
              "CG\_pc\_1",
              "CG\_pc\_2",
              "CG\_pc\_3",
              "CG\_pc\_4"]
    #'''

    '''
    labels = ["CG\_pc\_3 ($N=1$)",
              "CG\_pc\_3 ($N=2$)",
              "CG\_pc\_3 ($N=5$)",
              "CG\_pc\_3 ($N=10$)",
              "CG\_pc\_3 ($N=15$)",
              "CG\_pc\_3 ($N=20$)"]
    #'''

    for n, path in enumerate(args.files):
        h, kappa = np.loadtxt(path, comments="#", usecols=(0, 4), unpack=True)
        order = np.argsort(h)
        h, kappa = h[order], kappa[order]
 
        #label = os.path.basename(path).replace("_", r"\_")
        label = labels[n]
 
        if args.no_fit:
            ax.loglog(h, kappa, "+-", markersize=7, markeredgewidth=0.4, linewidth=0.4, label=label)
        else:
            if args.last_n is not None and args.last_n < len(h):
                h_fit, k_fit = h[:args.last_n], kappa[:args.last_n]
            else:
                h_fit, k_fit = h, kappa
 
            slope, intercept = np.polyfit(np.log(h_fit), np.log(k_fit), 1)
            fit = np.exp(intercept) * h ** slope
 
            line, = ax.loglog(h, kappa, "+-", markersize=7, markeredgewidth=0.4, linewidth=0.4,
                              label=f"{label}  ($\\sim h^{{{slope:.3f}}}$)")
            ax.loglog(h, fit, "--", linewidth=0.4, alpha=0.6, color=line.get_color())
 
            print(f"{label}: condition number slope = {slope:.4f}")
 
    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"condition number $\kappa$")
    leg = ax.legend(loc='best', facecolor='white', edgecolor='none', framealpha=1, fancybox=False)
    for line in leg.get_lines():
        line.set_linewidth(2.5)
    ax.invert_xaxis()

    ax.set_facecolor('#e6e6e6')

    #fmt = mticker.FuncFormatter(lambda y, _: rf'${y:g}^{{1}}$')
    #ax.yaxis.set_minor_formatter(fmt)
    #ax.yaxis.set_major_formatter(fmt)

    #ax.grid(True, which="both", color='black', linewidth=0.4, alpha=0.2)
    ax.grid(True, which="both",color='white', linestyle='-', linewidth=0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.tick_params(axis="y", pad=7.6)


 
    if args.linear:
        ax.set_xscale("linear")
        ax.set_yscale("linear")
    else:
        def sci_label(v):
            exp = int(np.floor(np.log10(v)))
            mant = round(v / 10**exp, 10)
            return rf"${mant:g} \times 10^{{{exp}}}$"
        
        ticks  = [0.5, 0.2, 0.1, 0.05]
        labels = [sci_label(t) for t in ticks]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        #yticks = [2, 3, 4, 5]
        #ylabel = [rf"$\phantom{{0}}{v}^{{1}}$" for v in yticks]  
        #pad_yticks_by_digits(ax, n=1) 
        
        #ax.set_yticks(yticks)
        #ax.set_yticklabels(ylabel)
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout()
 
    if args.output:
        fig.savefig(args.output, dpi=1000, bbox_inches='tight', pad_inches=0)
        print(f"saved figure to {args.output}")
    else:
        plt.show()
 
 
if __name__ == "__main__":
    main()
 