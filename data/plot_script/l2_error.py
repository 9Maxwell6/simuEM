import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
 
plt.rcParams['text.usetex'] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
 
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["DejaVu Serif"]
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["legend.fontsize"] = 9
 
def main():
    ap = argparse.ArgumentParser(description="Plot FEM convergence (h vs L2 error).")
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
    ax.set_prop_cycle(color=plt.cm.Set1.colors)


    color = ['palevioletred', 'royalblue']
    color = ['red']
    ax.set_prop_cycle(color=color)
    labels = ["$\mathbf{T}$-$\Omega$\_c solver"]

    '''
    labels = ["$-\mathrm{div}\,\mathbf{grad}\, \mathrm{u} = \mathrm{f}$",
              "$\mathbf{curl}\,\mathbf{curl}\, \mathrm{u} + \mathrm{u} = \mathrm{f}$"]
    #'''


    for n, path in enumerate(args.files):
        h, err = np.loadtxt(path, comments="#", usecols=(0, 2), unpack=True)
        order = np.argsort(h)
        h, err = h[order], err[order]
 
        #label = os.path.basename(path).replace("_", r"\_")
        label = labels[n]
 
        if args.no_fit:
            ax.loglog(h, err, "+-", markersize=7, markeredgewidth=0.4, linewidth=0.4, label=label)
        else:
            if args.last_n is not None and args.last_n < len(h):
                h_fit, err_fit = h[:args.last_n], err[:args.last_n]
            else:
                h_fit, err_fit = h, err
 
            slope, intercept = np.polyfit(np.log(h_fit), np.log(err_fit), 1)
            fit = np.exp(intercept) * h ** slope
 
            line, = ax.loglog(h, err, "+-", markersize=7, markeredgewidth=0.4, linewidth=0.4,
                              label=f"{label}  ($\\sim h^{{{slope:.3f}}}$)")
            ax.loglog(h, fit, "--", linewidth=0.7, alpha=0.6, color=line.get_color())
 
            print(f"{label}: convergence rate = {slope:.4f}")
 
    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"$\|u - u_h\|_{\text{\scriptsize 0}}$")
    #ax.grid(True, which="both", alpha=0.3)
    #ax.legend()

    leg = ax.legend(loc='best', facecolor='white', edgecolor='none', framealpha=1, fancybox=False)
    for line in leg.get_lines():
        line.set_linewidth(2.5)

    ax.invert_xaxis()

    ax.set_facecolor('#e6e6e6')
    ax.grid(True, which="both",color='white', linestyle='-', linewidth=0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
 
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
    fig.tight_layout()
 
    if args.output:
        fig.savefig(args.output, dpi=1000, bbox_inches='tight', pad_inches=0)
        print(f"saved figure to {args.output}")
    else:
        plt.show()
 
 
if __name__ == "__main__":
    main()
 