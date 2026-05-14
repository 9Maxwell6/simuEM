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
    ap.add_argument("files", nargs="+", help="one or more .dat files with columns: h, L2_error")
    ap.add_argument("-o", "--output", help="save figure to this file instead of showing it")
    ap.add_argument("-n", "--last-n", type=int, default=None, help="use only the last N refinements (smallest h) for the fit; "
                                                                   "default: use all points")
    args = ap.parse_args()
 
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.set_prop_cycle(color=plt.cm.Set1.colors)
    for path in args.files:
        h, err = np.loadtxt(path, comments="#", unpack=True)
        order = np.argsort(h)
        h, err = h[order], err[order]

        # subset used for the fit: the N smallest h values
        if args.last_n is not None and args.last_n < len(h):
            h_fit, err_fit = h[:args.last_n], err[:args.last_n]
        else:
            h_fit, err_fit = h, err
 
        # least-squares fit in log-log space
        slope, intercept = np.polyfit(np.log(h_fit), np.log(err_fit), 1)
        fit = np.exp(intercept) * h ** slope
 
        label = os.path.basename(path).replace("_", r"\_")
        line, = ax.loglog(h, err, "+-", markersize=7, linewidth=0.7, label=f"{label}  ($\sim h^{{{slope:.3f}}}$)")
        ax.loglog(h, fit, "--", linewidth=0.7, alpha=0.6, color=line.get_color())
 
        print(f"{label}: convergence rate = {slope:.4f}")
    
    
    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"$\|u - u_h\|_{\text{\scriptsize 0}}$")
    #ax.set_title(r"convergence rate")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    ax.invert_xaxis() 

    def sci_label(v):
        exp = int(np.floor(np.log10(v)))
        mant = round(v / 10**exp, 10)          # avoid 4.999... artifacts
        return rf"${mant:g} \times 10^{{{exp}}}$"

    ticks  = [0.5, 0.2, 0.1, 0.05]
    labels = [sci_label(t) for t in ticks]

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    fig.tight_layout()
 
    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"saved figure to {args.output}")
    else:
        plt.show()
 
 
if __name__ == "__main__":
    main()