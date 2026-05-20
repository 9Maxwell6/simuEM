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
    ap = argparse.ArgumentParser(description="Plot solve time vs h.")
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
 
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.set_prop_cycle(color=plt.cm.Set1.colors)
    for path in args.files:
        h, t_sol = np.loadtxt(path, comments="#", usecols=(0, 6), unpack=True)
        order = np.argsort(h)
        h, t_sol = h[order], t_sol[order]
 
        label = os.path.basename(path).replace("_", r"\_")
 
        if args.no_fit:
            ax.loglog(h, t_sol, "+-", markersize=7, linewidth=0.7, label=label)
        else:
            if args.last_n is not None and args.last_n < len(h):
                h_fit, t_fit = h[:args.last_n], t_sol[:args.last_n]
            else:
                h_fit, t_fit = h, t_sol
 
            slope, intercept = np.polyfit(np.log(h_fit), np.log(t_fit), 1)
            fit = np.exp(intercept) * h ** slope
 
            line, = ax.loglog(h, t_sol, "+-", markersize=7, linewidth=0.7,
                              label=f"{label}  ($\\sim h^{{{slope:.3f}}}$)")
            ax.loglog(h, fit, "--", linewidth=0.7, alpha=0.6, color=line.get_color())
 
            print(f"{label}: solve time slope = {slope:.4f}")
 
    ax.set_xlabel(r"$h$")
    ax.set_ylabel(r"$t_{\text{solve}}\ [\mathrm{s}]$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    ax.invert_xaxis()
 
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
        fig.savefig(args.output, dpi=150)
        print(f"saved figure to {args.output}")
    else:
        plt.show()
 
 
if __name__ == "__main__":
    main()
 