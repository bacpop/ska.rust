#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update({"font.size": 18})
import matplotlib.pyplot as plt

import argparse


def norm_series(series: list()):
    norm = max(series)
    series = [item / norm for item in series]
    return series


def main():
    parser = argparse.ArgumentParser(
        prog="plot_cov",
        description="Plot the results of `ska cov`",
        epilog="Requires matplotlib",
    )
    parser.add_argument("histfile", help="Input table (stdout from `ska cov`)")
    parser.add_argument("--output", help="Output prefix", default="coverage_histogram")
    args = parser.parse_args()

    cutoff = 0
    kmers = list()
    density = list()
    idx_xseries = list()
    with open(args.histfile, "r") as hfile:
        hfile.readline()
        for line in hfile:
            (idx, count, ll, comp) = line.rstrip().split("\t")
            kmers.append(int(count))
            density.append(float(ll))
            idx_xseries.append(int(idx))
            if comp == "Coverage" and cutoff == 0:
                cutoff = int(idx)

    k_norm = norm_series(kmers)
    # d_norm = norm_series(density)

    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle("Coverage histogram fit")
    fig.set_dpi(160)
    fig.set_facecolor("w")
    fig.set_edgecolor("k")
    fig.set_figwidth(11)
    fig.set_figheight(11)
    plt.tight_layout()

    ax1.set_xlabel("K-mer count")
    ax1.set_ylabel("Frequency")
    ax1.set_ylim(0, k_norm[1])
    ax1.plot(
        idx_xseries, k_norm, color="black", linewidth=2, label="K-mer count frequency"
    )
    ax1.plot(
        idx_xseries,
        density,
        color="red",
        linewidth=2,
        linestyle="--",
        label="Mixture model fit",
    )
    ax1.plot(
        [cutoff, cutoff],
        [0, 1],
        color="darkgray",
        linewidth=1,
        linestyle="-.",
        label=f"Count cutoff ({cutoff})",
    )
    ax1.legend(loc="upper right")

    ax2.set_yscale("log")
    ax2.set_xlabel("K-mer count")
    ax2.set_ylabel("log(Frequency)")
    ax2.set_ylim(min(k_norm), k_norm[1])
    ax2.plot(
        idx_xseries, k_norm, color="black", linewidth=2, label="K-mer count frequency"
    )
    ax2.plot(
        idx_xseries,
        density,
        color="red",
        linewidth=2,
        linestyle="--",
        label="Mixture model fit",
    )
    ax2.plot(
        [cutoff, cutoff],
        [0, 1],
        color="darkgray",
        linewidth=1,
        linestyle="-.",
        label=f"Count cutoff ({cutoff})",
    )

    plt.savefig(args.output + ".png", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
