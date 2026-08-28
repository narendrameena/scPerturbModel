#!/usr/bin/env python3
"""Name the reproducible context-response components by pathway enrichment.

Takes the split-half-validated components from test_rank1_viability.py and
tests, for each component and each loading direction, whether its extreme genes
are enriched for annotated gene sets.

Statistics: one-sided hypergeometric test against the correct background — the
responsive-gene universe the components were fitted in, NOT the whole genome
(using the genome would report the responsive-gene signature itself as if it
were a property of the component). Benjamini-Hochberg FDR is applied across all
(component, direction, gene set) tests within a library.

Only components that reproduced on held-out plates are named; the rest are
noise and are reported as such.

Outputs: results/tables/component_enrichment_<tag>.csv
         figure bundle results/figures/06_diagnostics/component_names_<tag>/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
GS = ROOT / "data" / "external" / "genesets"
FIG = ROOT / "results" / "figures" / "06_diagnostics"
TAB = ROOT / "results" / "tables"
TOP_N = 150          # genes per tail
MIN_OVERLAP = 3
MIN_SET, MAX_SET = 5, 500
BLUE, ORANGE = "#2a78d6", "#eb6834"


def read_gmt(path: Path, universe: set) -> dict:
    sets = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip().split("\t")
        if len(parts) < 3:
            continue
        genes = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
        genes &= universe
        if MIN_SET <= len(genes) <= MAX_SET:
            sets[parts[0]] = genes
    return sets


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - rank))
        q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="dev47")
    ap.add_argument("--libraries", nargs="*",
                    default=["MSigDB_Hallmark_2020", "Reactome_2022",
                             "GO_Biological_Process_2023"])
    args = ap.parse_args()

    load = pd.read_csv(TAB / f"rank1_component_loadings_{args.tag}.csv")
    comp_stats = pd.read_csv(TAB / f"rank1_viability_{args.tag}.csv")
    load["gene_symbol"] = load.gene_symbol.str.upper()
    load = load[~load.gene_symbol.str.startswith("ENSG")]
    universe = set(load.gene_symbol)
    good = comp_stats[comp_stats.split_half_r > 0.2].component.tolist()
    print(f"{len(universe)} background genes; reproducible components: {good}")

    recs = []
    for lib in args.libraries:
        path = GS / f"{lib}.gmt"
        if not path.exists():
            print(f"missing {lib}, skipping")
            continue
        sets = read_gmt(path, universe)
        M = len(universe)
        rows = []
        for c in good:
            col = f"comp{c}"
            order = load[col].to_numpy().argsort()
            tails = {"+": load.gene_symbol.to_numpy()[order[::-1][:TOP_N]],
                     "-": load.gene_symbol.to_numpy()[order[:TOP_N]]}
            for direction, gl in tails.items():
                q_set = set(gl)
                n = len(q_set)
                for name, gs in sets.items():
                    k = len(q_set & gs)
                    if k < MIN_OVERLAP:
                        continue
                    p = hypergeom.sf(k - 1, M, len(gs), n)
                    rows.append({"library": lib, "component": c,
                                 "direction": direction, "gene_set": name,
                                 "overlap": k, "set_size": len(gs),
                                 "p": p,
                                 "genes": ",".join(sorted(q_set & gs)[:8])})
        if not rows:
            continue
        df = pd.DataFrame(rows)
        df["q"] = bh(df.p.to_numpy())
        recs.append(df)
        print(f"{lib}: {len(df)} tests, {(df.q < 0.05).sum()} at FDR<0.05")

    if not recs:
        print("no enrichment results")
        return
    res = pd.concat(recs, ignore_index=True).sort_values(["component", "p"])
    res.to_csv(TAB / f"component_enrichment_{args.tag}.csv", index=False)

    print("\n=== component names (top enrichment per component/direction) ===")
    summary = []
    for c in good:
        sub = res[(res.component == c) & (res.q < 0.1)]
        line = [f"\ncomponent {c} "
                f"(split-half r={comp_stats.loc[comp_stats.component == c, 'split_half_r'].iloc[0]:.2f}, "
                f"{comp_stats.loc[comp_stats.component == c, 'share_of_reproducible'].iloc[0]:.0%} of "
                f"reproducible variance)"]
        for direction in ("+", "-"):
            d = sub[sub.direction == direction].nsmallest(3, "p")
            if len(d):
                for r in d.itertuples():
                    line.append(f"   {direction} {r.gene_set[:62]:64s} "
                                f"q={r.q:.1e} ({r.overlap}/{r.set_size})")
                    summary.append({"component": c, "direction": direction,
                                    "gene_set": r.gene_set, "q": r.q,
                                    "overlap": r.overlap})
            else:
                line.append(f"   {direction} (no enrichment at FDR<0.10)")
        print("\n".join(line))
    summ = pd.DataFrame(summary)

    # ---------------- figure ----------------
    if len(summ):
        plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                             "axes.spines.right": False, "axes.grid": True,
                             "grid.alpha": 0.25, "figure.facecolor": "white"})
        best = (summ.sort_values("q").groupby(["component", "direction"])
                .head(2).sort_values(["component", "q"]))
        h = max(4.0, 0.34 * len(best))
        fig, ax = plt.subplots(figsize=(11, h), constrained_layout=True)
        y = np.arange(len(best))[::-1]
        cols = [BLUE if d == "+" else ORANGE for d in best.direction]
        ax.barh(y, -np.log10(best.q.to_numpy()), color=cols, height=0.62)
        ax.set_yticks(y, [f"C{r.component}{r.direction}  {r.gene_set[:52]}"
                          for r in best.itertuples()], fontsize=7.5)
        ax.axvline(-np.log10(0.05), color="#888888", ls="--", lw=0.9)
        ax.set_xlabel("−log10 FDR q")
        ax.set_title("Reproducible context-response components, named by "
                     "pathway enrichment", loc="left", fontweight="bold",
                     fontsize=10)
        ax.text(0.99, 0.01, "blue = positive loading, orange = negative",
                transform=ax.transAxes, ha="right", fontsize=7.5,
                color="#666666")
        d = save_figure(fig, f"component_names_{args.tag}", FIG,
                        source_data={"enrichment": res, "summary": summ},
                        script=__file__)
        print(f"\nfigure bundle -> {d}")


if __name__ == "__main__":
    main()
