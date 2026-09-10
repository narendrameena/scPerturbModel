"""``perturbdesign`` -- the design calculation as a command anyone can run.

The calculation in :mod:`perturbmodel.design` is only useful to an atlas builder
if it can be run before the atlas exists. This is that interface.

Three subcommands:

``plan``    what a proposed design can detect, and what to change if it cannot.
``audit``   the same for a dataset that already exists, read from its metadata.
``budget``  given a fixed number of cells, where they should go.

The third is the one that changes decisions. Sequencing budget is the real
constraint on an atlas, and the calculation says something counter-intuitive
about how to spend it: cells enter the interaction estimate only through the
per-condition noise, which saturates, while replicates enter through the pair
count, which does not. ``budget`` maximises detectability over that trade-off
rather than assuming an answer.

Examples
--------
    perturbdesign plan --contexts 50 --perturbations 1100 --replicates 2
    perturbdesign --target 0.02 plan --contexts 50 --perturbations 1100
        (--target, --snr and --model-dim are GLOBAL and precede the subcommand)
    perturbdesign audit atlas.h5ad
    perturbdesign budget --cells 100000000 --contexts 50 --perturbations 1100
"""
from __future__ import annotations

import argparse
import json
import sys

from .design import (captured_fraction, interaction_se, min_detectable_share,
                     n_pairs, required_replicates)

SNR = 0.20   # the shared constant validated in RESULTS.md sec.46


def _fmt(x, nd=4):
    return "n/a" if x is None else (f"{x:.{nd}f}" if isinstance(x, float)
                                    else str(x))


def cmd_plan(a):
    # Only REPLICATED conditions contribute pairs. This is the single most
    # common way to overstate a design: Tahoe-100M profiles ~1100 compounds but
    # replicates 13.5% of its conditions, and it is the replicated ones the
    # estimator sees. Quoting the total inflates the pair count ~7x and turns an
    # underpowered atlas into an apparently adequate one.
    frac = getattr(a, "replicated_fraction", 1.0)
    eff = max(int(round(a.perturbations * frac)), 0)
    p = n_pairs(a.contexts, eff, a.replicates)
    mds = min_detectable_share(a.contexts, eff, a.replicates, a.snr)
    out = {"n_contexts": a.contexts, "n_perturbations": a.perturbations,
           "replicated_fraction": frac, "n_perturbations_replicated": eff,
           "n_replicates": a.replicates, "n_pairs": p,
           "min_detectable_share": mds,
           "se": interaction_se(a.contexts, eff, a.replicates, a.snr)}
    if a.json:
        need = required_replicates(a.target, a.contexts, eff, a.snr)
        out["target"] = a.target
        out["resolves_target"] = bool(mds <= a.target)
        out["replicates_for_target"] = need
        print(json.dumps(out, indent=2))
        return 0

    print(f"\n  {a.contexts} contexts x {a.perturbations} perturbations x "
          f"{a.replicates} replicates")
    if frac < 1.0:
        print(f"  {frac:.1%} of conditions replicated -> {eff} perturbations "
              f"actually contribute")
    print(f"  {p:,} independent cross-replicate pairs\n")
    if p == 0:
        print("  UNRESOLVABLE AT ANY EFFECT SIZE.")
        print("  With one replicate per condition there is no pair to correlate,")
        print("  so the interaction cannot be separated from noise however many")
        print("  cells are sequenced. This is the binding constraint; nothing")
        print("  else in the design matters until it is fixed.\n")
        print(f"  Two replicates would give "
              f"{n_pairs(a.contexts, eff, 2):,} pairs and a "
              f"detection limit of "
              f"{min_detectable_share(a.contexts, eff, 2, a.snr):.4f}.")
        return 1

    print(f"  Smallest detectable interaction share:  {mds:.4f}  "
          f"({mds:.2%} of reproducible variance)")
    print(f"  Effect you want to detect:              {a.target:.4f}")
    if mds <= a.target:
        print(f"\n  ADEQUATE. This design resolves that effect with "
              f"{a.target / mds:.1f}x margin.")
    else:
        need = required_replicates(a.target, a.contexts, eff, a.snr)
        print(f"\n  UNDERPOWERED by {mds / a.target:.1f}x.")
        if need:
            print(f"  {need} replicates per condition would reach it "
                  f"(you have {a.replicates}).")
            extra = (need - a.replicates) * a.contexts * eff
            print(f"  That is {extra:,} additional treated wells.")
        else:
            print("  No practical replicate count reaches it; that effect needs")
            print("  more contexts or perturbations as well.")
    print(f"\n  For scale: the pooled interaction share measured in Tahoe-100M "
          f"is 0.005\n  (sec.31) and in LINCS 0.57 (sec.34). A target of 0.05 "
          f"is generous for a\n  drug atlas and strict for a cytokine panel; "
          f"set --target from your own\n  pilot rather than from this default.")
    cap = captured_fraction(a.model_dim, a.contexts)
    print(f"\n  A model with a {a.model_dim}-dimensional context embedding "
          f"captures {cap:.0%} of the")
    print(f"  interaction at {a.contexts} contexts. Because the interaction's "
          f"dimensionality")
    print(f"  grows with contexts measured (sec.43), that same model would "
          f"capture")
    for nc in (100, 200, 500):
        if nc > a.contexts:
            print(f"    {captured_fraction(a.model_dim, nc):>4.0%} at {nc} "
                  f"contexts")
    return 0


def cmd_budget(a):
    """Where a fixed cell budget should go.

    Cells per condition and replicates trade against each other under a fixed
    total. More cells per condition lowers that condition's noise -- entering
    through ``snr`` with diminishing return, since biological variation does not
    average away -- while more replicates add pairs linearly. The optimum is
    therefore interior, and usually at fewer cells per condition than atlases
    actually use.
    """
    n_cond = a.contexts * a.perturbations
    print(f"\n  Budget {a.cells:,} cells over {n_cond:,} conditions "
          f"({a.contexts} x {a.perturbations})\n")
    rows, best = [], None
    for r in range(1, a.max_replicates + 1):
        per = int(a.cells / (n_cond * r))
        if per < a.min_cells:
            break
        # Noise falls with depth until biological variation floors it. The two
        # parameters of this saturating curve are ASSUMPTIONS, not measurements
        # -- see the note printed below.
        snr = a.snr_max / (1.0 + a.snr_half / max(per, 1))
        mds = min_detectable_share(a.contexts, a.perturbations, r, snr)
        rows.append((r, per, n_pairs(a.contexts, a.perturbations, r), snr, mds))
        if best is None or mds < best[4]:
            best = rows[-1]
    print(f"  {'reps':>5} {'cells/cond':>11} {'pairs':>12} "
          f"{'snr':>7} {'detects':>9}")
    print("  " + "-" * 50)
    for row in rows:
        print(f"  {row[0]:>5} {row[1]:>11,} {row[2]:>12,} {row[3]:>7.3f} "
              f"{row[4]:>9.4f}" + ("  <-- best" if row is best else ""))
    if best:
        r, per, _, _, mds = best
        print(f"\n  BEST: {r} replicates x {per:,} cells per condition, "
              f"detecting {mds:.4f}.")
        one = int(a.cells / n_cond)
        print(f"  The same budget as 1 replicate of {one:,} cells detects "
              f"NOTHING at any\n  effect size: the pair count is zero, so depth "
              f"cannot substitute for\n  replication. That is the qualitative "
              f"result, and it does not depend on\n  the curve below — one "
              f"replicate yields no pairs however it is modelled.")
        knee = next((x for x in rows if x[4] <= mds * 1.15), None)
        if knee and knee[0] < r:
            print(f"\n  Returns flatten early: {knee[0]} replicates already "
                  f"reaches {knee[4]:.4f}, within\n  15% of the best. Past that "
                  f"the budget is better spent on more contexts.")
    print(f"""
  ASSUMPTIONS in the snr column, which you should set from your own pilot:
    --snr-max  {a.snr_max}   signal-to-noise at infinite depth (biological floor)
    --snr-half {a.snr_half:g}   cells per condition giving half of that
  These two shape WHERE the optimum sits, not whether one exists. Swept over
  snr-half from 50 to 2000 the optimum is never 1 replicate, and is 3 or more
  unless noise saturates unusually slowly (snr-half > 1000); swept over snr-max
  from 0.1 to 0.8 it stays between 7 and 9. The direction of the advice is
  stable even though the exact number is not.""")
    return 0


def cmd_audit(a):
    try:
        from .atlas_meta import describe
    except ImportError:
        print("audit needs h5py and pandas installed", file=sys.stderr)
        return 2
    d = describe(a.path)
    if "error" in d:
        print(f"could not read {a.path}: {d['error']}", file=sys.stderr)
        return 1
    print(f"\n  {d['dataset']}")
    print(f"  {d['n_cells']:,} cells")
    print(f"  context      <- obs['{d['context_col']}']  "
          f"({d['n_contexts']} levels)")
    print(f"  perturbation <- obs['{d['pert_col']}']  "
          f"({d['n_perturbations']} levels)")
    if d["rep_col"]:
        print(f"  replicate    <- obs['{d['rep_col']}']  "
              f"({d['n_replicates']} per condition, median)")
    else:
        print("  replicate    <- NONE FOUND. Without an annotated replicate the "
              "interaction\n                  cannot be estimated at all; see "
              "the verdict below.")
    a.contexts, a.perturbations = d["n_contexts"], d["n_perturbations"]
    a.replicates = d["n_replicates"]
    return cmd_plan(a)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="perturbdesign",
        description="What a perturbation atlas must measure to resolve a "
                    "context x perturbation interaction.",
        epilog="Method: RESULTS.md sec.46. Validated on five published atlases "
               "with one shared noise constant and no per-atlas tuning.")
    ap.add_argument("--snr", type=float, default=SNR,
                    help=f"signal-to-noise per condition (default {SNR}, the "
                         "value validated across five atlases)")
    ap.add_argument("--target", type=float, default=0.05,
                    help="interaction share you want to detect (default 0.05)")
    ap.add_argument("--model-dim", type=int, default=10,
                    help="context-embedding width of the model you plan to fit")
    ap.add_argument("--json", action="store_true", help="machine-readable")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("plan", help="what a proposed design can detect")
    p.add_argument("--contexts", type=int, required=True)
    p.add_argument("--perturbations", type=int, required=True)
    p.add_argument("--replicates", type=int, default=2)
    p.add_argument("--replicated-fraction", type=float, default=1.0,
                   help="fraction of conditions you will actually replicate "
                        "(Tahoe-100M: 0.135). Only these contribute pairs.")
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("audit", help="the same, for an existing .h5ad")
    p.add_argument("path")
    p.set_defaults(func=cmd_audit)

    p = sub.add_parser("budget", help="where a fixed cell budget should go")
    p.add_argument("--cells", type=int, required=True)
    p.add_argument("--contexts", type=int, required=True)
    p.add_argument("--perturbations", type=int, required=True)
    p.add_argument("--max-replicates", type=int, default=8)
    p.add_argument("--min-cells", type=int, default=50,
                   help="fewest cells per condition worth pseudobulking")
    p.add_argument("--snr-max", type=float, default=0.35,
                   help="signal-to-noise at infinite depth (biological floor)")
    p.add_argument("--snr-half", type=float, default=300.0,
                   help="cells per condition at which snr is half its maximum")
    p.set_defaults(func=cmd_budget)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
