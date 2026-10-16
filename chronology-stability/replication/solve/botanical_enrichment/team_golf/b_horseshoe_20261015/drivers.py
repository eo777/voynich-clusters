"""Exact transition drivers of the frozen 84-unit Golf B-related page MDS.

Reads published CSVs, writes only drivers.json, retains no model arrays.
Peak array allocation is well below 500 MB; output is capped below 2 MB.
"""
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import rel_entr
from scipy.stats import pearsonr, spearmanr
from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DATA = ROOT / "apps/clustermatrix/research-data"
SPEC = importlib.util.spec_from_file_location("b_driver_reproduce", DATA / "reproduce.py")
REPRO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPRO)
K = 3


def read_counts(mode, ids, alphabet):
    positions = {uid: i for i, uid in enumerate(ids)}
    glyphs = {glyph: i for i, glyph in enumerate(alphabet)}
    counts = np.zeros((len(ids), len(alphabet), len(alphabet)))
    found = set()
    with (DATA / f"transitions-{mode}-counts.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["profile_id"] not in positions:
                continue
            key = (row["profile_id"], row["from_unit"])
            assert key not in found
            found.add(key)
            vector = np.array([int(row["to_" + glyph]) for glyph in alphabet])
            assert vector.sum() == int(row["row_total"]) and np.all(vector >= 0)
            counts[positions[row["profile_id"]], glyphs[row["from_unit"]]] = vector
    assert len(found) == len(ids) * len(alphabet)
    return counts


def decompose(counts, vectors, values):
    """Exact cell split of eigenvalues and coordinates, with unit weight."""
    n, u, _ = counts.shape
    totals = counts.sum(axis=2)
    shares = totals / np.maximum(totals.sum(axis=1)[:, None], 1)
    probabilities = counts / np.maximum(totals[:, :, None], 1)
    eigen_contributions = np.zeros((K, u, u))
    reconstructed_d = np.zeros((n, n))
    reconstructed_coords = np.zeros((n, K))
    v = vectors[:, :K]
    for context in range(u):
        p = probabilities[:, context, :]
        midpoint = (p[:, None, :] + p[None, :, :]) / 2
        cell_js = (rel_entr(p[:, None, :], midpoint) + rel_entr(p[None, :, :], midpoint)) / (2 * np.log(2))
        eligible = totals[:, context] >= 3
        weights = (shares[:, context, None] + shares[None, :, context]) / 2
        weights *= eligible[:, None] & eligible[None, :]
        contributions = cell_js * weights[:, :, None]
        reconstructed_d += contributions.sum(axis=2)
        eigen_contributions[:, context, :] = -0.5 * np.einsum("ia,ka,ikj->aj", v, v, contributions)
        projected = contributions.sum(axis=2) @ v
        projected -= projected.mean(axis=0)
        reconstructed_coords += -0.5 * projected / np.sqrt(values[:K])
    np.fill_diagonal(reconstructed_d, 0)
    return eigen_contributions, reconstructed_d, reconstructed_coords


def main():
    viewer = [r for r in REPRO.rows("mds/viewer-branches.csv")
              if r["version"] == "golf21" and r["branch"] == "b" and r["view"] == "pages"]
    assert len(viewer) == 84
    by_viewer = {r["shared_mds_id"]: r for r in REPRO.rows("profiles.csv") if r["shared_mds_id"]}
    ids = [by_viewer[r["viewer_id"]]["profile_id"] for r in viewer]
    assert len(set(ids)) == 84
    alphabet = [r["unit"] for r in REPRO.rows("alphabet.csv")]
    counts = {mode: read_counts(mode, ids, alphabet) for mode in ["word", "line"]}
    reconstructed = {mode: REPRO.weighted_js(x)[0] for mode, x in counts.items()}
    reconstructed["combined"] = (reconstructed["word"] + reconstructed["line"]) / 2
    checks = {}
    for mode in ["word", "line", "combined"]:
        source_ids, source_d = REPRO.square(f"pairwise/shared-{mode}-js.csv")
        ix = [source_ids.index(uid) for uid in ids]
        expected = source_d[np.ix_(ix, ix)]
        checks[f"{mode}_distance_max_error"] = float(abs(reconstructed[mode] - expected).max())
        np.testing.assert_allclose(reconstructed[mode], expected, atol=1e-12, rtol=0)
    d = reconstructed["combined"]
    gram = -0.5 * (d - d.mean(axis=0) - d.mean(axis=1)[:, None] + d.mean())
    with threadpool_limits(1):
        values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    reference = np.array([[float(r[f"MDS{k}"]) for k in range(1, 5)] for r in viewer])
    for k in range(4):
        if vectors[:, k] @ reference[:, k] < 0:
            vectors[:, k] *= -1
    xy = vectors[:, :4] * np.sqrt(values[:4])
    checks["saved_first4_max_error"] = float(abs(xy - reference).max())
    np.testing.assert_allclose(xy, reference, atol=1e-10, rtol=0)

    decomposed = {}
    for mode, x in counts.items():
        contribution, rebuilt_d, rebuilt_xy = decompose(x, vectors, values)
        checks[f"{mode}_cell_distance_max_error"] = float(abs(rebuilt_d - reconstructed[mode]).max())
        np.testing.assert_allclose(rebuilt_d, reconstructed[mode], atol=1e-12, rtol=0)
        decomposed[mode] = {"eigenvalues": contribution / 2, "coords": rebuilt_xy / 2}
    eigen_sums = sum(v["eigenvalues"].sum(axis=(1, 2)) for v in decomposed.values())
    coordinate_sums = sum(v["coords"] for v in decomposed.values())
    checks["eigenvalue_sum_max_error"] = float(abs(eigen_sums - values[:K]).max())
    checks["coordinate_sum_max_error"] = float(abs(coordinate_sums - xy[:, :K]).max())
    np.testing.assert_allclose(eigen_sums, values[:K], atol=1e-12, rtol=1e-10)
    np.testing.assert_allclose(coordinate_sums, xy[:, :K], atol=1e-10, rtol=0)

    codes = np.array([r["membership_in_view"] for r in viewer])
    groups = sorted(set(codes))

    def frequency(mode, context, target, group, axis):
        x = counts[mode]
        ix = np.flatnonzero(codes == group)
        total = int(x[ix, context, :].sum())
        count = int(x[ix, context, target].sum())
        return {"target_count": count, "context_count": total,
                "pooled_conditional_rate": count / total if total else None,
                "units": len(ix)}

    def cell_row(axis, context, target, modes):
        eigenvalue = sum(decomposed[mode]["eigenvalues"][axis, context, target] for mode in modes)
        out = {
            "from": alphabet[context], "to": alphabet[target], "views": modes,
            "share_of_axis_eigenvalue": float(eigenvalue / values[axis]),
            "view_shares": {mode: float(decomposed[mode]["eigenvalues"][axis, context, target] / values[axis]) for mode in modes},
            "conditional_rates_by_view_and_viewer_group": {},
            "rate_axis_correlations": {},
        }
        for mode in modes:
            out["conditional_rates_by_view_and_viewer_group"][mode] = {
                group: frequency(mode, context, target, group, axis) for group in groups}
            totals = counts[mode][:, context, :].sum(axis=1)
            eligible = totals >= 10
            rates = counts[mode][eligible, context, target] / totals[eligible]
            out["rate_axis_correlations"][mode] = {
                "units_with_at_least10_context_events": int(eligible.sum()),
                "spearman": float(spearmanr(rates, xy[eligible, axis]).statistic)
                if eligible.sum() >= 6 and np.ptp(rates) > 0 else None,
            }
        return out

    axes = {}
    for axis in range(K):
        candidates = []
        for context, origin in enumerate(alphabet):
            for target, destination in enumerate(alphabet):
                modesets = [["word"], ["line"]] if origin == "^" or destination == "$" else [["word", "line"]]
                for modes in modesets:
                    share = sum(decomposed[mode]["eigenvalues"][axis, context, target] for mode in modes) / values[axis]
                    if abs(share) > 1e-12:
                        candidates.append((float(share), context, target, modes))
        candidates.sort(key=lambda row: -row[0])
        positive_shares = np.array([row[0] for row in candidates if row[0] > 0])
        axes[str(axis + 1)] = {
            "positive_eigenvalue_mass_share": float(values[axis] / values[values > 0].sum()),
            "eigenvalue": float(values[axis]),
            "view_eigenvalue_shares": {mode: float(decomposed[mode]["eigenvalues"][axis].sum() / values[axis]) for mode in counts},
            "negative_cell_share_sum": sum(row[0] for row in candidates if row[0] < 0),
            "n_merged_positive_cells_to_reach_half_axis": int(np.searchsorted(np.cumsum(positive_shares), 0.5) + 1),
            "top12_merged_cells": [cell_row(axis, context, target, modes)
                                   for _, context, target, modes in candidates[:12]],
            "all_nonzero_cell_shares": [{"from": alphabet[context], "to": alphabet[target],
                                          "views": modes, "share": share}
                                         for share, context, target, modes in candidates],
            "viewer_group_coordinates": {
                group: {"n": int(sum(codes == group)), "mean": float(xy[codes == group, axis].mean()),
                        "median": float(np.median(xy[codes == group, axis]))} for group in groups},
        }

    shared = {r["profile_id"]: r for r in REPRO.rows("mds/shared-combined.csv")}
    global8 = np.array([[float(shared[uid][f"MDS{k}"]) for k in range(1, 9)] for uid in ids])
    comparisons = {}
    for axis in range(K):
        correlations = [{"global_shared_axis": k + 1,
                         "pearson": float(pearsonr(xy[:, axis], global8[:, k]).statistic),
                         "spearman": float(spearmanr(xy[:, axis], global8[:, k]).statistic)} for k in range(8)]
        fits = {}
        for label, selection in {"global1_2": [0, 1], "global3_8": list(range(2, 8)), "global1_8": list(range(8))}.items():
            design = np.column_stack([np.ones(len(ids)), global8[:, selection]])
            fitted = design @ np.linalg.lstsq(design, xy[:, axis], rcond=None)[0]
            fits[label] = float(1 - np.sum((xy[:, axis] - fitted) ** 2) / np.sum((xy[:, axis] - xy[:, axis].mean()) ** 2))
        comparisons[str(axis + 1)] = {"pairwise_correlations": correlations,
                                     "same_sample_linear_projection_r2": fits}

    inputs = ["profiles.csv", "alphabet.csv", "mds/viewer-branches.csv", "mds/shared-combined.csv",
              "transitions-word-counts.csv", "transitions-line-counts.csv", "pairwise/shared-word-js.csv",
              "pairwise/shared-line-js.csv", "pairwise/shared-combined-js.csv", "reproduce.py"]
    output = {
        "study": "Transition drivers of the exact 84-unit B-related individual-page plot",
        "date": "2026-10-15", "n": len(ids), "bifolia": len({r["bifolium"] for r in viewer}),
        "population_selector": "mds/viewer-branches.csv: version=golf21, branch=b, view=pages",
        "kind_counts": dict(Counter(r["kind"] for r in viewer)),
        "viewer_group_counts": dict(Counter(codes)),
        "points": [{"profile_id": uid, "viewer_id": r["viewer_id"], "kind": r["kind"],
                    "bifolium": r["bifolium"], "viewer_group": r["membership_in_view"],
                    "profile_export_group": by_viewer[r["viewer_id"]]["golf21_code"],
                    "words": int(r["readable_words"]), "local_first4": xy[i].tolist()}
                   for i, (uid, r) in enumerate(zip(ids, viewer))],
        "method": {
            "distance": "Equal-weight word and joined-line context-weighted conditional JS in bits; contexts need >=3 events in both profiles; weights not renormalized",
            "embedding": "B=-0.5*J*D*J, with D=JS (not JS squared again); signs aligned to frozen local viewer",
            "eigen_cell_contribution": "lambda[k,c]=-0.5*v[k].T*D[c]*v[k]; each view has coefficient 0.5",
            "coordinate_cell_check": "x[k,c]=-0.5*J*D[c]*v[k]/sqrt(lambda[k]); all cells reproduce local coordinates",
            "merged_cell_definition": "Sum word and joined-line contributions for non-boundary transitions; keep word/line ^-initial and $-final contributions separate",
            "conditional_rates": "Sum target-event counts divided by sum source-context counts over units carrying each viewer group; separate word and line views",
            "global_comparison": "Restrict original 273-profile shared embedding to these same 84 units, then correlate local axes with global axes1..8; joint R2 is descriptive in-sample regression",
        },
        "limits": [
            "Exact decomposition explains observed geometry; it is not causal, chronological, or independent evidence for derived group assignments.",
            "The viewer population includes 81 paragraphs and three non-paragraph profiles; it is neither all Currier B nor paragraph-only B.",
            "Group frequencies weight by context counts, whereas MDS observations have equal point weight.",
            "Positive and negative cell shares can coexist because of context weighting; all shares sum to the axis eigenvalue exactly.",
            "Correlations with global axes do not assign identical meaning to global and local coordinates; the local fit is newly centered and rotated.",
            "Local eigenvalues2 and3 are close; their exact axis orientation can be sensitive to perturbations even if their joint plane is stable.",
            "No uncertainty intervals, resampling significance tests, or historical direction are supplied.",
        ],
        "geometry": {"positive_dimensions": int((values > 1e-12).sum()),
                     "positive_eigenvalue_shares_first8": (values[:8] / values[values > 0].sum()).tolist(),
                     "lambda2_over_lambda3": float(values[1] / values[2]),
                     "negative_eigenvalue_mass": float(-values[values < 0].sum() / np.abs(values).sum())},
        "checks": checks, "axes": axes, "comparison_to_global_shared_axes": comparisons,
        "sha256": {name: hashlib.sha256((DATA / name).read_bytes()).hexdigest() for name in inputs},
    }
    output["sha256"][str(Path(__file__).relative_to(ROOT))] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(output, indent=2, allow_nan=False) + "\n"
    assert len(text.encode("utf-8")) < 2_000_000
    (HERE / "drivers.json").write_text(text, encoding="utf-8")
    print(json.dumps({"bytes": len(text.encode("utf-8")), "checks": checks,
                      "axes": {key: [{"from": row["from"], "to": row["to"], "views": row["views"],
                                      "share": row["share_of_axis_eigenvalue"]} for row in value["top12_merged_cells"][:6]]
                               for key, value in axes.items()},
                      "global_axis_correlations": comparisons}, indent=2))


if __name__ == "__main__":
    main()
