from __future__ import annotations

import argparse
import importlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np


REPO_ROOT = Path(r"E:\b2f10k")
AVERAGE_MAPS_ROOT = REPO_ROOT / "averagemaps"
AGE_GROUPS_SUMMARY_PATH = AVERAGE_MAPS_ROOT / "age_groups_summary.json"
OUTPUTS_ROOT = REPO_ROOT / "xialpha_hansen"
CACHE_ROOT = REPO_ROOT / ".cache"
REFERENCE_DIR = REPO_ROOT / "hansen_reference"
HANSEN_ROOT = REPO_ROOT / "hansen_receptors-main"

AGE_GROUPS = ["age0_20", "age20_40", "age40_60", "age60_80", "age80_100"]
FEATURE_ORDER = [
    "Xi_estimate_Power",
    "Xi_estimate_Width",
    "Xi_estimate_Exponent",
    "Alpha_estimate_Power",
    "Alpha_estimate_Width",
    "Alpha_estimate_Exponent",
    "Alpha_estimate_PAF",
]
DISPLAY_NAMES = {
    "Alpha_estimate_Exponent": "AE",
    "Alpha_estimate_PAF": "AF",
    "Alpha_estimate_Power": "AA",
    "Alpha_estimate_Width": "AB",
    "Xi_estimate_Exponent": "XE",
    "Xi_estimate_Power": "XA",
    "Xi_estimate_Width": "XB",
}
ANCHOR_THRESHOLD = 0.05


def import_toolkit() -> tuple[object, object, object, object, object]:
    toolkit_roots = [
        REPO_ROOT,
        Path(r"E:\matlabtopymaps\matlab_to_py_maps"),
        Path(r"E:\matlab_to_py_maps\matlab_to_py_maps"),
    ]
    tried: list[str] = []
    for root in toolkit_roots:
        if not root.exists():
            continue
        hansen_file = root / "brainstorm_fsaverage_toolkit" / "hansen_dominance.py"
        if not hansen_file.exists():
            tried.append(str(root))
            continue

        if str(root) not in sys.path:
            sys.path.insert(0, str(root))

        for mod_name in list(sys.modules):
            if mod_name == "brainstorm_fsaverage_toolkit" or mod_name.startswith("brainstorm_fsaverage_toolkit."):
                del sys.modules[mod_name]

        try:
            cfg = importlib.import_module("brainstorm_fsaverage_toolkit.config")
            hd = importlib.import_module("brainstorm_fsaverage_toolkit.hansen_dominance")
            return (
                cfg.DEFAULT_NEUROMAPS_DATA_DIR,
                hd.load_hansen_receptor_reference,
                hd.load_metric,
                hd.normalize_unit_interval,
                hd.run_hansen_dominance_analysis,
            )
        except Exception:
            tried.append(str(root))
            continue

    raise ImportError(
        "Cannot import brainstorm_fsaverage_toolkit.hansen_dominance. "
        "Checked toolkit roots: " + ", ".join(tried + [str(p) for p in toolkit_roots if str(p) not in tried])
    )


def ensure_reference_dir(reference_dir: Path, hansen_root: Path) -> None:
    reference_dir.mkdir(parents=True, exist_ok=True)
    required_map = {
        hansen_root / "results" / "receptor_data_scale100.csv": reference_dir / "receptor_data_scale100.csv",
        hansen_root / "data" / "receptor_names_pet.npy": reference_dir / "receptor_names_pet.npy",
        hansen_root / "data" / "schaefer" / "coordinates" / "Schaefer_100_centres.txt": reference_dir / "Schaefer_100_centres.txt",
        hansen_root / "data" / "colourmap.csv": reference_dir / "colourmap.csv",
    }
    missing_src = [str(src) for src in required_map if not src.exists()]
    if missing_src:
        raise FileNotFoundError("Missing Hansen source files:\n- " + "\n- ".join(missing_src))
    for src, dst in required_map.items():
        if not dst.exists():
            shutil.copy2(src, dst)


def feature_map_paths(age_group: str, feature_name: str) -> tuple[Path, Path]:
    source_name, measure_name = feature_name.rsplit("_", 1)
    maps_dir = AVERAGE_MAPS_ROOT / age_group / source_name / measure_name / "maps"
    left = maps_dir / f"{measure_name}_space-fsaverage10k_hemi-L.shape.gii"
    right = maps_dir / f"{measure_name}_space-fsaverage10k_hemi-R.shape.gii"
    if not left.exists() or not right.exists():
        raise FileNotFoundError(f"Could not find fsaverage10k maps for {feature_name} in {age_group}")
    return left, right


def load_subject_counts() -> dict[str, int]:
    if not AGE_GROUPS_SUMMARY_PATH.exists():
        raise FileNotFoundError(f"Missing summary file: {AGE_GROUPS_SUMMARY_PATH}")
    with AGE_GROUPS_SUMMARY_PATH.open("r", encoding="utf-8") as f_obj:
        summary = json.load(f_obj)
    counts = {age_group: int(summary["age_groups"][age_group]["subject_count"]) for age_group in AGE_GROUPS}
    if any(count <= 0 for count in counts.values()):
        raise ValueError(f"Invalid subject counts: {counts}")
    return counts


def load_age_bin_surface_maps(age_group: str, load_metric_fn) -> dict[str, dict[str, np.ndarray]]:
    surface_maps: dict[str, dict[str, np.ndarray]] = {}
    for feature_name in FEATURE_ORDER:
        left_path, right_path = feature_map_paths(age_group, feature_name)
        surface_maps[feature_name] = {
            "L": load_metric_fn(left_path),
            "R": load_metric_fn(right_path),
        }
    return surface_maps


def build_global_weighted_surface_maps(
    subject_counts: dict[str, int], load_metric_fn
) -> dict[str, dict[str, np.ndarray]]:
    total_weight = float(sum(subject_counts.values()))
    if total_weight <= 0.0:
        raise ValueError("Subject counts must sum to a positive value")

    surface_maps: dict[str, dict[str, np.ndarray]] = {}
    for feature_name in FEATURE_ORDER:
        left_metric: np.ndarray | None = None
        right_metric: np.ndarray | None = None
        for age_group in AGE_GROUPS:
            left_path, right_path = feature_map_paths(age_group, feature_name)
            age_left = load_metric_fn(left_path)
            age_right = load_metric_fn(right_path)
            weight = float(subject_counts[age_group]) / total_weight
            if left_metric is None:
                left_metric = np.zeros_like(age_left, dtype=np.float64)
                right_metric = np.zeros_like(age_right, dtype=np.float64)
            left_metric += weight * age_left
            right_metric += weight * age_right
        if left_metric is None or right_metric is None:
            raise RuntimeError(f"Could not build global weighted surface maps for {feature_name}")
        surface_maps[feature_name] = {"L": left_metric, "R": right_metric}
    return surface_maps


def signed_log1p(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return np.sign(values) * np.log1p(np.abs(values))


def apply_xialpha_surface_rules(
    surface_maps: dict[str, dict[str, np.ndarray]],
    normalize_unit_interval_fn,
) -> dict[str, dict[str, np.ndarray]]:
    transformed = {
        feature_name: {"L": hemi_maps["L"].copy(), "R": hemi_maps["R"].copy()}
        for feature_name, hemi_maps in surface_maps.items()
    }
    aa_name = "Alpha_estimate_Power"
    xa_name = "Xi_estimate_Power"

    aa_raw = np.concatenate([surface_maps[aa_name]["L"], surface_maps[aa_name]["R"]])
    xa_raw = np.concatenate([surface_maps[xa_name]["L"], surface_maps[xa_name]["R"]])
    aa_anchor = normalize_unit_interval_fn(aa_raw)
    aa_anchor[aa_anchor < ANCHOR_THRESHOLD] = 0.0
    xa_mask = np.where(xa_raw > 0.0, 1.0, 0.0)

    for feature_name in transformed:
        joined = np.concatenate([transformed[feature_name]["L"], transformed[feature_name]["R"]])
        joined = signed_log1p(joined)
        split = transformed[feature_name]["L"].shape[0]
        transformed[feature_name]["L"] = joined[:split]
        transformed[feature_name]["R"] = joined[split:]

    alpha_family = [
        "Alpha_estimate_Power",
        "Alpha_estimate_Width",
        "Alpha_estimate_Exponent",
        "Alpha_estimate_PAF",
    ]
    xi_family = ["Xi_estimate_Power", "Xi_estimate_Width", "Xi_estimate_Exponent"]

    for feature_name in alpha_family:
        joined = np.concatenate([transformed[feature_name]["L"], transformed[feature_name]["R"]])
        joined = joined * aa_anchor
        split = transformed[feature_name]["L"].shape[0]
        transformed[feature_name]["L"] = joined[:split]
        transformed[feature_name]["R"] = joined[split:]

    for feature_name in xi_family:
        joined = np.concatenate([transformed[feature_name]["L"], transformed[feature_name]["R"]])
        joined = joined * xa_mask
        split = transformed[feature_name]["L"].shape[0]
        transformed[feature_name]["L"] = joined[:split]
        transformed[feature_name]["R"] = joined[split:]

    return transformed


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Hansen-style dominance analysis on xialpha fsaverage10k maps.")
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument("--age-bin", default="age0_20", help="Single age-bin to analyze")
    mode_group.add_argument("--global-weighted", action="store_true", help="Use the subject-count-weighted average across age bins")
    parser.add_argument("--n-spins", type=int, default=10_000, help="Number of spin permutations")
    parser.add_argument("--out-root", type=Path, default=OUTPUTS_ROOT, help="Output root directory")
    args = parser.parse_args()

    (
        default_neuromaps_data_dir,
        load_hansen_receptor_reference,
        load_metric,
        normalize_unit_interval,
        run_hansen_dominance_analysis,
    ) = import_toolkit()

    if not AVERAGE_MAPS_ROOT.exists():
        raise FileNotFoundError(f"Average maps root not found: {AVERAGE_MAPS_ROOT}")

    ensure_reference_dir(REFERENCE_DIR, HANSEN_ROOT)
    subject_counts = load_subject_counts()
    args.out_root.mkdir(parents=True, exist_ok=True)
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)

    if args.global_weighted:
        dataset_key = "global_weighted"
        dataset_name = "Global weighted xialpha"
        output_stem = "global_weighted_hansen_dominance"
        output_dir = args.out_root / output_stem
        surface_maps = build_global_weighted_surface_maps(subject_counts, load_metric)
        metadata = {
            "analysis_scope": "global_weighted",
            "age_bins": ",".join(AGE_GROUPS),
            "subject_counts_json": json.dumps(subject_counts, sort_keys=True),
            "raw_maps_root": str(AVERAGE_MAPS_ROOT),
        }
    else:
        age_bin = args.age_bin
        if not (AVERAGE_MAPS_ROOT / age_bin).exists():
            raise FileNotFoundError(f"Age bin folder not found: {AVERAGE_MAPS_ROOT / age_bin}")
        dataset_key = age_bin
        dataset_name = f"{age_bin} xialpha"
        output_stem = f"{age_bin}_hansen_dominance"
        output_dir = args.out_root / output_stem
        surface_maps = load_age_bin_surface_maps(age_bin, load_metric)
        metadata = {
            "analysis_scope": "age_bin",
            "age_bin": age_bin,
            "subject_count": int(subject_counts.get(age_bin, 0)),
            "raw_maps_root": str(AVERAGE_MAPS_ROOT / age_bin),
        }

    surface_maps = apply_xialpha_surface_rules(surface_maps, normalize_unit_interval)

    notes = [
        "A signed log1p transform was applied to all xi/alpha maps before weighting, plotting, parcellation, and regression.",
        f"Alpha-family maps (AA, AB, AE, AF) were weighted by thresholded normalized raw AA anchors at {ANCHOR_THRESHOLD:.2f}.",
        "Xi-family maps (XA, XB, XE) were masked only by raw XA > 0 without normalization.",
        "All predictor and target columns are max-scaled before z-scoring and regression.",
        "Off-cortex vertices are excluded from both the figure and the parcel-wise analysis.",
        "The Schaefer100 atlas is projected to native fsaverage10k using neuromaps and cached locally.",
        f"Feature display order is fixed to: {', '.join(DISPLAY_NAMES[name] for name in FEATURE_ORDER)}.",
    ]

    reference = load_hansen_receptor_reference(REFERENCE_DIR)
    results = run_hansen_dominance_analysis(
        surface_maps=surface_maps,
        output_dir=output_dir,
        dataset_name=dataset_name,
        reference=reference,
        output_stem=output_stem,
        display_names=DISPLAY_NAMES,
        feature_order=FEATURE_ORDER,
        atlas_data_dir=default_neuromaps_data_dir,
        cache_dir=CACHE_ROOT,
        n_spins=args.n_spins,
        significance_metric="spin_p",
        metadata=metadata | {"dataset_key": dataset_key},
        notes=notes,
    )
    print(f"Saved figure: {results.output_png}")
    print(f"Output dir: {output_dir}")


if __name__ == "__main__":
    main()
