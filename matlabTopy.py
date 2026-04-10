from pathlib import Path
import os
import json
import re
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.io import loadmat

from brainstorm_fsaverage_toolkit.atlas import load_fsaverage10k_atlas
from brainstorm_fsaverage_toolkit.brainstorm import load_brainstorm_cortex
from brainstorm_fsaverage_toolkit.config import DEFAULT_NEUROMAPS_DATA_DIR
from brainstorm_fsaverage_toolkit.gifti_io import save_shape_gifti
from brainstorm_fsaverage_toolkit.morph import map_brainstorm_to_fsaverage10k
from brainstorm_fsaverage_toolkit.visualization import plot_source_vs_fsaverage


def sanitize_name(name: str) -> str:
    """Normalize a name so it is safe for folder and file names."""
    name = name.strip()
    name = re.sub(r"[\\/:*?\"<>|]", "-", name)
    name = re.sub(r"\s+", "_", name)
    return name.strip("._-")


def extract_feature_vectors(mat_path: Path, expected_length: int) -> dict[str, np.ndarray]:
    """Extract all numeric 1D vectors matching expected_length from a .mat file."""
    if not mat_path.exists():
        raise FileNotFoundError(f"Feature file not found: {mat_path}")

    mat = loadmat(mat_path)
    features: dict[str, np.ndarray] = {}

    for key, value in mat.items():
        if key.startswith("__"):
            continue

        arr = np.asarray(value).squeeze()
        if not np.issubdtype(arr.dtype, np.number):
            continue

        if arr.ndim == 1 and arr.shape[0] == expected_length:
            features[key] = arr.astype(np.float64)

    return features


def save_feature_outputs(
    *,
    subject: str,
    source_name: str,
    feature_name: str,
    feature_values: np.ndarray,
    cortex,
    atlas,
    output_root: Path,
    save_figures: bool,
) -> None:
    """
    Save mapping outputs for a single feature.
    Directory layout: output_root / subject / source_name / feature_name / {maps,reports,(figures)}
    """
    subject_name = sanitize_name(subject)
    source_dir_name = sanitize_name(source_name)
    feature_dir_name = sanitize_name(feature_name)

    feature_dir = output_root / subject_name / source_dir_name / feature_dir_name
    maps_dir = feature_dir / "maps"
    reports_dir = feature_dir / "reports"

    maps_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    results = map_brainstorm_to_fsaverage10k(cortex, feature_values, atlas)

    out_l = maps_dir / f"{feature_dir_name}_space-fsaverage10k_hemi-L.shape.gii"
    out_r = maps_dir / f"{feature_dir_name}_space-fsaverage10k_hemi-R.shape.gii"

    save_shape_gifti(results["L"].values, out_l, f"{subject_name}_{source_dir_name}_{feature_dir_name}_L")
    save_shape_gifti(results["R"].values, out_r, f"{subject_name}_{source_dir_name}_{feature_dir_name}_R")

    figure_path = None
    if save_figures:
        figures_dir = feature_dir / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        figure_path = figures_dir / f"{feature_dir_name}_source_vs_fsaverage10k.png"
        plot_source_vs_fsaverage(
            cortex=cortex,
            source_values={
                "L": feature_values[cortex.hemispheres["L"].original_indices],
                "R": feature_values[cortex.hemispheres["R"].original_indices],
            },
            atlas=atlas,
            target_values={
                "L": results["L"].values,
                "R": results["R"].values,
            },
            title=f"{subject_name} | {source_dir_name} | {feature_name}",
            source_caption=f"Original {feature_name}",
            target_caption=f"Estimated fsaverage10k {feature_name}",
            out_png=figure_path,
            display_align_source=True,
        )

    summary = {
        "subject": subject_name,
        "source": source_name,
        "feature": feature_name,
        "outputs": {
            "L": str(out_l),
            "R": str(out_r),
            "figure": str(figure_path) if figure_path is not None else None,
        },
        "mapping_diagnostics": {
            "L": results["L"].metrics,
            "R": results["R"].metrics,
        },
    }

    summary_path = reports_dir / f"{feature_dir_name}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    report_lines = [
        "# Mapping Report",
        "",
        f"- Subject: `{subject_name}`",
        f"- Source file: `{source_name}`",
        f"- Feature variable: `{feature_name}`",
        f"- Left output: `{out_l}`",
        f"- Right output: `{out_r}`",
    ]
    if figure_path is not None:
        report_lines.append(f"- Figure: `{figure_path}`")
    report_lines.extend(
        [
            "",
            "## Mapping Diagnostics",
            f"- Left hemisphere: `{json.dumps(results['L'].metrics)}`",
            f"- Right hemisphere: `{json.dumps(results['R'].metrics)}`",
        ]
    )

    report_path = reports_dir / f"{feature_dir_name}_report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")


def process_one_mat_file(
    *,
    mat_path: Path,
    source_name: str,
    subject: str,
    cortex,
    atlas,
    expected_length: int,
    output_root: Path,
    save_outputs: bool,
    save_figures: bool,
) -> None:
    """Process one feature .mat file, such as Alpha_estimate.mat or Xi_estimate.mat."""
    features = extract_feature_vectors(mat_path, expected_length)
    if not features:
        raise ValueError(
            f"No 1D numeric vector with length {expected_length} was found in {mat_path}."
        )

    if not save_outputs:
        return

    for feature_name, feature_values in features.items():
        save_feature_outputs(
            subject=subject,
            source_name=source_name,
            feature_name=feature_name,
            feature_values=feature_values,
            cortex=cortex,
            atlas=atlas,
            output_root=output_root,
            save_figures=save_figures,
        )


def process_subject(
    *,
    subject_dir: Path,
    cortex,
    atlas,
    expected_length: int,
    output_root: Path,
    save_outputs: bool,
    save_figures: bool,
) -> None:
    """Process Alpha and Xi files for one subject directory."""
    subject = subject_dir.name
    alpha_mat = subject_dir / "Alpha_estimate.mat"
    xi_mat = subject_dir / "Xi_estimate.mat"

    process_one_mat_file(
        mat_path=alpha_mat,
        source_name="Alpha_estimate",
        subject=subject,
        cortex=cortex,
        atlas=atlas,
        expected_length=expected_length,
        output_root=output_root,
        save_outputs=save_outputs,
        save_figures=save_figures,
    )
    process_one_mat_file(
        mat_path=xi_mat,
        source_name="Xi_estimate",
        subject=subject,
        cortex=cortex,
        atlas=atlas,
        expected_length=expected_length,
        output_root=output_root,
        save_outputs=save_outputs,
        save_figures=save_figures,
    )


def collect_subject_dirs(data_root: Path) -> tuple[list[Path], int, int]:
    """
    Return subject folders containing both Alpha_estimate.mat and Xi_estimate.mat.
    Also return how many root entries were skipped and total root entries.
    """
    subject_dirs: list[Path] = []
    skipped_entries = 0
    total_entries = 0

    for entry in sorted(data_root.iterdir(), key=lambda p: p.name):
        total_entries += 1
        if not entry.is_dir():
            skipped_entries += 1
            continue

        alpha_mat = entry / "Alpha_estimate.mat"
        xi_mat = entry / "Xi_estimate.mat"
        if alpha_mat.exists() and xi_mat.exists():
            subject_dirs.append(entry)
        else:
            skipped_entries += 1

    return subject_dirs, skipped_entries, total_entries


def run_subject_job(
    *,
    subject_dir: Path,
    index: int,
    total_subjects: int,
    cortex,
    atlas,
    expected_length: int,
    output_root: Path,
    save_outputs: bool,
    save_figures: bool,
) -> tuple[int, str, str | None]:
    """Run one subject and return (index, subject_name, error_message)."""
    subject = subject_dir.name
    print(f"[PROCESSING] ({index}/{total_subjects}) {subject}", flush=True)
    try:
        process_subject(
            subject_dir=subject_dir,
            cortex=cortex,
            atlas=atlas,
            expected_length=expected_length,
            output_root=output_root,
            save_outputs=save_outputs,
            save_figures=save_figures,
        )
        return index, subject, None
    except Exception as exc:
        return index, subject, f"{type(exc).__name__}: {exc}"


def main() -> None:
    cortex_mat = Path(r"E:\b2f10k\data\Cortex.mat")
    data_root = Path(r"E:\xialphanet_newresults22")
    output_root = Path(r"E:\b2f10k\resultfigur")
    save_outputs = True
    save_figures = False
    use_parallel = False
    max_workers = min(3, os.cpu_count() or 1)

    output_root.mkdir(parents=True, exist_ok=True)

    error_log_path = output_root / "batch_errors.txt"
    batch_summary_path = output_root / "batch_summary.json"

    start_time = time.time()
    started_at = datetime.now().isoformat(timespec="seconds")

    if not cortex_mat.exists():
        raise FileNotFoundError(f"Cortex.mat not found: {cortex_mat}")
    if not data_root.exists():
        raise FileNotFoundError(f"Data root not found: {data_root}")

    print("[INFO] Loading Cortex and atlas...", flush=True)
    cortex = load_brainstorm_cortex(cortex_mat)
    atlas = load_fsaverage10k_atlas(DEFAULT_NEUROMAPS_DATA_DIR)
    expected_length = cortex.total_vertices

    subject_dirs, skipped_entries, total_root_entries = collect_subject_dirs(data_root)

    successes: list[str] = []
    failures: list[dict[str, str]] = []
    error_lines: list[str] = []

    total_subjects = len(subject_dirs)
    print(
        f"[INFO] Candidate subjects: {total_subjects}, skipped entries: {skipped_entries}, "
        f"parallel: {use_parallel}, workers: {max_workers}",
        flush=True,
    )

    if use_parallel and max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_subject = {
                executor.submit(
                    run_subject_job,
                    subject_dir=subject_dir,
                    index=index,
                    total_subjects=total_subjects,
                    cortex=cortex,
                    atlas=atlas,
                    expected_length=expected_length,
                    output_root=output_root,
                    save_outputs=save_outputs,
                    save_figures=save_figures,
                ): subject_dir.name
                for index, subject_dir in enumerate(subject_dirs, start=1)
            }

            completed = 0
            for future in as_completed(future_to_subject):
                completed += 1
                index, subject, error_message = future.result()
                if error_message is None:
                    successes.append(subject)
                    print(f"[DONE] ({completed}/{total_subjects}) {subject}", flush=True)
                else:
                    failures.append({"subject": subject, "error": error_message})
                    error_lines.append(f"[{subject}] {error_message}")
                    print(f"[ERROR] ({completed}/{total_subjects}) {subject}: {error_message}", flush=True)
    else:
        for index, subject_dir in enumerate(subject_dirs, start=1):
            _, subject, error_message = run_subject_job(
                subject_dir=subject_dir,
                index=index,
                total_subjects=total_subjects,
                cortex=cortex,
                atlas=atlas,
                expected_length=expected_length,
                output_root=output_root,
                save_outputs=save_outputs,
                save_figures=save_figures,
            )
            if error_message is None:
                successes.append(subject)
                print(f"[DONE] ({index}/{total_subjects}) {subject}", flush=True)
            else:
                failures.append({"subject": subject, "error": error_message})
                error_lines.append(f"[{subject}] {error_message}")
                print(f"[ERROR] ({index}/{total_subjects}) {subject}: {error_message}", flush=True)

    if error_lines:
        error_log_path.write_text("\n".join(error_lines) + "\n", encoding="utf-8")
    else:
        error_log_path.write_text("No subject-level errors.\n", encoding="utf-8")

    finished_at = datetime.now().isoformat(timespec="seconds")
    elapsed_seconds = round(time.time() - start_time, 2)

    batch_summary = {
        "started_at": started_at,
        "finished_at": finished_at,
        "elapsed_seconds": elapsed_seconds,
        "data_root": str(data_root),
        "output_root": str(output_root),
        "save_outputs": save_outputs,
        "save_figures": save_figures,
        "use_parallel": use_parallel,
        "max_workers": max_workers,
        "expected_length": expected_length,
        "total_root_entries": total_root_entries,
        "candidate_subjects": len(subject_dirs),
        "skipped_entries": skipped_entries,
        "success_count": len(successes),
        "failure_count": len(failures),
        "success_subjects": successes,
        "failures": failures,
        "error_log_path": str(error_log_path),
    }

    batch_summary_path.write_text(
        json.dumps(batch_summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print(
        f"Finished. Success: {len(successes)}, Failed: {len(failures)}. "
        f"Summary: {batch_summary_path}",
        flush=True,
    )


if __name__ == "__main__":
    main()
