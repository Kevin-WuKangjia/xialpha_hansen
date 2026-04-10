from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import nibabel as nib
import numpy as np

from brainstorm_fsaverage_toolkit.gifti_io import save_shape_gifti


RESULT_ROOT = Path(r"E:\b2f10k\result")
META_JSON_PATH = Path(r"E:\xialphanet_newresults22\XIALPHANET.json")
OUTPUT_ROOT = Path(r"E:\b2f10k\averagemaps")

AGE_GROUPS = ("age0_20", "age20_40", "age40_60", "age60_80", "age80_100")


def age_group_name(age: float) -> str | None:
    if 0 <= age <= 20:
        return "age0_20"
    if 20 < age <= 40:
        return "age20_40"
    if 40 < age <= 60:
        return "age40_60"
    if 60 < age <= 80:
        return "age60_80"
    if 80 < age <= 100:
        return "age80_100"
    return None


def load_age_map(meta_json_path: Path) -> dict[str, float]:
    if not meta_json_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_json_path}")

    data = json.loads(meta_json_path.read_text(encoding="utf-8"))
    participants = data.get("Participants")
    if not isinstance(participants, list):
        raise ValueError("Invalid JSON format: key 'Participants' is missing or not a list.")

    out: dict[str, float] = {}
    for item in participants:
        subid = str(item.get("SubID", "")).strip()
        if not subid:
            continue
        age_raw = item.get("Age")
        try:
            age = float(age_raw)
        except (TypeError, ValueError):
            continue
        out[subid] = age
    return out


def load_gii_values(path: Path) -> np.ndarray:
    img = nib.load(str(path))
    return np.asarray(img.darrays[0].data, dtype=np.float64).ravel()


def pick_hemi_file(maps_dir: Path, feature_name: str, hemi: str) -> Path | None:
    expected = maps_dir / f"{feature_name}_space-fsaverage10k_hemi-{hemi}.shape.gii"
    if expected.exists():
        return expected
    candidates = sorted(maps_dir.glob(f"*_hemi-{hemi}.shape.gii"))
    if not candidates:
        return None
    return candidates[0]


@dataclass
class FeatureAccumulator:
    sum_l: np.ndarray | None = None
    sum_r: np.ndarray | None = None
    count: int = 0
    subjects: set[str] = field(default_factory=set)
    ages: list[float] = field(default_factory=list)

    def add(self, *, left_values: np.ndarray, right_values: np.ndarray, subject: str, age: float) -> None:
        l_arr = np.asarray(left_values, dtype=np.float64).ravel()
        r_arr = np.asarray(right_values, dtype=np.float64).ravel()

        if self.sum_l is None:
            self.sum_l = np.zeros_like(l_arr, dtype=np.float64)
        elif self.sum_l.shape != l_arr.shape:
            raise ValueError(
                f"L hemisphere length mismatch: got {l_arr.shape[0]}, expected {self.sum_l.shape[0]}"
            )

        if self.sum_r is None:
            self.sum_r = np.zeros_like(r_arr, dtype=np.float64)
        elif self.sum_r.shape != r_arr.shape:
            raise ValueError(
                f"R hemisphere length mismatch: got {r_arr.shape[0]}, expected {self.sum_r.shape[0]}"
            )

        self.sum_l += l_arr
        self.sum_r += r_arr
        self.count += 1
        self.subjects.add(subject)
        self.ages.append(age)

    def mean_l(self) -> np.ndarray:
        if self.sum_l is None or self.count <= 0:
            raise ValueError("No L hemisphere data to average.")
        return (self.sum_l / self.count).astype(np.float32)

    def mean_r(self) -> np.ndarray:
        if self.sum_r is None or self.count <= 0:
            raise ValueError("No R hemisphere data to average.")
        return (self.sum_r / self.count).astype(np.float32)

    def mean_age(self) -> float:
        if not self.ages:
            return 0.0
        return float(np.mean(self.ages))


def main() -> None:
    if not RESULT_ROOT.exists():
        raise FileNotFoundError(f"Result root not found: {RESULT_ROOT}")

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    age_map = load_age_map(META_JSON_PATH)

    accumulators: dict[tuple[str, str, str], FeatureAccumulator] = {}
    group_subject_ages: dict[str, dict[str, float]] = {g: {} for g in AGE_GROUPS}

    scanned_subject_dirs = 0
    skipped_no_age = 0
    skipped_age_out_of_range = 0
    skipped_missing_maps = 0
    skipped_missing_hemi = 0
    skipped_non_subject_entries = 0

    subject_dirs: list[Path] = []
    for entry in sorted(RESULT_ROOT.iterdir(), key=lambda p: p.name):
        if entry.is_dir():
            subject_dirs.append(entry)
        else:
            skipped_non_subject_entries += 1
    total_subject_dirs = len(subject_dirs)
    print(f"[INFO] Subject folders to scan: {total_subject_dirs}", flush=True)

    for idx, subject_dir in enumerate(subject_dirs, start=1):
        scanned_subject_dirs += 1

        if idx % 200 == 0 or idx == total_subject_dirs:
            print(f"[SCAN] {idx}/{total_subject_dirs} subject folders", flush=True)

        subject = subject_dir.name
        age = age_map.get(subject)
        if age is None:
            skipped_no_age += 1
            continue

        group = age_group_name(age)
        if group is None:
            skipped_age_out_of_range += 1
            continue

        subject_used = False

        for source_dir in sorted(subject_dir.iterdir(), key=lambda p: p.name):
            if not source_dir.is_dir():
                continue
            source_name = source_dir.name

            for feature_dir in sorted(source_dir.iterdir(), key=lambda p: p.name):
                if not feature_dir.is_dir():
                    continue
                feature_name = feature_dir.name
                maps_dir = feature_dir / "maps"
                if not maps_dir.is_dir():
                    skipped_missing_maps += 1
                    continue

                left_file = pick_hemi_file(maps_dir, feature_name, "L")
                right_file = pick_hemi_file(maps_dir, feature_name, "R")
                if left_file is None or right_file is None:
                    skipped_missing_hemi += 1
                    continue

                left_values = load_gii_values(left_file)
                right_values = load_gii_values(right_file)

                key = (group, source_name, feature_name)
                if key not in accumulators:
                    accumulators[key] = FeatureAccumulator()
                accumulators[key].add(
                    left_values=left_values,
                    right_values=right_values,
                    subject=subject,
                    age=age,
                )
                subject_used = True

        if subject_used:
            group_subject_ages[group][subject] = age

    feature_summaries: list[dict[str, object]] = []
    written_files = 0

    age_group_feature_items: dict[str, list[tuple[str, str, FeatureAccumulator]]] = {g: [] for g in AGE_GROUPS}
    for (group, source_name, feature_name), acc in sorted(accumulators.items()):
        age_group_feature_items[group].append((source_name, feature_name, acc))

    completed_age_groups = 0
    total_age_groups = len(AGE_GROUPS)

    for group in AGE_GROUPS:
        items = age_group_feature_items[group]
        if not items:
            completed_age_groups += 1
            print(
                f"[GROUP] ({completed_age_groups}/{total_age_groups}) {group} finished: no valid feature maps",
                flush=True,
            )
            continue

        group_feature_written = 0
        group_gii_written = 0

        for source_name, feature_name, acc in items:
            out_feature_dir = OUTPUT_ROOT / group / source_name / feature_name
            out_maps_dir = out_feature_dir / "maps"
            out_maps_dir.mkdir(parents=True, exist_ok=True)

            out_l = out_maps_dir / f"{feature_name}_space-fsaverage10k_hemi-L.shape.gii"
            out_r = out_maps_dir / f"{feature_name}_space-fsaverage10k_hemi-R.shape.gii"

            save_shape_gifti(acc.mean_l(), out_l, f"{group}_{source_name}_{feature_name}_L_mean")
            save_shape_gifti(acc.mean_r(), out_r, f"{group}_{source_name}_{feature_name}_R_mean")
            written_files += 2
            group_gii_written += 2

            feature_summary = {
                "age_group": group,
                "source": source_name,
                "feature": feature_name,
                "subject_count": acc.count,
                "subjects": sorted(acc.subjects),
                "mean_age": round(acc.mean_age(), 4),
                "outputs": {
                    "L": str(out_l),
                    "R": str(out_r),
                },
            }
            feature_summary_path = out_feature_dir / "summary.json"
            feature_summary_path.write_text(
                json.dumps(feature_summary, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            feature_summaries.append(feature_summary)
            group_feature_written += 1

        completed_age_groups += 1
        print(
            f"[GROUP] ({completed_age_groups}/{total_age_groups}) {group} finished: "
            f"feature groups={group_feature_written}, gii files={group_gii_written}",
            flush=True,
        )

    group_summaries: dict[str, dict[str, object]] = {}
    for group in AGE_GROUPS:
        subject_ages = group_subject_ages[group]
        ages = list(subject_ages.values())
        group_summaries[group] = {
            "subject_count": len(subject_ages),
            "subjects": sorted(subject_ages.keys()),
            "mean_age": round(float(np.mean(ages)), 4) if ages else None,
            "min_age": round(float(np.min(ages)), 4) if ages else None,
            "max_age": round(float(np.max(ages)), 4) if ages else None,
        }

    overall_summary = {
        "result_root": str(RESULT_ROOT),
        "meta_json": str(META_JSON_PATH),
        "output_root": str(OUTPUT_ROOT),
        "scanned_subject_dirs": scanned_subject_dirs,
        "written_feature_groups": len(feature_summaries),
        "written_gii_files": written_files,
        "skipped_non_subject_entries": skipped_non_subject_entries,
        "skipped_no_age": skipped_no_age,
        "skipped_age_out_of_range": skipped_age_out_of_range,
        "skipped_missing_maps": skipped_missing_maps,
        "skipped_missing_hemi": skipped_missing_hemi,
        "age_groups": group_summaries,
    }

    overall_summary_path = OUTPUT_ROOT / "age_groups_summary.json"
    overall_summary_path.write_text(
        json.dumps(overall_summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    detail_summary_path = OUTPUT_ROOT / "averagemaps_summary.json"
    detail_summary_path.write_text(
        json.dumps(feature_summaries, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("Average map generation finished.", flush=True)
    print(f"Written feature groups: {len(feature_summaries)}", flush=True)
    print(f"Written gii files: {written_files}", flush=True)
    print(f"Group summary: {overall_summary_path}", flush=True)


if __name__ == "__main__":
    main()
