from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd


@dataclass
class ModalityBlock:
    name: str
    values: np.ndarray
    presence: np.ndarray
    feature_names: Sequence[str]


def align_modalities(
    clinical: pd.DataFrame,
    modality_tables: dict[str, pd.DataFrame],
    patient_id: str = "patient_id",
    label_col: str = "label",
) -> tuple[pd.DataFrame, dict[str, ModalityBlock]]:
    """Align modality tables by patient_id without inventing rows.

    The returned clinical frame defines the patient universe. For a modality
    that is absent for a patient, the values are NaN and presence=0. Downstream
    models must explicitly use the presence mask so missing data is not treated
    as biological zero.
    """
    if patient_id not in clinical.columns:
        raise ValueError(f"Clinical table must contain {patient_id!r}")
    if clinical[patient_id].duplicated().any():
        raise ValueError("Clinical table contains duplicated patient IDs")

    master = clinical.copy()
    blocks: dict[str, ModalityBlock] = {}

    for name, df in modality_tables.items():
        if patient_id not in df.columns:
            raise ValueError(f"Modality {name!r} lacks {patient_id!r}")
        if df[patient_id].duplicated().any():
            raise ValueError(f"Modality {name!r} contains duplicated patient IDs")

        feature_cols = [c for c in df.columns if c != patient_id]
        numeric = df[feature_cols].apply(pd.to_numeric, errors="coerce")
        aligned = numeric.set_index(df[patient_id]).reindex(master[patient_id])
        values = aligned.to_numpy(dtype=np.float32)
        presence = np.isfinite(values).any(axis=1).astype(np.float32)
        values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
        blocks[name] = ModalityBlock(
            name=name,
            values=values,
            presence=presence,
            feature_names=[str(c) for c in feature_cols],
        )

    return master, blocks
