#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import torch

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "src"))
from models.gated_fusion import GatedMultimodalClassifier  # noqa: E402


def main() -> None:
    torch.manual_seed(42)
    model = GatedMultimodalClassifier({
        "clinical": 12,
        "mirna": 80,
        "methylation": 493,
        "cnv": 40,
        "mutations": 100,
    }, hidden_dim=64)
    batch = 8
    x = {
        "clinical": torch.randn(batch, 12),
        "mirna": torch.randn(batch, 80),
        "methylation": torch.randn(batch, 493),
        "cnv": torch.randn(batch, 40),
        "mutations": torch.randn(batch, 100),
    }
    presence = {
        "clinical": torch.ones(batch),
        "mirna": torch.ones(batch),
        "methylation": torch.tensor([1,1,1,0,1,0,1,1], dtype=torch.float32),
        "cnv": torch.ones(batch),
        "mutations": torch.tensor([1,0,1,1,0,1,1,0], dtype=torch.float32),
    }
    logits, gates = model(x, presence)
    assert logits.shape == (batch,)
    assert gates["methylation"][3].item() == 0.0 or presence["methylation"][3].item() == 0.0
    print("Smoke test passed.")
    print("Logits shape:", tuple(logits.shape))
    print("Gate sums:", torch.stack(list(gates.values()), dim=1).sum(dim=1).detach().numpy())


if __name__ == "__main__":
    main()
