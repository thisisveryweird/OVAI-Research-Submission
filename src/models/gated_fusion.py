from __future__ import annotations

from typing import Mapping

import torch
from torch import nn


class ModalityEncoder(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int = 128, dropout: float = 0.20):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class GatedMultimodalClassifier(nn.Module):
    """Late-fusion model with explicit missing-modality masks.

    Each modality first becomes a patient embedding. A learned gate then
    weights the embeddings, while the presence mask prevents a missing assay
    from being interpreted as zero biological signal.
    """

    def __init__(
        self,
        modality_dims: Mapping[str, int],
        hidden_dim: int = 128,
        dropout: float = 0.20,
    ) -> None:
        super().__init__()
        self.names = list(modality_dims)
        self.encoders = nn.ModuleDict({
            name: ModalityEncoder(dim, hidden_dim, dropout)
            for name, dim in modality_dims.items()
        })
        self.gates = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim // 2),
                nn.GELU(),
                nn.Linear(hidden_dim // 2, 1),
            )
            for name in self.names
        })
        self.classifier = nn.Sequential(
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(
        self,
        x: Mapping[str, torch.Tensor],
        presence: Mapping[str, torch.Tensor],
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        embeddings = {name: self.encoders[name](x[name]) for name in self.names}
        gate_logits = []
        gate_values = {}
        for name in self.names:
            g = self.gates[name](embeddings[name]).squeeze(-1)
            # Hard mask before softmax: a missing modality gets no weight.
            g = g.masked_fill(presence[name] <= 0, -1e9)
            gate_logits.append(g)

        logits = torch.stack(gate_logits, dim=1)
        gates = torch.softmax(logits, dim=1)
        for i, name in enumerate(self.names):
            gate_values[name] = gates[:, i]

        stacked = torch.stack([embeddings[name] for name in self.names], dim=1)
        fused = torch.sum(stacked * gates.unsqueeze(-1), dim=1)
        out = self.classifier(fused).squeeze(-1)
        return out, gate_values
