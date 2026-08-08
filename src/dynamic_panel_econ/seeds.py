"""Worker-order-invariant random seed derivation."""

from __future__ import annotations

import hashlib

import numpy as np


def stable_uint32(value: object) -> int:
    digest = hashlib.sha256(str(value).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little", signed=False)


def seed_sequence(master_seed: int, *keys: object) -> np.random.SeedSequence:
    """Derive a deterministic seed from semantic identifiers, never worker order."""

    return np.random.SeedSequence([int(master_seed), *(stable_uint32(key) for key in keys)])


def rng_for(master_seed: int, *keys: object) -> np.random.Generator:
    return np.random.default_rng(seed_sequence(master_seed, *keys))
