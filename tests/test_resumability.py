from __future__ import annotations

import json
from pathlib import Path

import pytest

from dynamic_panel_econ.resumability import (
    TaskCheckpointStore,
    scientific_fingerprint,
    semantic_task_specification,
)


def _fingerprint(marker: str = "same") -> dict[str, object]:
    return scientific_fingerprint(
        code_commit="abc123",
        source_tree_hash=marker,
        config_hash="config",
        calibration_hash="calibration",
        selector_method="revision10_ridge_ratio",
        master_seed=8675309,
        scientific_configuration={"B": 10, "caps": [3, 3, 3]},
    )


def _specifications(count: int = 6) -> list[dict[str, object]]:
    return [
        semantic_task_specification(
            dgp=99,
            n=8,
            t=7,
            true_rank=(index % 3, 0, 1),
            replication=index,
            master_seed=8675309,
            selector_method="revision10_ridge_ratio",
        )
        for index in range(count)
    ]


def _save(store: TaskCheckpointStore, specification: dict[str, object], state: str) -> None:
    task_id = str(specification["semantic_task_id"])
    store.mark_running(task_id, worker_hint="test")
    store.save_terminal(
        task_id,
        terminal_state=state,
        seed_metadata={"master_seed": 8675309, "replication": specification["replication"]},
        rows=[
            {
                "semantic_task_id": task_id,
                "status": state,
                "objective": float(specification["replication"]) / 10.0,
                "rank": [int(specification["replication"]) % 3, 0, 1],
            }
        ],
        metrics={"worker_pid": 1},
    )


def _scientific_payloads(store: TaskCheckpointStore) -> list[dict[str, object]]:
    payloads = []
    for task_id in sorted(store.specifications):
        payload = store.load_terminal(task_id)
        assert payload is not None
        payloads.append(
            {
                "task": payload["task_specification"],
                "seed": payload["seed_metadata"],
                "state": payload["terminal_state"],
                "rows": payload["rows"],
            }
        )
    return payloads


def test_uninterrupted_equals_interrupted_and_resumed(tmp_path: Path) -> None:
    specifications = _specifications()
    uninterrupted = TaskCheckpointStore(
        tmp_path / "uninterrupted", _fingerprint(), specifications, resume=False
    )
    states = ["completed", "completed", "failed", "unresolved", "completed", "completed"]
    for specification, state in zip(specifications, states, strict=True):
        _save(uninterrupted, specification, state)

    resumed_root = tmp_path / "resumed"
    interrupted = TaskCheckpointStore(resumed_root, _fingerprint(), specifications, resume=False)
    for specification, state in zip(specifications[:3], states[:3], strict=True):
        _save(interrupted, specification, state)
    interrupted.record_interruption("ctrl_c", "intentional deterministic interruption")
    stale_id = str(specifications[3]["semantic_task_id"])
    interrupted.mark_running(stale_id, worker_hint="dead-worker")

    resumed = TaskCheckpointStore(resumed_root, _fingerprint(), specifications, resume=True)
    assert resumed.manifest["tasks"][stale_id]["state"] == "expected"
    for specification, state in zip(specifications[3:], states[3:], strict=True):
        _save(resumed, specification, state)
    assert _scientific_payloads(resumed) == _scientific_payloads(uninterrupted)
    assert resumed.states() == {
        "completed": 4,
        "failed": 1,
        "unresolved": 1,
        "corrupt": 0,
    }


def test_resume_refuses_fingerprint_or_task_set_mismatch(tmp_path: Path) -> None:
    root = tmp_path / "run"
    specifications = _specifications(2)
    TaskCheckpointStore(root, _fingerprint(), specifications, resume=False)
    with pytest.raises(ValueError, match="fingerprint mismatch"):
        TaskCheckpointStore(root, _fingerprint("different"), specifications, resume=True)
    with pytest.raises(ValueError, match="expected-task set mismatch"):
        TaskCheckpointStore(root, _fingerprint(), _specifications(3), resume=True)


def test_corrupt_and_incomplete_bundles_are_quarantined(tmp_path: Path) -> None:
    specifications = _specifications(2)
    root = tmp_path / "run"
    store = TaskCheckpointStore(root, _fingerprint(), specifications, resume=False)
    first_id = str(specifications[0]["semantic_task_id"])
    second_id = str(specifications[1]["semantic_task_id"])
    first_path = store._bundle_path(first_id)
    first_path.parent.mkdir(parents=True, exist_ok=True)
    first_path.write_text("not-json", encoding="utf-8")
    partial = store._bundle_path(second_id).with_name("orphan.partial")
    partial.parent.mkdir(parents=True, exist_ok=True)
    partial.write_text("partial", encoding="utf-8")
    resumed = TaskCheckpointStore(root, _fingerprint(), specifications, resume=True)
    assert resumed.manifest["tasks"][first_id]["state"] == "expected"
    assert resumed.states()["corrupt"] == 2
    assert len(list((root / "task_quarantine").glob("*.corrupt"))) == 2


def test_manifest_exposes_pending_running_and_terminal_sets(tmp_path: Path) -> None:
    specifications = _specifications(3)
    store = TaskCheckpointStore(tmp_path / "run", _fingerprint(), specifications, resume=False)
    first_id = str(specifications[0]["semantic_task_id"])
    second_id = str(specifications[1]["semantic_task_id"])
    store.mark_running(first_id)
    _save(store, specifications[1], "unresolved")
    manifest = json.loads(store.manifest_path.read_text(encoding="utf-8"))
    assert first_id in manifest["running"]
    assert second_id in manifest["unresolved"]
    assert len(manifest["pending"]) == 1
    assert sorted(manifest["expected_task_ids"]) == sorted(store.specifications)
