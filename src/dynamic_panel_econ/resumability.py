"""Atomic task checkpoints and deterministic resume state for long simulations."""

from __future__ import annotations

import hashlib
import json
import math
import os
import time
import uuid
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np

BUNDLE_SCHEMA = "dynamic-panel-task-bundle-v1"
MANIFEST_SCHEMA = "dynamic-panel-live-manifest-v1"
TERMINAL_STATES = frozenset({"completed", "failed", "unresolved"})


def _safe_json(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [_safe_json(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _safe_json(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return {"__nonfinite_float__": str(value)}
    if isinstance(value, Mapping):
        return {str(key): _safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe_json(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _restore_json(value: Any) -> Any:
    if isinstance(value, dict):
        if set(value) == {"__nonfinite_float__"}:
            return float(value["__nonfinite_float__"])
        return {key: _restore_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_json(item) for item in value]
    return value


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        _safe_json(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def sha256_value(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def semantic_task_specification(
    *,
    dgp: int,
    n: int,
    t: int,
    true_rank: tuple[int, ...],
    replication: int,
    master_seed: int,
    selector_method: str,
) -> dict[str, Any]:
    specification = {
        "dgp": int(dgp),
        "N": int(n),
        "T": int(t),
        "true_rank": [int(value) for value in true_rank],
        "replication": int(replication),
        "master_seed": int(master_seed),
        "selector_method": str(selector_method),
    }
    specification["semantic_task_id"] = sha256_value(specification)
    return specification


def scientific_fingerprint(
    *,
    code_commit: str,
    source_tree_hash: str,
    config_hash: str,
    calibration_hash: str,
    selector_method: str,
    master_seed: int,
    scientific_configuration: Mapping[str, Any],
) -> dict[str, Any]:
    fingerprint = {
        "code_commit": code_commit,
        "source_tree_hash": source_tree_hash,
        "config_hash": config_hash,
        "calibration_hash": calibration_hash,
        "selector_method": selector_method,
        "master_seed": int(master_seed),
        "scientific_configuration": _safe_json(scientific_configuration),
        "bundle_schema": BUNDLE_SCHEMA,
    }
    fingerprint["fingerprint_hash"] = sha256_value(fingerprint)
    return fingerprint


def _atomic_write(destination: Path, content: bytes) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.partial")
    with temporary.open("wb") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, destination)
    if os.name != "nt":
        descriptor = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


class TaskCheckpointStore:
    """Sole-writer durable manifest and one immutable bundle per semantic task."""

    def __init__(
        self,
        root: Path,
        fingerprint: Mapping[str, Any],
        task_specifications: Iterable[Mapping[str, Any]],
        *,
        resume: bool,
    ) -> None:
        self.root = Path(root)
        self.task_root = self.root / "task_bundles"
        self.corrupt_root = self.root / "task_quarantine"
        self.manifest_path = self.root / "task_manifest.json"
        self.fingerprint = dict(fingerprint)
        self.fingerprint_hash = str(self.fingerprint["fingerprint_hash"])
        specification_list = [dict(item) for item in task_specifications]
        self.specifications = {
            str(item["semantic_task_id"]): item for item in specification_list
        }
        if len(self.specifications) != len(specification_list):
            raise ValueError("duplicate semantic task IDs")
        if self.manifest_path.exists():
            if not resume:
                raise FileExistsError("task manifest exists; resume is required")
            self.manifest = self._read_manifest()
            if self.manifest["fingerprint_hash"] != self.fingerprint_hash:
                raise ValueError("resume fingerprint mismatch")
            if set(self.manifest["tasks"]) != set(self.specifications):
                raise ValueError("resume expected-task set mismatch")
        else:
            self.manifest = {
                "schema": MANIFEST_SCHEMA,
                "fingerprint": self.fingerprint,
                "fingerprint_hash": self.fingerprint_hash,
                "tasks": {
                    task_id: {
                        "state": "expected",
                        "bundle_hash": None,
                        "attempt": 0,
                    }
                    for task_id in sorted(self.specifications)
                },
                "interruptions": [],
                "corrupt": [],
                "updated_unix": time.time(),
            }
            self._write_manifest()
        self.reconcile()

    def _bundle_path(self, task_id: str) -> Path:
        return self.task_root / task_id[:2] / f"{task_id}.bundle.json"

    def _read_manifest(self) -> dict[str, Any]:
        data = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        if data.get("schema") != MANIFEST_SCHEMA:
            raise ValueError("unsupported task-manifest schema")
        return data

    def _write_manifest(self) -> None:
        self.manifest["expected_task_ids"] = sorted(self.manifest["tasks"])
        for state in ("expected", "running", "completed", "failed", "unresolved"):
            label = "pending" if state == "expected" else state
            self.manifest[label] = sorted(
                task_id
                for task_id, record in self.manifest["tasks"].items()
                if record["state"] == state
            )
        self.manifest["updated_unix"] = time.time()
        _atomic_write(self.manifest_path, json.dumps(self.manifest, indent=2, sort_keys=True).encode("utf-8") + b"\n")

    def _validate_bundle(self, path: Path, task_id: str) -> tuple[dict[str, Any], str]:
        raw = path.read_bytes()
        file_hash = hashlib.sha256(raw).hexdigest()
        wrapper = json.loads(raw.decode("utf-8"))
        if wrapper.get("schema") != BUNDLE_SCHEMA:
            raise ValueError("unsupported task-bundle schema")
        payload = wrapper.get("payload")
        if wrapper.get("payload_hash") != sha256_value(payload):
            raise ValueError("task-bundle payload hash mismatch")
        if payload.get("semantic_task_id") != task_id:
            raise ValueError("task-bundle semantic ID mismatch")
        if payload.get("fingerprint_hash") != self.fingerprint_hash:
            raise ValueError("task-bundle fingerprint mismatch")
        if payload.get("terminal_state") not in TERMINAL_STATES:
            raise ValueError("task bundle is not terminal")
        if payload.get("task_specification") != _safe_json(self.specifications[task_id]):
            raise ValueError("task-bundle specification mismatch")
        if not isinstance(payload.get("rows"), list):
            raise ValueError("task-bundle rows are malformed")
        return _restore_json(payload), file_hash

    def _quarantine(self, path: Path, reason: str) -> None:
        self.corrupt_root.mkdir(parents=True, exist_ok=True)
        destination = self.corrupt_root / f"{path.name}.{uuid.uuid4().hex}.corrupt"
        os.replace(path, destination)
        self.manifest["corrupt"].append(
            {"source": str(path), "quarantine": str(destination), "reason": reason, "time": time.time()}
        )

    def reconcile(self) -> None:
        for temporary in self.task_root.glob("**/*.partial"):
            self._quarantine(temporary, "incomplete_temporary_bundle")
        for task_id in sorted(self.specifications):
            path = self._bundle_path(task_id)
            record = self.manifest["tasks"][task_id]
            if path.exists():
                try:
                    payload, file_hash = self._validate_bundle(path, task_id)
                except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
                    self._quarantine(path, f"{type(exc).__name__}: {exc}")
                    record.update({"state": "expected", "bundle_hash": None})
                else:
                    record.update(
                        {"state": payload["terminal_state"], "bundle_hash": file_hash}
                    )
            elif record["state"] == "running":
                record.update({"state": "expected", "bundle_hash": None})
        self._write_manifest()

    def mark_running(self, task_id: str, *, worker_hint: str | None = None) -> None:
        record = self.manifest["tasks"][task_id]
        if record["state"] in TERMINAL_STATES:
            return
        record.update(
            {
                "state": "running",
                "attempt": int(record.get("attempt", 0)) + 1,
                "worker_hint": worker_hint,
                "started_unix": time.time(),
            }
        )
        self._write_manifest()

    def save_terminal(
        self,
        task_id: str,
        *,
        terminal_state: str,
        seed_metadata: Mapping[str, Any],
        rows: list[dict[str, Any]],
        metrics: Mapping[str, Any] | None = None,
    ) -> dict[str, float | str]:
        if terminal_state not in TERMINAL_STATES:
            raise ValueError("invalid terminal state")
        serialization_started = time.perf_counter()
        payload = {
            "semantic_task_id": task_id,
            "task_specification": self.specifications[task_id],
            "fingerprint_hash": self.fingerprint_hash,
            "terminal_state": terminal_state,
            "seed_metadata": dict(seed_metadata),
            "rows": rows,
            "metrics": dict(metrics or {}),
        }
        wrapper = {
            "schema": BUNDLE_SCHEMA,
            "payload_hash": sha256_value(payload),
            "payload": _safe_json(payload),
        }
        content = json.dumps(wrapper, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        serialization_seconds = time.perf_counter() - serialization_started
        write_started = time.perf_counter()
        destination = self._bundle_path(task_id)
        _atomic_write(destination, content)
        write_seconds = time.perf_counter() - write_started
        _, file_hash = self._validate_bundle(destination, task_id)
        self.manifest["tasks"][task_id].update(
            {
                "state": terminal_state,
                "bundle_hash": file_hash,
                "completed_unix": time.time(),
            }
        )
        self._write_manifest()
        return {
            "bundle_hash": file_hash,
            "serialization_seconds": serialization_seconds,
            "write_seconds": write_seconds,
        }

    def load_terminal(self, task_id: str) -> dict[str, Any] | None:
        if self.manifest["tasks"][task_id]["state"] not in TERMINAL_STATES:
            return None
        payload, _ = self._validate_bundle(self._bundle_path(task_id), task_id)
        return payload

    def record_interruption(self, kind: str, detail: str) -> None:
        self.manifest["interruptions"].append(
            {"kind": kind, "detail": detail, "time": time.time()}
        )
        self._write_manifest()

    def states(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for record in self.manifest["tasks"].values():
            state = str(record["state"])
            counts[state] = counts.get(state, 0) + 1
        counts["corrupt"] = len(self.manifest["corrupt"])
        return counts
