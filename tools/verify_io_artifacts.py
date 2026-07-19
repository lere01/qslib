#!/usr/bin/env python3
"""Independently inspect qslib checkpoint and Parquet artifacts.

This verifier intentionally does not import qslib or reproduce its Rust data
structures. It checks durable JSON envelopes, little-endian NPY headers and
payload shapes, Parquet framing, and the exact completion marker. If `pyarrow`
is installed it also reads every trajectory part and checks the column names.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import struct
from pathlib import Path

try:
    import blake3  # type: ignore
except ImportError:
    blake3 = None


def read_npy(path: Path) -> tuple[list[float], tuple[int, ...]]:
    data = path.read_bytes()
    if len(data) < 12 or data[:6] != b"\x93NUMPY" or data[6:8] != b"\x01\x00":
        raise ValueError(f"{path}: unsupported NPY header")
    header_size = struct.unpack("<H", data[8:10])[0]
    end = 10 + header_size
    if end > len(data):
        raise ValueError(f"{path}: truncated NPY header")
    header = ast.literal_eval(data[10:end].decode("ascii"))
    if header.get("descr") != "<f8" or header.get("fortran_order") is not False:
        raise ValueError(f"{path}: expected little-endian C-order f64")
    shape = tuple(header.get("shape", ()))
    if not shape or any(not isinstance(size, int) or size <= 0 for size in shape):
        raise ValueError(f"{path}: invalid shape")
    count = 1
    for size in shape:
        count *= size
    payload = data[end:]
    if len(payload) != count * 8:
        raise ValueError(f"{path}: payload length does not match shape")
    values = list(struct.unpack(f"<{count}d", payload))
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"{path}: non-finite value")
    return values, shape


def verify_checkpoint(directory: Path) -> None:
    envelope = json.loads((directory / "checkpoint.json").read_text())
    required = {
        "schema_version",
        "convention_schema",
        "config_checksum",
        "run_id",
        "accepted_step",
        "parameter_layout_fingerprint",
        "rng_state_checksum",
        "rng_state_path",
        "payload_checksum",
        "payload_path",
        "accepted_state",
        "rng_state",
        "arrays",
    }
    if set(envelope) != required:
        raise ValueError("checkpoint envelope fields differ from qslib-checkpoint-v1")
    accepted = envelope["accepted_state"]
    rng = envelope["rng_state"]
    if accepted["schema_version"] != "qslib-accepted-state-v1" or accepted["dtype"] != "f64" or accepted["order"] != "C":
        raise ValueError("unsupported accepted-state metadata")
    evolution = accepted["evolution"]
    if evolution["method"] not in ("Euler", "Heun") or evolution["error_metric"] not in ("Euclidean", "Qgt"):
        raise ValueError("unsupported evolution controls")
    if evolution["seed_algorithm_version"] != 1 or evolution["dt_min"] > evolution["dt_max"]:
        raise ValueError("invalid evolution controls")
    if not (evolution["dt_min"] <= accepted["next_step"] <= evolution["dt_max"]):
        raise ValueError("accepted step is outside evolution bounds")
    if rng["schema_version"] != "qslib-rng-state-v1" or rng["algorithm"] != "chacha20":
        raise ValueError("unsupported RNG-state metadata")
    if not isinstance(rng["position_blocks"], int) or rng["position_blocks"] < 0:
        raise ValueError("invalid ChaCha20 block position")
    for field in ("payload_path", "rng_state_path"):
        if "/" in envelope[field] or "\\" in envelope[field] or envelope[field].startswith("."):
            raise ValueError(f"unsafe checkpoint path: {envelope[field]}")
        if not (directory / envelope[field]).is_file():
            raise ValueError(f"missing checkpoint payload: {envelope[field]}")
        if blake3 is not None:
            actual = blake3.blake3((directory / envelope[field]).read_bytes()).hexdigest()
            if actual != envelope[field.replace("_path", "_checksum")]:
                raise ValueError(f"checksum mismatch: {envelope[field]}")
    payload = json.loads((directory / envelope["payload_path"]).read_text())
    if payload != accepted:
        raise ValueError("accepted-state payload does not match its envelope metadata")
    rng_payload = json.loads((directory / envelope["rng_state_path"]).read_text())
    if rng_payload != rng:
        raise ValueError("RNG-state payload does not match its envelope metadata")
    for array in envelope["arrays"]:
        array_path = directory / array["path"]
        if "/" in array["path"] or "\\" in array["path"] or array["path"].startswith("."):
            raise ValueError(f"unsafe checkpoint array path: {array['path']}")
        values, shape = read_npy(array_path)
        if blake3 is not None and blake3.blake3(array_path.read_bytes()).hexdigest() != array["checksum"]:
            raise ValueError(f"checksum mismatch: {array['name']}")
        if tuple(array["shape"]) != shape or array["dtype"] != "f64" or array["order"] != "C":
            raise ValueError(f"array metadata mismatch: {array['name']}")
        if len(values) == 0:
            raise ValueError(f"empty checkpoint array: {array['name']}")
    print(f"verified checkpoint envelope and {len(envelope['arrays'])} NPY arrays: {directory}")


def verify_dataset(directory: Path) -> None:
    manifest = json.loads((directory / "manifest.json").read_text())
    if manifest.get("schema_version") != "qslib-parquet-dataset-v1":
        raise ValueError("unsupported Parquet dataset schema")
    if manifest.get("convention_schema") != "qslib-conventions-v1":
        raise ValueError("unsupported convention schema")
    marker = (directory / "COMPLETE").read_bytes()
    if marker != b"qslib-dataset-complete-v1\n" or not manifest.get("complete"):
        raise ValueError("dataset completion marker is missing or invalid")
    for part in manifest["parts"]:
        path = directory / part["path"]
        raw = path.read_bytes()
        if raw[:4] != b"PAR1" or raw[-4:] != b"PAR1":
            raise ValueError(f"invalid Parquet framing: {path}")
        if blake3 is not None and blake3.blake3(raw).hexdigest() != part["checksum"]:
            raise ValueError(f"checksum mismatch: {path}")
        try:
            import pyarrow.parquet as parquet  # type: ignore
        except ImportError:
            continue
        table = parquet.read_table(path)
        if table.column_names != ["step", "time", "energy"]:
            raise ValueError(f"unexpected trajectory columns: {path}")
    print(f"verified completed Parquet dataset with {len(manifest['parts'])} parts: {directory}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--dataset", type=Path)
    args = parser.parse_args()
    if args.checkpoint is None and args.dataset is None:
        parser.error("provide --checkpoint and/or --dataset")
    if args.checkpoint is not None:
        verify_checkpoint(args.checkpoint)
    if args.dataset is not None:
        verify_dataset(args.dataset)


if __name__ == "__main__":
    main()
