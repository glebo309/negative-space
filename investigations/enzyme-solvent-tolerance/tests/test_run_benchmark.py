from pathlib import Path

import pytest

from enzsolv.run_benchmark import resolve_mmseqs_binary


def test_resolve_mmseqs_binary_prefers_environment_override(monkeypatch, tmp_path):
    binary = tmp_path / "mmseqs-custom"
    binary.write_text("")
    monkeypatch.setenv("MMSEQS_BIN", str(binary))

    assert resolve_mmseqs_binary(tmp_path) == binary


def test_resolve_mmseqs_binary_uses_path(monkeypatch, tmp_path):
    monkeypatch.delenv("MMSEQS_BIN", raising=False)
    binary = tmp_path / "mmseqs"
    binary.write_text("")
    monkeypatch.setattr("enzsolv.run_benchmark.shutil.which", lambda name: str(binary))

    assert resolve_mmseqs_binary(tmp_path) == binary


def test_resolve_mmseqs_binary_reports_installation_options(monkeypatch, tmp_path):
    monkeypatch.delenv("MMSEQS_BIN", raising=False)
    monkeypatch.setattr("enzsolv.run_benchmark.shutil.which", lambda name: None)

    with pytest.raises(FileNotFoundError, match="MMSEQS_BIN"):
        resolve_mmseqs_binary(tmp_path)
