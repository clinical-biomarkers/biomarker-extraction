"""
log_utils.py
Logging utilities for the SalivaDB extraction pipeline.
"""

from __future__ import annotations

import logging
import sys
from datetime import datetime
from pathlib import Path


class DeduplicatingFilter(logging.Filter):
    """
    A logging filter that passes each unique message exactly once.

    Uniqueness is determined by the fully-rendered message string
    (i.e. after %-interpolation), so records that differ only in
    their timestamp or level are still considered duplicates if the
    message text is identical.

    Usage
    -----
    Attach to the root logger to deduplicate across all handlers::

        import logging
        from log_utils import setup_logging

        setup_logging(Path("logs"))

    The deduplicating filter is applied automatically by setup_logging.
    To attach it manually to a single handler instead::

        handler.addFilter(DeduplicatingFilter())
    """

    def __init__(self) -> None:
        super().__init__()
        self._seen: set[str] = set()

    def filter(self, record: logging.LogRecord) -> bool:
        msg = record.getMessage()
        if msg in self._seen:
            return False
        self._seen.add(msg)
        return True

def setup_logging(log_dir: Path, script_stem: str = "generate_salivadb") -> None:
    """
    Configure root-level logging to a timestamped file and stdout,
    with duplicate lines suppressed.

    Parameters
    ----------
    log_dir:
        Directory in which to create the log file (created if absent).
    script_stem:
        Base name used in the log filename, e.g. ``"generate_salivadb"``
        produces ``generate_salivadb_20260101_120000.log``.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{script_stem}_{timestamp}.log"

    dedup = DeduplicatingFilter()

    file_handler = logging.FileHandler(log_file)
    file_handler.addFilter(dedup)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.addFilter(dedup)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[file_handler, stream_handler],
    )
