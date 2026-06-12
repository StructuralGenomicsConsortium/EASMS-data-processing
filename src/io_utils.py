# -*- coding: utf-8 -*-
"""Filesystem helpers that work transparently on the local disk and on remote
object stores (e.g. Google Cloud Storage ``gs://`` paths).

The pipeline was written for local paths using ``os``/``open`` directly. To let
the exact same code run against a GCP bucket, every path is classified by a
simple rule: if it contains ``"://"`` it is treated as a remote URL and handled
via ``fsspec`` (which uses ``gcsfs`` for ``gs://`` URLs and your Application
Default Credentials); otherwise the original local ``os`` behaviour is kept
byte-for-byte.

Only the filesystem *primitives* live here. Pandas already reads/writes
``gs://`` URLs natively through fsspec, so ``pd.read_csv`` / ``df.to_csv`` /
``read_parquet`` / ``to_parquet`` / ``read_excel`` are left untouched at the call
sites -- only the path handed to them is built with :func:`pjoin`.
"""

import os
import posixpath
import contextlib


def _is_remote(path):
    """True if ``path`` is a remote URL (e.g. ``gs://bucket/folder``)."""
    return isinstance(path, str) and "://" in path


def pjoin(base, *parts):
    """Join paths, choosing the separator from ``base``.

    Remote URLs must use forward slashes, so on Windows ``os.path.join`` (which
    inserts ``\\``) would corrupt a ``gs://`` URL. For remote bases we use
    ``posixpath.join``; for local bases the behaviour is identical to
    ``os.path.join``.
    """
    if _is_remote(base):
        return posixpath.join(base, *parts)
    return os.path.join(base, *parts)


def dirname(path):
    """Parent directory of ``path``. Uses forward slashes for remote URLs so a
    ``gs://`` path stays a valid URL on Windows."""
    if _is_remote(path):
        return posixpath.dirname(path)
    return os.path.dirname(path)


def get_fs(path):
    """Return an fsspec filesystem object for ``path`` (gcsfs for ``gs://``)."""
    import fsspec

    fs, _ = fsspec.core.url_to_fs(path)
    return fs


def listdir(path):
    """List the *names* of entries in ``path`` (like ``os.listdir``).

    For remote paths, fsspec's ``ls`` returns full keys (without the URL
    scheme), so we strip them down to basenames. This keeps the existing
    idiom -- ``for f in listdir(d): ... pjoin(d, f)`` -- working unchanged.
    """
    if not _is_remote(path):
        return os.listdir(path)
    fs = get_fs(path)
    entries = fs.ls(path, detail=False)
    return [e.rstrip("/").split("/")[-1] for e in entries]


def makedirs(path, exist_ok=True):
    """Create a directory tree. No-op on object stores (they have no real
    directories -- a prefix exists once an object is written under it)."""
    if not _is_remote(path):
        os.makedirs(path, exist_ok=exist_ok)
        return
    try:
        get_fs(path).makedirs(path, exist_ok=exist_ok)
    except Exception:
        # Object stores have no directories; nothing to create.
        pass


def exists(path):
    """True if ``path`` exists (local file/dir or remote object/prefix)."""
    if not _is_remote(path):
        return os.path.exists(path)
    return get_fs(path).exists(path)


def isdir(path):
    """True if ``path`` is a directory (local) or an existing prefix (remote)."""
    if not _is_remote(path):
        return os.path.isdir(path)
    return get_fs(path).isdir(path)


def getsize(path):
    """Size of ``path`` in bytes (local file or remote object)."""
    if not _is_remote(path):
        return os.path.getsize(path)
    return get_fs(path).size(path)


@contextlib.contextmanager
def open_file(path, mode="r", encoding=None):
    """Open ``path`` for reading/writing on either backend.

    Supports text ("r"/"w"/"a") and binary ("wb"/"rb") modes. Note: true
    append is not supported by object stores, so an "a" mode on a remote path
    may overwrite rather than append; the pipeline only uses "a" for rare
    best-effort error notes, so this is acceptable.
    """
    if not _is_remote(path):
        with open(path, mode, encoding=encoding) as fh:
            yield fh
        return
    import fsspec

    with fsspec.open(path, mode, encoding=encoding) as fh:
        yield fh
