import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
from pathlib import Path
from zlib import crc32

log = logging.getLogger("CAT_pack")

MAX_WORKERS = 16
MIN_TO_RUN_PARALLEL = 128 * 1024**2
CHUNK_BYTES = 32 * 1024**2
FILTER_BITS_PER_HIT = 32
MAX_FILTER_BYTES = 256 * 1024**2

_hit_filter = None
_hit_mask = 0


def read_fastaid2LCAtaxid(fastaid2LCAtaxid_file: Path, all_hits: set[str], workers: int = 1) -> dict[str, str]:
    log.info(f"Loading file {fastaid2LCAtaxid_file}.")
    workers = min(workers, MAX_WORKERS)
    if workers == 1 or fastaid2LCAtaxid_file.stat().st_size < MIN_TO_RUN_PARALLEL:
        return read_serial(fastaid2LCAtaxid_file, all_hits)
    else:
        log.info(f"Reading {fastaid2LCAtaxid_file.name} with {workers} workers.")
        return read_parallel(fastaid2LCAtaxid_file, all_hits, workers)


def read_serial(fastaid2LCAtaxid_file: Path, all_hits: set[str]) -> dict[str, str]:
    fastaid2LCAtaxid = {}
    with fastaid2LCAtaxid_file.open(encoding="utf-8", buffering=1024 * 1024) as source:
        for row in source:
            fastaid, LCAtaxid = row.rstrip().split("\t")
            if fastaid in all_hits:
                fastaid2LCAtaxid[fastaid] = LCAtaxid
    return fastaid2LCAtaxid


def read_parallel(fastaid2LCAtaxid_file: Path, all_hits: set[str], workers: int) -> dict[str, str]:
    size = fastaid2LCAtaxid_file.stat().st_size
    path = str(fastaid2LCAtaxid_file)
    hit_filter, hit_mask = _make_hit_filter(all_hits)
    fastaid2LCAtaxid = {}
    max_pending_chunks = workers * 2
    try:
        with ProcessPoolExecutor(
            workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=open_hit_filter,
            initargs=(hit_filter.name, hit_mask),
        ) as pool:
            pending_chunks = {}
            chunk_to_submit = 0
            chunk_to_collect = 0
            total_chunks = (size + CHUNK_BYTES - 1) // CHUNK_BYTES
            while chunk_to_collect < total_chunks:
                while chunk_to_submit < total_chunks and len(pending_chunks) < max_pending_chunks:
                    chunk_start = chunk_to_submit * CHUNK_BYTES
                    pending_chunks[chunk_to_submit] = pool.submit(scan_chunk, path, chunk_start, chunk_start + CHUNK_BYTES)
                    chunk_to_submit += 1
                for row in pending_chunks.pop(chunk_to_collect).result().splitlines():
                    fastaid, LCAtaxid = row.decode("utf-8").split("\t")
                    if fastaid in all_hits:
                        fastaid2LCAtaxid[fastaid] = LCAtaxid
                chunk_to_collect += 1
    finally:
        hit_filter.close()
        hit_filter.unlink()
    return fastaid2LCAtaxid


def _make_hit_filter(all_hits: set[str]) -> tuple[shared_memory.SharedMemory, int]:
    """One CRC32 bit per hit in a power of two size filter, it is however capped at MAX_FILTER_BYTES"""
    bits = 1 << (len(all_hits) * FILTER_BITS_PER_HIT - 1).bit_length()
    bits = min(bits, MAX_FILTER_BYTES * 8)
    hit_mask = bits - 1
    hit_filter = shared_memory.SharedMemory(create=True, size=bits // 8)
    for hit in all_hits:
        bit = crc32(hit.encode("utf-8")) & hit_mask
        hit_filter.buf[bit >> 3] |= 1 << (bit & 7)
    return hit_filter, hit_mask

# Yay, snared memory between workers :)
def open_hit_filter(name: str, hit_mask: int) -> None:
    global _hit_filter, _hit_mask
    _hit_filter = shared_memory.SharedMemory(name=name)
    _hit_mask = hit_mask


def scan_chunk(path: str, start: int, stop: int) -> bytes:
    hit_filter, hit_mask = _hit_filter.buf, _hit_mask
    rows = []
    with open(path, "rb") as file:
        if start:
            file.seek(start - 1)
            file.readline()
        position = file.tell()
        while position < stop:
            row = file.readline()
            if not row:
                break
            position += len(row)
            bit = crc32(row.split(b"\t", 1)[0]) & hit_mask
            if hit_filter[bit >> 3] & (1 << (bit & 7)):
                rows.append(row)
    return b"".join(rows)
