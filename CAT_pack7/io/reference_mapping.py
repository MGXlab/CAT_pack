import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from multiprocessing import shared_memory
from pathlib import Path
from zlib import crc32

log = logging.getLogger("CAT_pack")

MAX_WORKERS = 4
CHUNK_BYTES = 8 * 1024**2
FILTER_BITS_PER_HIT = 32
MAX_FILTER_BYTES = 64 * 1024**2

_hit_filter = None
_hit_mask = 0


def read_fastaid2LCAtaxid(fastaid2LCAtaxid_file: Path, all_hits: set[str], workers: int = 1) -> dict[str, str]:
    workers = min(workers, MAX_WORKERS)
    if workers == 1:
        return read_serial(fastaid2LCAtaxid_file, all_hits)
    else:
        log.info(f"Reading {fastaid2LCAtaxid_file.name} with {workers} workers.")
        return read_parallel(fastaid2LCAtaxid_file, all_hits, workers)


def read_serial(fastaid2LCAtaxid_file: Path, all_hits: set[str]) -> dict[str, str]:
    fastaid2LCAtaxid = {}
    with fastaid2LCAtaxid_file.open(encoding="utf-8") as source:
        for row in source:
            fastaid, LCAtaxid = row.rstrip().split("\t")
            if fastaid in all_hits:
                fastaid2LCAtaxid[fastaid] = LCAtaxid
    return fastaid2LCAtaxid


def read_parallel(fastaid2LCAtaxid_file: Path, all_hits: set[str], workers: int) -> dict[str, str]:
    starts = range(0, fastaid2LCAtaxid_file.stat().st_size, CHUNK_BYTES)
    stops = [start + CHUNK_BYTES for start in starts]
    hit_filter, hit_mask = _make_hit_filter(all_hits)
    fastaid2LCAtaxid = {}
    try:
        with ProcessPoolExecutor(
            workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=open_hit_filter,
            initargs=(hit_filter.name, hit_mask),
        ) as pool:
            # map() yields chunks in file order, so a repeated fastaid keeps its last row.
            for rows in pool.map(scan_chunk, repeat(str(fastaid2LCAtaxid_file)), starts, stops):
                for row in rows:
                    fastaid, LCAtaxid = row.decode("utf-8").rstrip().split("\t")
                    if fastaid in all_hits:
                        fastaid2LCAtaxid[fastaid] = LCAtaxid
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


def scan_chunk(path: str, start: int, stop: int) -> list[bytes]:
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
    return rows
