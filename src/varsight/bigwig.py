"""Read single values from a bigWig file, such as phyloP conservation scores.

The phyloP files are 9 GB each. bigWig is built for random access: a small index says which compressed
block holds a position, so one lookup reads a few kilobytes, with nothing to convert or load first.
Format: Kent et al., Bioinformatics 2010.
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

from .errors import DataError

_BIGWIG_MAGIC = 0x888FFC26
_CHROM_TREE_MAGIC = 0x78CA8C91
_RTREE_MAGIC = 0x2468ACE0
_BLOCK_CACHE = 64


class BigWig:
    def __init__(self, source):
        """source is a path, or any object with seek and read."""
        self._owns = not hasattr(source, "read")
        self._file = open(Path(source), "rb") if self._owns else source
        header = self._read(0, 64)
        magic, _, _, tree_at, _, index_at = struct.unpack_from("<IHHQQQ", header)
        if magic != _BIGWIG_MAGIC:
            raise DataError("this file is not a little-endian bigWig")
        self._compressed = struct.unpack_from("<I", header, 52)[0] > 0
        self._index_at = index_at
        self._chroms = self._read_chroms(tree_at)
        self._blocks: dict[int, bytes] = {}

    def _read(self, offset: int, size: int) -> bytes:
        self._file.seek(offset)
        return self._file.read(size)

    def _read_chroms(self, offset: int) -> dict[str, int]:
        head = self._read(offset, 32)
        magic, _, key_size, value_size = struct.unpack_from("<IIII", head)
        if magic != _CHROM_TREE_MAGIC or value_size != 8:
            raise DataError("this bigWig has an unexpected chromosome index")
        chroms: dict[str, int] = {}
        pending = [offset + 32]
        while pending:
            at = pending.pop()
            is_leaf, _, count = struct.unpack("<BBH", self._read(at, 4))
            at += 4
            item_size = key_size + 8
            items = self._read(at, item_size * count)
            for number in range(count):
                start = number * item_size
                if is_leaf:
                    name = items[start : start + key_size].rstrip(b"\0").decode("ascii")
                    chroms[name] = struct.unpack_from("<I", items, start + key_size)[0]
                else:
                    pending.append(struct.unpack_from("<Q", items, start + key_size)[0])
        return chroms

    def _find_block(self, chrom: int, pos: int) -> tuple[int, int] | None:
        head = self._read(self._index_at, 48)
        if struct.unpack_from("<I", head)[0] != _RTREE_MAGIC:
            raise DataError("this bigWig has an unexpected data index")
        pending = [self._index_at + 48]
        while pending:
            at = pending.pop()
            is_leaf, _, count = struct.unpack("<BBH", self._read(at, 4))
            at += 4
            item_size = 32 if is_leaf else 24
            items = self._read(at, item_size * count)
            for number in range(count):
                start_chrom, start_base, end_chrom, end_base = struct.unpack_from("<IIII", items, number * item_size)
                inside = (start_chrom, start_base) <= (chrom, pos) and (chrom, pos) < (end_chrom, end_base)
                if not inside:
                    continue
                if is_leaf:
                    return struct.unpack_from("<QQ", items, number * item_size + 16)
                pending.append(struct.unpack_from("<Q", items, number * item_size + 16)[0])
        return None

    def _block(self, offset: int, size: int) -> bytes:
        cached = self._blocks.get(offset)
        if cached is not None:
            return cached
        raw = self._read(offset, size)
        data = zlib.decompress(raw) if self._compressed else raw
        if len(self._blocks) >= _BLOCK_CACHE:
            self._blocks.pop(next(iter(self._blocks)))
        self._blocks[offset] = data
        return data

    def has(self, contig: str) -> bool:
        return contig in self._chroms

    def value(self, contig: str, pos: int) -> float | None:
        """The score at a 0-based position, or None where the file has none."""
        chrom = self._chroms.get(contig)
        if chrom is None:
            return None
        found = self._find_block(chrom, pos)
        if found is None:
            return None
        data = self._block(*found)
        _, start, _, step, span, kind, _, count = struct.unpack_from("<IIIIIBBH", data)
        body = 24
        if kind == 3:  # fixed step: every value is at start + step * i
            index = (pos - start) // step
            if 0 <= index < count and (pos - start) % step < span:
                return struct.unpack_from("<f", data, body + index * 4)[0]
            return None
        if kind == 2:  # variable step: each value has its own start
            for number in range(count):
                item_start, value = struct.unpack_from("<If", data, body + number * 8)
                if item_start <= pos < item_start + span:
                    return value
            return None
        for number in range(count):  # bedGraph: each value has a start and an end
            item_start, item_end, value = struct.unpack_from("<IIf", data, body + number * 12)
            if item_start <= pos < item_end:
                return value
        return None

    def close(self) -> None:
        if self._owns:
            self._file.close()
