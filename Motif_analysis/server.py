"""Local DNA / BED / bigWig workbench, with cached public reference genomes."""
from __future__ import annotations

import argparse
import copy
from collections import OrderedDict
from contextlib import contextmanager
import gzip
import hashlib
import heapq
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from http.cookies import SimpleCookie, CookieError
import json
import math
from pathlib import Path
import re
import secrets
import shutil
import tempfile
import threading
import time
import uuid
from urllib.parse import urlsplit
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

ROOT = Path(__file__).resolve().parent
MAX_REGION = 2_000_000
MAX_ROWS = 100_000
MAX_BED_BASES = 50_000_000
MAX_UPLOAD = 20 * 1024**3
MAX_JSON = 20 * 1024**2
SESSION_COOKIE = "motif_session"
SESSION_TTL = 24 * 60 * 60
GENOME_CATALOG = {
    "mm10": {"name": "マウス GRCm38（mm10）", "source": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/"},
    "mm39": {"name": "マウス GRCm39（mm39）", "source": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/"},
}
CHIP_ATLAS_METADATA = "https://chip-atlas.org/data/exp_metadata.json"
CHIP_ATLAS_DATA = "https://chip-atlas.dbcls.jp/data"
CHIP_ATLAS_HOSTS = {"chip-atlas.org", "chip-atlas.dbcls.jp"}
IUPAC = dict(A="A", C="C", G="G", T="T", R="AG", Y="CT", S="CG", W="AT", K="GT", M="AC", B="CGT", D="AGT", H="ACT", V="ACG", N="ACGT")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def integer(value, name, minimum=0, maximum=2**53 - 1):
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError(f"{name}: {minimum}〜{maximum}の整数を指定してください。")
    return value


def chip_atlas_accession(value):
    accession = str(value or "").strip().upper()
    if not re.fullmatch(r"[SED]RX\d{6,12}", accession):
        raise ValueError("ChIP-Atlas IDはSRX・ERX・DRXで始まる実験IDを指定してください（例: SRX000410）。")
    return accession


def trusted_response(response, hosts):
    final = urlsplit(response.geturl())
    if final.scheme != "https" or final.hostname not in hosts or final.username or final.password:
        raise ValueError("ChIP-Atlasから予期しない転送先が返されました。")


def chip_atlas_metadata(accession):
    accession = chip_atlas_accession(accession)
    url = CHIP_ATLAS_METADATA + "?" + urlencode({"expid": accession})
    try:
        with urlopen(url, timeout=60) as response:
            trusted_response(response, CHIP_ATLAS_HOSTS)
            raw = response.read(2 * 1024**2 + 1)
    except HTTPError as exc:
        raise ValueError(f"ChIP-Atlasで{accession}を確認できませんでした（HTTP {exc.code}）。") from exc
    except URLError as exc:
        raise ValueError("ChIP-Atlasへ接続できませんでした。しばらく待って再実行してください。") from exc
    if len(raw) > 2 * 1024**2:
        raise ValueError("ChIP-Atlasのメタデータ応答が大きすぎます。")
    try:
        rows = json.loads(raw)
    except (json.JSONDecodeError, UnicodeError) as exc:
        raise ValueError("ChIP-Atlasのメタデータを読み取れませんでした。") from exc
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise ValueError("ChIP-Atlasのメタデータ形式が不正です。")
    rows = [row for row in rows if row.get("expid") == accession]
    if not rows:
        raise ValueError(f"{accession}はChIP-Atlasに登録されていません。")
    return rows


def import_chip_atlas_bigwig(store, genomes, accession, genome):
    accession = chip_atlas_accession(accession)
    if genome not in GENOME_CATALOG:
        raise ValueError("ChIP-Atlasと対応する参照ゲノム版を選択してください。")
    rows = chip_atlas_metadata(accession)
    available = sorted({str(row.get("genome")) for row in rows if row.get("genome")})
    metadata = next((row for row in rows if row.get("genome") == genome), None)
    if metadata is None:
        versions = "、".join(available) or "なし"
        raise ValueError(f"{accession}の{genome}版bigWigはありません。ChIP-Atlas上の利用可能な版: {versions}")
    if metadata.get("agClass") == "Bisulfite-Seq":
        raise ValueError("Bisulfite-seqはmethylation値とcoverage値が別ファイルです。現在のSRX取得はChIP/ATAC/DNase-seqに対応しています。")
    key = (accession, genome)
    with store.lock:
        cached = store.chip_atlas.get(key)
        if cached and cached["path"].is_file():
            return dict(cached["result"])
    token = uuid.uuid4().hex
    path = store.directory / (token + ".bw")
    partial = store.directory / (token + ".bw.part")
    url = f"{CHIP_ATLAS_DATA}/{genome}/eachData/bw/{accession}.bw"
    try:
        try:
            with urlopen(url, timeout=120) as response, partial.open("wb") as out:
                trusted_response(response, {"chip-atlas.dbcls.jp"})
                total = int(response.headers.get("Content-Length", 0)) or None
                if total is not None and total > MAX_UPLOAD:
                    raise ValueError("ChIP-AtlasのbigWigが20 GiBを超えています。")
                downloaded = 0
                while chunk := response.read(1024**2):
                    downloaded += len(chunk)
                    if downloaded > MAX_UPLOAD:
                        raise ValueError("ChIP-AtlasのbigWigが20 GiBを超えています。")
                    out.write(chunk)
                if total is not None and downloaded != total:
                    raise ValueError("ChIP-AtlasのbigWigダウンロードが中断されました。")
        except HTTPError as exc:
            if exc.code == 404:
                raise ValueError(f"ChIP-Atlasに{accession}の{genome}版bigWigが見つかりません。") from exc
            raise ValueError(f"ChIP-AtlasからbigWigを取得できませんでした（HTTP {exc.code}）。") from exc
        except URLError as exc:
            raise ValueError("ChIP-Atlasへ接続できませんでした。しばらく待って再実行してください。") from exc
        import pyBigWig
        with pyBigWig.open(str(partial)) as bw:
            if not bw.isBigWig():
                raise ValueError("ChIP-Atlasから取得したファイルがbigWig形式ではありません。")
            sizes = bw.chroms()
        reference = genomes.references.get(genome)
        if reference is not None:
            mismatched = [chrom for chrom, size in sizes.items()
                          if chrom not in reference.chroms or reference.chroms[chrom] != size]
            if mismatched:
                raise ValueError(f"ChIP-AtlasのbigWigと参照ゲノムの染色体長が一致しません: {mismatched[0]}")
        partial.replace(path)
        result = dict(id=token, chroms=sizes, accession=accession, genome=genome,
                      title=str(metadata.get("title") or ""), antigen=str(metadata.get("agSubClass") or ""),
                      cell_type=str(metadata.get("clSubClass") or ""), experiment_type=str(metadata.get("agClass") or ""),
                      source_url=url, bytes=path.stat().st_size)
        with store.lock:
            store.bigwigs[token] = path
            store.chip_atlas[key] = dict(path=path, result=result)
        return result
    except Exception:
        partial.unlink(missing_ok=True)
        path.unlink(missing_ok=True)
        raise


def motif_pattern(value):
    if not isinstance(value, str):
        raise ValueError("モチーフを文字列で指定してください。")
    motif = re.sub(r"\s", "", value).upper()
    if not motif or len(motif) > 100 or set(motif) - IUPAC.keys():
        raise ValueError("モチーフはIUPAC塩基記号で1〜100文字にしてください。")
    return motif


def iter_motif_starts(sequence, motif, both=False):
    """Merge strands without retaining millions of matches in a Python set."""
    patterns = [motif]
    reverse = motif.translate(COMPLEMENT)[::-1]
    if both and reverse != motif:
        patterns.append(reverse)
    streams = []
    for pattern in patterns:
        expression = "(?=(" + "".join("[" + IUPAC[b] + "]" for b in pattern) + "))"
        streams.append((match.start() for match in re.finditer(expression, sequence)))
    previous = None
    for start in heapq.merge(*streams):
        if start != previous:
            yield start
            previous = start


def motif_starts(sequence, motif, both=False):
    return list(iter_motif_starts(sequence, motif, both))


def parse_bed(text):
    rows = []
    for number, raw in enumerate(str(text).splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "track ", "browser ")):
            continue
        fields = line.split()
        if len(fields) < 3 or not re.fullmatch(r"\d+", fields[1]) or not re.fullmatch(r"\d+", fields[2]):
            raise ValueError(f"BED {number}行目: 染色体名・開始・終了の3列が必要です。")
        start, end = int(fields[1]), int(fields[2])
        if not 0 <= start < end <= 2**53 - 1:
            raise ValueError(f"BED {number}行目: 0 ≤ start < end を満たしていません。")
        rows.append(dict(chrom=fields[0], start=start, end=end, name=fields[3] if len(fields) > 3 else f"region_{len(rows)+1}"))
        if len(rows) > MAX_ROWS:
            raise ValueError(f"BEDは{MAX_ROWS:,}領域までです。")
    if not rows:
        raise ValueError("BEDに有効な領域がありません。")
    return rows


class Reference:
    """Normalize FASTA into a disk-backed, seekable sequence (one pass)."""
    def __init__(self, source, destination):
        self.path = Path(destination)
        self.chroms = OrderedDict()
        self.offsets = {}
        opener = gzip.open if Path(source).suffix.lower() == ".gz" else open
        current = None
        offset = 0
        with opener(source, "rb") as inp, self.path.open("wb") as out:
            while True:
                raw = inp.readline(1_048_577)
                if not raw:
                    break
                if len(raw) > 1_048_576:
                    raise ValueError("FASTAは1行1,048,576文字以下に改行してください。")
                if raw.startswith(b">"):
                    if current is not None and not self.chroms[current]:
                        raise ValueError(f"空の配列です: {current}")
                    header = raw[1:].strip().split()
                    name = header[0].decode("ascii") if header else ""
                    if not name or name in self.chroms:
                        raise ValueError("FASTAの配列名が空か、重複しています。")
                    current = name
                    self.chroms[current] = 0
                    self.offsets[current] = offset
                else:
                    chunk = raw.strip().upper()
                    if not chunk:
                        continue
                    if current is None:
                        raise ValueError("FASTAヘッダー（>染色体名）が必要です。")
                    if re.search(rb"[^ACGTRYSWKMBDHVN]", chunk):
                        raise ValueError(f"FASTAに不正な塩基文字があります: {current}")
                    out.write(chunk)
                    self.chroms[current] += len(chunk)
                    offset += len(chunk)
                    if offset > MAX_UPLOAD:
                        raise ValueError("参照配列は展開後20 GiB以下にしてください。")
        if current is None or not self.chroms[current]:
            raise ValueError("FASTAに空でない配列が必要です。")

    def read(self, chrom, start, end):
        if chrom not in self.chroms or not 0 <= start < end <= self.chroms[chrom]:
            raise ValueError("領域が参照配列の範囲外です。染色体名と座標を確認してください。")
        with self.path.open("rb") as inp:
            inp.seek(self.offsets[chrom] + start)
            return inp.read(end - start).decode("ascii")

    @classmethod
    def restore(cls, path, chroms):
        if not chroms or any(not isinstance(name, str) or not isinstance(size, int) or size < 1 for name, size in chroms.items()):
            raise ValueError("参照ゲノムの索引が不正です。")
        reference = cls.__new__(cls)
        reference.path = Path(path)
        reference.chroms = OrderedDict(chroms)
        reference.offsets = {}
        offset = 0
        for name, size in reference.chroms.items():
            reference.offsets[name] = offset
            offset += size
        if reference.path.stat().st_size != offset:
            raise ValueError("参照ゲノムのサイズが索引と一致しません。")
        return reference


class Genomes:
    """One persistent, read-only reference per assembly, shared across sessions."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.lock = threading.Lock()
        self.download_lock = threading.Lock()
        self.references = {}
        self.entries = {}
        for assembly, info in GENOME_CATALOG.items():
            self.entries[assembly] = dict(assembly=assembly, **info, status="available", downloaded=0, total=None)
            try:
                folder = self.directory / assembly
                metadata = json.loads((folder / "metadata.json").read_text())
                if metadata.get("assembly") != assembly or metadata.get("schema") != 1:
                    raise ValueError("参照ゲノムの索引が不正です。")
                self.references[assembly] = Reference.restore(folder / "reference.seq", metadata["chroms"])
                self.entries[assembly]["status"] = "ready"
            except (OSError, ValueError, KeyError, TypeError):
                pass

    def catalog(self):
        with self.lock:
            return [dict(entry) for entry in self.entries.values()]

    def update(self, assembly, **values):
        with self.lock:
            self.entries[assembly].update(values)

    def start(self, assembly):
        if not isinstance(assembly, str) or assembly not in self.entries:
            raise ValueError("対応するゲノム版を選択してください。")
        with self.lock:
            entry = self.entries[assembly]
            if entry["status"] in ("available", "error"):
                entry.update(status="queued", downloaded=0, total=None, error=None)
                threading.Thread(target=self.worker, args=(assembly,), daemon=True).start()
            return dict(entry)

    def worker(self, assembly):
        try:
            with self.download_lock:
                self.prepare(assembly)
        except Exception as exc:
            self.update(assembly, status="error", error=str(exc))
            print(f"参照ゲノム {assembly} の準備に失敗: {exc}", flush=True)

    def prepare(self, assembly):
        if assembly not in GENOME_CATALOG:
            raise ValueError("対応するゲノム版を選択してください。")
        if assembly in self.references:
            return
        info = GENOME_CATALOG[assembly]
        folder = self.directory / assembly
        folder.mkdir(exist_ok=True)
        filename = assembly + ".fa.gz"
        archive = folder / filename
        part = folder / (filename + ".part")
        self.update(assembly, status="downloading")
        # Verify against the publisher's checksum before normalizing the FASTA.
        with urlopen(info["source"] + "md5sum.txt", timeout=60) as response:
            checksums = response.read(1024 * 1024).decode("ascii")
        expected = next((fields[0].lower() for line in checksums.splitlines()
                         if len(fields := line.split()) == 2 and fields[1].lstrip("*") == filename), None)
        if expected is None or not re.fullmatch(r"[0-9a-f]{32}", expected):
            raise ValueError("配布元のチェックサムを確認できませんでした。")
        try:
            if not archive.exists():
                with urlopen(info["source"] + filename, timeout=120) as response, part.open("wb") as out:
                    total = int(response.headers.get("Content-Length", 0)) or None
                    downloaded = 0
                    self.update(assembly, total=total)
                    while chunk := response.read(1024**2):
                        downloaded += len(chunk)
                        if downloaded > MAX_UPLOAD:
                            raise ValueError("参照ゲノムが20 GiBを超えています。")
                        out.write(chunk)
                        self.update(assembly, downloaded=downloaded)
                    if total is not None and downloaded != total:
                        raise ValueError("参照ゲノムのダウンロードが中断されました。")
                part.replace(archive)
            self.update(assembly, status="verifying")
            digest = hashlib.md5()
            with archive.open("rb") as inp:
                for chunk in iter(lambda: inp.read(1024**2), b""):
                    digest.update(chunk)
            if digest.hexdigest() != expected:
                archive.unlink()
                raise ValueError("参照ゲノムのチェックサムが一致しません。もう一度実行してください。")
            with urlopen(info["source"] + assembly + ".chrom.sizes", timeout=60) as response:
                sizes = dict((name, int(size)) for name, size in
                             (line.split() for line in response.read(1024 * 1024).decode("ascii").splitlines() if line.strip()))
            self.update(assembly, status="indexing")
            temporary = folder / "reference.seq.part"
            reference = Reference(archive, temporary)
            if reference.chroms != sizes:
                raise ValueError("参照ゲノムの染色体長が配布元の情報と一致しません。")
            temporary.replace(folder / "reference.seq")
            reference.path = folder / "reference.seq"
            manifest = folder / "metadata.json.part"
            manifest.write_text(json.dumps(dict(schema=1, assembly=assembly, source=info["source"], md5=expected, chroms=reference.chroms)))
            manifest.replace(folder / "metadata.json")
            with self.lock:
                self.references[assembly] = reference
                self.entries[assembly].update(status="ready", error=None)
        finally:
            part.unlink(missing_ok=True)
            (folder / "reference.seq.part").unlink(missing_ok=True)
            (folder / "metadata.json.part").unlink(missing_ok=True)


def sequence_bed(reference, rows, motif, both):
    motif = motif_pattern(motif)
    if sum(row["end"] - row["start"] for row in rows) > MAX_BED_BASES:
        raise ValueError("一度に解析するBED領域の合計長は50 Mb以下にしてください。")
    result = []
    for row in rows:
        row = dict(row, count=None, density=None, status="ok")
        chrom, start, end = row["chrom"], row["start"], row["end"]
        if chrom not in reference.chroms:
            row["status"] = "chromosome_not_found"
        elif end > reference.chroms[chrom]:
            row["status"] = "outside_reference"
        else:
            # Match starts inside BED; allow following bases beyond BED end.
            seq = reference.read(chrom, start, min(reference.chroms[chrom], end + len(motif) - 1))
            if re.search("[^ACGT]", seq[:end-start]):
                row["status"] = "ambiguous_bases"
            else:
                row["count"] = sum(1 for position in iter_motif_starts(seq, motif, both) if position < end-start)
                row["density"] = row["count"] / (end-start)
        result.append(row)
    return result


def search_reference(reference, query, limit=50, chunk_size=1024**2, max_edits=0):
    """Find exact query/RC matches on disk, including across chunk boundaries.

    Full-genome record names have origin zero. A coordinate-qualified FASTA
    record, e.g. chr12:101-200, maps its local bases back to genomic coordinates.
    Never return a partial list as a unique match when the hit limit is reached.
    """
    if not isinstance(query, str):
        raise ValueError("検索配列を文字列で指定してください。")
    query = re.sub(r"\s", "", query).upper()
    if not 20 <= len(query) <= MAX_REGION or re.search("[^ACGT]", query):
        raise ValueError("位置検索にはACGTのみの20 bp〜2 Mbの配列を指定してください。")
    targets = [(query, "+")]
    reverse = query.translate(COMPLEMENT)[::-1]
    if reverse != query:
        targets.append((reverse, "-"))
    targets = [(value.encode("ascii"), strand) for value, strand in targets]
    matches = []
    for record, length in reference.chroms.items():
        qualified = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", record)
        chrom, origin = record, 0
        if qualified:
            chrom = qualified[1]
            origin, stop = int(qualified[2])-1, int(qualified[3])
            if origin < 0 or stop-origin != length:
                raise ValueError(f"座標付きFASTAの長さがヘッダーと一致しません: {record}")
        if length < len(query):
            continue
        with reference.path.open("rb") as inp:
            inp.seek(reference.offsets[record])
            consumed, tail = 0, b""
            while consumed < length:
                chunk = inp.read(min(chunk_size, length-consumed))
                if not chunk:
                    raise ValueError("検索用参照の読み込みが中断されました。")
                segment = tail+chunk
                local_origin = consumed-len(tail)
                for target, strand in targets:
                    position = segment.find(target)
                    while position >= 0:
                        local = local_origin+position
                        if local+len(query) > consumed:
                            matches.append(dict(chrom=chrom, start=origin+local, end=origin+local+len(query), strand=strand, record=record, local_start=local, local_end=local+len(query), edits=0, difference="exact"))
                            if len(matches) > limit:
                                return dict(matches=matches[:limit], truncated=True)
                        position = segment.find(target, position+1)
                consumed += len(chunk)
                tail = segment[-(len(query)-1):]
    if not matches and max_edits:
        return search_one_edit(reference, query, limit, chunk_size)
    matches.sort(key=lambda row: (row["chrom"], row["start"], row["strand"]))
    return dict(matches=matches, truncated=False)


def one_edit_kind(query, sequence):
    """Return the single edit relative to the forward genomic reference."""
    if abs(len(query)-len(sequence)) > 1:
        return None
    if len(query) == len(sequence):
        return "substitution" if sum(a != b for a, b in zip(query, sequence)) == 1 else None
    index = 0
    while index < min(len(query), len(sequence)) and query[index] == sequence[index]:
        index += 1
    if len(sequence) > len(query):
        return "query_deletion" if query[index:] == sequence[index+1:] else None
    return "query_insertion" if query[index+1:] == sequence[index:] else None


def search_one_edit(reference, query, limit=50, chunk_size=1024**2):
    """Seed and verify exactly one edit, without a genome-sized in-memory copy."""
    targets = [(query, "+")]
    reverse = query.translate(COMPLEMENT)[::-1]
    if reverse != query:
        targets.append((reverse, "-"))
    width = min(64, len(query)//2)
    matches, seen = [], set()
    for record, length in reference.chroms.items():
        qualified = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", record)
        chrom, origin = (qualified[1], int(qualified[2])-1) if qualified else (record, 0)
        candidates = set()
        seeds = [(target[:width].encode(), 0, strand) for target, strand in targets]
        seeds += [(target[-width:].encode(), len(target)-width, strand) for target, strand in targets]
        with reference.path.open("rb") as inp:
            inp.seek(reference.offsets[record])
            consumed, tail = 0, b""
            while consumed < length:
                chunk = inp.read(min(chunk_size, length-consumed))
                if not chunk:
                    raise ValueError("検索用参照を読み込めません。")
                segment, local_origin = tail+chunk, consumed-len(tail)
                for seed, offset, strand in seeds:
                    position = segment.find(seed)
                    while position >= 0:
                        if local_origin+position+len(seed) > consumed:
                            predicted = local_origin+position-offset
                            for shift in (-1, 0, 1):
                                if 0 <= predicted+shift < length:
                                    candidates.add((predicted+shift, strand))
                        if len(candidates) > 5000:
                            raise ValueError("検索配列の端が繰り返し配列のため候補が多すぎます。特徴的な領域を含む配列を使ってください。")
                        position = segment.find(seed, position+1)
                consumed += len(chunk)
                tail = segment[-(width-1):]
        for start, strand in sorted(candidates):
            target = query if strand == "+" else reverse
            for size in (len(query), len(query)-1, len(query)+1):
                end = start+size
                if end > length:
                    continue
                kind = one_edit_kind(target, reference.read(record, start, end))
                key = (chrom, origin+start, origin+end, strand)
                if kind is not None and key not in seen:
                    seen.add(key)
                    matches.append(dict(chrom=chrom, start=origin+start, end=origin+end, strand=strand, record=record, local_start=start, local_end=end, edits=1, difference=kind))
                    if len(matches) > limit:
                        return dict(matches=matches[:limit], truncated=True)
    matches.sort(key=lambda row: (row["chrom"], row["start"], row["strand"]))
    return dict(matches=matches, truncated=False)


def window_rows(chrom, start, end, flank, step):
    count = max(0, (end - start - 2*flank - 1) // step + 1)
    if count < 1:
        raise ValueError("選択領域がウィンドウ幅より短すぎます。")
    if count > MAX_ROWS:
        raise ValueError("ウィンドウが100,000件を超えます。ステップ幅を大きくしてください。")
    return [dict(chrom=chrom, start=i-flank, end=i+flank+1, name=f"window_{j+1}", center=i)
            for j, i in enumerate(range(start+flank, end-flank, step))]


def bin_rows(chrom, start, end, bin_size):
    if bin_size < 1 or end <= start:
        raise ValueError("Bin幅と対象領域を確認してください。")
    if (end - start + bin_size - 1) // bin_size > MAX_ROWS:
        raise ValueError("Binが100,000件を超えます。Bin幅を大きくしてください。")
    return [dict(chrom=chrom, start=i, end=min(i+bin_size, end),
                 name=f"bin_{j+1}", center=(i+min(i+bin_size, end)-1)//2)
            for j, i in enumerate(range(start, end, bin_size))]


def bigwig_summary(bw, rows):
    sizes = bw.chroms()
    result = []
    for row in rows:
        row = dict(row, mean=None, minimum=None, maximum=None, coverage=None, status="ok")
        chrom, start, end = row["chrom"], row["start"], row["end"]
        if chrom not in sizes:
            row["status"] = "chromosome_not_found"
        elif end > sizes[chrom]:
            row["status"] = "outside_reference"
        else:
            values = {}
            for column, kind in (("mean", "mean"), ("minimum", "min"), ("maximum", "max"), ("coverage", "coverage")):
                value = bw.stats(chrom, start, end, type=kind, exact=True)[0]
                values[column] = value if value is not None and math.isfinite(value) else None
            row.update(values)
            if row["mean"] is None:
                row["status"] = "no_data"
        result.append(row)
    return result


def bigwig_stored_values(bw, chrom, start, end):
    """Return stored interval values, without smoothing; gaps remain missing."""
    if end-start > MAX_REGION:
        raise ValueError("元のbigWig値の表示領域は2 Mb以下にしてください。")
    sizes = bw.chroms()
    if chrom not in sizes or end > sizes[chrom]:
        return [dict(chrom=chrom,start=start,end=end,name="stored_1",mean=None,minimum=None,maximum=None,coverage=None,status="chromosome_not_found" if chrom not in sizes else "outside_reference")]
    rows = []
    def append(left, right, value):
        if left >= right:
            return
        value = value if value is not None and math.isfinite(value) else None
        if rows and rows[-1]["end"] == left and rows[-1]["mean"] == value:
            rows[-1]["end"] = right
            return
        if len(rows) >= MAX_ROWS:
            raise ValueError("bigWigの区間が100,000件を超えます。表示領域を狭くするかSliding windowを選んでください。")
        rows.append(dict(chrom=chrom,start=left,end=right,name=f"stored_{len(rows)+1}",mean=value,minimum=value,maximum=value,coverage=1.0 if value is not None else None,status="ok" if value is not None else "no_data"))
    for chunk_start in range(start,end,100_000):
        chunk_end = min(end,chunk_start+100_000)
        current = chunk_start
        for left,right,value in bw.intervals(chrom,chunk_start,chunk_end) or ():
            left,right = max(left,chunk_start),min(right,chunk_end)
            append(current,left,None)
            append(left,right,value)
            current = right
        append(current,chunk_end,None)
    return rows


class Store:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.references = {}
        self.bigwigs = {}
        self.lock = threading.Lock()
        self.startup = {}
        self.last_used = time.monotonic()
        self.active = 0
        self.genomes = {}
        self.chip_atlas = {}


class Sessions:
    """Browser-scoped file stores; active requests keep their files alive."""
    def __init__(self, directory, shared_store=None, startup=None):
        self.directory = Path(directory)
        self.shared_store = shared_store
        self.startup = startup or {}
        self.stores = {}
        self.lock = threading.Lock()

    @contextmanager
    def acquire(self, cookie_header):
        cookie = SimpleCookie()
        try:
            if len(cookie_header) <= 4096:
                cookie.load(cookie_header)
        except CookieError:
            cookie.clear()
        value = cookie.get(SESSION_COOKIE)
        token = value.value if value else None
        with self.lock:
            now = time.monotonic()
            for expired, store in list(self.stores.items()):
                if not store.active and now - store.last_used >= SESSION_TTL:
                    shutil.rmtree(store.directory)
                    del self.stores[expired]
            store = self.stores.get(token)
            if store is None:
                token = secrets.token_hex(32)
                directory = self.directory / token
                directory.mkdir()
                store = Store(directory)
                # CLI-loaded inputs are intentional shared, read-only templates.
                # Each browser receives its own file IDs, never another session's.
                store.startup = copy.deepcopy(self.startup)
                for key, collection in (("genome", "references"), ("bw", "bigwigs")):
                    if key in store.startup:
                        entry = store.startup[key]
                        new_id = uuid.uuid4().hex
                        getattr(store, collection)[new_id] = getattr(self.shared_store, collection)[entry["id"]]
                        entry["id"] = new_id
                self.stores[token] = store
            store.active += 1
        try:
            yield token, store
        finally:
            with self.lock:
                store.active -= 1
                store.last_used = time.monotonic()


class Handler(BaseHTTPRequestHandler):
    server_version = "MotifAnalysis/1.0"

    def setup(self):
        super().setup()
        self.connection.settimeout(120)

    def end_headers(self):
        if getattr(self, "session_token", None):
            self.send_header("Set-Cookie", f"{SESSION_COOKIE}={self.session_token}; Path=/; HttpOnly; SameSite=Strict; Max-Age={SESSION_TTL}")
        super().end_headers()

    def allowed_host(self):
        """Reject DNS rebinding while keeping localhost / explicit LAN use available."""
        try:
            host = urlsplit("//" + self.headers.get("Host", ""))
            return (host.hostname in self.server.allowed_hosts and
                    host.username is None and host.password is None and
                    not host.path and not host.query and not host.fragment and
                    (host.port is None or host.port > 0))
        except ValueError:
            return False

    def send_json(self, value, status=200):
        body = json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        if status == 503:
            self.send_header("Retry-After", "3")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if not self.allowed_host():
            return self.send_json(dict(error="このホスト名ではアクセスできません。localhostを使用してください。"), 403)
        with self.server.sessions.acquire(self.headers.get("Cookie", "")) as (token, store):
            self.session_token, self.store = token, store
            return self.get()

    def get(self):
        path = urlsplit(self.path).path
        if path == "/api/genomes":
            return self.send_json(dict(genomes=self.server.genomes.catalog(), default="mm10"))
        if path == "/api/health":
            try:
                import pyBigWig
                available = True
            except ImportError:
                available = False
            health = dict(bigwig=available)
            if self.store.startup:
                health["startup"] = self.store.startup
            return self.send_json(health)
        files = {"/": "index.html", "/index.html": "index.html", "/app.js": "app.js", "/core.js": "core.js", "/style.css": "style.css"}
        if path not in files:
            return self.send_json(dict(error="Not found"), 404)
        file = ROOT / files[path]
        body = file.read_bytes()
        mime = {".html": "text/html", ".js": "text/javascript", ".css": "text/css"}[file.suffix]
        self.send_response(200)
        self.send_header("Content-Type", mime + "; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; img-src 'self' data: blob:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        # No cross-origin writes or requests from third-party web pages.
        if not self.allowed_host():
            return self.send_json(dict(error="このホスト名ではアクセスできません。localhostを使用してください。"), 403)
        origin = self.headers.get("Origin")
        if (origin and origin != "http://" + self.headers.get("Host", "")) or self.headers.get("Sec-Fetch-Site") == "cross-site":
            return self.send_json(dict(error="別のサイトからのリクエストは受け付けません。"), 403)
        if not self.server.jobs.acquire(blocking=False):
            return self.send_json(dict(error="サーバーが混み合っています。少し待ってから再実行してください。"), 503)
        try:
            with self.server.sessions.acquire(self.headers.get("Cookie", "")) as (token, store):
                self.session_token, self.store = token, store
                return self.post()
        finally:
            self.server.jobs.release()

    def post(self):
        try:
            path = urlsplit(self.path).path
            size = int(self.headers.get("Content-Length", "0"))
            if size <= 0:
                raise ValueError("データが空です。")
            if self.headers.get("Transfer-Encoding"):
                raise ValueError("Content-Lengthを指定して送信してください。")
            if path in ("/api/reference", "/api/bigwig"):
                if size > MAX_UPLOAD:
                    raise ValueError("ファイルは20 GiB以下にしてください。")
                return self.upload(path, size)
            if size > MAX_JSON:
                raise ValueError("リクエストが大きすぎます。")
            raw = self.rfile.read(size)
            if len(raw) != size:
                raise ValueError("リクエストが中断されました。")
            data = json.loads(raw)
            if not isinstance(data, dict):
                raise ValueError("JSONオブジェクトを送信してください。")
            if "id" in data and not isinstance(data["id"], str):
                raise ValueError("読み込んだファイルのIDが不正です。")
            if "both" in data and not isinstance(data["both"], bool):
                raise ValueError("両鎖検索はtrueまたはfalseで指定してください。")
            return self.action(path, data)
        except (ValueError, KeyError, UnicodeError, OSError, EOFError, RuntimeError) as exc:
            self.send_json(dict(error=str(exc)), 400)
        except Exception as exc:
            self.log_error("%s", exc)
            self.send_json(dict(error="処理に失敗しました。サーバーのログを確認してください。"), 500)

    def upload(self, route, size):
        token = uuid.uuid4().hex
        suffix = ".gz" if self.headers.get("X-File-Gzip") == "1" else ".input"
        path = self.store.directory / (token + suffix)
        normalized = self.store.directory / (token + ".seq")
        try:
            with path.open("wb") as out:
                remaining = size
                while remaining:
                    chunk = self.rfile.read(min(1024**2, remaining))
                    if not chunk:
                        raise ValueError("アップロードが中断されました。")
                    out.write(chunk)
                    remaining -= len(chunk)
            if route == "/api/reference":
                reference = Reference(path, normalized)
                with self.store.lock:
                    self.store.references[token] = reference
                path.unlink()
                return self.send_json(dict(id=token, chroms=reference.chroms))
            try:
                import pyBigWig
            except ImportError:
                raise ValueError("bigWig用依存関係がありません。pip install -r requirements.txt を実行してください。")
            with pyBigWig.open(str(path)) as bw:
                if not bw.isBigWig():
                    raise ValueError("bigWig形式のファイルを選択してください。")
                sizes = bw.chroms()
            with self.store.lock:
                self.store.bigwigs[token] = path
            return self.send_json(dict(id=token, chroms=sizes))
        except Exception:
            with self.store.lock:
                self.store.references.pop(token, None)
                self.store.bigwigs.pop(token, None)
            path.unlink(missing_ok=True)
            normalized.unlink(missing_ok=True)
            raise

    def action(self, route, data):
        store = self.store
        if route == "/api/chip-atlas/import":
            result = import_chip_atlas_bigwig(store, self.server.genomes, data.get("accession"), data.get("genome"))
            return self.send_json(result)
        if route == "/api/genomes/select":
            assembly = data.get("assembly")
            entry = self.server.genomes.start(assembly)
            if entry["status"] != "ready":
                return self.send_json(entry, 202)
            reference = self.server.genomes.references[assembly]
            with store.lock:
                token = store.genomes.get(assembly)
                if token is None:
                    token = uuid.uuid4().hex
                    store.genomes[assembly] = token
                    store.references[token] = reference
            return self.send_json(dict(entry, id=token, chroms=reference.chroms))
        if route == "/api/reference/search":
            ref = store.references.get(data.get("id"))
            if ref is None:
                raise ValueError("検索用の参照ゲノムFASTAを読み込んでください。")
            edits = integer(data.get("max_edits", 0), "許容する塩基差", 0, 1)
            return self.send_json(search_reference(ref, data.get("query"), max_edits=edits))
        if route == "/api/reference/region":
            ref = store.references.get(data.get("id"))
            if ref is None:
                raise ValueError("参照FASTAを読み込み直してください。")
            start = integer(data.get("start"), "開始位置")
            end = integer(data.get("end"), "終了位置", start+1)
            if end-start > MAX_REGION:
                raise ValueError("表示領域は2 Mb以下にしてください。BED集計は全参照を対象にできます。")
            return self.send_json(dict(sequence=ref.read(str(data.get("chrom")), start, end)))
        if route == "/api/reference/bed":
            ref = store.references.get(data.get("id"))
            if ref is None:
                raise ValueError("参照FASTAを読み込み直してください。")
            rows = sequence_bed(ref, parse_bed(data.get("bed", "")), data.get("motif", "CG"), bool(data.get("both")))
            return self.send_json(dict(rows=rows))
        if route == "/api/bigwig/analyze":
            path = store.bigwigs.get(data.get("id"))
            if path is None:
                raise ValueError("bigWigを読み込み直してください。")
            import pyBigWig
            mode = data.get("mode")
            if mode == "raw":
                start = integer(data.get("start"), "開始位置")
                end = integer(data.get("end"), "終了位置", start+1)
                with pyBigWig.open(str(path)) as bw:
                    return self.send_json(dict(rows=bigwig_stored_values(bw,str(data.get("chrom")),start,end)))
            if mode == "bed":
                rows = parse_bed(data.get("bed", ""))
            elif mode == "bin":
                start = integer(data.get("start"), "開始位置")
                end = integer(data.get("end"), "終了位置", start+1)
                bin_size = integer(data.get("bin_size"), "Bin幅", 1, MAX_REGION)
                rows = bin_rows(str(data.get("chrom")), start, end, bin_size)
            elif mode == "sliding":
                start = integer(data.get("start"), "開始位置")
                end = integer(data.get("end"), "終了位置", start+1)
                flank = integer(data.get("flank"), "片側幅", 1, MAX_REGION)
                step = integer(data.get("step"), "ステップ幅", 1, MAX_REGION)
                rows = window_rows(str(data.get("chrom")), start, end, flank, step)
            else:
                raise ValueError("解析方法が不正です。")
            with pyBigWig.open(str(path)) as bw:
                rows = bigwig_summary(bw, rows)
            return self.send_json(dict(rows=rows))
        raise ValueError("不明なAPIです。")


def main():
    parser = argparse.ArgumentParser(description="Motif Analysis — ローカル配列・bigWig解析ツール")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--allowed-host", action="append", default=[], help="接続を許可する追加のホスト名・IP（Docker等で使用、複数指定可）")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--max-jobs", type=int, default=4, help="サーバー側で同時処理する最大リクエスト数（既定: 4）")
    parser.add_argument("--search-genome", type=Path, help="起動時に読み込む検索用参照FASTA")
    parser.add_argument("--query", type=Path, help="起動時に画面へ入力する単一配列FASTA")
    parser.add_argument("--bigwig", type=Path, help="起動時に読み込むbigWig")
    parser.add_argument("--genome-cache", type=Path, default=ROOT / "genome_cache", help="選択式参照ゲノムの保存先")
    parser.add_argument("--prepare-genome", action="append", choices=GENOME_CATALOG, help="参照ゲノムを事前取得して終了（複数指定可）")
    args = parser.parse_args()
    if not 1 <= args.max_jobs <= 32:
        parser.error("--max-jobsは1〜32で指定してください。")
    genomes = Genomes(args.genome_cache)
    if args.prepare_genome:
        for assembly in args.prepare_genome:
            genomes.prepare(assembly)
            print(f"{GENOME_CATALOG[assembly]['name']} の準備ができました。", flush=True)
        return
    with tempfile.TemporaryDirectory(prefix="motif-analysis-") as directory:
        server = ThreadingHTTPServer((args.host, args.port), Handler)
        server.allowed_hosts = {"localhost", "127.0.0.1", "::1"}
        if args.host not in ("0.0.0.0", "::"):
            server.allowed_hosts.add(args.host.lower())
        server.allowed_hosts.update(host.lower() for host in args.allowed_host)
        server.store = Store(directory)
        server.genomes = genomes
        server.startup = {}
        if args.search_genome:
            token = uuid.uuid4().hex
            reference = Reference(args.search_genome, Path(directory)/(token+".seq"))
            server.store.references[token] = reference
            server.startup["genome"] = dict(id=token, chroms=reference.chroms, name=args.search_genome.name)
        if args.query:
            if args.query.stat().st_size > 5*1024**2:
                parser.error("起動時の検索配列は5 MB以下にしてください。")
            server.startup["query"] = args.query.read_text(encoding="utf-8-sig")
        if args.bigwig:
            import pyBigWig
            token = uuid.uuid4().hex
            with pyBigWig.open(str(args.bigwig)) as bw:
                if not bw.isBigWig():
                    parser.error("--bigwigにはbigWigを指定してください。")
                sizes = bw.chroms()
            server.store.bigwigs[token] = args.bigwig
            server.startup["bw"] = dict(id=token, chroms=sizes, name=args.bigwig.name)
        server.sessions = Sessions(directory, server.store, server.startup)
        server.jobs = threading.BoundedSemaphore(args.max_jobs)
        display_host = args.host if args.host not in ("0.0.0.0", "::") else "localhost"
        print(f"Motif Analysis: ブラウザで http://{display_host}:{args.port} を開いてください。終了: Ctrl+C", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    main()
