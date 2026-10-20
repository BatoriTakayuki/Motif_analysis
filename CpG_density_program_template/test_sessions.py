"""Run with: python3 -m unittest -v test_sessions.py"""
from concurrent.futures import ThreadPoolExecutor
import http.client
import importlib.util
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
from urllib.error import HTTPError

from server import Handler, Reference, Sessions, Store, Genomes, GENOME_CATALOG, SESSION_TTL, ThreadingHTTPServer


class QuietHandler(Handler):
    def log_message(self, *args):
        pass


class FakeResponse(io.BytesIO):
    def __init__(self, payload, url):
        super().__init__(payload)
        self.url = url
        self.headers = {"Content-Length": str(len(payload))}

    def geturl(self):
        return self.url

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class SessionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), QuietHandler)
        self.server.allowed_hosts = {"127.0.0.1"}
        self.server.sessions = Sessions(self.temp.name)
        self.server.genomes = Genomes(Path(self.temp.name) / "cache")
        self.server.jobs = threading.BoundedSemaphore(4)
        self.thread = threading.Thread(target=self.server.serve_forever)
        self.thread.start()

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temp.cleanup()

    def request(self, method, path, body=None, cookie=None, headers=None):
        headers = dict(headers or {})
        if cookie:
            headers["Cookie"] = cookie
        if isinstance(body, dict):
            body = json.dumps(body)
            headers["Content-Type"] = "application/json"
        connection = http.client.HTTPConnection(*self.server.server_address, timeout=10)
        try:
            connection.request(method, path, body, headers)
            response = connection.getresponse()
            payload = response.read()
            content = json.loads(payload) if "application/json" in response.getheader("Content-Type", "") else payload
            return response.status, content, dict(response.getheaders())
        finally:
            connection.close()

    def browser(self):
        status, _, headers = self.request("GET", "/")
        self.assertEqual(status, 200)
        self.assertIn("HttpOnly", headers["Set-Cookie"])
        self.assertIn("SameSite=Strict", headers["Set-Cookie"])
        return headers["Set-Cookie"].split(";", 1)[0]

    def upload(self, cookie, sequence=b"ACGCGT"):
        status, result, _ = self.request("POST", "/api/reference", b">chr1\n" + sequence + b"\n", cookie)
        self.assertEqual(status, 200, result)
        return result["id"]

    def test_reference_isolation_and_same_browser_reuse(self):
        a, b = self.browser(), self.browser()
        file_id = self.upload(a)
        data = dict(id=file_id, chrom="chr1", start=0, end=6)
        self.assertEqual(self.request("POST", "/api/reference/region", data, a)[:2], (200, {"sequence": "ACGCGT"}))
        for cookie in (b, None, "motif_session=made-up"):
            self.assertEqual(self.request("POST", "/api/reference/region", data, cookie)[0], 400)
        for route, payload in (("search", dict(id=file_id, query="A" * 20)), ("bed", dict(id=file_id, bed="chr1 0 6"))):
            self.assertEqual(self.request("POST", "/api/reference/" + route, payload, b)[0], 400)

    def test_bigwig_ids_are_isolated_before_file_open(self):
        a, b = self.browser(), self.browser()
        with self.server.sessions.acquire(a) as (_, store):
            store.bigwigs["private-bw"] = Path(self.temp.name) / "private.bw"
        status, result, _ = self.request("POST", "/api/bigwig/analyze", dict(id="private-bw", mode="raw", chrom="chr1", start=0, end=6), b)
        self.assertEqual(status, 400)
        self.assertIn("bigWig", result["error"])

    @unittest.skipUnless(importlib.util.find_spec("pyBigWig"), "pyBigWig is not installed")
    def test_real_bigwig_upload_summary_and_isolation(self):
        import pyBigWig
        path = Path(self.temp.name) / "sample.bw"
        with pyBigWig.open(str(path), "w") as bw:
            bw.addHeader([("chr1", 30)])
            bw.addEntries(["chr1", "chr1"], [0, 10], ends=[10, 20], values=[0.0, 2.0])
        a, b = self.browser(), self.browser()
        status, result, _ = self.request("POST", "/api/bigwig", path.read_bytes(), a)
        self.assertEqual(status, 200, result)
        payload = dict(id=result["id"], chrom="chr1", start=0, end=30)
        for config in ({"mode": "raw"}, {"mode": "bin", "bin_size": 10},
                       {"mode": "sliding", "flank": 4, "step": 5},
                       {"mode": "bed", "bed": "chr1 0 30"}):
            query = dict(payload, **config)
            status, result, _ = self.request("POST", "/api/bigwig/analyze", query, a)
            self.assertEqual(status, 200, result)
            rows = result["rows"]
            if config["mode"] in ("raw", "bin"):
                self.assertEqual([row["mean"] for row in rows], [0.0, 2.0, None])
                self.assertEqual([row["status"] for row in rows], ["ok", "ok", "no_data"])
            elif config["mode"] == "bed":
                self.assertEqual(rows[0]["mean"], 1.0)
                self.assertEqual(rows[0]["minimum"], 0.0)
                self.assertEqual(rows[0]["maximum"], 2.0)
                self.assertAlmostEqual(rows[0]["coverage"], 2 / 3)
            else:
                self.assertEqual(rows[0]["mean"], 0.0)
            self.assertEqual(self.request("POST", "/api/bigwig/analyze", query, b)[0], 400)

    @unittest.skipUnless(importlib.util.find_spec("pyBigWig"), "pyBigWig is not installed")
    def test_chip_atlas_import_validates_metadata_bigwig_and_session(self):
        import pyBigWig
        source = Path(self.temp.name) / "chip-atlas.bw"
        with pyBigWig.open(str(source), "w") as bw:
            bw.addHeader([("chr1", 30)])
            bw.addEntries(["chr1"], [0], ends=[30], values=[2.0])
        metadata = json.dumps([dict(expid="SRX000410", genome="mm10", agClass="Histone",
                                    agSubClass="H3K27ac", clSubClass="C2C12", title="Mouse test")]).encode()
        bigwig = source.read_bytes()
        calls = []
        def fetch(url, timeout):
            calls.append(url)
            if "exp_metadata" in url:
                return FakeResponse(metadata, "https://chip-atlas.org/data/exp_metadata.json?expid=SRX000410")
            return FakeResponse(bigwig, "https://chip-atlas.dbcls.jp/data/mm10/eachData/bw/SRX000410.bw")
        reference = Path(self.temp.name) / "mm10.fa"
        reference.write_text(">chr1\n" + "A" * 30 + "\n")
        self.server.genomes.references["mm10"] = Reference(reference, Path(self.temp.name) / "mm10.seq")
        self.server.genomes.update("mm10", status="ready")
        a, b = self.browser(), self.browser()
        payload = {"accession": "srx000410", "genome": "mm10"}
        with patch("server.urlopen", side_effect=fetch):
            status, result, _ = self.request("POST", "/api/chip-atlas/import", payload, a)
            self.assertEqual(status, 200, result)
            self.assertEqual(result["accession"], "SRX000410")
            self.assertEqual(result["antigen"], "H3K27ac")
            self.assertEqual(len(calls), 2)
            cached = self.request("POST", "/api/chip-atlas/import", payload, a)
            self.assertEqual(cached[0], 200)
            self.assertEqual(cached[1]["id"], result["id"])
            self.assertEqual(len(calls), 3)  # metadata is rechecked; bigWig is reused
        analyze = dict(id=result["id"], mode="bin", chrom="chr1", start=0, end=30, bin_size=30)
        self.assertEqual(self.request("POST", "/api/bigwig/analyze", analyze, a)[1]["rows"][0]["mean"], 2.0)
        self.assertEqual(self.request("POST", "/api/bigwig/analyze", analyze, b)[0], 400)

    def test_chip_atlas_rejects_bad_ids_mismatched_genomes_and_bisulfite(self):
        cookie = self.browser()
        for accession in ("../../etc/passwd", "https://example.com/a.bw", "SRR000001", "SRX1", [], None):
            with patch("server.urlopen", side_effect=AssertionError("network must not be used")):
                self.assertEqual(self.request("POST", "/api/chip-atlas/import", {"accession": accession, "genome": "mm10"}, cookie)[0], 400)
        cases = [
            ([dict(expid="SRX000410", genome="hg38", agClass="Histone")], "利用可能な版: hg38"),
            ([dict(expid="SRX000410", genome="mm10", agClass="Bisulfite-Seq")], "Bisulfite-seq"),
        ]
        for rows, message in cases:
            response = FakeResponse(json.dumps(rows).encode(), "https://chip-atlas.org/data/exp_metadata.json?expid=SRX000410")
            with patch("server.urlopen", return_value=response):
                status, result, _ = self.request("POST", "/api/chip-atlas/import", {"accession": "SRX000410", "genome": "mm10"}, cookie)
            self.assertEqual(status, 400)
            self.assertIn(message, result["error"])

    def test_chip_atlas_rejects_redirect_to_untrusted_host(self):
        cookie = self.browser()
        response = FakeResponse(b"[]", "https://evil.example/metadata.json")
        with patch("server.urlopen", return_value=response):
            status, result, _ = self.request("POST", "/api/chip-atlas/import", {"accession": "SRX000410", "genome": "mm10"}, cookie)
        self.assertEqual(status, 400)
        self.assertIn("転送先", result["error"])

    def test_parallel_browsers_keep_independent_sequences(self):
        cookies = [self.browser() for _ in range(4)]
        def analyze(item):
            index, cookie = item
            sequence = ("ACGT"[index] * 30).encode()
            file_id = self.upload(cookie, sequence)
            status, result, _ = self.request("POST", "/api/reference/region", dict(id=file_id, chrom="chr1", start=0, end=30), cookie)
            return status, result["sequence"], sequence.decode()
        with ThreadPoolExecutor(max_workers=4) as pool:
            for status, actual, expected in pool.map(analyze, enumerate(cookies)):
                self.assertEqual(status, 200)
                self.assertEqual(actual, expected)

    def test_busy_server_rejects_extra_work_and_recovers(self):
        cookie = self.browser()
        for _ in range(4):
            self.server.jobs.acquire()
        try:
            status, _, headers = self.request("POST", "/api/reference", b">chr1\nACGT\n", cookie)
            self.assertEqual(status, 503)
            self.assertEqual(headers["Retry-After"], "3")
            self.assertEqual(self.request("GET", "/api/health", cookie=cookie)[0], 200)
        finally:
            for _ in range(4):
                self.server.jobs.release()
        self.upload(cookie)

    def test_cross_site_post_is_rejected(self):
        cookie = self.browser()
        self.assertEqual(self.request("POST", "/api/reference", b">chr1\nACGT\n", cookie, {"Origin": "http://elsewhere.example"})[0], 403)

    def test_shared_genome_selection_has_session_specific_ids(self):
        source = Path(self.temp.name) / "public.fa"
        sequence = "ACGTGCTAGCTTAGGCTACGATCG"
        source.write_text(">chr1\n" + sequence + "\n")
        genomes = self.server.genomes
        genomes.references["mm10"] = Reference(source, Path(self.temp.name) / "public.seq")
        genomes.update("mm10", status="ready")
        a, b = self.browser(), self.browser()
        self.assertEqual(self.request("GET", "/api/genomes", cookie=a)[1]["default"], "mm10")
        selected = []
        for cookie in (a, b):
            status, result, _ = self.request("POST", "/api/genomes/select", {"assembly": "mm10"}, cookie)
            self.assertEqual(status, 200)
            selected.append(result["id"])
            self.assertEqual(self.request("POST", "/api/genomes/select", {"assembly": "mm10"}, cookie)[1]["id"], result["id"])
            status, matches, _ = self.request("POST", "/api/reference/search", {"id": result["id"], "query": sequence}, cookie)
            self.assertEqual(status, 200)
            self.assertEqual(matches["matches"][0]["start"], 0)
        self.assertNotEqual(*selected)
        self.assertEqual(self.request("POST", "/api/reference/region", {"id": selected[0], "chrom": "chr1", "start": 0, "end": 5}, b)[0], 400)
        for assembly in ("../../etc/passwd", "https://example.com/file.fa.gz", [], None):
            self.assertEqual(self.request("POST", "/api/genomes/select", {"assembly": assembly}, a)[0], 400)

    def test_pending_genome_does_not_block_other_users(self):
        self.server.genomes.update("mm10", status="indexing")
        a = self.browser()
        status, result, _ = self.request("POST", "/api/genomes/select", {"assembly": "mm10"}, a)
        self.assertEqual(status, 202)
        self.assertEqual(result["status"], "indexing")
        self.assertEqual(self.request("GET", "/api/health", cookie=a)[0], 200)
        self.upload(self.browser())

    def test_expiry_removes_files_and_invalidates_old_ids(self):
        cookie = self.browser()
        file_id = self.upload(cookie)
        with self.server.sessions.acquire(cookie) as (_, store):
            directory = store.directory
        store.last_used = time.monotonic() - SESSION_TTL - 1
        status, _, headers = self.request("GET", "/api/health", cookie=cookie)
        self.assertEqual(status, 200)
        self.assertNotEqual(headers["Set-Cookie"].split(";", 1)[0], cookie)
        self.assertFalse(directory.exists())
        self.assertEqual(self.request("POST", "/api/reference/region", dict(id=file_id, chrom="chr1", start=0, end=6), cookie)[0], 400)

    def test_active_request_prevents_expiry(self):
        cookie = self.browser()
        with self.server.sessions.acquire(cookie) as (_, store):
            store.last_used = time.monotonic() - SESSION_TTL - 1
            self.browser()
            self.assertTrue(store.directory.exists())
        self.assertLess(time.monotonic() - store.last_used, 1)

    def test_startup_inputs_get_browser_specific_ids(self):
        shared = Store(self.temp.name)
        source = Path(self.temp.name) / "shared.fa"
        source.write_text(">chr1\nACGCGT\n")
        shared.references["shared-reference"] = Reference(source, Path(self.temp.name) / "shared.seq")
        shared.bigwigs["shared-bw"] = source
        self.server.sessions = Sessions(self.temp.name, shared, {
            "genome": {"id": "shared-reference", "chroms": {"chr1": 6}},
            "bw": {"id": "shared-bw", "chroms": {"chr1": 6}},
        })
        a, b = self.browser(), self.browser()
        data_a = self.request("GET", "/api/health", cookie=a)[1]["startup"]
        data_b = self.request("GET", "/api/health", cookie=b)[1]["startup"]
        for key in ("genome", "bw"):
            self.assertNotEqual(data_a[key]["id"], data_b[key]["id"])
        payload = dict(id=data_a["genome"]["id"], chrom="chr1", start=0, end=6)
        self.assertEqual(self.request("POST", "/api/reference/region", payload, a)[0], 200)
        self.assertEqual(self.request("POST", "/api/reference/region", payload, b)[0], 400)


class GenomeCacheTests(unittest.TestCase):
    def fixture(self, checksum=None, sizes=b"chr1\t24\n"):
        compressed = gzip.compress(b">chr1\nACGTGCTAGCTTAGGCTACGATCG\n")
        checksum = checksum or hashlib.md5(compressed).hexdigest()
        files = {"mm10.fa.gz": compressed, "md5sum.txt": (checksum + "  mm10.fa.gz\n").encode(), "mm10.chrom.sizes": sizes}
        def download(url, timeout):
            self.assertTrue(url.startswith(GENOME_CATALOG["mm10"]["source"]))
            payload = files[url.rsplit("/", 1)[1]]
            response = io.BytesIO(payload)
            response.headers = {"Content-Length": str(len(payload))}
            return response
        return download

    def test_download_verify_prepare_and_restore_after_restart(self):
        with tempfile.TemporaryDirectory() as directory:
            genomes = Genomes(directory)
            with patch("server.urlopen", side_effect=self.fixture()) as fetch:
                genomes.prepare("mm10")
                self.assertEqual(fetch.call_count, 3)
            self.assertEqual(genomes.entries["mm10"]["status"], "ready")
            self.assertEqual(genomes.references["mm10"].read("chr1", 0, 5), "ACGTG")
            restarted = Genomes(directory)
            with patch("server.urlopen", side_effect=AssertionError("cache should avoid network")):
                restarted.prepare("mm10")
            self.assertEqual(restarted.references["mm10"].read("chr1", 0, 5), "ACGTG")

    def test_checksum_mismatch_never_publishes_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            genomes = Genomes(directory)
            with patch("server.urlopen", side_effect=self.fixture("0" * 32)):
                with self.assertRaisesRegex(ValueError, "チェックサム"):
                    genomes.prepare("mm10")
            self.assertNotIn("mm10", genomes.references)
            self.assertFalse((Path(directory) / "mm10" / "mm10.fa.gz").exists())
            self.assertFalse((Path(directory) / "mm10" / "metadata.json").exists())

    def test_wrong_chromosome_sizes_never_publish_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            genomes = Genomes(directory)
            with patch("server.urlopen", side_effect=self.fixture(sizes=b"chr1\t25\n")):
                with self.assertRaisesRegex(ValueError, "染色体長"):
                    genomes.prepare("mm10")
            self.assertNotIn("mm10", genomes.references)
            self.assertFalse((Path(directory) / "mm10" / "reference.seq.part").exists())

    def test_truncated_cache_is_prepared_again(self):
        with tempfile.TemporaryDirectory() as directory:
            genomes = Genomes(directory)
            with patch("server.urlopen", side_effect=self.fixture()):
                genomes.prepare("mm10")
            (Path(directory) / "mm10" / "reference.seq").write_bytes(b"A")
            restarted = Genomes(directory)
            self.assertNotIn("mm10", restarted.references)
            self.assertEqual(restarted.entries["mm10"]["status"], "available")


if __name__ == "__main__":
    unittest.main()
