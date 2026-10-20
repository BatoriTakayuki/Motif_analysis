/* Shared browser/Node sequence analysis. Genomic coordinates are 0-based, half-open. */
(function (root, factory) {
  "use strict";
  var api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  if (root) root.MotifCore = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  var LIMITS = Object.freeze({ sequence: 2000000, motif: 100, rows: 100000, guides: 100, guideLength: 200, guideHits: 20000 });
  var MASK = Object.freeze({ A: 1, C: 2, G: 4, T: 8, R: 5, Y: 10, S: 6, W: 9, K: 12, M: 3, B: 14, D: 13, H: 11, V: 7, N: 15 });
  var COMPLEMENT = Object.freeze({ A: "T", C: "G", G: "C", T: "A", R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K", B: "V", D: "H", H: "D", V: "B", N: "N" });

  function fail(message) { throw new Error(message); }

  function validateSequence(sequence) {
    if (!sequence.length) fail("参照配列を入力してください。");
    if (sequence.length > LIMITS.sequence) fail("参照配列は2,000,000 bp以内にしてください。全ゲノムbigWigは.bw入力をご利用ください。");
    if (!/^[ACGTRYSWKMBDHVN]+$/.test(sequence)) fail("参照配列にはDNAのIUPAC塩基（ACGT、Nなど）のみ使用できます。");
    return sequence;
  }

  function parseSequence(text) {
    if (typeof text !== "string") fail("参照配列はテキストまたはFASTAで入力してください。");
    var lines = text.replace(/^\uFEFF/, "").split(/\r?\n/);
    var headers = 0;
    var sequence = [];
    for (var i = 0; i < lines.length; i++) {
      var line = lines[i].trim();
      if (!line) continue;
      if (line[0] === ">") {
        headers++;
        if (headers > 1 || sequence.length) fail("参照配列には1レコードのFASTAを指定してください。染色体と開始座標を別途指定できます。");
      } else {
        sequence.push(line.replace(/\s/g, ""));
      }
    }
    return validateSequence(sequence.join("").toUpperCase());
  }

  // Returns the normalized IUPAC motif; only canonical bases in the reference match it.
  function motifPattern(text) {
    if (typeof text !== "string") fail("モチーフ配列を入力してください。");
    var motif = text.replace(/\s/g, "").toUpperCase();
    if (!motif.length || motif.length > LIMITS.motif) fail("モチーフは1～100 bpで入力してください（CpGの場合はCG）。");
    if (!/^[ACGTRYSWKMBDHVN]+$/.test(motif)) fail("モチーフにはDNAのIUPAC塩基のみ使用できます（CpGの場合はCG）。");
    return motif;
  }

  function reverseComplement(text) {
    if (typeof text !== "string") fail("逆相補鎖に変換する配列を指定してください。");
    var sequence = text.replace(/\s/g, "").toUpperCase();
    if (!/^[ACGTRYSWKMBDHVN]*$/.test(sequence)) fail("逆相補鎖にはDNAのIUPAC塩基のみ使用できます。");
    var result = new Array(sequence.length);
    for (var i = 0; i < sequence.length; i++) result[sequence.length - 1 - i] = COMPLEMENT[sequence[i]];
    return result.join("");
  }

  function markPattern(sequence, pattern, flags, bit, bases) {
    if (/^[ACGT]+$/.test(pattern)) {
      var position = sequence.indexOf(pattern);
      while (position !== -1) {
        flags[position] |= bit;
        position = sequence.indexOf(pattern, position + 1);
      }
      return;
    }
    var masks = new Uint8Array(pattern.length);
    var anchor = 0;
    var smallest = 5;
    for (var i = 0; i < pattern.length; i++) {
      var mask = MASK[pattern[i]];
      masks[i] = mask;
      var members = (mask & 1) + ((mask >> 1) & 1) + ((mask >> 2) & 1) + ((mask >> 3) & 1);
      if (members < smallest) { smallest = members; anchor = i; }
    }
    for (var start = 0; start <= sequence.length - pattern.length; start++) {
      if (!(bases[start + anchor] & masks[anchor])) continue;
      var matches = true;
      for (var j = 0; j < masks.length; j++) {
        if (!(bases[start + j] & masks[j])) { matches = false; break; }
      }
      if (matches) flags[start] |= bit;
    }
  }

  function scanMotif(sequence, motif, both) {
    sequence = parseSequence(sequence);
    motif = motifPattern(motif);
    var flags = new Uint8Array(sequence.length);
    var prefix = new Uint32Array(sequence.length + 1);
    var unknown = new Uint32Array(sequence.length + 1);
    var bases = new Uint8Array(sequence.length);
    for (var i = 0; i < sequence.length; i++) {
      var mask = MASK[sequence[i]];
      bases[i] = mask === 1 || mask === 2 || mask === 4 || mask === 8 ? mask : 0;
      unknown[i + 1] = unknown[i] + (bases[i] ? 0 : 1);
    }
    markPattern(sequence, motif, flags, 1, bases);
    if (both) markPattern(sequence, reverseComplement(motif), flags, 2, bases);
    for (var p = 0; p < flags.length; p++) prefix[p + 1] = prefix[p] + (flags[p] ? 1 : 0);
    return { flags: flags, prefix: prefix, unknown: unknown };
  }

  function coordinate(value, label, minimum) {
    if (value === "" || value === null || typeof value === "boolean") fail(label + "は" + minimum + "以上の整数にしてください。");
    var number = Number(value);
    if (!Number.isSafeInteger(number) || number < minimum) fail(label + "は" + minimum + "以上の整数にしてください。");
    return number;
  }

  function parseBed(text) {
    if (typeof text !== "string") fail("BEDファイルを読み込んでください。");
    var lines = text.replace(/^\uFEFF/, "").split(/\r?\n/);
    var rows = [];
    for (var i = 0; i < lines.length; i++) {
      var line = lines[i].trim();
      if (!line || line[0] === "#" || /^(track|browser)(\s|$)/.test(line)) continue;
      var fields = line.split(/\s+/);
      if (fields.length < 3 || !/^\d+$/.test(fields[1]) || !/^\d+$/.test(fields[2])) fail("BEDの" + (i + 1) + "行目はchrom・start・endが必要です。座標には0-basedの整数を指定してください。");
      var start = coordinate(fields[1], "BED start（" + (i + 1) + "行目）", 0);
      var end = coordinate(fields[2], "BED end（" + (i + 1) + "行目）", 0);
      if (end <= start) fail("BEDの" + (i + 1) + "行目はendをstartより大きくしてください。");
      rows.push({ chrom: fields[0], start: start, end: end, name: fields[3] || "region_" + (rows.length + 1) });
      if (rows.length > LIMITS.rows) fail("BED領域は100,000件以内にしてください。");
    }
    if (!rows.length) fail("BEDファイルに集計対象の領域がありません。");
    return rows;
  }

  function analyzeSequence(sequence, options) {
    options = options || {};
    sequence = parseSequence(sequence);
    var chrom = String(options.chrom || "reference").trim();
    if (!chrom || /\s/.test(chrom)) fail("染色体名に空白は使用できません。");
    var referenceStart = coordinate(options.start === undefined ? 0 : options.start, "参照配列の開始座標", 0);
    var referenceEnd = referenceStart + sequence.length;
    if (!Number.isSafeInteger(referenceEnd)) fail("参照配列の終端座標が大きすぎます。");
    var scan = scanMotif(sequence, options.motif, options.both);
    var rows = [];
    var mode = options.mode || "sliding";

    function measure(row) {
      row.count = null;
      row.density = null;
      if (row.chrom !== chrom) row.status = "chromosome_not_found";
      else if (row.start < referenceStart || row.end > referenceEnd) row.status = "outside_reference";
      else {
        var localStart = row.start - referenceStart;
        var localEnd = row.end - referenceStart;
        if (scan.unknown[localEnd] !== scan.unknown[localStart]) row.status = "ambiguous_bases";
        else {
          row.count = scan.prefix[localEnd] - scan.prefix[localStart];
          row.density = row.count / (row.end - row.start);
          row.status = "ok";
        }
      }
      rows.push(row);
    }

    if (mode === "sliding") {
      var flank = coordinate(options.flank === undefined ? 100 : options.flank, "Windowの片側幅", 0);
      var step = coordinate(options.step === undefined ? 1 : options.step, "Step", 1);
      var width = flank * 2 + 1;
      if (!Number.isSafeInteger(width)) fail("Window幅が大きすぎます。");
      if (width > sequence.length) fail("参照配列がWindow幅（" + width + " bp）より短いため、解析できません。");
      var total = Math.floor((sequence.length - width) / step) + 1;
      if (total > LIMITS.rows) fail("解析結果が100,000行を超えます。Stepを大きくするか、参照配列を短くしてください。");
      for (var local = 0; local + width <= sequence.length; local += step) {
        var center = referenceStart + local + flank;
        measure({ chrom: chrom, start: referenceStart + local, end: referenceStart + local + width, name: "window_" + (rows.length + 1), center: center });
      }
    } else if (mode === "bin") {
      var binSize = coordinate(options.binSize, "Bin幅", 1);
      if (Math.ceil(sequence.length / binSize) > LIMITS.rows) fail("Binが100,000件を超えます。Bin幅を大きくしてください。");
      for (var offset = 0; offset < sequence.length; offset += binSize) {
        var binStart = referenceStart + offset;
        var binEnd = referenceStart + Math.min(offset + binSize, sequence.length);
        measure({ chrom: chrom, start: binStart, end: binEnd, name: "bin_" + (rows.length + 1), center: Math.floor((binStart + binEnd - 1) / 2) });
      }
    } else if (mode === "bed") {
      var bedRows = parseBed(options.bed);
      for (var r = 0; r < bedRows.length; r++) measure(bedRows[r]);
    } else fail("解析方法はsliding、bin、bedのいずれかを指定してください。");
    return { rows: rows, scan: scan };
  }

  function parseGuides(text) {
    if (typeof text !== "string") fail("gRNA配列を入力してください。");
    var lines = text.replace(/^\uFEFF/, "").split(/\r?\n/).map(function (line) { return line.trim(); }).filter(function (line) { return line && line[0] !== "#"; });
    if (!lines.length) fail("gRNA配列を入力してください。");
    var guides = [];
    function add(name, sequence) {
      sequence = sequence.replace(/\s/g, "").toUpperCase();
      if (!/^[ACGT]+$/.test(sequence)) fail("gRNA「" + name + "」はACGTのみで入力してください。");
      if (sequence.length > LIMITS.guideLength) fail("gRNA「" + name + "」は200 bp以内にしてください。");
      guides.push({ name: name || "gRNA_" + (guides.length + 1), sequence: sequence, matches: 0 });
      if (guides.length > LIMITS.guides) fail("gRNAは100本以内にしてください。");
    }
    if (lines[0][0] === ">") {
      var name = null;
      var fragments = [];
      for (var i = 0; i < lines.length; i++) {
        if (lines[i][0] === ">") {
          if (name !== null) add(name, fragments.join(""));
          name = lines[i].slice(1).trim() || "gRNA_" + (guides.length + 1);
          fragments = [];
        } else fragments.push(lines[i]);
      }
      if (name !== null) add(name, fragments.join(""));
    } else {
      for (var j = 0; j < lines.length; j++) {
        var parts = lines[j].split(/\s+/);
        if (parts.length === 1) add("gRNA_" + (guides.length + 1), parts[0]);
        else if (parts.length === 2) add(parts[0], parts[1]);
        else fail("gRNAの各行は「配列」または「名前 配列」の形式にしてください。FASTAにも対応しています。");
      }
    }
    return guides;
  }

  function findGuides(sequence, text, options) {
    options = options || {};
    sequence = parseSequence(sequence);
    var guides = parseGuides(text);
    var referenceStart = coordinate(options.start === undefined ? 0 : options.start, "参照配列の開始座標", 0);
    if (!Number.isSafeInteger(referenceStart + sequence.length)) fail("参照配列の終端座標が大きすぎます。");
    var hits = [];
    var both = options.both !== false;
    function find(guide, target, strand) {
      var position = sequence.indexOf(target);
      while (position !== -1) {
        if (hits.length >= LIMITS.guideHits) fail("gRNAの一致が20,000件を超えます。短い配列や繰り返し配列を避け、gRNAの指定数を減らしてください。");
        hits.push({ name: guide.name, start: referenceStart + position, end: referenceStart + position + target.length, strand: strand, sequence: guide.sequence });
        guide.matches++;
        position = sequence.indexOf(target, position + 1);
      }
    }
    for (var i = 0; i < guides.length; i++) {
      find(guides[i], guides[i].sequence, "+");
      if (both) find(guides[i], reverseComplement(guides[i].sequence), "-");
    }
    hits.sort(function (a, b) { return a.start - b.start || a.end - b.end || a.name.localeCompare(b.name) || a.strand.localeCompare(b.strand); });
    return { guides: guides, hits: hits, unmatched: guides.filter(function (guide) { return !guide.matches; }).map(function (guide) { return guide.name; }) };
  }

  return Object.freeze({ LIMITS: LIMITS, parseSequence: parseSequence, motifPattern: motifPattern, reverseComplement: reverseComplement, scanMotif: scanMotif, parseBed: parseBed, analyzeSequence: analyzeSequence, findGuides: findGuides });
});
