#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scaling benchmark: load 10^2, 10^3 and 10^4 synthetic reports with the real
loaders and time the three use-case queries.

    python3 tools/scale_benchmark.py                       # sizes 100 1000 10000
    python3 tools/scale_benchmark.py --sizes 100 1000 --json results.json

For every size N the script writes N vendor files in the native formats
(Foundation Medicine and GenMineTOP XML, Guardant XLSX) into a temporary
directory, loads them with load_foundation.py, load_genminetop.py and
load_guardant.py plus the official gene lists (panels/), and reports the
wall-clock load time, database size, row counts and the median latency of the
recruitment, longitudinal and institutional queries, in SQL and end to end
through the CGI scripts (tools/qc_report.py --benchmark).

The reports are synthetic.  Their mix and density follow an institutional
registry of 345 reports (73% FoundationOne CDx, 7% FoundationOne Liquid CDx,
7% GenMineTOP, 13% Guardant360 CDx; about 15 variant rows per report and about
280 unreported candidate rows per Guardant workbook).  Variants are drawn from
the public hotspot variants of tests/make_fixtures.py; one report in five
belongs to a patient who already has a report, so that the longitudinal and
per-patient queries have work to do.  The random seed is fixed.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import random
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(ROOT / "tools"))
import make_fixtures as mf  # noqa: E402
import qc_report  # noqa: E402

MIX = (("F1", 253), ("F1L", 25), ("GM", 23), ("G360", 44))       # institutional registry
DISEASES = [  # FMI disease, disease label, OncoTree code, tissue
    ("COLON", "Colon adenocarcinoma (CRC)", "COAD", "Colon"),
    ("LUNG", "Lung adenocarcinoma", "LUAD", "Lung"),
    ("PROSTATE", "Prostate acinar adenocarcinoma", "PRAD", "Prostate"),
    ("BREAST", "Breast invasive ductal carcinoma (IDC)", "IDC", "Breast"),
    ("PANCREAS", "Pancreas ductal adenocarcinoma", "PAAD", "Pancreas"),
    ("STOMACH", "Stomach adenocarcinoma", "STAD", "Stomach"),
]

# Foundation Medicine pool (GRCh37), as in tests/make_fixtures.py
F1_SV = [
    mf.sv("KRAS", "chr12:25398284", "35G>A", "G12D", "missense", 0.3120, 812, strand="-", transcript="NM_004985"),
    mf.sv("APC", "chr5:112175639", "4348C>T", "R1450*", "nonsense", 0.4410, 690, transcript="NM_000038"),
    mf.sv("TP53", "chr17:7578406", "524G>A", "R175H", "missense", 0.5230, 745, strand="-", transcript="NM_000546"),
    mf.sv("SMAD4", "chr18:48586290", "904+1G>A", "splice site 904+1G>A", "splice", 0.2870, 603, status="likely",
          transcript="NM_005359"),
    mf.sv("TERT", "chr5:1295228", "-124C>T", "promoter -124C>T", "promoter", 0.2210, 402, status="likely",
          strand="-", transcript="NM_198253"),
    mf.sv("FBXW7", "chr4:153249384", "1394G>A", "R465H", "missense", 0.1850, 588, strand="-", transcript="NM_033632"),
    mf.sv("EGFR", "chr7:55242465", "2235_2249delGGAATTAAGAGAAGC", "E746_A750del", "nonframeshift", 0.3380, 921,
          transcript="NM_005228"),
    mf.sv("TP53", "chr17:7577120", "818G>A", "R273H", "missense", 0.4050, 733, strand="-", transcript="NM_000546"),
    mf.sv("STK11", "chr19:1220676", "580G>T", "D194Y", "missense", 0.1120, 510, status="unknown",
          transcript="NM_000455"),
    mf.sv("WHSC1L1", "chr8:38174537", "1780C>T", "R594*", "nonsense", 0.2010, 604, status="likely", strand="-",
          transcript="NM_023034"),
    mf.sv("BRCA2", "chr13:32914438", "5946delT", "S1982fs*22", "frameshift", 0.0480, 4120, status="likely",
          transcript="NM_000059"),
    mf.sv("ARID1A", "chr1:27100181", "4636C>T", "R1546*", "nonsense", 0.0310, 3980, status="likely",
          transcript="NM_006015"),
    mf.sv("TP53", "chr17:7578176", "672+1G>A", "splice site 672+1G>A", "splice", 0.0270, 4400, status="likely",
          strand="-", transcript="NM_000546"),
    mf.sv("EGFR", "chr7:55259515", "2573T>G", "L858R", "missense", 0.2510, 880, transcript="NM_005228"),
    mf.sv("PIK3CA", "chr3:178952085", "3140A>G", "H1047R", "missense", 0.2230, 760, transcript="NM_006218"),
    mf.sv("BRAF", "chr7:140453136", "1799T>A", "V600E", "missense", 0.2810, 690, strand="-", transcript="NM_004333"),
]
F1_CNV = [
    {"copy-number": "11", "equivocal": "false", "gene": "MYC", "number-of-exons": "3 of 3",
     "position": "chr8:128748315-128753680", "ratio": "2.41", "status": "known", "type": "amplification"},
    {"copy-number": "0", "equivocal": "false", "gene": "CDKN2A", "number-of-exons": "3 of 3",
     "position": "chr9:21967751-21995300", "ratio": "0.12", "status": "known", "type": "loss"},
    {"copy-number": "14", "equivocal": "false", "gene": "AR", "number-of-exons": "8 of 8",
     "position": "chrX:66764465-66950461", "ratio": "3.02", "status": "known", "type": "amplification"},
    {"copy-number": "18", "equivocal": "false", "gene": "ERBB2", "number-of-exons": "27 of 27",
     "position": "chr17:37844167-37884915", "ratio": "3.40", "status": "known", "type": "amplification"},
    {"copy-number": "9", "equivocal": "false", "gene": "MET", "number-of-exons": "21 of 21",
     "position": "chr7:116312459-116436178", "ratio": "2.20", "status": "known", "type": "amplification"},
]
F1_REARR = {"allele-fraction": "0.2140", "description": "EML4(NM_019063)-ALK(NM_004304) fusion (E13;A20)",
            "equivocal": "false", "in-frame": "Yes", "other-gene": "ALK", "percent-reads": "21.40",
            "pos1": "chr2:42522656", "pos2": "chr2:29446394", "status": "known", "supporting-read-pairs": "118",
            "targeted-gene": "EML4", "type": "fusion"}

# GenMineTOP pool (GRCh38)
GM_SHORT = [
    ("KRAS", "NM_004985", "chr12:25245350", "C", "T", "12p12.1", "snv", "c.35G>A", "p.G12D"),
    ("PIK3CA", "NM_006218", "chr3:179234297", "A", "G", "3q26.32", "snv", "c.3140A>G", "p.H1047R"),
    ("ARID1A", "NM_006015", "chr1:26772513", "G", "GA", "1p36.11", "insertion", "c.4001_4002insA", "p.G1335Rfs*20"),
    ("SMAD4", "NM_005359", "chr18:51065540", "TCTTTAAAAG", "T", "18q21.2", "deletion", "c.905-9_905del", ""),
    ("TP53", "NM_000546", "chr17:7675088", "C", "T", "17p13.1", "snv", "c.524G>A", "p.Arg175His"),
    ("EGFR", "NM_005228", "chr7:55191822", "T", "G", "7p11.2", "snv", "c.2573T>G", "p.L858R"),
    ("TP53", "NM_000546", "chr17:7673802", "C", "T", "17p13.1", "snv", "c.818G>A", "p.R273H"),
    ("BRAF", "NM_004333", "chr7:140753336", "A", "T", "7q34", "snv", "c.1799T>A", "p.V600E"),
    ("APC", "NM_000038", "chr5:112839942", "C", "T", "5q22.2", "snv", "c.4348C>T", "p.R1450*"),
    ("FBXW7", "NM_033632", "chr4:152328233", "C", "T", "4q31.3", "snv", "c.1394G>A", "p.R465H"),
]

# Guardant pool (GRCh37): reported rows
G_SNV = [
    dict(gene="KRAS", chrom="12", position=25398284, mut_nt="C>T", mut_aa="G12D", cdna="c.35G>A",
         transcript_id="NM_004985.5", exon="2.0"),
    dict(gene="APC", chrom="5", position=112175639, mut_nt="C>T", mut_aa="R1450*", cdna="c.4348C>T",
         transcript_id="NM_000038.6", exon="16.0"),
    dict(gene="TP53", chrom="17", position=7578406, mut_nt="C>T", mut_aa="R175H", cdna="c.524G>A",
         transcript_id="NM_000546.6", exon="5.0"),
    dict(gene="EGFR", chrom="7", position=55259515, mut_nt="T>G", mut_aa="L858R", cdna="c.2573T>G",
         transcript_id="NM_005228.5", exon="21.0"),
    dict(gene="EGFR", chrom="7", position=55249071, mut_nt="C>T", mut_aa="T790M", cdna="c.2369C>T",
         transcript_id="NM_005228.5", exon="20.0"),
    dict(gene="PIK3CA", chrom="3", position=178952085, mut_nt="A>G", mut_aa="H1047R", cdna="c.3140A>G",
         transcript_id="NM_006218.4", exon="21.0"),
]
G_INDEL = [
    dict(gene="APC", chrom="5", position=112174757, mut_nt="CA>C", mut_aa=None, cdna="c.3466del", length=1,
         exon="16.0", type="Deletion", transcript_id="NM_000038.6", mut_aa_short=None),
    dict(gene="EGFR", chrom="7", position=55242464, mut_nt="AGGAATTAAGAGAAGC>A", mut_aa="E746_A750del",
         cdna="c.2235_2249del", length=15, exon="19.0", type="Deletion", transcript_id="NM_005228.5",
         mut_aa_short="E746_A750del"),
]
G_UNCALLED_GENES = [("TP53", "17", 7579472), ("EGFR", "7", 55249063), ("ATM", "11", 108138003),
                    ("BRCA2", "13", 32929232), ("APC", "5", 112176325), ("KRAS", "12", 25362845)]


def expression_genes():
    f = ROOT / "panels" / "GenMineTOP_expression_genes.txt"
    return [g for g in f.read_text(encoding="utf-8").split() if g and not g.startswith("#")][:27]


def generate(n: int, out: Path, rng: random.Random) -> dict:
    """Write n synthetic reports under out/ and return counts per vendor."""
    for sub in ("foundation", "genminetop", "guardant"):
        (out / sub).mkdir(parents=True)
    total = sum(k for _, k in MIX)
    kinds = [kind for kind, k in MIX for _ in range(round(n * k / total))]
    while len(kinds) < n:
        kinds.append("F1")
    kinds = kinds[:n]
    rng.shuffle(kinds)
    patients, sidecar, counts = [], ["panel_name,report_id,patient_id,disease,oncotree_code,tissue_of_origin"], {}
    expr = expression_genes()
    for i, kind in enumerate(kinds, 1):
        counts[kind] = counts.get(kind, 0) + 1
        if patients and rng.random() < 0.2:
            pid, dis = rng.choice(patients)
        else:
            pid, dis = f"SCALEPT-{len(patients) + 1:06d}", rng.choice(DISEASES)
            patients.append((pid, dis))
        fmi_dis, label, oncotree, tissue = dis
        rid = f"SCALE-{kind}-{i:06d}"
        if kind in ("F1", "F1L"):
            k = rng.randint(7, 14) if kind == "F1" else rng.randint(9, 15)
            svs = [dict(v) for v in rng.sample(F1_SV, min(k, len(F1_SV)))]
            cnvs = rng.sample(F1_CNV, rng.randint(2, 5)) if kind == "F1" else rng.sample(F1_CNV, rng.randint(0, 2))
            rearr = [F1_REARR] if rng.random() < 0.4 else []
            attrs = {"disease": fmi_dis, "disease-ontology": label, "gender": rng.choice(["female", "male"]),
                     "pathology-diagnosis": label, "pipeline-version": "v3.29.0", "study": "Synthetic scaling",
                     "test-type": "FoundationOneDx" if kind == "F1" else "FoundationOneLiquidDx",
                     "tissue-of-origin": tissue}
            (out / "foundation" / f"{rid}.xml").write_text(mf.foundation_xml(
                rid, pid, attrs, short_variants=svs, cnvs=cnvs, rearrangements=rearr, msi="MSS",
                tmb={"score": f"{rng.uniform(0.5, 25):.2f}", "status": "low", "unit": "mutations-per-megabase"}),
                encoding="utf-8")
        elif kind == "GM":
            items = []
            for g in rng.sample(GM_SHORT, rng.randint(7, 10)):
                af = f"{rng.randint(20, 400)}/{rng.randint(401, 900)}"
                items.append(mf.gm_short(*g, af))
            if rng.random() < 0.5:
                items.append(mf.gm_item([
                    "<gene>ERBB2</gene>", "<transcript>NM_004448</transcript>",
                    "<locus>chr17:39687914-39730426</locus>", "<cytoband>17q12</cytoband>", "<origin>somatic</origin>",
                    "<type>cnv-amplification</type>", "<num-copy>9.412300</num-copy>", "<ratio>3.918200</ratio>",
                    "<status>finding</status>", "<clinical-relevance/>"]))
            for g in rng.sample(expr, 25):
                items.append(mf.gm_expression(g, "NM_000000", "1p1", f"{rng.uniform(1, 2000):.2f}",
                                              str(rng.randint(100, 50000)), "112", "88.10", "41.20"))
            (out / "genminetop" / f"{rid}.xml").write_text(mf.genminetop_xml(
                rid, "2025-06-18", pid, rng.choice(["female", "male"]), str(rng.randint(30, 85)), f"SP-{i:06d}",
                label, rng.randint(1, 30), f"{rng.uniform(0.2, 10):.6f}", "45", "50", items), encoding="utf-8")
        else:
            snv = [dict(r, percentage=round(rng.uniform(0.2, 15), 2), call=1, reporting_category="protein_coding",
                        rm_reportable=1) for r in rng.sample(G_SNV, rng.randint(1, 4))]
            indels = [dict(r, percentage=round(rng.uniform(0.2, 10), 2), call=1, reporting_category="protein_coding",
                           rm_reportable=1) for r in rng.sample(G_INDEL, rng.randint(0, 1))]
            for j in range(280):                       # unreported candidates (germline polymorphisms etc.)
                gene, chrom, pos = rng.choice(G_UNCALLED_GENES)
                snv.append(dict(gene=gene, chrom=chrom, position=pos + j, mut_nt="G>A", mut_aa=None,
                                cdna=f"c.{100 + j}+14G>A", percentage=rng.choice([49.8, 50.3, 99.7]), call=0,
                                transcript_id="NM_000000.0", exon=None, reporting_category="non_coding",
                                rm_reportable=0))
            cnas = [dict(chrom="17", gene="ERBB2", copy_number=3.46, call=1)] if rng.random() < 0.2 else []
            fus = ([dict(gene_a="EML4", gene_b="ALK", downstream_gene="ALK", chrom_a="2", chrom_b="2", pos_a=42522656,
                         pos_b=29446394, percentage=0.61, call=1)] if rng.random() < 0.05 else [])
            msi = [dict(runid=f"RUN{i}", run_sample_id=f"S{i}", msi_score=rng.randint(0, 30), msi_status="MSS/MSI-L",
                        max_maf=0.05, mean_maf=None)]
            mf.write_book(out / "guardant" / f"Interim_{pid}_{rid}.xlsx", {
                "SNV": (mf.SNV_COLS, snv), "Indels": (mf.INDEL_COLS, indels), "CNAs": (mf.CNA_COLS, cnas),
                "Fusions": (mf.FUSION_COLS, fus), "MSI": (mf.MSI_COLS, msi)})
            sidecar.append(f"Guardant,{rid},,{label},{oncotree},{tissue}")
    (out / "case_metadata.csv").write_text("\n".join(sidecar) + "\n", encoding="utf-8")
    return {"reports_by_kind": counts, "patients": len(patients)}


def load(out: Path, db: Path) -> float:
    t0 = time.perf_counter()
    for loader, sub in (("load_foundation.py", "foundation"), ("load_genminetop.py", "genminetop"),
                        ("load_guardant.py", "guardant")):
        if any((out / sub).iterdir()):
            subprocess.run([sys.executable, str(ROOT / loader), str(db), str(out / sub), "--case-metadata",
                            str(out / "case_metadata.csv"), "--quiet"], check=True, capture_output=True)
    subprocess.run([sys.executable, str(ROOT / "load_panel_genes.py"), str(db), str(ROOT / "panels")],
                   check=True, capture_output=True)
    return time.perf_counter() - t0


def run_size(n: int, workdir: Path, seed: int) -> dict:
    out = workdir / f"n{n}"
    if out.exists():
        shutil.rmtree(out)
    rng = random.Random(seed + n)
    info = generate(n, out, rng)
    db = out / "panels.db"
    load_s = load(out, db)
    conn = sqlite3.connect(db)
    rows = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
            for t in ("cases", "variants", "biomarkers", "panel_genes")}
    multi = conn.execute("SELECT COUNT(*) FROM (SELECT patient_id FROM cases WHERE patient_id IS NOT NULL "
                         "GROUP BY patient_id HAVING COUNT(*) >= 2)").fetchone()[0]
    conn.close()
    b = qc_report.benchmark(str(db))
    return {"reports": n, **info, "patients_with_2_or_more_reports": multi, "load_seconds": round(load_s, 1),
            "database_mb": round(os.path.getsize(db) / 1e6, 1), "rows": rows,
            "sql_median_ms": b["sql_median_ms"], "cgi_median_ms": b.get("cgi_end_to_end_median_ms", {})}


def machine() -> str:
    if sys.platform == "darwin":
        cpu = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True).stdout.strip()
        mem = int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout) / 2**30
        return f"{cpu}, {mem:.0f} GB, macOS {platform.mac_ver()[0]}"
    cpu = next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")), "?")
    return f"{cpu}, {platform.system()} {platform.release()}"


def markdown(results: dict) -> str:
    sizes = results["sizes"]
    head = "| | " + " | ".join(f"{r['reports']:,} reports" for r in sizes) + " |"
    rows = [head, "|---" * (len(sizes) + 1) + "|"]

    def line(label, fn):
        rows.append(f"| {label} | " + " | ".join(fn(r) for r in sizes) + " |")

    line("variant rows", lambda r: f"{r['rows']['variants']:,}")
    line("biomarker rows", lambda r: f"{r['rows']['biomarkers']:,}")
    line("patients with ≥ 2 reports", lambda r: f"{r['patients_with_2_or_more_reports']:,}")
    line("database size (MB)", lambda r: f"{r['database_mb']}")
    line("loading, all loaders and gene lists (s)", lambda r: f"{r['load_seconds']}")
    for key in sizes[0]["sql_median_ms"]:
        name = key.split(" ")[0]
        line(f"{name} query, SQL (ms)", lambda r, k=key: f"{r['sql_median_ms'][k]:.2f}")
        line(f"{name} query, web interface (ms)", lambda r, n=name: f"{r['cgi_median_ms'].get(n, float('nan')):.1f}")
    return "\n".join([
        "# Scaling benchmark",
        "",
        "Generated by `python3 tools/scale_benchmark.py --json docs/scale_benchmark.json --markdown docs/BENCHMARK.md`.",
        "Synthetic reports in the native vendor formats (vendor mix and variant density of a 345-report "
        "institutional registry) were loaded with the loaders and the official gene lists; query times are "
        "medians (SQL: 21 runs; web interface: 3 runs of the CGI script, including Perl start-up).  Queries: "
        "recruitment = KRAS p.G12D with a diagnosis containing 'colon', per patient; longitudinal = all variants "
        "of the patient with most reports; institutional = panel-aware gene frequencies (n_mutated/n_tested) over "
        "all reports.  The SQL figures are those of the queries in tools/qc_report.py; the web interface runs its "
        "own SQL (panel_stats.cgi aggregates the tested reports per gene before joining the variants), so the two "
        "institutional figures are not directly comparable.",
        "",
        f"Machine: {results['machine']}; Python {results['python']}; SQLite {results['sqlite']}; "
        f"seed {results['seed']}.",
        "",
        *rows,
        ""])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sizes", type=int, nargs="+", default=[100, 1000, 10000])
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--workdir", help="where to write the synthetic reports (default: a temporary directory)")
    ap.add_argument("--keep", action="store_true", help="keep the generated reports and databases")
    ap.add_argument("--json", metavar="FILE", help="also write the results as JSON")
    ap.add_argument("--markdown", metavar="FILE", help="also write the results as a Markdown table")
    ap.add_argument("--from-json", metavar="FILE", help="only render --markdown from an earlier --json file")
    args = ap.parse_args(argv)
    if args.from_json:
        Path(args.markdown).write_text(markdown(json.loads(Path(args.from_json).read_text())), encoding="utf-8")
        return 0

    workdir = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="oncounify-scale-"))
    workdir.mkdir(parents=True, exist_ok=True)
    results = {"machine": machine(), "python": platform.python_version(), "sqlite": sqlite3.sqlite_version,
               "seed": args.seed, "sizes": []}
    try:
        for n in args.sizes:
            r = run_size(n, workdir, args.seed)
            results["sizes"].append(r)
            s, c = r["sql_median_ms"], r["cgi_median_ms"]
            print(f"n={n:>6}  load {r['load_seconds']:>7.1f} s  db {r['database_mb']:>6.1f} MB  "
                  f"variants {r['rows']['variants']:>7}  patients>=2 reports {r['patients_with_2_or_more_reports']:>5}",
                  flush=True)
            for k in s:
                name = k.split(" ")[0]
                cgi = "n/a" if c.get(name) is None else f"{c[name]:.1f}"
                print(f"          {name:<14} SQL {s[k]:>9.2f} ms   CGI {cgi:>8} ms", flush=True)
    finally:
        if not args.keep and not args.workdir:
            shutil.rmtree(workdir, ignore_errors=True)
    print(f"machine: {results['machine']}; Python {results['python']}; SQLite {results['sqlite']}")
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        Path(args.markdown).write_text(markdown(results), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
