#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate the synthetic vendor fixtures under tests/data/.

Every identifier, date, and patient attribute in these files is invented.
Variants are well-known, publicly documented cancer hotspots (e.g. KRAS
p.G12D, EGFR p.L858R) placed on the assembly each vendor reports on
(Foundation Medicine and Guardant: GRCh37; GenMineTOP: GRCh38).  Only the
file *structure* (element, attribute, sheet and column names) mirrors the
vendor deliverables; no value was copied from a real report.

The fixtures are designed to exercise every code path of the loaders:
each variant class, vendor-supplied vs. inferred consequences, splice-site
and splice-region coding changes, in-frame indels, promoter variants,
uncalled Guardant rows, both Guardant fusion-sheet spellings, germline
findings, fusions and RNA exon skipping, TMB/MSI in every vendor flavour,
and one patient tested by two vendors (SYNPT-0001).

Run from the repository root:
    python3 tests/make_fixtures.py
"""

from __future__ import annotations

from pathlib import Path
from xml.sax.saxutils import escape, quoteattr

from openpyxl import Workbook

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
INVALID = HERE / "data_invalid"


# ---------------------------------------------------------------------------
# Foundation Medicine XML
# ---------------------------------------------------------------------------

def _attrs(d):
    return " ".join(f"{k}={quoteattr(str(v))}" for k, v in d.items() if v is not None)


def foundation_xml(ref_id, mrn, report_attrs, short_variants=(), cnvs=(), rearrangements=(),
                   msi=None, tmb=None, non_human=(), sample="SQ-SYN0000001-1"):
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<rr:ResultsReport xmlns:rr="http://integration.foundationmedicine.com/reporting" '
           'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
           'xsi:schemaLocation="http://integration.foundationmedicine.com/reporting '
           'http://integration.foundationmedicine.com/reporting/ResultsReport.2.1.xsd">',
           '    <rr:CustomerInformation>',
           f'        <rr:ReferenceID>{escape(ref_id)}</rr:ReferenceID>',
           '        <rr:CSN/>',
           f'        <rr:TRF>ORD-SYN-{escape(ref_id[-4:])}-01</rr:TRF>',
           f'        <rr:MRN>{escape(mrn)}</rr:MRN>' if mrn else '        <rr:MRN/>',
           '        <rr:PhysicianId/>',
           '        <rr:NPI/>',
           '    </rr:CustomerInformation>',
           '    <rr:ResultsPayload>',
           '        <variant-report xmlns="http://foundationmedicine.com/compbio/variant-report-external" '
           'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
           'xsi:schemaLocation="http://foundationmedicine.com/compbio/variant-report-external '
           'http://integration.foundationmedicine.com/reporting/variant-report-external-2.2.xsd" '
           + _attrs(report_attrs) + '>',
           '            <samples>',
           f'                <sample bait-set="DX2" mean-exon-depth="812.40" name="{sample}" nucleic-acid-type="DNA"/>',
           '            </samples>',
           '            <quality-control status="Pass"/>',
           '            <short-variants>']
    for sv in short_variants:
        out.append(f'                <short-variant {_attrs(sv)}>')
        out.append(f'                    <dna-evidence sample="{sample}"/>')
        out.append('                </short-variant>')
    out.append('            </short-variants>')
    out.append('            <copy-number-alterations>')
    for c in cnvs:
        out.append(f'                <copy-number-alteration {_attrs(c)}>')
        out.append(f'                    <dna-evidence sample="{sample}"/>')
        out.append('                </copy-number-alteration>')
    out.append('            </copy-number-alterations>')
    out.append('            <rearrangements>')
    for r in rearrangements:
        out.append(f'                <rearrangement {_attrs(r)}>')
        out.append(f'                    <dna-evidence sample="{sample}"/>')
        out.append('                </rearrangement>')
    out.append('            </rearrangements>')
    out.append('            <biomarkers>')
    if msi:
        out.append(f'                <microsatellite-instability status="{msi}"/>')
    if tmb:
        out.append(f'                <tumor-mutation-burden {_attrs(tmb)}/>')
    out.append('            </biomarkers>')
    if non_human:
        out.append('            <non-human-content>')
        for nh in non_human:
            out.append(f'                <non-human {_attrs(nh)}>')
            out.append(f'                    <dna-evidence sample="{sample}"/>')
            out.append('                </non-human>')
        out.append('            </non-human-content>')
    else:
        out.append('            <non-human-content/>')
    out += ['        </variant-report>', '    </rr:ResultsPayload>', '</rr:ResultsReport>', '']
    return "\n".join(out)


def sv(gene, pos, cds, protein, fe, af, depth, status="known", strand="+", transcript="NM_000000", **kw):
    d = {"allele-fraction": f"{af:.4f}", "cds-effect": cds, "depth": str(depth), "equivocal": "false",
         "functional-effect": fe, "gene": gene, "percent-reads": f"{af * 100:.2f}", "position": pos,
         "protein-effect": protein, "status": status, "strand": strand, "transcript": transcript}
    d.update(kw)
    return d


def write_foundation():
    d = DATA / "foundation"
    d.mkdir(parents=True, exist_ok=True)

    # 1. Colon adenocarcinoma, FoundationOne CDx; MRN left empty (linked via sidecar)
    (d / "SYN-F1-0001.xml").write_text(foundation_xml(
        "SYN-F1-0001", None,
        {"disease": "COLON", "disease-ontology": "Colon adenocarcinoma (CRC)", "flowcell-analysis": "9000000001",
         "gender": "female", "pathology-diagnosis": "Colon adenocarcinoma (CRC)", "percent-tumor-nuclei": "40",
         "pipeline-version": "v3.29.0", "purity-assessment": "38.50", "specimen": "ORD-SYN-0001-01*SYN0000001.01",
         "study": "Synthetic F1CDx", "test-request": "ORD-SYN-0001-01", "test-type": "FoundationOneDx",
         "tissue-of-origin": "Colon"},
        short_variants=[
            sv("KRAS", "chr12:25398284", "35G>A", "G12D", "missense", 0.3120, 812, strand="-", transcript="NM_004985"),
            sv("APC", "chr5:112175639", "4348C>T", "R1450*", "nonsense", 0.4410, 690, transcript="NM_000038"),
            sv("TP53", "chr17:7578406", "524G>A", "R175H", "missense", 0.5230, 745, strand="-", transcript="NM_000546"),
            sv("SMAD4", "chr18:48586290", "904+1G>A", "splice site 904+1G>A", "splice", 0.2870, 603,
               status="likely", transcript="NM_005359"),
            sv("TERT", "chr5:1295228", "-124C>T", "promoter -124C>T", "promoter", 0.2210, 402,
               status="likely", strand="-", transcript="NM_198253"),
            sv("FBXW7", "chr4:153249384", "1394G>A", "R465H", "missense", 0.1850, 588, strand="-",
               transcript="NM_033632"),
        ],
        cnvs=[{"copy-number": "11", "equivocal": "false", "gene": "MYC", "number-of-exons": "3 of 3",
               "position": "chr8:128748315-128753680", "ratio": "2.41", "status": "known", "type": "amplification"}],
        msi="MSS", tmb={"score": "6.30", "status": "low", "unit": "mutations-per-megabase"},
        sample="SQ-SYN0000001-1"), encoding="utf-8")

    # 2. Lung adenocarcinoma, FoundationOne CDx: in-frame EGFR deletion, EML4-ALK
    (d / "SYN-F1-0002.xml").write_text(foundation_xml(
        "SYN-F1-0002", "SYNPT-0003",
        {"disease": "LUNG", "disease-ontology": "Lung adenocarcinoma", "flowcell-analysis": "9000000002",
         "gender": "male", "pathology-diagnosis": "Lung adenocarcinoma", "percent-tumor-nuclei": "60",
         "pipeline-version": "v3.29.0", "purity-assessment": "55.10", "specimen": "ORD-SYN-0002-01*SYN0000002.01",
         "study": "Synthetic F1CDx", "test-request": "ORD-SYN-0002-01", "test-type": "FoundationOneDx",
         "tissue-of-origin": "Lung"},
        short_variants=[
            sv("EGFR", "chr7:55242465", "2235_2249delGGAATTAAGAGAAGC", "E746_A750del", "nonframeshift",
               0.3380, 921, transcript="NM_005228"),
            sv("TP53", "chr17:7577120", "818G>A", "R273H", "missense", 0.4050, 733, strand="-",
               transcript="NM_000546"),
            sv("STK11", "chr19:1220676", "580G>T", "D194Y", "missense", 0.1120, 510, status="unknown",
               transcript="NM_000455"),
            # previous HGNC symbol as used in Foundation Medicine reports (current: NSD3)
            sv("WHSC1L1", "chr8:38174537", "1780C>T", "R594*", "nonsense", 0.2010, 604, status="likely",
               strand="-", transcript="NM_023034"),
        ],
        cnvs=[{"copy-number": "0", "equivocal": "false", "gene": "CDKN2A", "number-of-exons": "3 of 3",
               "position": "chr9:21967751-21995300", "ratio": "0.12", "status": "known", "type": "loss"}],
        rearrangements=[{"allele-fraction": "0.2140", "description": "EML4(NM_019063)-ALK(NM_004304) fusion (E13;A20)",
                         "equivocal": "false", "in-frame": "Yes", "other-gene": "ALK", "percent-reads": "21.40",
                         "pos1": "chr2:42522656", "pos2": "chr2:29446394", "status": "known",
                         "supporting-read-pairs": "118", "targeted-gene": "EML4", "type": "fusion"}],
        msi="MSS", tmb={"score": "4.10", "status": "low", "unit": "mutations-per-megabase"},
        sample="SQ-SYN0000002-1"), encoding="utf-8")

    # 3. Prostate, FoundationOne Liquid CDx: AR amplification (vs ARID1A), non-human content
    (d / "SYN-F1L-0003.xml").write_text(foundation_xml(
        "SYN-F1L-0003", "SYNPT-0005",
        {"disease": "PROSTATE", "disease-ontology": "Prostate acinar adenocarcinoma", "flowcell-analysis": "9000000003",
         "gender": "male", "non-human-content": "0.0035", "pathology-diagnosis": "Prostate acinar adenocarcinoma",
         "pipeline-version": "v3.8.2", "specimen": "ORD-SYN-0003-01*SYN0000003.01", "study": "Synthetic F1LCDx",
         "test-request": "ORD-SYN-0003-01", "test-type": "FoundationOneLiquidDx", "tissue-of-origin": "Prostate"},
        short_variants=[
            sv("BRCA2", "chr13:32914438", "5946delT", "S1982fs*22", "frameshift", 0.0480, 4120,
               status="likely", transcript="NM_000059"),
            sv("ARID1A", "chr1:27100181", "4636C>T", "R1546*", "nonsense", 0.0310, 3980, status="likely",
               transcript="NM_006015"),
            sv("TP53", "chr17:7578176", "672+1G>A", "splice site 672+1G>A", "splice", 0.0270, 4400,
               status="likely", strand="-", transcript="NM_000546"),
        ],
        cnvs=[{"copy-number": "14", "equivocal": "false", "gene": "AR", "number-of-exons": "8 of 8",
               "position": "chrX:66764465-66950461", "ratio": "3.02", "status": "known", "type": "amplification"}],
        msi="MSS", tmb={"score": "7.58", "status": "low", "unit": "mutations-per-megabase"},
        non_human=[{"organism": "HPV-16", "reads-per-million": "31", "status": "unknown"},
                   {"organism": "HHV-4", "reads-per-million": "12", "status": "unknown"}],
        sample="SQ-SYN0000003-1"), encoding="utf-8")


# ---------------------------------------------------------------------------
# GenMineTOP XML
# ---------------------------------------------------------------------------

def gm_item(inner):
    return "      <item>\n" + "\n".join("        " + line for line in inner) + "\n      </item>"


def gm_short(gene, transcript, locus, ref, alt, cytoband, typ, cds, protein, af, origin="somatic",
             status="finding", clinvar=None, maf=None, ag_class=None):
    x = [f"<gene>{gene}</gene>", f"<transcript>{transcript}</transcript>", f"<locus>{locus}</locus>",
         f"<ref>{ref}</ref>", f"<alt>{alt}</alt>", f"<cytoband>{cytoband}</cytoband>", f"<origin>{origin}</origin>"]
    if ag_class:
        x.append(f"<ag-class>{ag_class}</ag-class>")
    x += [f"<type>{typ}</type>", f"<coding-dna-alteration>{escape(cds)}</coding-dna-alteration>",
          f"<protein-alteration>{escape(protein)}</protein-alteration>", f"<allele-frequency>{af}</allele-frequency>",
          f"<status>{status}</status>"]
    if clinvar:
        x += ["<dbs>", f"  <clinvar><item>https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar[0]}/</item></clinvar>",
              "</dbs>"]
    x.append("<clinical-relevance>")
    if maf:
        x += ["  <minor-allele-frequency>", f"    <key name=\"1000genomes\">{maf[0]}</key>",
              f"    <hgvd>{maf[1]}</hgvd>", f"    <tommo-8p3kjpn>{maf[2]}</tommo-8p3kjpn>", "  </minor-allele-frequency>"]
    if clinvar:
        x += ["  <dbs>", "    <clinvar>", f"      <id><item>{clinvar[0]}</item></id>",
              f"      <clinical-significance><item>{clinvar[1]}</item></clinical-significance>",
              "      <match><item>exact</item></match>", "    </clinvar>", "  </dbs>",
              "  <content>", f"    <pathogenic>{clinvar[2]}</pathogenic>", "  </content>"]
    x.append("</clinical-relevance>")
    return gm_item(x)


def gm_expression(gene, transcript, cytoband, tpm, reads, n, mean, sd):
    return gm_item([
        f"<gene>{gene}</gene>", f"<transcript>{transcript}</transcript>", f"<cytoband>{cytoband}</cytoband>",
        "<origin>somatic</origin>", "<type>expression</type>", f"<num-reads>{reads}</num-reads>", f"<tpm>{tpm}</tpm>",
        "<reference>", "  <normal-expression>", "    <tpm>", f"      <n>{n}</n>", f"      <mean>{mean}</mean>",
        f"      <sd>{sd}</sd>", "    </tpm>", "  </normal-expression>", "</reference>",
        "<status>finding</status>", "<clinical-relevance/>"])


def gm_breakpoint_item(genes, transcripts, loci, cytobands, typ, reads, wt, frame, vid, bps, sv_type=None):
    x = []
    x.append("<gene>" + "".join(f"<item>{g}</item>" for g in genes) + "</gene>" if len(genes) > 1
             else f"<gene>{genes[0]}</gene>")
    x.append("<transcript>" + "".join(f"<item>{t}</item>" for t in transcripts) + "</transcript>")
    x.append("<locus>" + "".join(f"<item>{loc}</item>" for loc in loci) + "</locus>")
    x.append(f"<id><vendor>{vid}</vendor></id>")
    x.append("<cytoband>" + "".join(f"<item>{c}</item>" for c in cytobands) + "</cytoband>")
    x += ["<origin>somatic</origin>", f"<type>{typ}</type>"]
    if sv_type:
        x.append(f"<sv-type>{sv_type}</sv-type>")
    x += [f"<num-reads>{reads}</num-reads>", "<num-wt-reads>" + "".join(f"<item>{w}</item>" for w in wt) + "</num-wt-reads>",
          f"<frame>{frame}</frame>", "<breakpoint>"]
    for b in bps:
        x.append("  <item>" + "".join(f"<{k}>{v}</{k}>" for k, v in b.items()) + "</item>")
    x += ["</breakpoint>", "<status>finding</status>", "<clinical-relevance/>"]
    return gm_item(x)


def genminetop_xml(report_id, accepted, patient_id, sex, age, specimen_id, pathology, tmb_n, tmb_f,
                   est, nuclei, items):
    return "\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>',
        "<root>",
        "  <version>2</version>",
        "  <vendor-id>todai-oncopanel</vendor-id>",
        "  <target/>",
        "  <report>",
        f"    <id>{report_id}</id>",
        "    <owner><hospital>Synthetic Hospital</hospital><doctor>Synthetic Doctor</doctor></owner>",
        f"    <date><accepted>{accepted}</accepted><reported/><updated/></date>",
        f"    <patient><id>{patient_id}</id><sex>{sex}</sex><age>{age}</age><c-cat-id>SYN-CCAT-{report_id[-4:]}</c-cat-id></patient>",
        f"    <specimen><id>{specimen_id}</id><pathology>{escape(pathology)}</pathology></specimen>",
        "    <preference><germline-disclosure>true</germline-disclosure><delivery/></preference>",
        "    <test><type>todai-oncopanel</type><panel><item>TDv6</item><item>TRv6</item></panel></test>",
        "    <reference>",
        "      <db><item><name>ClinVar</name><version>202401</version></item>"
        "<item><name>COSMIC</name><version>95</version></item></db>",
        "      <snp-db><item><name>1000 Genomes</name><version>Phase_3(20130502)</version><genome>hg38</genome></item>"
        "<item><name>ToMMo</name><version>8p3kjpn-20200831</version><genome>hg38</genome></item>"
        "<item><name>HGVD</name><version>v2.30</version><genome>hg38</genome></item></snp-db>",
        "    </reference>",
        "    <result>",
        "      <reference-genome>hg38</reference-genome>",
        f"      <marker><tmb><exon><num-non-synonymous-alterations>{tmb_n}</num-non-synonymous-alterations>"
        f"<frequency-non-synonymous-alterations>{tmb_f}</frequency-non-synonymous-alterations></exon></tmb></marker>",
        "      <alterations>",
        *items,
        "      </alterations>",
        "    </result>",
        "    <qc><tumor-content>"
        f"<estimated><value>{est}</value><unit>percentage</unit></estimated>"
        f"<nuclei><value>{nuclei}</value><unit>percentage</unit></nuclei></tumor-content></qc>",
        "  </report>",
        "  <input><options><genome-version>hg38</genome-version></options></input>",
        "</root>",
        ""])


def write_genminetop():
    d = DATA / "genminetop"
    d.mkdir(parents=True, exist_ok=True)

    items1 = [
        gm_short("KRAS", "NM_004985", "chr12:25245350", "C", "T", "12p12.1", "snv", "c.35G>A", "p.G12D", "201/634",
                 clinvar=("12582", "Pathogenic", "12"), maf=("0", "0", "0"), ag_class="1"),
        gm_short("PIK3CA", "NM_006218", "chr3:179234297", "A", "G", "3q26.32", "snv", "c.3140A>G", "p.H1047R",
                 "143/520", clinvar=("13652", "Pathogenic", "10")),
        gm_short("ARID1A", "NM_006015", "chr1:26772513", "G", "GA", "1p36.11", "insertion", "c.4001_4002insA",
                 "p.G1335Rfs*20", "88/401"),
        gm_short("SMAD4", "NM_005359", "chr18:51065540", "TCTTTAAAAG", "T", "18q21.2", "deletion",
                 "c.905-9_905del", "", "97/455"),
        gm_short("BRCA1", "NM_007294", "chr17:43057063", "G", "GC", "17q21.31", "insertion", "c.5266dup",
                 "p.Gln1756ProfsTer74", "231/470", origin="germline", status="notice",
                 clinvar=("17677", "Pathogenic", "40"), maf=("0.0002", "0", "0.0001")),
        gm_short("TP53", "NM_000546", "chr17:7675088", "C", "T", "17p13.1", "snv", "c.524G>A", "p.Arg175His",
                 "312/590", clinvar=("12374", "Pathogenic", "25"), ag_class="1"),
        gm_item(["<gene>ERBB2</gene>", "<transcript>NM_004448</transcript>",
                 "<locus>chr17:39687914-39730426</locus>", "<cytoband>17q12</cytoband>", "<origin>somatic</origin>",
                 "<type>cnv-amplification</type>", "<num-copy>9.412300</num-copy>", "<ratio>3.918200</ratio>",
                 "<status>finding</status>", "<clinical-relevance/>"]),
        gm_expression("ERBB2", "NM_004448", "17q12", "1523.41", "40211", "112", "88.10", "41.20"),
        gm_expression("MYC", "NM_002467", "8q24.21", "412.77", "9834", "112", "150.30", "70.90"),
    ]
    (d / "SYN-GM-0001.xml").write_text(genminetop_xml(
        "SYN-GM-0001", "2025-06-18", "SYNPT-0002", "male", "63", "SYN-SP-0001",
        "Adenocarcinoma of the colon", 5, "3.412884", "45", "50", items1), encoding="utf-8")

    items2 = [
        gm_short("EGFR", "NM_005228", "chr7:55191822", "T", "G", "7p11.2", "snv", "c.2573T>G", "p.L858R",
                 "256/701", clinvar=("16609", "Pathogenic", "18"), ag_class="1"),
        gm_breakpoint_item(["MET"], ["NM_001127500"], ["chr7:116771976", "chr7:116774880"], ["7q31.2", "7q31.2"],
                           "splicing-variant", "1812", ["5120", "4987"], "in-frame", "S000001",
                           [{"region": "exon", "index": "13", "length": "36", "gene": "MET",
                             "transcript": "NM_001127500", "chr": "chr7", "pos": "116771976"},
                            {"region": "exon", "index": "15", "length": "74", "gene": "MET",
                             "transcript": "NM_001127500", "chr": "chr7", "pos": "116774880"}],
                           sv_type="skipping"),
        gm_breakpoint_item(["EML4", "ALK"], ["NM_019063", "NM_004304"], ["chr2:42295517", "chr2:29223528"],
                           ["2p21", "2p23.2"], "fusion", "642", ["10201", "85"], "in-frame", "F0001",
                           [{"region": "exon", "index": "13", "length": "112", "gene": "EML4",
                             "transcript": "NM_019063", "chr": "chr2", "pos": "42295517"},
                            {"region": "exon", "index": "20", "length": "86", "gene": "ALK",
                             "transcript": "NM_004304", "chr": "chr2", "pos": "29223528"}]),
        gm_expression("MET", "NM_001127500", "7q31.2", "845.20", "19002", "112", "60.40", "25.80"),
    ]
    (d / "SYN-GM-0002.xml").write_text(genminetop_xml(
        "SYN-GM-0002", "2025-08-02", "SYNPT-0004", "female", "71", "SYN-SP-0002",
        "Adenocarcinoma of the lung", 2, "1.365526", "60", "70", items2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Guardant XLSX
# ---------------------------------------------------------------------------

SNV_COLS = ["gene", "chrom", "position", "mut_nt", "mut_aa", "cdna", "percentage", "call", "transcript_id",
            "exon", "reporting_category", "rm_reportable"]
INDEL_COLS = ["gene", "chrom", "position", "mut_nt", "mut_aa", "cdna", "length", "exon", "type", "percentage",
              "call", "transcript_id", "reporting_category", "mut_aa_short", "rm_reportable"]
CNA_COLS = ["chrom", "gene", "copy_number", "call"]
FUSION_COLS = ["gene_a", "gene_b", "downstream_gene", "chrom_a", "chrom_b", "pos_a", "pos_b", "percentage", "call"]
MSI_COLS = ["runid", "run_sample_id", "msi_score", "msi_status", "max_maf", "mean_maf"]
QC_COLS = ["runid", "run_sample_id", "category", "metric", "verbose_name", "unit", "value", "status"]


def write_book(path, sheets):
    wb = Workbook()
    wb.remove(wb.active)
    for name, (cols, rows) in sheets.items():
        ws = wb.create_sheet(name)
        ws.append(cols)
        for r in rows:
            ws.append([r.get(c) for c in cols])
    wb.save(path)


def write_guardant():
    d = DATA / "guardant"
    d.mkdir(parents=True, exist_ok=True)

    # Colon, same patient as SYN-F1-0001 (longitudinal liquid biopsy)
    snv = [
        dict(gene="KRAS", chrom="12", position=25398284, mut_nt="C>T", mut_aa="G12D", cdna="c.35G>A",
             percentage=4.21, call=1, transcript_id="NM_004985.5", exon="2.0", reporting_category="protein_coding",
             rm_reportable=1),
        dict(gene="APC", chrom="5", position=112175639, mut_nt="C>T", mut_aa="R1450*", cdna="c.4348C>T",
             percentage=3.87, call=1, transcript_id="NM_000038.6", exon="16.0", reporting_category="protein_coding",
             rm_reportable=1),
        dict(gene="TP53", chrom="17", position=7578406, mut_nt="C>T", mut_aa="R175H", cdna="c.524G>A",
             percentage=5.02, call=1, transcript_id="NM_000546.6", exon="5.0", reporting_category="protein_coding",
             rm_reportable=1),
        # uncalled: heterozygous / homozygous germline polymorphisms and a synonymous change
        dict(gene="TP53", chrom="17", position=7579472, mut_nt="G>C", mut_aa="P72R", cdna="c.215C>G",
             percentage=49.8, call=0, transcript_id="NM_000546.6", exon="4.0", reporting_category="protein_coding",
             rm_reportable=0),
        dict(gene="EGFR", chrom="7", position=55249063, mut_nt="G>A", mut_aa="Q787Q", cdna="c.2361G>A",
             percentage=99.7, call=0, transcript_id="NM_005228.5", exon="20.0", reporting_category="protein_coding",
             rm_reportable=0),
        dict(gene="ATM", chrom="11", position=108138003, mut_nt="T>C", mut_aa=None, cdna="c.2921+14T>C",
             percentage=51.2, call=0, transcript_id="NM_000051.4", exon="19.0", reporting_category="non_coding",
             rm_reportable=0),
    ]
    indels = [
        dict(gene="APC", chrom="5", position=112174757, mut_nt="CA>C", mut_aa=None, cdna="c.3466del", length=1,
             exon="16.0", type="Deletion", percentage=2.12, call=1, transcript_id="NM_000038.6",
             reporting_category="protein_coding", mut_aa_short=None, rm_reportable=1),
        dict(gene="BRCA2", chrom="13", position=32929232, mut_nt="AAAG>A", mut_aa=None, cdna="c.7008-15_7008-13del",
             length=3, exon=None, type="Deletion", percentage=48.9, call=0, transcript_id="NM_000059.4",
             reporting_category="non_coding", mut_aa_short=None, rm_reportable=0),
    ]
    cnas = [dict(chrom="17", gene="ERBB2", copy_number=3.46, call=1),
            dict(chrom="7", gene="MET", copy_number=2.08, call=0)]
    fusions = [dict(gene_a="NTRK1", gene_b="TPM3", downstream_gene="-", chrom_a=0, chrom_b=0, pos_a=0, pos_b=0,
                    percentage=0, call=0)]
    msi = [dict(runid="SYNRUN_0001", run_sample_id="SYNSAMPLE01", msi_score=1, msi_status="MSS/MSI-L",
                max_maf=0.0502, mean_maf=None)]
    qc = [dict(runid="SYNRUN_0001", run_sample_id="SYNSAMPLE01", category="sample", metric="coverage",
               verbose_name="Coverage", unit="families", value=2100, status="PASS")]
    write_book(d / "Interim_SYNPT-0001_SYN-G360-0001.xlsx", {
        "SNV": (SNV_COLS, snv), "Indels": (INDEL_COLS, indels), "CNAs": (CNA_COLS, cnas),
        "Fusions": (FUSION_COLS, fusions), "MSI": (MSI_COLS, msi), "QC": (QC_COLS, qc)})

    # Lung; older layout with a singular 'Fusion' sheet and an MSI-High call
    snv2 = [
        dict(gene="EGFR", chrom="7", position=55259515, mut_nt="T>G", mut_aa="L858R", cdna="c.2573T>G",
             percentage=8.40, call=1, transcript_id="NM_005228.5", exon="21.0", reporting_category="protein_coding",
             rm_reportable=1),
        dict(gene="EGFR", chrom="7", position=55249071, mut_nt="C>T", mut_aa="T790M", cdna="c.2369C>T",
             percentage=1.35, call=1, transcript_id="NM_005228.5", exon="20.0", reporting_category="protein_coding",
             rm_reportable=1),
    ]
    indels2 = [
        dict(gene="EGFR", chrom="7", position=55242464, mut_nt="AGGAATTAAGAGAAGC>A", mut_aa="E746_A750del",
             cdna="c.2235_2249del", length=15, exon="19.0", type="Deletion", percentage=0.92, call=1,
             transcript_id="NM_005228.5", reporting_category="protein_coding", mut_aa_short="E746_A750del",
             rm_reportable=1),
    ]
    cnas2 = [dict(chrom="7", gene="MET", copy_number=4.12, call=1)]
    fusions2 = [dict(gene_a="EML4", gene_b="ALK", downstream_gene="ALK", chrom_a="2", chrom_b="2",
                     pos_a=42522656, pos_b=29446394, percentage=0.61, call=1)]
    msi2 = [dict(runid="SYNRUN_0002", run_sample_id="SYNSAMPLE02", msi_score=24, msi_status="MSI-High",
                 max_maf=0.084, mean_maf=None)]
    write_book(d / "Interim_SYNPT-0004_SYN-G360-0002.xlsx", {
        "SNV": (SNV_COLS, snv2), "Indels": (INDEL_COLS, indels2), "CNAs": (CNA_COLS, cnas2),
        "Fusion": (FUSION_COLS, fusions2), "MSI": (MSI_COLS, msi2)})


# ---------------------------------------------------------------------------
# Case-metadata sidecar and deliberately invalid inputs
# ---------------------------------------------------------------------------

SIDECAR = """panel_name,report_id,patient_id,disease,oncotree_code,tissue_of_origin,note
FoundationOne,SYN-F1-0001,SYNPT-0001,,COAD,,MRN absent in the vendor file; linked by the curator
FoundationOne,SYN-F1-0002,,,LUAD,,
FoundationOneLiquid,SYN-F1L-0003,,,PRAD,,
GenMineTOP,SYN-GM-0001,,Colon adenocarcinoma,COAD,Colon,
GenMineTOP,SYN-GM-0002,,Lung adenocarcinoma,LUAD,Lung,
Guardant,SYN-G360-0001,,Colon adenocarcinoma,COAD,Colon,disease is not part of the Guardant deliverable
Guardant,SYN-G360-0002,,Lung adenocarcinoma,LUAD,Lung,
"""


SYN_PANELS = {
    # file stem: (panel_name, panel_version, genes with short-variant coverage, CNA genes, fusion genes)
    "synthetic_FoundationOneDx": (
        "FoundationOne", "FoundationOneDx",
        "ALK APC AR ARID1A BRAF BRCA1 BRCA2 CDKN2A EGFR ERBB2 FBXW7 KRAS MET MYC NTRK1 PIK3CA PTEN SMAD4 STK11 TERT TP53 "
        "WHSC1L1",
        "AR CDKN2A EGFR ERBB2 MET MYC", "ALK EGFR MET NTRK1"),
    "synthetic_FoundationOneLiquidDx": (
        "FoundationOneLiquid", "FoundationOneLiquidDx",
        "ALK APC AR ARID1A BRAF BRCA1 BRCA2 EGFR ERBB2 KRAS MET NTRK1 PIK3CA PTEN TP53",
        "AR EGFR ERBB2 MET", "ALK NTRK1"),
    "synthetic_GenMineTOP_TDv6+TRv6": (
        "GenMineTOP", "TDv6+TRv6",
        "ALK APC AR ARID1A BRAF BRCA1 BRCA2 CDKN2A EGFR ERBB2 FBXW7 KRAS MET MYC NTRK1 PIK3CA PTEN SMAD4 STK11 TERT TP53",
        "AR CDKN2A EGFR ERBB2 MET MYC", "ALK EML4 MET NTRK1"),
    "synthetic_Guardant360CDx": (
        "Guardant", "Guardant360 CDx",
        "ALK APC AR ARID1A BRAF BRCA1 BRCA2 EGFR ERBB2 KRAS MET NTRK1 PIK3CA PTEN TP53",
        "ERBB2 MET", "ALK NTRK1"),
}


def write_panels():
    d = DATA / "panels"
    d.mkdir(parents=True, exist_ok=True)
    for stem, (name, version, sv_genes, cna_genes, fus_genes) in SYN_PANELS.items():
        sv_set, cna_set, fus_set = set(sv_genes.split()), set(cna_genes.split()), set(fus_genes.split())
        lines = [f"# panel_name: {name}", f"# panel_version: {version}",
                 "# description: SYNTHETIC gene list for tests - NOT the vendor's official assay content",
                 "# source: tests/make_fixtures.py", "gene\tshort_variants\tcopy_number\trearrangements"]
        for g in sorted(sv_set | cna_set | fus_set):
            lines.append(f"{g}\t{int(g in sv_set)}\t{int(g in cna_set)}\t{int(g in fus_set)}")
        (d / f"{stem}.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_misc():
    (DATA / "case_metadata.csv").write_text(SIDECAR, encoding="utf-8")
    (INVALID / "foundation").mkdir(parents=True, exist_ok=True)
    (INVALID / "guardant").mkdir(parents=True, exist_ok=True)
    (INVALID / "foundation" / "not_a_foundation_report.xml").write_text(
        '<?xml version="1.0"?>\n<report><id>X</id></report>\n', encoding="utf-8")
    (INVALID / "foundation" / "truncated.xml").write_text(
        '<?xml version="1.0"?>\n<rr:ResultsReport xmlns:rr="http://integration.foundationmedicine.com/reporting">\n',
        encoding="utf-8")
    # a valid workbook whose file name carries no identifiers
    write_book(INVALID / "guardant" / "results.xlsx", {
        "SNV": (SNV_COLS, []), "Indels": (INDEL_COLS, []), "CNAs": (CNA_COLS, [])})
    # a workbook with the right name but a missing required sheet
    write_book(INVALID / "guardant" / "Interim_SYNPT-9999_SYN-BAD-0001.xlsx", {"SNV": (SNV_COLS, [])})


if __name__ == "__main__":
    write_foundation()
    write_genminetop()
    write_guardant()
    write_panels()
    write_misc()
    print(f"fixtures written under {DATA} and {INVALID}")
