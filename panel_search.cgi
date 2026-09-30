#!/usr/bin/perl
# panel_search.cgi — cross-vendor search over an OncoUnify database.
#
# Reference implementation of a read-only query interface.  Every query is a
# parameterized SELECT over the canonical schema; see docs/SCHEMA.md.
#
# Configuration (environment, e.g. Apache SetEnv):
#   ONCOUNIFY_DB        path to panels.db             (default /var/www/data/panels.db)
#   ONCOUNIFY_CGI_URL   URL prefix of the CGI scripts (default /cgi-bin)
#   ONCOUNIFY_HTML_URL  URL prefix of search.html     (default /panel)
#   ONCOUNIFY_MAX_ROWS  hard cap on returned rows      (default 100000)
use strict;
use warnings;
use utf8;
use CGI qw(:standard);
use DBI;
use DBD::SQLite ();

binmode(STDOUT, ':encoding(UTF-8)');

my $DB_FILE  = $ENV{ONCOUNIFY_DB}       || '/var/www/data/panels.db';
my $CGI_URL  = $ENV{ONCOUNIFY_CGI_URL}  || '/cgi-bin';
my $HTML_URL = $ENV{ONCOUNIFY_HTML_URL} || '/panel';
my $MAX_ROWS = ($ENV{ONCOUNIFY_MAX_ROWS} // '') =~ /^\d+$/ ? $ENV{ONCOUNIFY_MAX_ROWS} : 100000;

my $q = CGI->new;

sub _p { my $v = $q->param($_[0]); $v = '' unless defined $v; $v =~ s/^\s+|\s+$//g; return $v; }

my $gene             = _p('gene');
my $gene_match       = lc(_p('gene_match')) eq 'substring' ? 'substring' : 'exact';
my $protein_effect   = _p('protein_effect');
my $patient_id       = _p('patient_id');
my $variant_type     = _p('variant_type');
my $consequence      = lc(_p('consequence'));
my $disease          = _p('disease');
my $panel_name       = _p('panel_name');
my $include_uncalled = _p('include_uncalled') eq '1' ? 1 : 0;
my $view_mode        = lc(_p('view_mode')) || 'variant';
my $limit            = _p('limit');
my $want_tsv         = lc(_p('download')) eq 'tsv';

$view_mode = 'variant' unless $view_mode =~ /^(variant|case|patient)$/;
$variant_type = '' unless $variant_type =~ /^(short_variant|cnv|rearrangement|expression)$/;
$consequence = '' unless $consequence =~ /^(protein_altering|missense|nonsense|frameshift|splice|inframe|other|silent|noncoding)$/;
$limit = 500 if $limit !~ /^\d+$/ || $limit < 1;
$limit = $MAX_ROWS if $limit > $MAX_ROWS;

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
sub _safe { defined $_[0] ? $_[0] : '' }

sub _h {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/&/&amp;/g; $s =~ s/</&lt;/g; $s =~ s/>/&gt;/g; $s =~ s/"/&quot;/g; $s =~ s/'/&#39;/g;
    return $s;
}

# TSV cell: tabs and line breaks inside a value would shift the columns
sub _t { my $s = _safe($_[0]); $s =~ s/[\t\r\n]+/ /g; return $s; }

# Escape LIKE metacharacters; used with "LIKE ? ESCAPE '\'"
sub _like { my $s = shift; $s =~ s/([\\%_])/\\$1/g; return $s; }

sub _fail {
    my ($status, $msg) = @_;
    print $q->header(-type => 'text/plain', -charset => 'utf-8', -status => $status);
    print "$msg\n";
    exit;
}

my %AA3 = (Ala=>'A',Arg=>'R',Asn=>'N',Asp=>'D',Cys=>'C',Gln=>'Q',Glu=>'E',Gly=>'G',His=>'H',Ile=>'I',
           Leu=>'L',Lys=>'K',Met=>'M',Phe=>'F',Pro=>'P',Ser=>'S',Thr=>'T',Trp=>'W',Tyr=>'Y',Val=>'V',
           Ter=>'*',Sec=>'U',Pyl=>'O',Xaa=>'X');
my $AA3_RE = join('|', sort { length($b) <=> length($a) } keys %AA3);

# Same normalization as oncounify_core.canonical_protein(): p.Gly12Asp, pG12D, G12D -> p.G12D
sub canonical_protein {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/\s+//g;
    return '' if $s eq '';
    $s =~ s/^\(?p\.//;
    $s =~ s/^p(?=[A-Z*(])//;
    $s =~ s/^\(+|\)+$//g;
    $s =~ s/($AA3_RE)/$AA3{$1}/g;
    $s =~ s/(?<=\d)X$/*/;
    $s =~ s/fsX/fs*/g;
    if ($s =~ /^([A-Z*])(\d+)([A-Z*])$/ && $1 eq $3) { $s = "$1$2=" }
    return "p.$s";
}

sub _split_tokens { return grep { length } map { s/^\s+|\s+$//gr } split /,/, ($_[0] // '') }

sub _view_label { $_[0] eq 'patient' ? 'per patient' : $_[0] eq 'case' ? 'per case' : 'per variant' }

sub _type_badge {
    my ($t) = @_;
    my %m = (short_variant => ['SNV/indel', 'tag-snv'], cnv => ['CNV', 'tag-cnv'],
             rearrangement => ['Rearrangement', 'tag-rearr'], expression => ['Expression', 'tag-exp']);
    my ($label, $cls) = @{ $m{_safe($t)} || [_safe($t), 'tag-default'] };
    return qq{<span class="tag $cls">} . _h($label) . qq{</span>};
}

# ---------------------------------------------------------------------------
# WHERE clause (static fragments + bound parameters only)
# ---------------------------------------------------------------------------
my (@where, @bind);

if ($gene ne '') {
    my @tok = _split_tokens($gene);
    if ($gene_match eq 'exact') {
        my %seen; my @vals = grep { !$seen{$_}++ } map { ($_, uc $_) } @tok;
        my $ph = join(',', ('?') x @vals);
        # previous HGNC symbols (e.g. WHSC1L1) resolve to the current symbol (NSD3)
        my $cur = "SELECT symbol FROM gene_symbol_map WHERE previous_symbol IN ($ph)";
        push @where, "(variants.gene IN ($ph) OR variants.other_gene IN ($ph) "
                   . "OR variants.gene IN ($cur) OR variants.other_gene IN ($cur))";
        push @bind, @vals, @vals, @vals, @vals;
    } else {
        push @where, '(' . join(' OR ', map { "variants.gene LIKE ? ESCAPE '\\' OR variants.other_gene LIKE ? ESCAPE '\\'" } @tok) . ')';
        push @bind, map { ('%' . _like($_) . '%') x 2 } @tok;
    }
}

if ($protein_effect ne '') {
    my $cp = canonical_protein($protein_effect);
    if ($cp =~ /^p\.([A-Z*])(\d+)$/) {                       # residue only, e.g. G12 -> any change at G12
        push @where, "variants.hgvs_p GLOB ?";
        push @bind, "p.$1$2\[^0-9\]*";                       # SQLite GLOB negates with '^'
    } elsif ($cp =~ /^p\.([A-Z*])(\d+)[A-Z*]?fs/) {           # any frameshift starting at the residue
        push @where, "(variants.hgvs_p GLOB ? OR variants.hgvs_p GLOB ?)";
        push @bind, "p.$1$2fs*", "p.$1$2\[A-Z*\]fs*";
    } else {
        push @where, 'variants.hgvs_p = ?';
        push @bind, $cp;
    }
}

if ($patient_id ne '') {
    push @where, "cases.patient_id LIKE ? ESCAPE '\\'";
    push @bind, '%' . _like($patient_id) . '%';
}

if ($variant_type ne '') {
    push @where, 'variants.variant_type = ?';
    push @bind, $variant_type;
}

if ($consequence ne '') {
    if ($consequence eq 'protein_altering') {
        push @where, "variants.functional_effect IN (SELECT term FROM so_terms WHERE display_group NOT IN ('silent','noncoding'))";
    } else {
        push @where, 'variants.functional_effect IN (SELECT term FROM so_terms WHERE display_group = ?)';
        push @bind, $consequence;
    }
}

if ($disease ne '') {
    my $pat = '%' . _like($disease) . '%';
    push @where, '(' . join(' OR ', map { "cases.$_ LIKE ? ESCAPE '\\'" }
                               qw(disease disease_ontology oncotree_code tissue_of_origin pathology_diagnosis)) . ')';
    push @bind, ($pat) x 5;
}

if ($panel_name ne '') {
    push @where, 'cases.panel_name = ?';
    push @bind, $panel_name;
}

push @where, "(variants.status IS NULL OR variants.status != 'not_called')" unless $include_uncalled;

my $where_sql = @where ? join("\n  AND ", @where) : '1=1';

# ---------------------------------------------------------------------------
# queries
# ---------------------------------------------------------------------------
my $NH_SUMMARY = <<'SQL';
(SELECT group_concat(nh.organism || '(' || printf('%.0f', COALESCE(nh.reads_per_million, 0.0)) || ')', ', ')
   FROM non_human_contents nh WHERE nh.case_id = cases.case_id)
SQL
my $BM_SUMMARY = <<'SQL';
(SELECT group_concat(b.name || COALESCE(' ' || b.call, '') ||
                     CASE WHEN b.value IS NOT NULL THEN ' (' || b.value || COALESCE(' ' || b.unit, '') || ')' ELSE '' END, '; ')
   FROM biomarkers b WHERE b.case_id = cases.case_id)
SQL

my $sql;
if ($view_mode eq 'variant') {
    $sql = <<"SQL";
SELECT cases.case_id, cases.panel_name, cases.report_id, cases.patient_id, cases.date, cases.genome_build,
       cases.disease, cases.oncotree_code, cases.tissue_of_origin, cases.pathology_diagnosis,
       variants.gene, variants.other_gene, variants.variant_type, variants.variant_subtype,
       variants.chrom, variants.pos, variants.chrom2, variants.pos2, variants.ref, variants.alt,
       variants.transcript, variants.hgvs_c, variants.hgvs_p, variants.protein_effect,
       variants.functional_effect, variants.functional_effect_so, variants.functional_effect_source,
       so_terms.display_group,
       variants.status, variants.origin, variants.classification,
       variants.allele_fraction, variants.depth, variants.copy_number, variants.cnv_type, variants.cnv_ratio,
       variants.in_frame, variants.clinvar_id, variants.clinvar_sig, variants.clinvar_match,
       variants.maf_1kg, variants.maf_hgvd, variants.maf_tommo,
       variants.tpm, variants.tpm_normal_mean, variants.tpm_normal_sd,
       $NH_SUMMARY AS non_human_summary
FROM variants
JOIN cases ON variants.case_id = cases.case_id
LEFT JOIN so_terms ON so_terms.term = variants.functional_effect
WHERE $where_sql
ORDER BY cases.case_id, variants.gene, variants.pos
LIMIT ?
SQL
} elsif ($view_mode eq 'case') {
    $sql = <<"SQL";
SELECT cases.case_id, cases.panel_name, cases.report_id, cases.patient_id, cases.date, cases.genome_build,
       cases.disease, cases.oncotree_code, cases.tissue_of_origin, cases.pathology_diagnosis,
       COUNT(*) AS matched_variant_count,
       COUNT(DISTINCT variants.gene) AS matched_gene_count,
       group_concat(DISTINCT variants.gene) AS matched_genes,
       group_concat(DISTINCT variants.variant_type) AS matched_variant_types,
       $BM_SUMMARY AS biomarker_summary,
       $NH_SUMMARY AS non_human_summary
FROM variants
JOIN cases ON variants.case_id = cases.case_id
WHERE $where_sql
GROUP BY cases.case_id
ORDER BY COALESCE(cases.date, '') DESC, cases.case_id DESC
LIMIT ?
SQL
} else {
    $sql = <<"SQL";
WITH f AS (
    SELECT cases.case_id, cases.panel_name, cases.report_id, cases.patient_id, cases.date,
           cases.disease, cases.oncotree_code, cases.tissue_of_origin, cases.pathology_diagnosis,
           variants.gene
    FROM variants
    JOIN cases ON variants.case_id = cases.case_id
    WHERE $where_sql
)
SELECT CASE WHEN f.patient_id IS NOT NULL AND TRIM(f.patient_id) != '' THEN f.patient_id
            ELSE '[case:' || f.case_id || ']' END AS patient_group,
       COUNT(DISTINCT f.case_id) AS case_count,
       MAX(COALESCE(f.date, '')) AS latest_date,
       group_concat(DISTINCT f.panel_name) AS panel_names,
       group_concat(DISTINCT f.disease) AS diseases,
       group_concat(DISTINCT f.oncotree_code) AS oncotree_codes,
       group_concat(DISTINCT f.tissue_of_origin) AS tissues,
       group_concat(DISTINCT f.pathology_diagnosis) AS pathologies,
       COUNT(*) AS matched_variant_count,
       COUNT(DISTINCT f.gene) AS matched_gene_count,
       group_concat(DISTINCT f.gene) AS matched_genes,
       group_concat(DISTINCT f.case_id || '|' || COALESCE(f.report_id, '') || '|' || f.panel_name) AS case_links
FROM f
GROUP BY patient_group
ORDER BY latest_date DESC, matched_variant_count DESC, patient_group
LIMIT ?
SQL
}

_fail('500 Internal Server Error', "OncoUnify database not found or not readable: set ONCOUNIFY_DB") unless -r $DB_FILE;
my $dbh = eval {
    DBI->connect("dbi:SQLite:dbname=$DB_FILE", '', '', {
        RaiseError => 1, PrintError => 0, sqlite_unicode => 1, AutoCommit => 1, ReadOnly => 1,
        sqlite_open_flags => DBD::SQLite::OPEN_READONLY(),
    });
} or _fail('500 Internal Server Error', 'Cannot open the OncoUnify database.');

my $sth = $dbh->prepare($sql);
$sth->execute(@bind, $limit);
my @rows;
while (my $r = $sth->fetchrow_hashref) { push @rows, $r; }
$sth->finish;
$dbh->disconnect;

sub _loc {
    my ($r) = @_;
    return '' unless _safe($r->{chrom}) ne '';
    my $s = $r->{chrom} . (defined $r->{pos} ? ':' . $r->{pos} : '');
    if (_safe($r->{chrom2}) ne '' || defined $r->{pos2}) {
        $s .= ' / ' . (_safe($r->{chrom2}) ne '' ? $r->{chrom2} . ':' : '') . _safe($r->{pos2});
    }
    return $s;
}
sub _num { my ($v, $f) = @_; defined $v ? sprintf($f, $v) : '' }

# ---------------------------------------------------------------------------
# TSV
# ---------------------------------------------------------------------------
if ($want_tsv) {
    print "Content-Type: text/tab-separated-values; charset=utf-8\r\n";
    print qq{Content-Disposition: attachment; filename="oncounify_${view_mode}.tsv"\r\n\r\n};
    my @cols;
    if ($view_mode eq 'variant') {
        @cols = qw(case_id panel_name report_id patient_id date genome_build disease oncotree_code tissue_of_origin
                   pathology_diagnosis gene other_gene variant_type variant_subtype chrom pos chrom2 pos2 ref alt
                   transcript hgvs_c hgvs_p protein_effect functional_effect functional_effect_so
                   functional_effect_source status origin classification allele_fraction depth copy_number cnv_type
                   cnv_ratio in_frame clinvar_id clinvar_sig clinvar_match maf_1kg maf_hgvd maf_tommo tpm
                   tpm_normal_mean tpm_normal_sd non_human_summary);
    } elsif ($view_mode eq 'case') {
        @cols = qw(case_id panel_name report_id patient_id date genome_build disease oncotree_code tissue_of_origin
                   pathology_diagnosis matched_variant_count matched_gene_count matched_genes matched_variant_types
                   biomarker_summary non_human_summary);
    } else {
        @cols = qw(patient_group case_count latest_date panel_names diseases oncotree_codes tissues pathologies
                   matched_variant_count matched_gene_count matched_genes case_links);
    }
    print join("\t", @cols), "\n";
    for my $r (@rows) { print join("\t", map { _t($r->{$_}) } @cols), "\n"; }
    exit;
}

# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
print $q->header(-type => 'text/html', -charset => 'utf-8');
my $row_count = scalar @rows;

print <<'HTML';
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>OncoUnify search results</title>
  <style>
    * { box-sizing: border-box; }
    body { margin: 0; padding: 1.2rem; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           background: radial-gradient(circle at top, #eff6ff 0, #e0f2fe 35%, #f3f4f6 100%); color: #111827; }
    a { color: #2563eb; text-decoration: none; } a:hover { text-decoration: underline; }
    .card { background: #ffffffee; border-radius: 1rem; padding: 1.2rem 1.4rem 1rem;
            box-shadow: 0 18px 45px rgba(15,23,42,.16), 0 0 0 1px rgba(255,255,255,.85); }
    .card-header { display: flex; justify-content: space-between; align-items: baseline; gap: 1rem; margin-bottom: .75rem; }
    h1 { margin: 0; font-size: 1.35rem; font-weight: 650; letter-spacing: .03em; }
    .subtitle { font-size: .8rem; color: #6b7280; }
    .meta-block { text-align: right; font-size: .78rem; color: #6b7280; }
    .filter-summary { margin: .4rem 0 .5rem; padding: .45rem .75rem; border-radius: 999px; background: #f9fafb;
                      font-size: .78rem; color: #4b5563; display: flex; flex-wrap: wrap; gap: .4rem; align-items: center; }
    .filter-pill { display: inline-flex; gap: .25rem; padding: .1rem .5rem; border-radius: 999px; background: #e5f3ff; color: #1d4ed8; }
    .filter-pill .key { font-weight: 600; }
    .download-form { margin-left: auto; }
    .download-button { border: none; border-radius: 999px; padding: .25rem .8rem; font-size: .75rem; font-weight: 600;
                       text-transform: uppercase; cursor: pointer; background: linear-gradient(135deg,#0ea5e9,#2563eb); color: #fff; }
    .table-wrapper { margin-top: .4rem; border-radius: .7rem; border: 1px solid #e5e7eb; background: #f9fafb;
                     max-height: 72vh; overflow: auto; }
    table { width: max-content; min-width: 100%; border-collapse: collapse; font-size: .8rem; }
    th, td { padding: .35rem .45rem; border-bottom: 1px solid #e5e7eb; border-right: 1px solid #e5e7eb; vertical-align: top; }
    th { position: sticky; top: 0; z-index: 2; font-weight: 600; color: #374151; text-align: left; white-space: nowrap;
         background: linear-gradient(120deg,#eff6ff,#e0f2fe); }
    tbody tr:nth-child(even) td { background: #f3f4f6; } tbody tr:hover td { background: #e0f2fe; }
    .num { text-align: right; white-space: nowrap; }
    .mono, .gene { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
    .gene { font-weight: 600; }
    .small { font-size: .72rem; color: #6b7280; }
    .wrap { max-width: 360px; min-width: 160px; white-space: normal; word-break: break-word; }
    .tag { display: inline-flex; justify-content: center; min-width: 76px; padding: .1rem .4rem; border-radius: 999px;
           font-size: .7rem; font-weight: 600; text-transform: uppercase; border: 1px solid transparent; }
    .tag-snv { background: #eef2ff; color: #4f46e5; border-color: #c7d2fe; }
    .tag-cnv { background: #fef3c7; color: #92400e; border-color: #fde68a; }
    .tag-rearr { background: #fee2e2; color: #b91c1c; border-color: #fecaca; }
    .tag-exp { background: #ecfdf5; color: #047857; border-color: #a7f3d0; }
    .tag-default { background: #f3f4f6; color: #4b5563; border-color: #e5e7eb; }
    .cq { display: inline-block; width: .6rem; height: .6rem; border-radius: 999px; margin-right: .3rem; }
    .cq-missense { background: #4f46e5; } .cq-nonsense { background: #ef4444; } .cq-frameshift { background: #ec4899; }
    .cq-splice { background: #0ea5e9; } .cq-inframe { background: #f59e0b; } .cq-other { background: #9ca3af; }
    .cq-silent { background: #d1d5db; } .cq-noncoding { background: #e5e7eb; }
    .footer { margin-top: .6rem; font-size: .78rem; color: #6b7280; display: flex; justify-content: space-between; }
    .empty-note { padding: 1rem; color: #6b7280; }
  </style>
</head>
<body>
  <div class="card">
    <div class="card-header">
      <div>
        <h1>OncoUnify search results</h1>
        <div class="subtitle">Cross-vendor results over the canonical OncoUnify schema.</div>
      </div>
HTML
print qq{      <div class="meta-block"><a href="} . _h("$HTML_URL/search.html") . qq{">Modify search</a></div>\n    </div>\n};

my @pills = (['View', _view_label($view_mode)]);
push @pills, ['Gene', $gene . ($gene_match eq 'substring' ? ' (substring)' : '')] if $gene ne '';
push @pills, ['Protein', $protein_effect . ' → ' . canonical_protein($protein_effect)] if $protein_effect ne '';
push @pills, ['Consequence', $consequence] if $consequence ne '';
push @pills, ['Patient ID', $patient_id] if $patient_id ne '';
push @pills, ['Type', $variant_type] if $variant_type ne '';
push @pills, ['Disease/Tissue', $disease] if $disease ne '';
push @pills, ['Panel', $panel_name] if $panel_name ne '';
push @pills, ['Uncalled rows', 'included'] if $include_uncalled;
push @pills, ['Limit', "$limit rows"];

print qq{    <div class="filter-summary">\n};
print qq{      <span class="filter-pill"><span class="key">} . _h($_->[0]) . qq{</span><span>} . _h($_->[1]) . qq{</span></span>\n} for @pills;
print qq{      <form class="download-form" method="get" action="} . _h("$CGI_URL/panel_search.cgi") . qq{">\n};
for my $k (qw(gene gene_match protein_effect consequence patient_id variant_type disease panel_name include_uncalled view_mode limit)) {
    my $v = _p($k);
    next if $v eq '';
    print qq{        <input type="hidden" name="} . _h($k) . qq{" value="} . _h($v) . qq{">\n};
}
print qq{        <input type="hidden" name="download" value="tsv"><button type="submit" class="download-button">Download TSV</button>\n};
print qq{      </form>\n    </div>\n    <div class="table-wrapper">\n};

sub _case_link {
    my ($cid, $label) = @_;
    return qq{<a href="} . _h("$CGI_URL/case_detail.cgi?case_id=$cid") . qq{">} . _h($label) . qq{</a>};
}

if ($view_mode eq 'variant') {
    print qq{<table><thead><tr>} . join('', map { "<th>$_</th>" }
        ('case', 'panel', 'report', 'patient', 'disease / OncoTree', 'gene', 'partner', 'type', 'subtype',
         'location (build)', 'cDNA', 'protein', 'consequence (SO)', 'status', 'origin', 'AF', 'depth',
         'copy#', 'ClinVar', 'MAF (1KG / HGVD / ToMMo)', 'TPM tumor', 'TPM normal (mean±SD)', 'non-human'))
        . qq{</tr></thead><tbody>\n};
    for my $r (@rows) {
        my $dis = join(' / ', grep { length } map { _safe($_) } ($r->{disease} // $r->{pathology_diagnosis}, $r->{oncotree_code}));
        my $loc = _loc($r);
        $loc .= ' (' . _safe($r->{genome_build}) . ')' if $loc ne '';
        my $protein = _safe($r->{hgvs_p});
        my $prot_cell = $protein ne '' ? qq{<span title="vendor: } . _h($r->{protein_effect}) . qq{">} . _h($protein) . '</span>'
                                       : '<span class="small">' . _h($r->{protein_effect}) . '</span>';
        my $cq = '';
        if (_safe($r->{functional_effect}) ne '') {
            my $g = _safe($r->{display_group}) || 'other';
            $cq = qq{<span class="cq cq-$g"></span>} . _h($r->{functional_effect})
                . qq{<br><span class="small">} . _h($r->{functional_effect_so}) . ' · ' . _h($r->{functional_effect_source}) . '</span>';
        }
        my $clin = '';
        if (_safe($r->{clinvar_id}) ne '') {
            $clin = qq{<a href="} . _h('https://www.ncbi.nlm.nih.gov/clinvar/variation/' . $r->{clinvar_id}) . qq{" target="_blank" rel="noopener">}
                  . _h($r->{clinvar_id}) . '</a>';
            $clin .= '<br><span class="small">' . _h($r->{clinvar_sig}) . '</span>' if _safe($r->{clinvar_sig}) ne '';
        }
        my $maf = join('<br>', map { _h($_) } grep { defined }
            (defined $r->{maf_1kg} ? '1KG ' . sprintf('%.4g', $r->{maf_1kg}) : undef,
             defined $r->{maf_hgvd} ? 'HGVD ' . sprintf('%.4g', $r->{maf_hgvd}) : undef,
             defined $r->{maf_tommo} ? 'ToMMo ' . sprintf('%.4g', $r->{maf_tommo}) : undef));
        my $tpmn = defined $r->{tpm_normal_mean}
            ? sprintf('%.2f', $r->{tpm_normal_mean}) . (defined $r->{tpm_normal_sd} ? ' ± ' . sprintf('%.2f', $r->{tpm_normal_sd}) : '') : '';
        my $cn = _num($r->{copy_number}, '%.2f');
        $cn .= ' ' . $r->{cnv_type} if _safe($r->{cnv_type}) ne '';
        print '<tr>',
            '<td class="num">', _h($r->{case_id}), '</td>',
            '<td>', _h($r->{panel_name}), '</td>',
            '<td>', _case_link($r->{case_id}, _safe($r->{report_id})), '</td>',
            '<td class="mono">', _h($r->{patient_id}), '</td>',
            '<td class="wrap">', _h($dis), '</td>',
            '<td class="gene">', _h($r->{gene}), '</td>',
            '<td class="gene">', _h($r->{other_gene}), '</td>',
            '<td>', _type_badge($r->{variant_type}), '</td>',
            '<td>', _h($r->{variant_subtype}), '</td>',
            '<td class="mono">', _h($loc), '</td>',
            '<td class="mono">', _h($r->{hgvs_c}), '<br><span class="small">', _h($r->{transcript}), '</span></td>',
            '<td class="mono">', $prot_cell, '</td>',
            '<td>', $cq, '</td>',
            '<td>', _h($r->{status}), '</td>',
            '<td>', _h($r->{origin}), '</td>',
            '<td class="num">', _h(_num($r->{allele_fraction}, '%.4f')), '</td>',
            '<td class="num">', _h($r->{depth}), '</td>',
            '<td class="num">', _h($cn), '</td>',
            '<td>', $clin, '</td>',
            '<td class="wrap">', $maf, '</td>',
            '<td class="num">', _h(_num($r->{tpm}, '%.2f')), '</td>',
            '<td class="num">', _h($tpmn), '</td>',
            '<td class="wrap">', _h($r->{non_human_summary}), '</td>',
            "</tr>\n";
    }
    print "</tbody></table>\n";
} elsif ($view_mode eq 'case') {
    print qq{<table><thead><tr>} . join('', map { "<th>$_</th>" }
        ('case', 'panel', 'report', 'patient', 'date', 'build', 'disease', 'OncoTree', 'tissue', 'pathology',
         'matched variants', 'matched genes', 'gene list', 'variant types', 'TMB / MSI', 'non-human'))
        . qq{</tr></thead><tbody>\n};
    for my $r (@rows) {
        print '<tr>',
            '<td class="num">', _h($r->{case_id}), '</td>',
            '<td>', _h($r->{panel_name}), '</td>',
            '<td>', _case_link($r->{case_id}, _safe($r->{report_id})), '</td>',
            '<td class="mono">', _h($r->{patient_id}), '</td>',
            '<td class="mono">', _h($r->{date}), '</td>',
            '<td>', _h($r->{genome_build}), '</td>',
            '<td class="wrap">', _h($r->{disease}), '</td>',
            '<td>', _h($r->{oncotree_code}), '</td>',
            '<td>', _h($r->{tissue_of_origin}), '</td>',
            '<td class="wrap">', _h($r->{pathology_diagnosis}), '</td>',
            '<td class="num">', _h($r->{matched_variant_count}), '</td>',
            '<td class="num">', _h($r->{matched_gene_count}), '</td>',
            '<td class="wrap gene">', _h($r->{matched_genes}), '</td>',
            '<td>', _h($r->{matched_variant_types}), '</td>',
            '<td class="wrap">', _h($r->{biomarker_summary}), '</td>',
            '<td class="wrap">', _h($r->{non_human_summary}), '</td>',
            "</tr>\n";
    }
    print "</tbody></table>\n";
} else {
    print qq{<table><thead><tr>} . join('', map { "<th>$_</th>" }
        ('patient', 'cases', 'latest date', 'reports', 'panels', 'disease', 'OncoTree', 'tissue', 'pathology',
         'matched variants', 'matched genes', 'gene list'))
        . qq{</tr></thead><tbody>\n};
    for my $r (@rows) {
        my @links;
        for my $item (split /,/, _safe($r->{case_links})) {
            my ($cid, $rid, $pn) = split /\|/, $item, 3;
            next unless defined $cid && $cid =~ /^\d+$/;
            push @links, _case_link($cid, (_safe($rid) ne '' ? $rid : "case:$cid")) . ' <span class="small">' . _h($pn) . '</span>';
        }
        print '<tr>',
            '<td class="mono">', _h($r->{patient_group}), '</td>',
            '<td class="num">', _h($r->{case_count}), '</td>',
            '<td class="mono">', _h($r->{latest_date}), '</td>',
            '<td class="wrap">', join('<br>', @links), '</td>',
            '<td class="wrap">', _h($r->{panel_names}), '</td>',
            '<td class="wrap">', _h($r->{diseases}), '</td>',
            '<td>', _h($r->{oncotree_codes}), '</td>',
            '<td class="wrap">', _h($r->{tissues}), '</td>',
            '<td class="wrap">', _h($r->{pathologies}), '</td>',
            '<td class="num">', _h($r->{matched_variant_count}), '</td>',
            '<td class="num">', _h($r->{matched_gene_count}), '</td>',
            '<td class="wrap gene">', _h($r->{matched_genes}), '</td>',
            "</tr>\n";
    }
    print "</tbody></table>\n";
}

print qq{<div class="empty-note">No matching results.</div>} if $row_count == 0;
print qq{    </div>\n    <div class="footer"><div><strong>$row_count</strong> rows shown (limit } . _h($limit)
    . qq{).</div><div>} . _h(_view_label($view_mode)) . qq{ view</div></div>\n  </div>\n</body>\n</html>\n};
