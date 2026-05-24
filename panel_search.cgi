#!/usr/bin/perl
use strict;
use warnings;
use utf8;
use CGI qw(:standard);
use DBI;
use DBD::SQLite ();

binmode(STDOUT, ':encoding(UTF-8)');

my $DB_FILE = '/var/www/data/panels.db';

my $q = CGI->new;

my $gene           = $q->param('gene')           // '';
my $protein_effect = $q->param('protein_effect') // '';
my $patient_id     = $q->param('patient_id')     // '';
my $variant_type   = $q->param('variant_type')   // '';
my $disease        = $q->param('disease')        // '';
my $panel_name     = $q->param('panel_name')     // '';
my $view_mode      = lc($q->param('view_mode')   // 'variant');
my $limit          = $q->param('limit')          // 100000;
my $download       = $q->param('download')       // '';
my $want_tsv       = (defined $download && lc($download) eq 'tsv');

for ($gene, $protein_effect, $patient_id, $variant_type, $disease, $panel_name, $view_mode) {
    $_ //= '';
    s/^\s+|\s+$//g;
}

$view_mode = 'variant' unless $view_mode =~ /^(variant|case|patient)$/;
$limit = 200 if !$limit || $limit !~ /^\d+$/;
$limit = 100000 if $limit > 100000;

my $dbh = DBI->connect(
    "dbi:SQLite:dbname=$DB_FILE",
    "",
    "",
    {
        RaiseError        => 1,
        sqlite_unicode    => 1,
        AutoCommit        => 1,
        ReadOnly          => 1,
        sqlite_open_flags => DBD::SQLite::OPEN_READONLY(),
    }
) or die $DBI::errstr;

sub _safe {
    my ($v) = @_;
    return defined $v ? $v : '';
}

sub _h {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/&/&amp;/g;
    $s =~ s/</&lt;/g;
    $s =~ s/>/&gt;/g;
    $s =~ s/"/&quot;/g;
    $s =~ s/'/&#39;/g;
    return $s;
}

sub _view_label {
    my ($v) = @_;
    return 'per patient' if $v eq 'patient';
    return 'per case'    if $v eq 'case';
    return 'per variant';
}

sub _type_badge {
    my ($t) = @_;
    $t ||= '';
    my $label = $t;
    my $cls   = 'tag-default';
    if ($t eq 'short_variant') {
        $label = 'SNV/indel';
        $cls   = 'tag-snv';
    } elsif ($t eq 'cnv') {
        $label = 'CNV';
        $cls   = 'tag-cnv';
    } elsif ($t eq 'rearrangement') {
        $label = 'Rearrangement';
        $cls   = 'tag-rearr';
    } elsif ($t eq 'expression') {
        $label = 'Expression';
        $cls   = 'tag-exp';
    } elsif ($t eq 'biomarker') {
        $label = 'Biomarker';
        $cls   = 'tag-bm';
    }
    return qq{<span class="tag $cls">} . _h($label) . qq{</span>};
}

sub _split_csv_tokens {
    my ($s) = @_;
    return grep { length($_) } map { s/^\s+|\s+$//gr } split /,/, ($s // '');
}

sub _build_where_and_bind {
    my @where;
    my @bind;

    if ($gene) {
        my @genes = _split_csv_tokens($gene);
        if (@genes) {
            push @where, '(' . join(' OR ', map { 'variants.gene LIKE ?' } @genes) . ')';
            push @bind, map { '%' . $_ . '%' } @genes;
        }
    }

    if ($protein_effect) {
        my $raw = $protein_effect;
        my $norm = $raw;
        $norm =~ s/^\s*p\.?\s*//i;
        $norm =~ s/\s+//g;

        my @patterns;
        push @patterns, '%' . $raw  . '%';
        push @patterns, '%' . $norm . '%';
        push @patterns, '%p.' . $norm . '%';
        push @patterns, '%p'  . $norm . '%';
        my %seen;
        @patterns = grep { !$seen{$_}++ } @patterns;

        push @where, '(' . join(' OR ', map { 'variants.protein_effect LIKE ?' } @patterns) . ')';
        push @bind, @patterns;
    }

    if ($patient_id) {
        push @where, 'cases.patient_id LIKE ?';
        push @bind, '%' . $patient_id . '%';
    }

    if ($variant_type) {
        push @where, 'variants.variant_type = ?';
        push @bind, $variant_type;
    }

    if ($disease) {
        my $pat = '%' . $disease . '%';
        push @where, '(cases.disease LIKE ? OR cases.tissue_of_origin LIKE ? OR cases.pathology_diagnosis LIKE ?)';
        push @bind, ($pat, $pat, $pat);
    }

    if ($panel_name) {
        push @where, 'cases.panel_name = ?';
        push @bind, $panel_name;
    }

    my $where_sql = @where ? join(' AND ', @where) : '1=1';
    return ($where_sql, @bind);
}

sub _render_case_links {
    my ($raw) = @_;
    return '' unless defined $raw && length $raw;
    my @items;
    for my $pair (split /,/, $raw) {
        next unless length $pair;
        my ($cid, $rid) = split /\|/, $pair, 2;
        next unless defined $cid && $cid ne '';
        my $label = (defined $rid && $rid ne '') ? $rid : ('case:' . $cid);
        push @items, qq{<a href="/cgi-bin/case_detail.cgi?case_id=} . _h($cid) . qq{">} . _h($label) . qq{</a>};
    }
    return join('<br>', @items);
}

my ($where_sql, @bind) = _build_where_and_bind();

my ($sql, @exec_bind);

if ($view_mode eq 'variant') {
    $sql = <<"SQL";
SELECT
    cases.case_id,
    cases.panel_name,
    cases.report_id,
    cases.patient_id,
    cases.disease,
    cases.tissue_of_origin,
    cases.pathology_diagnosis,
    variants.gene,
    variants.variant_type,
    variants.variant_subtype,
    variants.cds_effect,
    variants.protein_effect,
    variants.strand,
    variants.transcript,
    variants.functional_effect,
    variants.status,
    variants.origin,
    variants.classification,
    variants.allele_fraction,
    variants.depth,
    variants.copy_number,
    variants.cnv_ratio,
    variants.clinvar_id,
    variants.clinvar_sig,
    variants.clinvar_match,
    variants.maf_1kg,
    variants.maf_hgvd,
    variants.maf_tommo,
    variants.tpm,
    variants.tpm_normal_mean,
    variants.tpm_normal_sd,
    (
      SELECT group_concat(
               nh.organism || '(' || printf('%.0f', COALESCE(nh.reads_per_million, 0.0)) || ')',
               ', '
             )
      FROM non_human_contents nh
      WHERE nh.case_id = cases.case_id
    ) AS non_human_summary
FROM variants
JOIN cases ON variants.case_id = cases.case_id
WHERE $where_sql
ORDER BY cases.case_id, variants.gene
LIMIT ?
SQL
    @exec_bind = (@bind, $limit);
}
elsif ($view_mode eq 'case') {
    $sql = <<"SQL";
WITH filtered AS (
    SELECT
        cases.case_id,
        cases.panel_name,
        cases.report_id,
        cases.patient_id,
        cases.date,
        cases.disease,
        cases.tissue_of_origin,
        cases.pathology_diagnosis,
        variants.gene,
        variants.variant_type,
        variants.variant_subtype,
        variants.protein_effect,
        variants.functional_effect,
        variants.status
    FROM variants
    JOIN cases ON variants.case_id = cases.case_id
    WHERE $where_sql
)
SELECT
    f.case_id,
    f.panel_name,
    f.report_id,
    f.patient_id,
    f.date,
    f.disease,
    f.tissue_of_origin,
    f.pathology_diagnosis,
    COUNT(*) AS matched_variant_count,
    COUNT(DISTINCT f.gene) AS matched_gene_count,
    group_concat(DISTINCT f.gene) AS matched_genes,
    group_concat(DISTINCT f.variant_type) AS matched_variant_types,
    (
      SELECT group_concat(
               nh.organism || '(' || printf('%.0f', COALESCE(nh.reads_per_million, 0.0)) || ')',
               ', '
             )
      FROM non_human_contents nh
      WHERE nh.case_id = f.case_id
    ) AS non_human_summary
FROM filtered f
GROUP BY
    f.case_id, f.panel_name, f.report_id, f.patient_id, f.date,
    f.disease, f.tissue_of_origin, f.pathology_diagnosis
ORDER BY COALESCE(f.date, '') DESC, f.case_id DESC
LIMIT ?
SQL
    @exec_bind = (@bind, $limit);
}
else {
    $sql = <<"SQL";
WITH filtered AS (
    SELECT
        cases.case_id,
        cases.panel_name,
        cases.report_id,
        cases.patient_id,
        cases.date,
        cases.disease,
        cases.tissue_of_origin,
        cases.pathology_diagnosis,
        variants.gene,
        variants.variant_type,
        variants.variant_subtype,
        variants.protein_effect,
        variants.functional_effect,
        variants.status
    FROM variants
    JOIN cases ON variants.case_id = cases.case_id
    WHERE $where_sql
)
SELECT
    CASE
      WHEN f.patient_id IS NOT NULL AND TRIM(f.patient_id) != '' THEN f.patient_id
      ELSE '[case:' || f.case_id || ']'
    END AS patient_group,
    CASE
      WHEN f.patient_id IS NOT NULL AND TRIM(f.patient_id) != '' THEN f.patient_id
      ELSE ''
    END AS patient_id_display,
    COUNT(DISTINCT f.case_id) AS case_count,
    MAX(COALESCE(f.date, '')) AS latest_date,
    group_concat(DISTINCT f.report_id) AS report_ids,
    group_concat(DISTINCT f.panel_name) AS panel_names,
    group_concat(DISTINCT f.disease) AS diseases,
    group_concat(DISTINCT f.tissue_of_origin) AS tissues,
    group_concat(DISTINCT f.pathology_diagnosis) AS pathologies,
    COUNT(*) AS matched_variant_count,
    COUNT(DISTINCT f.gene) AS matched_gene_count,
    group_concat(DISTINCT f.gene) AS matched_genes,
    group_concat(DISTINCT f.case_id || '|' || COALESCE(f.report_id, '')) AS case_links
FROM filtered f
GROUP BY patient_group
ORDER BY latest_date DESC, matched_variant_count DESC, patient_group
LIMIT ?
SQL
    @exec_bind = (@bind, $limit);
}

my $sth = $dbh->prepare($sql);
$sth->execute(@exec_bind);

if ($want_tsv) {
    my $fname = "panel_search_${view_mode}.tsv";
    print "Content-Type: text/tab-separated-values; charset=utf-8\r\n";
    print qq{Content-Disposition: attachment; filename="$fname"\r\n\r\n};

    if ($view_mode eq 'variant') {
        print join("\t", qw(case_id panel report_id patient_id disease tissue_of_origin pathology_diagnosis gene variant_type variant_subtype cds_effect protein_effect strand transcript functional_effect status origin classification allele_fraction depth copy_number cnv_ratio clinvar_id clinvar_sig clinvar_match maf_1kg maf_hgvd maf_tommo tpm_tumor tpm_normal_mean tpm_normal_sd non_human_summary)), "\n";
        while (my $row = $sth->fetchrow_hashref) {
            my @cols = (
                _safe($row->{case_id}),
                _safe($row->{panel_name}),
                _safe($row->{report_id}),
                _safe($row->{patient_id}),
                _safe($row->{disease}),
                _safe($row->{tissue_of_origin}),
                _safe($row->{pathology_diagnosis}),
                _safe($row->{gene}),
                _safe($row->{variant_type}),
                _safe($row->{variant_subtype}),
                _safe($row->{cds_effect}),
                _safe($row->{protein_effect}),
                _safe($row->{strand}),
                _safe($row->{transcript}),
                _safe($row->{functional_effect}),
                _safe($row->{status}),
                _safe($row->{origin}),
                _safe($row->{classification}),
                defined $row->{allele_fraction} ? sprintf('%.4f', $row->{allele_fraction}) : '',
                defined $row->{depth} ? $row->{depth} : '',
                defined $row->{copy_number} ? sprintf('%.2f', $row->{copy_number}) : '',
                defined $row->{cnv_ratio} ? sprintf('%.2f', $row->{cnv_ratio}) : '',
                _safe($row->{clinvar_id}),
                _safe($row->{clinvar_sig}),
                _safe($row->{clinvar_match}),
                defined $row->{maf_1kg} ? sprintf('%.4g', $row->{maf_1kg}) : '',
                defined $row->{maf_hgvd} ? sprintf('%.4g', $row->{maf_hgvd}) : '',
                defined $row->{maf_tommo} ? sprintf('%.4g', $row->{maf_tommo}) : '',
                defined $row->{tpm} ? sprintf('%.2f', $row->{tpm}) : '',
                defined $row->{tpm_normal_mean} ? sprintf('%.2f', $row->{tpm_normal_mean}) : '',
                defined $row->{tpm_normal_sd} ? sprintf('%.2f', $row->{tpm_normal_sd}) : '',
                _safe($row->{non_human_summary}),
            );
            print join("\t", @cols), "\n";
        }
    }
    elsif ($view_mode eq 'case') {
        print join("\t", qw(case_id panel report_id patient_id date disease tissue_of_origin pathology_diagnosis matched_variant_count matched_gene_count matched_genes matched_variant_types non_human_summary)), "\n";
        while (my $row = $sth->fetchrow_hashref) {
            my @cols = (
                _safe($row->{case_id}), _safe($row->{panel_name}), _safe($row->{report_id}), _safe($row->{patient_id}), _safe($row->{date}),
                _safe($row->{disease}), _safe($row->{tissue_of_origin}), _safe($row->{pathology_diagnosis}),
                _safe($row->{matched_variant_count}), _safe($row->{matched_gene_count}), _safe($row->{matched_genes}), _safe($row->{matched_variant_types}),
                _safe($row->{non_human_summary}),
            );
            print join("\t", @cols), "\n";
        }
    }
    else {
        print join("\t", qw(patient_group patient_id case_count latest_date report_ids panel_names diseases tissues pathologies matched_variant_count matched_gene_count matched_genes case_links)), "\n";
        while (my $row = $sth->fetchrow_hashref) {
            my @cols = (
                _safe($row->{patient_group}), _safe($row->{patient_id_display}), _safe($row->{case_count}), _safe($row->{latest_date}),
                _safe($row->{report_ids}), _safe($row->{panel_names}), _safe($row->{diseases}), _safe($row->{tissues}), _safe($row->{pathologies}),
                _safe($row->{matched_variant_count}), _safe($row->{matched_gene_count}), _safe($row->{matched_genes}), _safe($row->{case_links}),
            );
            print join("\t", @cols), "\n";
        }
    }

    $sth->finish;
    $dbh->disconnect;
    exit;
}

print $q->header(-type => 'text/html', -charset => 'utf-8');

my $row_count = 0;

print <<'HTML';
<!DOCTYPE html>
<html lang="ja">
<head>
  <meta charset="UTF-8">
  <title>OncoUnify search results</title>
  <style>
    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 1.2rem;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: radial-gradient(circle at top, #eff6ff 0, #e0f2fe 35%, #f3f4f6 100%);
      color: #111827;
    }
    a { color: #2563eb; text-decoration: none; }
    a:hover { text-decoration: underline; }

    .shell {
      width: 100%;
      max-width: none;
      margin: 0;
    }

    .card {
      width: 100%;
      background: #ffffffee;
      border-radius: 1rem;
      padding: 1.2rem 1.4rem 1rem;
      box-shadow:
        0 18px 45px rgba(15, 23, 42, 0.16),
        0 0 0 1px rgba(255, 255, 255, 0.85);
      backdrop-filter: blur(4px);
    }

    .card-header {
      display: flex;
      justify-content: space-between;
      align-items: baseline;
      gap: 1rem;
      margin-bottom: 0.75rem;
    }

    .title-block {
      display: flex;
      flex-direction: column;
      gap: 0.25rem;
    }

    h1 {
      margin: 0;
      font-size: 1.35rem;
      font-weight: 650;
      color: #111827;
      letter-spacing: 0.03em;
    }

    .subtitle {
      font-size: 0.8rem;
      color: #6b7280;
    }

    .meta-block {
      text-align: right;
      font-size: 0.78rem;
      color: #6b7280;
      line-height: 1.7;
    }

    .meta-link {
      font-weight: 500;
    }

    .meta-chip {
      display: inline-flex;
      align-items: center;
      gap: 0.35rem;
      padding: 0.15rem 0.6rem;
      border-radius: 999px;
      background: #eff6ff;
      color: #1d4ed8;
      border: 1px solid #bfdbfe;
      margin-top: 0.25rem;
    }
    .meta-dot {
      width: 7px;
      height: 7px;
      border-radius: 999px;
      background: #22c55e;
      box-shadow: 0 0 0 3px rgba(34, 197, 94, 0.25);
    }
    .logout-link {
      color: #b91c1c;
      font-weight: 600;
      margin-left: 0.5rem;
    }

    .filter-summary {
      margin: 0.4rem 0 0.5rem;
      padding: 0.45rem 0.75rem;
      border-radius: 999px;
      background: #f9fafb;
      font-size: 0.78rem;
      color: #4b5563;
      display: inline-flex;
      flex-wrap: wrap;
      gap: 0.4rem;
      align-items: center;
      width: 100%;
    }
    .filter-pill {
      display: inline-flex;
      align-items: center;
      gap: 0.25rem;
      padding: 0.1rem 0.5rem;
      border-radius: 999px;
      background: #e5f3ff;
      color: #1d4ed8;
    }
    .filter-pill span.key { font-weight: 600; }

    .download-form {
      margin-left: auto;
      display: inline-flex;
      align-items: center;
      gap: 0.5rem;
    }
    .download-button {
      border: none;
      border-radius: 999px;
      padding: 0.25rem 0.8rem;
      font-size: 0.75rem;
      font-weight: 600;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      cursor: pointer;
      background: linear-gradient(135deg, #0ea5e9, #2563eb);
      color: #ffffff;
      box-shadow: 0 4px 10px rgba(37, 99, 235, 0.35);
    }

    .table-wrapper {
      margin-top: 0.4rem;
      border-radius: 0.7rem;
      border: 1px solid #e5e7eb;
      box-shadow: inset 0 1px 0 rgba(255,255,255,0.6);
      background: #f9fafb;
      max-height: 72vh;
      overflow-x: auto;
      overflow-y: auto;
    }

    table {
      width: max-content;
      min-width: 100%;
      border-collapse: collapse;
      font-size: 0.8rem;
    }

    thead { background: #f9fafb; }

    th, td {
      padding: 0.35rem 0.45rem;
      border-bottom: 1px solid #e5e7eb;
      border-right: 1px solid #e5e7eb;
      vertical-align: top;
    }
    th:last-child, td:last-child { border-right: none; }

    th {
      position: sticky;
      top: 0;
      z-index: 2;
      font-weight: 600;
      color: #374151;
      text-align: left;
      white-space: nowrap;
      background: linear-gradient(120deg, #eff6ff, #e0f2fe);
      background-clip: padding-box;
    }

    tbody tr:nth-child(even) td { background: #f3f4f6; }
    tbody tr:hover td { background: #e0f2fe; }

    .num { text-align: right; white-space: nowrap; }
    .mono {
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
    }
    .gene-cell {
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
      font-size: 0.78rem;
      font-weight: 600;
      color: #111827;
    }
    .panel-cell {
      font-size: 0.76rem;
      color: #4b5563;
    }
    .wrap-cell {
      max-width: 420px;
      min-width: 220px;
      white-space: normal;
      word-break: break-word;
      overflow: visible;
      text-overflow: clip;
    }
    .narrow-wrap {
      max-width: 360px;
      min-width: 180px;
      white-space: normal;
      word-break: break-word;
    }

    .tag {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      min-width: 76px;
      padding: 0.1rem 0.4rem;
      border-radius: 999px;
      font-size: 0.7rem;
      font-weight: 600;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      border: 1px solid transparent;
    }
    .tag-snv { background: #eef2ff; color: #4f46e5; border-color: #c7d2fe; }
    .tag-cnv { background: #fef3c7; color: #92400e; border-color: #fde68a; }
    .tag-rearr { background: #fee2e2; color: #b91c1c; border-color: #fecaca; }
    .tag-exp { background: #ecfdf5; color: #047857; border-color: #a7f3d0; }
    .tag-bm { background: #e0f2fe; color: #0369a1; border-color: #bae6fd; }
    .tag-default { background: #f3f4f6; color: #4b5563; border-color: #e5e7eb; }

    .clin-sig { font-size: 0.72rem; color: #6b21a8; }
    .clin-match { font-size: 0.7rem; color: #4b5563; }

    .footer {
      margin-top: 0.6rem;
      font-size: 0.78rem;
      color: #6b7280;
      display: flex;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 0.4rem;
    }
    .footer strong { color: #111827; }
    .empty-note {
      padding: 1rem;
      color: #6b7280;
    }

    @media (max-width: 768px) {
      body { padding: 0.75rem; }
      .card { padding: 1rem 1rem 0.9rem; }
      .card-header { flex-direction: column; align-items: flex-start; }
      .meta-block { text-align: left; }
      .wrap-cell { min-width: 180px; }
    }
  </style>
</head>
<body>
  <div class="shell">
    <div class="card">
      <div class="card-header">
        <div class="title-block">
          <h1>OncoUnify search results</h1>
          <div class="subtitle">Cross-vendor results from FoundationOne, GenMineTOP, and Guardant.</div>
        </div>
        <div class="meta-block">
          <div class="meta-link">
            <a href="/panel/search.html">Modify search</a>
            <a class="logout-link" href="/cgi-bin/logout.cgi">Logout</a>
          </div>
          <div class="meta-chip">
            <span class="meta-dot"></span>
            Result viewer
          </div>
        </div>
      </div>
HTML

print qq{      <div class="filter-summary">\n};
print qq{        <span class="filter-pill"><span class="key">View</span><span>} . _h(_view_label($view_mode)) . qq{</span></span>\n};
if ($gene ne '') {
    print qq{        <span class="filter-pill"><span class="key">Gene</span><span>} . _h($gene) . qq{</span></span>\n};
}
if ($protein_effect ne '') {
    print qq{        <span class="filter-pill"><span class="key">Protein</span><span>} . _h($protein_effect) . qq{</span></span>\n};
}
if ($patient_id ne '') {
    print qq{        <span class="filter-pill"><span class="key">Patient ID</span><span>} . _h($patient_id) . qq{</span></span>\n};
}
if ($variant_type ne '') {
    print qq{        <span class="filter-pill"><span class="key">Type</span><span>} . _h($variant_type) . qq{</span></span>\n};
}
if ($disease ne '') {
    print qq{        <span class="filter-pill"><span class="key">Disease/Tissue</span><span>} . _h($disease) . qq{</span></span>\n};
}
if ($panel_name ne '') {
    print qq{        <span class="filter-pill"><span class="key">Panel</span><span>} . _h($panel_name) . qq{</span></span>\n};
}
print qq{        <span class="filter-pill"><span class="key">Limit</span><span>} . _h($limit) . qq{ rows</span></span>\n};

print qq{
        <form class="download-form" method="get" action="/cgi-bin/panel_search.cgi">
          <input type="hidden" name="gene" value="} . _h($gene) . qq{" />
          <input type="hidden" name="protein_effect" value="} . _h($protein_effect) . qq{" />
          <input type="hidden" name="patient_id" value="} . _h($patient_id) . qq{" />
          <input type="hidden" name="variant_type" value="} . _h($variant_type) . qq{" />
          <input type="hidden" name="disease" value="} . _h($disease) . qq{" />
          <input type="hidden" name="panel_name" value="} . _h($panel_name) . qq{" />
          <input type="hidden" name="view_mode" value="} . _h($view_mode) . qq{" />
          <input type="hidden" name="limit" value="} . _h($limit) . qq{" />
          <input type="hidden" name="download" value="tsv" />
          <button type="submit" class="download-button">Download TSV</button>
        </form>
};
print qq{      </div>\n};

print qq{      <div class="table-wrapper">\n};

if ($view_mode eq 'variant') {
    print <<'HTML';
        <table>
          <thead>
            <tr>
              <th>case_id</th>
              <th>panel</th>
              <th>report_id</th>
              <th>patient_id</th>
              <th>disease / tissue</th>
              <th>pathology</th>
              <th>gene</th>
              <th>variant_type</th>
              <th>subtype</th>
              <th>cDNA</th>
              <th>protein</th>
              <th>strand</th>
              <th>transcript</th>
              <th>func_effect</th>
              <th>status</th>
              <th>origin</th>
              <th>class</th>
              <th>AF</th>
              <th>depth</th>
              <th>copy#</th>
              <th>CNV ratio</th>
              <th>ClinVar</th>
              <th>MAF (1KG / HGVD / ToMMo)</th>
              <th>TPM (tumor)</th>
              <th>TPM normal (mean±SD)</th>
              <th>Non-human</th>
            </tr>
          </thead>
          <tbody>
HTML

    while (my $row = $sth->fetchrow_hashref) {
        $row_count++;
        my $af  = defined $row->{allele_fraction} ? sprintf('%.4f', $row->{allele_fraction}) : '';
        my $cn  = defined $row->{copy_number} ? sprintf('%.2f', $row->{copy_number}) : '';
        my $cnr = defined $row->{cnv_ratio} ? sprintf('%.2f', $row->{cnv_ratio}) : '';

        my $disease_text = join(' / ', grep { length $_ } map { _safe($_) } ($row->{disease}, $row->{tissue_of_origin}));
        my $badge_html = _type_badge($row->{variant_type});

        my $clinvar_cell = '';
        if (_safe($row->{clinvar_id}) ne '') {
            my $url = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' . _safe($row->{clinvar_id});
            $clinvar_cell = qq{<a href="} . _h($url) . qq{" target="_blank" rel="noopener">} . _h($row->{clinvar_id}) . qq{</a>};
            $clinvar_cell .= qq{<br><span class="clin-sig">} . _h($row->{clinvar_sig}) . qq{</span>} if _safe($row->{clinvar_sig}) ne '';
            $clinvar_cell .= qq{<br><span class="clin-match">} . _h($row->{clinvar_match}) . qq{</span>} if _safe($row->{clinvar_match}) ne '';
        }

        my @maf_parts;
        push @maf_parts, '1KG: '   . sprintf('%.4g', $row->{maf_1kg})   if defined $row->{maf_1kg};
        push @maf_parts, 'HGVD: '  . sprintf('%.4g', $row->{maf_hgvd})  if defined $row->{maf_hgvd};
        push @maf_parts, 'ToMMo: ' . sprintf('%.4g', $row->{maf_tommo}) if defined $row->{maf_tommo};
        my $maf_cell = join('<br>', map { _h($_) } @maf_parts);

        my $tpm_tumor = defined $row->{tpm} ? sprintf('%.2f', $row->{tpm}) : '';
        my $tpm_norm_cell = '';
        if (defined $row->{tpm_normal_mean}) {
            $tpm_norm_cell = sprintf('%.2f', $row->{tpm_normal_mean});
            $tpm_norm_cell .= ' ± ' . sprintf('%.2f', $row->{tpm_normal_sd}) if defined $row->{tpm_normal_sd};
        }

        my $cid = _safe($row->{case_id});
        my $rid = _safe($row->{report_id});
        my $detail_url = "/cgi-bin/case_detail.cgi?case_id=$cid";

        print "<tr>";
        print qq{<td class="num">} . _h($cid) . qq{</td>};
        print qq{<td class="panel-cell">} . _h($row->{panel_name}) . qq{</td>};
        print $rid ne '' ? qq{<td><a href="} . _h($detail_url) . qq{">} . _h($rid) . qq{</a></td>} : qq{<td></td>};
        print qq{<td>} . _h($row->{patient_id}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($disease_text) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{pathology_diagnosis}) . qq{</td>};
        print qq{<td class="gene-cell">} . _h($row->{gene}) . qq{</td>};
        print qq{<td>$badge_html</td>};
        print qq{<td>} . _h($row->{variant_subtype}) . qq{</td>};
        print qq{<td class="mono">} . _h($row->{cds_effect}) . qq{</td>};
        print qq{<td class="mono">} . _h($row->{protein_effect}) . qq{</td>};
        print qq{<td>} . _h($row->{strand}) . qq{</td>};
        print qq{<td class="mono">} . _h($row->{transcript}) . qq{</td>};
        print qq{<td>} . _h($row->{functional_effect}) . qq{</td>};
        print qq{<td>} . _h($row->{status}) . qq{</td>};
        print qq{<td>} . _h($row->{origin}) . qq{</td>};
        print qq{<td>} . _h($row->{classification}) . qq{</td>};
        print qq{<td class="num">} . _h($af) . qq{</td>};
        print qq{<td class="num">} . _h($row->{depth}) . qq{</td>};
        print qq{<td class="num">} . _h($cn) . qq{</td>};
        print qq{<td class="num">} . _h($cnr) . qq{</td>};
        print qq{<td>} . $clinvar_cell . qq{</td>};
        print qq{<td class="narrow-wrap">} . $maf_cell . qq{</td>};
        print qq{<td class="num">} . _h($tpm_tumor) . qq{</td>};
        print qq{<td class="num">} . _h($tpm_norm_cell) . qq{</td>};
        print qq{<td class="narrow-wrap">} . _h($row->{non_human_summary}) . qq{</td>};
        print "</tr>\n";
    }
    print "</tbody></table>\n";
}
elsif ($view_mode eq 'case') {
    print <<'HTML';
        <table>
          <thead>
            <tr>
              <th>case_id</th>
              <th>panel</th>
              <th>report_id</th>
              <th>patient_id</th>
              <th>date</th>
              <th>disease / tissue</th>
              <th>pathology</th>
              <th>matched variants</th>
              <th>matched genes</th>
              <th>gene list</th>
              <th>variant types</th>
              <th>Non-human</th>
            </tr>
          </thead>
          <tbody>
HTML
    while (my $row = $sth->fetchrow_hashref) {
        $row_count++;
        my $cid = _safe($row->{case_id});
        my $rid = _safe($row->{report_id});
        my $detail_url = "/cgi-bin/case_detail.cgi?case_id=$cid";
        my $disease_text = join(' / ', grep { length $_ } map { _safe($_) } ($row->{disease}, $row->{tissue_of_origin}));

        print "<tr>";
        print qq{<td class="num">} . _h($cid) . qq{</td>};
        print qq{<td class="panel-cell">} . _h($row->{panel_name}) . qq{</td>};
        print $rid ne '' ? qq{<td><a href="} . _h($detail_url) . qq{">} . _h($rid) . qq{</a></td>} : qq{<td></td>};
        print qq{<td>} . _h($row->{patient_id}) . qq{</td>};
        print qq{<td class="mono">} . _h($row->{date}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($disease_text) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{pathology_diagnosis}) . qq{</td>};
        print qq{<td class="num">} . _h($row->{matched_variant_count}) . qq{</td>};
        print qq{<td class="num">} . _h($row->{matched_gene_count}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{matched_genes}) . qq{</td>};
        print qq{<td>} . _h($row->{matched_variant_types}) . qq{</td>};
        print qq{<td class="narrow-wrap">} . _h($row->{non_human_summary}) . qq{</td>};
        print "</tr>\n";
    }
    print "</tbody></table>\n";
}
else {
    print <<'HTML';
        <table>
          <thead>
            <tr>
              <th>patient</th>
              <th>case count</th>
              <th>latest date</th>
              <th>reports</th>
              <th>panels</th>
              <th>disease</th>
              <th>tissue</th>
              <th>pathology</th>
              <th>matched variants</th>
              <th>matched genes</th>
              <th>gene list</th>
            </tr>
          </thead>
          <tbody>
HTML
    while (my $row = $sth->fetchrow_hashref) {
        $row_count++;
        my $patient_disp = _safe($row->{patient_id_display});
        $patient_disp = _safe($row->{patient_group}) if $patient_disp eq '';

        print "<tr>";
        print qq{<td class="mono">} . _h($patient_disp) . qq{</td>};
        print qq{<td class="num">} . _h($row->{case_count}) . qq{</td>};
        print qq{<td class="mono">} . _h($row->{latest_date}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _render_case_links($row->{case_links}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{panel_names}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{diseases}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{tissues}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{pathologies}) . qq{</td>};
        print qq{<td class="num">} . _h($row->{matched_variant_count}) . qq{</td>};
        print qq{<td class="num">} . _h($row->{matched_gene_count}) . qq{</td>};
        print qq{<td class="wrap-cell">} . _h($row->{matched_genes}) . qq{</td>};
        print "</tr>\n";
    }
    print "</tbody></table>\n";
}

if ($row_count == 0) {
    print qq{<div class="empty-note">No matching results.</div>};
}

print qq{      </div>\n};
print qq{      <div class="footer"><div><strong>$row_count</strong> rows shown (limit: } . _h($limit) . qq{).</div><div>} . _h(_view_label($view_mode)) . qq{ view</div></div>\n};
print qq{    </div>\n  </div>\n</body>\n</html>\n};

$sth->finish;
$dbh->disconnect;
