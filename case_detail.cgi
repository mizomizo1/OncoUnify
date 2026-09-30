#!/usr/bin/perl
# case_detail.cgi — per-report drill-down (read-only).
# Configuration: ONCOUNIFY_DB, ONCOUNIFY_CGI_URL, ONCOUNIFY_HTML_URL (see panel_search.cgi).
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

my $q = CGI->new;
my $case_id = $q->param('case_id') // '';

sub _safe { defined $_[0] ? $_[0] : '' }
sub _h {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/&/&amp;/g; $s =~ s/</&lt;/g; $s =~ s/>/&gt;/g; $s =~ s/"/&quot;/g; $s =~ s/'/&#39;/g;
    return $s;
}
sub _num { my ($v, $f) = @_; defined $v ? sprintf($f, $v) : '' }
sub _fail {
    my ($status, $msg) = @_;
    print $q->header(-type => 'text/html', -charset => 'utf-8', -status => $status);
    print '<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><title>OncoUnify</title></head><body><p>'
        . _h($msg) . '</p><p><a href="' . _h("$HTML_URL/search.html") . '">Back to search</a></p></body></html>';
    exit;
}

_fail('400 Bad Request', 'Invalid case_id.') unless $case_id =~ /^\d+$/;
_fail('500 Internal Server Error', 'OncoUnify database not found or not readable.') unless -r $DB_FILE;

my $dbh = eval {
    DBI->connect("dbi:SQLite:dbname=$DB_FILE", '', '', {
        RaiseError => 1, PrintError => 0, sqlite_unicode => 1, AutoCommit => 1, ReadOnly => 1,
        sqlite_open_flags => DBD::SQLite::OPEN_READONLY(),
    });
} or _fail('500 Internal Server Error', 'Cannot open the OncoUnify database.');

my $case = $dbh->selectrow_hashref('SELECT * FROM cases WHERE case_id = ?', undef, $case_id);
_fail('404 Not Found', "Case not found (case_id = $case_id).") unless $case;

my $biomarkers = $dbh->selectall_arrayref(
    'SELECT name, value, unit, call, call_raw, assay, source_field FROM biomarkers WHERE case_id = ? ORDER BY name',
    { Slice => {} }, $case_id);
my $non_humans = $dbh->selectall_arrayref(
    'SELECT organism, reads_per_million, status, sample FROM non_human_contents WHERE case_id = ?
     ORDER BY reads_per_million DESC, organism', { Slice => {} }, $case_id);
my $variants = $dbh->selectall_arrayref(<<'SQL', { Slice => {} }, $case_id);
SELECT v.*, t.display_group
FROM variants v LEFT JOIN so_terms t ON t.term = v.functional_effect
WHERE v.case_id = ?
ORDER BY CASE v.variant_type WHEN 'short_variant' THEN 1 WHEN 'cnv' THEN 2 WHEN 'rearrangement' THEN 3 ELSE 4 END,
         v.gene, v.pos
SQL
my ($n_tested) = $dbh->selectrow_array(
    'SELECT COUNT(*) FROM v_case_genes_tested WHERE case_id = ? AND short_variants = 1', undef, $case_id);
$dbh->disconnect;

print $q->header(-type => 'text/html', -charset => 'utf-8');
print <<'HTML';
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>OncoUnify case detail</title>
  <style>
    * { box-sizing: border-box; }
    body { margin: 0; padding: 1.5rem; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           background: radial-gradient(circle at top, #f5f3ff 0, #eff6ff 40%, #f3f4f6 100%); color: #111827; }
    a { color: #2563eb; text-decoration: none; } a:hover { text-decoration: underline; }
    .card { max-width: 1400px; margin: 0 auto; background: #ffffffee; border-radius: 1rem; padding: 1.4rem 1.6rem;
            box-shadow: 0 18px 45px rgba(15,23,42,.16), 0 0 0 1px rgba(255,255,255,.85); }
    h1 { margin: 0 0 .4rem; font-size: 1.3rem; font-weight: 650; }
    h2 { margin: 1.2rem 0 .4rem; font-size: 1rem; font-weight: 600; color: #1f2937; }
    .subtitle { font-size: .8rem; color: #6b7280; margin-bottom: .6rem; }
    table { width: 100%; border-collapse: collapse; font-size: .8rem; margin-top: .2rem; }
    th, td { border: 1px solid #e5e7eb; padding: .3rem .4rem; vertical-align: top; overflow-wrap: anywhere; }
    th { background: #f9fafb; font-weight: 600; color: #374151; white-space: nowrap; text-align: left; }
    td.label { background: #f9fafb; font-weight: 600; color: #374151; }
    .num { text-align: right; white-space: nowrap; }
    .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
    .small { font-size: .72rem; color: #6b7280; }
    .tag-panel { display: inline-flex; padding: .1rem .5rem; border-radius: 999px; background: #eff6ff; color: #1d4ed8;
                 border: 1px solid #bfdbfe; font-size: .75rem; }
    .scroll-box { max-height: 60vh; overflow: auto; border-radius: .5rem; border: 1px solid #e5e7eb; background: #f9fafb; }
    .scroll-box table { width: max-content; min-width: 100%; }
    .cq { display: inline-block; width: .6rem; height: .6rem; border-radius: 999px; margin-right: .3rem; }
    .cq-missense { background: #4f46e5; } .cq-nonsense { background: #ef4444; } .cq-frameshift { background: #ec4899; }
    .cq-splice { background: #0ea5e9; } .cq-inframe { background: #f59e0b; } .cq-other { background: #9ca3af; }
    .cq-silent { background: #d1d5db; } .cq-noncoding { background: #e5e7eb; }
    .extra { max-width: 420px; white-space: normal; word-break: break-all; font-size: .7rem; color: #4b5563; }
  </style>
</head>
<body>
  <div class="card">
HTML

my $c = $case;
print '<div class="small"><a href="', _h("$HTML_URL/search.html"), '">&laquo; Back to search</a></div>';
print '<h1>Case detail</h1><div class="subtitle">case_id ', _h($c->{case_id}), ' &middot; report ', _h($c->{report_id}),
      ' &middot; <span class="tag-panel">', _h($c->{panel_name}), '</span></div>';

my $curated = _safe($c->{curated_fields});
my %cur = map { $_ => 1 } split /,/, $curated;
sub _cv { my ($k) = @_; _h($c->{$k}) . ($cur{$k} ? ' <span class="small">(curated)</span>' : '') }

my @info = (
    ['Report ID', _h($c->{report_id}), 'Panel / version', _h($c->{panel_name}) . ' / ' . _h($c->{panel_version})],
    ['Vendor', _h($c->{vendor}), 'Patient ID', _cv('patient_id')],
    ['Report date', _cv('date'), 'Specimen ID', _h($c->{specimen_id})],
    ['Sex / Age', _cv('sex') . ' / ' . _cv('age'), 'Test type', _h($c->{test_type})],
    ['Genome build', _h($c->{genome_build}) . ' <span class="small">(' . _h($c->{genome_build_source}) . ')</span>',
     'Genes tested (SNV/indel)', ($n_tested ? _h($n_tested) : '<span class="small">no gene list loaded for this panel version</span>')],
    ['Disease', _cv('disease'), 'OncoTree', _cv('oncotree_code')],
    ['Vendor disease classification', _h($c->{disease_ontology}), 'Tissue of origin', _cv('tissue_of_origin')],
    ['Pathology diagnosis', _cv('pathology_diagnosis'), 'Tumor nuclei % / purity %',
     _h($c->{percent_tumor_nuclei}) . ' / ' . _h($c->{purity})],
    ['Non-human content (scalar)', _h($c->{non_human_content}), 'Vendor extras', '<span class="extra">' . _h($c->{other_info}) . '</span>'],
    ['Source file', '<span class="mono">' . _h($c->{source_file}) . '</span><br><span class="small mono">sha256 '
     . _h($c->{source_sha256}) . '</span>', 'Loader / format',
     _h($c->{loader}) . '<br><span class="small">' . _h($c->{format_version}) . ' &middot; loaded ' . _h($c->{loaded_at}) . '</span>'],
);
print '<h2>Report</h2><table style="table-layout:fixed"><colgroup><col style="width:16%"><col style="width:34%">'
    . '<col style="width:16%"><col style="width:34%"></colgroup>';
for my $r (@info) {
    print '<tr><td class="label">', $r->[0], '</td><td>', $r->[1], '</td><td class="label">', $r->[2], '</td><td>', $r->[3], "</td></tr>\n";
}
print '</table>';
print '<div class="small">Fields marked (curated) were supplied through the case-metadata sidecar.</div>' if $curated ne '';

print '<h2>Biomarkers</h2>';
if (@$biomarkers) {
    print '<table><thead><tr><th>Biomarker</th><th>Call</th><th>Value</th><th>Unit</th><th>Vendor call</th><th>Assay</th><th>Source field</th></tr></thead><tbody>';
    for my $b (@$biomarkers) {
        print '<tr><td>', _h($b->{name}), '</td><td>', _h($b->{call}), '</td><td class="num">', _h($b->{value}),
              '</td><td>', _h($b->{unit}), '</td><td>', _h($b->{call_raw}), '</td><td>', _h($b->{assay}),
              '</td><td class="small">', _h($b->{source_field}), "</td></tr>\n";
    }
    print '</tbody></table><div class="small">Values are comparable only within the same assay.</div>';
} else {
    print '<div class="small">No biomarker results in this report.</div>';
}

print '<h2>Non-human contents</h2>';
if (@$non_humans) {
    print '<table><thead><tr><th>Organism</th><th>Reads per million</th><th>Status</th><th>Sample</th></tr></thead><tbody>';
    for my $n (@$non_humans) {
        print '<tr><td class="mono">', _h($n->{organism}), '</td><td class="num">', _h(_num($n->{reads_per_million}, '%.0f')),
              '</td><td>', _h($n->{status}), '</td><td class="mono">', _h($n->{sample}), "</td></tr>\n";
    }
    print '</tbody></table>';
} else {
    print '<div class="small">No non-human-content records for this report.</div>';
}

print '<h2>Variants (', scalar(@$variants), ' rows)</h2><div class="scroll-box"><table><thead><tr>';
print map { "<th>$_</th>" } ('Gene', 'Partner', 'Type', 'Subtype', 'Location (' . _h($c->{genome_build}) . ')', 'Ref/Alt',
    'Transcript', 'cDNA', 'Protein (canonical)', 'Protein (vendor)', 'Consequence (SO)', 'Status', 'Origin', 'Class',
    'AF', 'Depth', 'Copy#', 'CNV type', 'In-frame', 'Reads', 'TPM', 'Description', 'Vendor extras');
print "</tr></thead><tbody>\n";
for my $v (@$variants) {
    my $loc = _safe($v->{chrom}) ne '' ? $v->{chrom} . (defined $v->{pos} ? ':' . $v->{pos} : '') : '';
    if (_safe($v->{chrom2}) ne '' || defined $v->{pos2}) {
        $loc .= ($v->{variant_type} eq 'cnv' ? '-' : ' / ') . (_safe($v->{chrom2}) ne '' ? $v->{chrom2} . ':' : '') . _safe($v->{pos2});
    }
    my $ra = (_safe($v->{ref}) ne '' || _safe($v->{alt}) ne '') ? _safe($v->{ref}) . '>' . _safe($v->{alt}) : '';
    my $cq = '';
    if (_safe($v->{functional_effect}) ne '') {
        my $g = _safe($v->{display_group}) || 'other';
        $cq = qq{<span class="cq cq-$g"></span>} . _h($v->{functional_effect}) . '<br><span class="small">'
            . _h($v->{functional_effect_so}) . ' &middot; ' . _h($v->{functional_effect_source})
            . (_safe($v->{functional_effect_raw}) ne '' ? ' &middot; vendor: ' . _h($v->{functional_effect_raw}) : '') . '</span>';
    }
    print '<tr>',
        '<td class="mono">', _h($v->{gene}), '</td>',
        '<td class="mono">', _h($v->{other_gene}), '</td>',
        '<td>', _h($v->{variant_type}), '</td>',
        '<td>', _h($v->{variant_subtype}), '</td>',
        '<td class="mono">', _h($loc), '</td>',
        '<td class="mono">', _h($ra), '</td>',
        '<td class="mono">', _h($v->{transcript}), '</td>',
        '<td class="mono">', _h($v->{hgvs_c}), '</td>',
        '<td class="mono">', _h($v->{hgvs_p}), '</td>',
        '<td class="mono small">', _h($v->{protein_effect}), '</td>',
        '<td>', $cq, '</td>',
        '<td>', _h($v->{status}), '</td>',
        '<td>', _h($v->{origin}), '</td>',
        '<td>', _h($v->{classification}), '</td>',
        '<td class="num">', _h(_num($v->{allele_fraction}, '%.4f')), '</td>',
        '<td class="num">', _h($v->{depth}), '</td>',
        '<td class="num">', _h(_num($v->{copy_number}, '%.2f')), '</td>',
        '<td>', _h($v->{cnv_type}), '</td>',
        '<td>', _h($v->{in_frame}), '</td>',
        '<td class="num">', _h($v->{read_count} // $v->{supporting_read_pairs}), '</td>',
        '<td class="num">', _h(_num($v->{tpm}, '%.2f')), '</td>',
        '<td class="extra">', _h($v->{effect}), '</td>',
        '<td class="extra">', _h($v->{extra}), '</td>',
        "</tr>\n";
}
print "</tbody></table></div>\n  </div>\n</body>\n</html>\n";
