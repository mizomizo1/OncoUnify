#!/usr/bin/perl
# panel_stats.cgi — registry overview with panel-aware mutation frequencies.
#
# A gene's frequency is n_mutated / n_tested, where n_tested counts only the
# cases whose assay (panel_name + panel_version) interrogates the gene for
# short variants according to the loaded gene lists (panel_genes).  Cases
# whose assay has no gene list are excluded from both numerator and
# denominator and reported separately, so that genes covered by one vendor
# are never ranked against genes covered by all of them.
#
# Parameters: src=disease|oncotree|tissue|pathology, organ=<value>, genes=all|common
# Configuration: ONCOUNIFY_DB (see panel_search.cgi).
use strict;
use warnings;
use utf8;
use CGI qw(:standard);
use DBI;
use DBD::SQLite ();

binmode(STDOUT, ':encoding(UTF-8)');

my $DB_FILE = $ENV{ONCOUNIFY_DB} || '/var/www/data/panels.db';

my $q = CGI->new;
my $src   = lc($q->param('src') // 'oncotree');
my $organ = $q->param('organ') // '';
my $genes = lc($q->param('genes') // 'all') eq 'common' ? 'common' : 'all';
$organ =~ s/^\s+|\s+$//g;

my %SRC = (disease => 'disease', oncotree => 'oncotree_code', tissue => 'tissue_of_origin', pathology => 'pathology_diagnosis');
$src = 'oncotree' unless exists $SRC{$src};
my $col = $SRC{$src};   # from a fixed whitelist, never from user input

sub _safe { defined $_[0] ? $_[0] : '' }
sub _h {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/&/&amp;/g; $s =~ s/</&lt;/g; $s =~ s/>/&gt;/g; $s =~ s/"/&quot;/g; $s =~ s/'/&#39;/g;
    return $s;
}

unless (-r $DB_FILE) {
    print $q->header(-type => 'text/plain', -charset => 'utf-8', -status => '500 Internal Server Error');
    print "OncoUnify database not found or not readable.\n";
    exit;
}
my $dbh = DBI->connect("dbi:SQLite:dbname=$DB_FILE", '', '', {
    RaiseError => 1, PrintError => 0, sqlite_unicode => 1, AutoCommit => 1, ReadOnly => 1,
    sqlite_open_flags => DBD::SQLite::OPEN_READONLY(),
});

# Qualifying alterations: short variants (and RNA exon skipping) whose
# consequence alters the protein or splicing; silent and non-coding excluded.
my $QUALIFY = "v.functional_effect IN (SELECT term FROM so_terms WHERE display_group NOT IN ('silent','noncoding'))"
            . " AND (v.status IS NULL OR v.status != 'not_called')";

my $panels = $dbh->selectall_arrayref(<<'SQL', { Slice => {} });
SELECT c.panel_name, COALESCE(c.panel_version, '') AS panel_version, COUNT(*) AS n_cases,
       (SELECT COUNT(*) FROM panel_versions pv JOIN panel_genes pg ON pg.panel_version_id = pv.panel_version_id
         WHERE pv.panel_name = c.panel_name AND pv.panel_version = c.panel_version AND pg.short_variants = 1) AS n_genes
FROM cases c
GROUP BY c.panel_name, c.panel_version
ORDER BY c.panel_name, c.panel_version
SQL

my $organs = $dbh->selectall_arrayref(<<"SQL", { Slice => {} });
SELECT $col AS organ, COUNT(*) AS n_cases FROM cases
WHERE $col IS NOT NULL AND TRIM($col) != ''
GROUP BY $col ORDER BY n_cases DESC, organ LIMIT 300
SQL

my @scope_bind;
my $scope = 'SELECT c.case_id, c.panel_name, c.panel_version FROM cases c';
if ($organ ne '') { $scope .= " WHERE c.$col = ?"; push @scope_bind, $organ; }

my ($n_scope, $n_with_list) = $dbh->selectrow_array(<<"SQL", undef, @scope_bind);
WITH scope AS ($scope)
SELECT COUNT(*),
       SUM(CASE WHEN EXISTS (SELECT 1 FROM panel_versions pv
                             WHERE pv.panel_name = scope.panel_name AND pv.panel_version = scope.panel_version)
                THEN 1 ELSE 0 END)
FROM scope
SQL
$n_with_list ||= 0;

# number of distinct assays (with gene lists) in scope, for the 'common genes' option
my ($n_assays) = $dbh->selectrow_array(<<"SQL", undef, @scope_bind);
WITH scope AS ($scope)
SELECT COUNT(DISTINCT pv.panel_version_id) FROM scope
JOIN panel_versions pv ON pv.panel_name = scope.panel_name AND pv.panel_version = scope.panel_version
SQL

$n_assays = int($n_assays || 0);
my $common_filter = $genes eq 'common'
    ? "HAVING COUNT(DISTINCT t.panel_version_id) = $n_assays"   # integer computed above, not user input
    : '';

my $freq = $dbh->selectall_arrayref(<<"SQL", { Slice => {} }, @scope_bind);
WITH scope AS ($scope),
tested AS (
    SELECT s.case_id, pg.gene, pv.panel_version_id
    FROM scope s
    JOIN panel_versions pv ON pv.panel_name = s.panel_name AND pv.panel_version = s.panel_version
    JOIN panel_genes pg    ON pg.panel_version_id = pv.panel_version_id AND pg.short_variants = 1
),
gene_tested AS (
    SELECT t.gene, COUNT(DISTINCT t.case_id) AS n_tested
    FROM tested t GROUP BY t.gene $common_filter
),
hits AS (
    SELECT DISTINCT v.case_id, v.gene
    FROM variants v JOIN tested t ON t.case_id = v.case_id AND t.gene = v.gene
    WHERE $QUALIFY
)
SELECT g.gene, g.n_tested, COUNT(h.case_id) AS n_mutated
FROM gene_tested g LEFT JOIN hits h ON h.gene = g.gene
GROUP BY g.gene, g.n_tested
HAVING n_mutated > 0
ORDER BY n_mutated DESC, CAST(n_mutated AS REAL) / g.n_tested DESC, g.gene
LIMIT 20
SQL

# consequence composition of the qualifying variants in the tested cases
my %comp;
if (@$freq) {
    my @g = map { $_->{gene} } @$freq;
    my $ph = join(',', ('?') x @g);
    my $rows = $dbh->selectall_arrayref(<<"SQL", { Slice => {} }, @scope_bind, @g);
WITH scope AS ($scope),
tested AS (
    SELECT s.case_id, pg.gene
    FROM scope s
    JOIN panel_versions pv ON pv.panel_name = s.panel_name AND pv.panel_version = s.panel_version
    JOIN panel_genes pg    ON pg.panel_version_id = pv.panel_version_id AND pg.short_variants = 1
)
SELECT v.gene, st.display_group AS grp, COUNT(*) AS n
FROM variants v
JOIN tested t ON t.case_id = v.case_id AND t.gene = v.gene
JOIN so_terms st ON st.term = v.functional_effect
WHERE $QUALIFY AND v.gene IN ($ph)
GROUP BY v.gene, st.display_group
SQL
    $comp{$_->{gene}}{$_->{grp}} = $_->{n} for @$rows;
}

# variants reported in genes outside the declared assay content (data-quality signal)
my ($n_outside) = $dbh->selectrow_array(<<"SQL", undef, @scope_bind);
WITH scope AS ($scope)
SELECT COUNT(*) FROM variants v JOIN scope s ON s.case_id = v.case_id
JOIN panel_versions pv ON pv.panel_name = s.panel_name AND pv.panel_version = s.panel_version
WHERE $QUALIFY
  AND NOT EXISTS (SELECT 1 FROM panel_genes pg WHERE pg.panel_version_id = pv.panel_version_id
                  AND pg.gene = v.gene AND pg.short_variants = 1)
SQL
$dbh->disconnect;

my @GROUPS = (['missense', '#4f46e5'], ['nonsense', '#ef4444'], ['frameshift', '#ec4899'],
              ['splice', '#0ea5e9'], ['inframe', '#f59e0b'], ['other', '#9ca3af']);

print $q->header(-type => 'text/html', -charset => 'utf-8');
print <<'HTML';
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>OncoUnify registry overview</title>
  <style>
    body { margin: 0; padding: .6rem .8rem; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background: #f9fafb; color: #111827; }
    h3 { margin: .5rem 0 .25rem; font-size: .82rem; font-weight: 600; color: #374151; }
    table { border-collapse: collapse; width: 100%; font-size: .75rem; }
    th, td { border: 1px solid #e5e7eb; padding: .15rem .35rem; text-align: left; }
    th { background: #f3f4f6; font-weight: 600; }
    .num { text-align: right; }
    .note { font-size: .7rem; color: #6b7280; margin-top: .25rem; }
    .warn { font-size: .72rem; color: #92400e; background: #fffbeb; border: 1px solid #fde68a; border-radius: .4rem; padding: .3rem .5rem; margin: .3rem 0; }
    .toolbar { display: flex; flex-wrap: wrap; align-items: center; gap: .5rem; padding: .35rem .5rem; border: 1px solid #e5e7eb;
               border-radius: .7rem; background: #fff; }
    .toolbar label { font-size: .75rem; font-weight: 600; color: #374151; }
    .toolbar select { font-size: .75rem; padding: .2rem .45rem; border-radius: .6rem; border: 1px solid #d1d5db; background: #f9fafb; }
    .bar-row { display: flex; align-items: center; gap: .35rem; margin: .15rem 0; }
    .bar-label { width: 5.8rem; font-size: .72rem; font-family: ui-monospace, Menlo, monospace; font-weight: 600; }
    .bar-track { flex: 1; height: 11px; border-radius: 999px; background: #e5e7eb; overflow: hidden; display: flex; }
    .bar-total { width: 8.5rem; font-size: .7rem; text-align: right; color: #4b5563; white-space: nowrap; }
    .legend { display: flex; flex-wrap: wrap; gap: .25rem .6rem; margin-top: .35rem; font-size: .7rem; color: #4b5563; }
    .dot { display: inline-block; width: 10px; height: 10px; border-radius: 999px; margin-right: .2rem; vertical-align: middle; }
  </style>
</head>
<body>
HTML

print '<h3>Registered reports per assay</h3><table><thead><tr><th>Panel</th><th>Version</th><th class="num">Reports</th><th class="num">Genes in list (SNV/indel)</th></tr></thead><tbody>';
for my $p (@$panels) {
    print '<tr><td>', _h($p->{panel_name}), '</td><td>', _h($p->{panel_version}), '</td><td class="num">', _h($p->{n_cases}),
          '</td><td class="num">', ($p->{n_genes} ? _h($p->{n_genes}) : '<span style="color:#b45309">not loaded</span>'), "</td></tr>\n";
}
print '</tbody></table>';

my %sel = map { $_ => ($_ eq $src ? ' selected' : '') } keys %SRC;
print '<h3>Mutation frequency in the selected scope</h3><div class="toolbar">';
print '<label for="srcSel">Category</label><select id="srcSel">';
print qq{<option value="$_"$sel{$_}>} . _h($SRC{$_}) . '</option>' for qw(oncotree disease tissue pathology);
print '</select><label for="organSel">Value</label><select id="organSel"><option value="">(all reports)</option>';
for my $o (@$organs) {
    my $v = _safe($o->{organ});
    print '<option value="', _h($v), '"', ($v eq $organ ? ' selected' : ''), '>', _h($v), ' (', _h($o->{n_cases}), ')</option>';
}
print '</select><label for="genesSel">Genes</label><select id="genesSel">';
print '<option value="all"', ($genes eq 'all' ? ' selected' : ''), '>all genes tested in scope</option>';
print '<option value="common"', ($genes eq 'common' ? ' selected' : ''), '>only genes common to every assay in scope</option>';
print '</select></div>';

my $excluded = $n_scope - $n_with_list;
print '<div class="note">Scope: ', _h($organ eq '' ? 'all reports' : "$SRC{$src} = $organ"), ' &middot; ', _h($n_scope),
      ' reports, ', _h($n_with_list), ' with a gene list.</div>';
print '<div class="warn">', _h($excluded), ' report(s) in scope belong to an assay without a loaded gene list and are excluded ',
      'from both numerator and denominator (load it with load_panel_genes.py).</div>' if $excluded > 0;
print '<div class="warn">', _h($n_outside), ' qualifying variant(s) fall in genes outside the declared content of their assay ',
      'and are not counted; the gene list may be outdated.</div>' if $n_outside;

if (@$freq) {
    for my $r (@$freq) {
        my $pct = $r->{n_tested} ? 100 * $r->{n_mutated} / $r->{n_tested} : 0;
        my $c = $comp{$r->{gene}} || {};
        my $tot = 0; $tot += $_ for values %$c; $tot ||= 1;
        print '<div class="bar-row"><div class="bar-label">', _h($r->{gene}), '</div><div class="bar-track">';
        for my $g (@GROUPS) {
            my $n = $c->{$g->[0]} or next;
            printf '<span title="%s: %d" style="height:100%%;width:%.2f%%;background:%s"></span>',
                   _h($g->[0]), $n, $pct * $n / $tot, $g->[1];
        }
        printf '</div><div class="bar-total">%d / %d (%.1f%%)</div></div>' . "\n", $r->{n_mutated}, $r->{n_tested}, $pct;
    }
    print '<div class="legend">';
    print '<span><span class="dot" style="background:', $_->[1], '"></span>', _h($_->[0]), '</span>' for @GROUPS;
    print '</div><div class="note">Bar length = n_mutated / n_tested (reports carrying at least one protein-altering or splice ',
          'short variant, among reports whose assay covers the gene). Colours show the consequence mix (Sequence Ontology groups). ',
          'Silent, non-coding, and uncalled rows are excluded.</div>';
} elsif ($n_with_list) {
    print '<div class="note">No qualifying variants in this scope.</div>';
} else {
    print '<div class="note">Frequencies need assay gene lists. Load them with <code>load_panel_genes.py</code>.</div>';
}

print <<'HTML';
<script>
(function () {
  function go(reset) {
    const u = new URL(window.location.href);
    u.searchParams.set('src', document.getElementById('srcSel').value);
    u.searchParams.set('genes', document.getElementById('genesSel').value);
    const organ = reset ? '' : document.getElementById('organSel').value;
    if (organ) u.searchParams.set('organ', organ); else u.searchParams.delete('organ');
    window.location.href = u.toString();
  }
  document.getElementById('srcSel').addEventListener('change', () => go(true));
  document.getElementById('organSel').addEventListener('change', () => go(false));
  document.getElementById('genesSel').addEventListener('change', () => go(false));
})();
</script>
</body>
</html>
HTML
