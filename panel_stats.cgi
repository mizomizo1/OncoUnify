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
my $src   = $q->param('src')   // 'disease';   # disease | tissue | pathology
my $organ = $q->param('organ') // '';          # selected tissue/disease (empty = no filter)

# Sanitise the src parameter
$src = lc($src);
$src = 'disease'   unless $src eq 'disease' || $src eq 'tissue' || $src eq 'pathology';

my %src_to_col = (
    disease   => 'cases.disease',
    tissue    => 'cases.tissue_of_origin',
    pathology => 'cases.pathology_diagnosis',
);

my $col_expr = $src_to_col{$src};

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

sub _safe { defined $_[0] ? $_[0] : '' }

# -----------------------------
# Case count per panel
# -----------------------------
my $sth_panel = $dbh->prepare(<<'SQL');
SELECT panel_name, COUNT(DISTINCT case_id) AS n_cases
FROM cases
GROUP BY panel_name
ORDER BY panel_name
SQL
$sth_panel->execute();
my @panel_rows;
while (my $r = $sth_panel->fetchrow_hashref) {
    push @panel_rows, $r;
}
$sth_panel->finish;

# -----------------------------
# Tissue / disease candidate list (per src)
#  - exclude NULL / empty
#  - sorted by descending count
# -----------------------------
my $sth_org = $dbh->prepare(<<"SQL");
SELECT $col_expr AS organ, COUNT(DISTINCT cases.case_id) AS n_cases
FROM cases
WHERE $col_expr IS NOT NULL AND TRIM($col_expr) != ''
GROUP BY organ
ORDER BY n_cases DESC, organ ASC
LIMIT 200
SQL
$sth_org->execute();
my @org_list;
while (my $r = $sth_org->fetchrow_hashref) {
    push @org_list, $r;  # {organ, n_cases}
}
$sth_org->finish;

# -----------------------------
# Counts of gene x functional_effect (all data, or filtered by tissue)
# -----------------------------
my @bind = ();
#my $where_sql = "WHERE variants.gene IS NOT NULL AND variants.functional_effect IS NOT NULL";
my $where_sql = "WHERE variants.gene IS NOT NULL AND variants.functional_effect IS NOT NULL AND LOWER(variants.functional_effect) != 'synonymous'";

if ($organ ne '') {
    $where_sql .= " AND $col_expr = ?";
    push @bind, $organ;
}

my $sth_gene = $dbh->prepare(<<"SQL");
SELECT variants.gene AS gene, variants.functional_effect AS functional_effect, COUNT(*) AS n
FROM variants
JOIN cases ON cases.case_id = variants.case_id
$where_sql
GROUP BY gene, functional_effect
SQL
$sth_gene->execute(@bind);

my %gene_func_counts;
my %gene_totals;
while (my $r = $sth_gene->fetchrow_hashref) {
    my $gene = _safe($r->{gene});
    my $fe   = _safe($r->{functional_effect});
    my $n    = $r->{n} || 0;
    $gene_func_counts{$gene}{$fe} = $n;
    $gene_totals{$gene} += $n;
}
$sth_gene->finish;

$dbh->disconnect;

# Top20
my @top_genes = sort { ($gene_totals{$b}||0) <=> ($gene_totals{$a}||0) } keys %gene_totals;
splice @top_genes, 20 if @top_genes > 20;

my @known_effects = qw(missense nonsense frameshift splice);

print $q->header(-type => 'text/html', -charset => 'utf-8');

print <<'HTML';
<!DOCTYPE html>
<html lang="ja">
<head>
  <meta charset="UTF-8">
  <title>Panel stats</title>
  <style>
    :root {
      --c-missense:   #4f46e5;
      --c-nonsense:   #ef4444;
      --c-frameshift: #ec4899;
      --c-splice:     #0ea5e9;
      --c-other:      #9ca3af;
    }

    body {
      margin: 0;
      padding: 0.6rem 0.8rem;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f9fafb;
      color: #111827;
    }

    h2 {
      margin: 0.1rem 0 0.3rem;
      font-size: 0.9rem;
      font-weight: 600;
      color: #111827;
    }

    h3 {
      margin: 0.4rem 0 0.2rem;
      font-size: 0.8rem;
      font-weight: 600;
      color: #374151;
    }

    .section { margin-bottom: 0.7rem; }

    table {
      border-collapse: collapse;
      width: 100%;
      font-size: 0.75rem;
    }
    th, td {
      border: 1px solid #e5e7eb;
      padding: 0.15rem 0.3rem;
      text-align: left;
    }
    th {
      background: #f3f4f6;
      font-weight: 600;
    }
    .num { text-align: right; }

    .small-note {
      font-size: 0.7rem;
      color: #6b7280;
      margin-top: 0.25rem;
    }

    .toolbar {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 0.5rem;
      padding: 0.35rem 0.5rem;
      border: 1px solid #e5e7eb;
      border-radius: 0.7rem;
      background: #ffffff;
      box-shadow: 0 4px 12px rgba(15, 23, 42, 0.06);
    }
    .toolbar label {
      font-size: 0.75rem;
      color: #374151;
      font-weight: 600;
    }
    .toolbar select {
      font-size: 0.75rem;
      padding: 0.22rem 0.45rem;
      border-radius: 0.6rem;
      border: 1px solid #d1d5db;
      background: #f9fafb;
      outline: none;
    }
    .toolbar select:focus {
      border-color: #0ea5e9;
      box-shadow: 0 0 0 2px rgba(14,165,233,0.25);
      background: #fff;
    }

    .bar-row {
      display: flex;
      align-items: center;
      gap: 0.35rem;
      margin: 0.15rem 0;
    }
    .bar-label {
      width: 5.8rem;
      font-size: 0.7rem;
      color: #374151;
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .bar-track {
      flex: 1;
      height: 10px;
      border-radius: 999px;
      background: #e5e7eb;
      overflow: hidden;
      display: flex;
    }
    .bar-seg { height: 100%; }
    .bar-seg.missense   { background: var(--c-missense); }
    .bar-seg.nonsense   { background: var(--c-nonsense); }
    .bar-seg.frameshift { background: var(--c-frameshift); }
    .bar-seg.splice     { background: var(--c-splice); }
    .bar-seg.other      { background: var(--c-other); }

    .bar-total {
      width: 2.8rem;
      font-size: 0.7rem;
      text-align: right;
      color: #4b5563;
    }

    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 0.25rem 0.6rem;
      margin-top: 0.35rem;
    }
    .legend-item {
      display: inline-flex;
      align-items: center;
      gap: 0.25rem;
      font-size: 0.7rem;
      color: #4b5563;
    }
    .legend-color {
      width: 10px;
      height: 10px;
      border-radius: 999px;
      border: 1px solid #e5e7eb;
    }
  </style>
</head>
<body>
  <h2>Registry overview</h2>
HTML

# -----------------------------
# Panel counts
# -----------------------------
print qq{  <div class="section">\n};
print qq{    <h3>Registered cases per panel</h3>\n};
if (@panel_rows) {
print qq{    <table>\n      <thead><tr><th>Panel</th><th>Cases</th></tr></thead>\n      <tbody>\n};
    for my $r (@panel_rows) {
        print "<tr><td>" . _safe($r->{panel_name}) . "</td><td class=\"num\">" . _safe($r->{n_cases}) . "</td></tr>\n";
    }
    print qq{      </tbody>\n    </table>\n};
} else {
    print qq{    <div class="small-note">No cases are currently registered.</div>\n};
}
print qq{  </div>\n};

# -----------------------------
# Toolbar (organ selector)
# -----------------------------
my %sel_src = (disease => '', tissue => '', pathology => '');
$sel_src{$src} = 'selected';

my $title_src = ($src eq 'disease') ? 'disease' : ($src eq 'tissue' ? 'tissue_of_origin' : 'pathology_diagnosis');

print qq{  <div class="section">\n};
print qq{    <h3>Select a tissue category to show the Top 20</h3>\n};
print qq{    <div class="toolbar">\n};
print qq{      <label>Category:</label>\n};
print qq{      <select id="srcSel">\n};
print qq{        <option value="disease" $sel_src{disease}>disease</option>\n};
print qq{        <option value="tissue" $sel_src{tissue}>tissue_of_origin</option>\n};
print qq{        <option value="pathology" $sel_src{pathology}>pathology_diagnosis</option>\n};
print qq{      </select>\n};
print qq{      <label>Tissue:</label>\n};
print qq{      <select id="organSel">\n};
print qq{        <option value="">(all)</option>\n};

for my $o (@org_list) {
    my $val = _safe($o->{organ});
    my $ncs = _safe($o->{n_cases});
    my $sel = ($organ ne '' && $val eq $organ) ? 'selected' : '';
    # Display labels as e.g. "colon (12)"
    print qq{        <option value="} . _html_escape($val) . qq{" $sel>} . _html_escape($val) . qq{ ($ncs)</option>\n};
}

print qq{      </select>\n};
print qq{    </div>\n};

# NOTE: HTML escape helper
sub _html_escape {
    my ($s) = @_;
    $s = '' unless defined $s;
    $s =~ s/&/&amp;/g;
    $s =~ s/</&lt;/g;
    $s =~ s/>/&gt;/g;
    $s =~ s/"/&quot;/g;
    $s =~ s/'/&#39;/g;
    return $s;
}

# -----------------------------
# Bars
# -----------------------------
my $scope_label = ($organ eq '') ? "all (switch $title_src above)" : "$title_src = " . $organ;

print qq{    <div class="small-note">Scope: } . _html_escape($scope_label) . qq{</div>\n};
print qq{    <h3>Top 20 mutated genes</h3>\n};

if (@top_genes) {
    my $max_total = 0;
    for my $g (@top_genes) {
        $max_total = $gene_totals{$g} if ($gene_totals{$g}||0) > $max_total;
    }
    $max_total ||= 1;

    for my $g (@top_genes) {
        my $total = $gene_totals{$g} || 0;
        next if $total == 0;

        print qq{    <div class="bar-row">\n};
        print qq{      <div class="bar-label">} . _html_escape($g) . qq{</div>\n};
        print qq{      <div class="bar-track">\n};

        my $func_map = $gene_func_counts{$g} || {};
        my $other_bucket = 0;

        for my $fe (@known_effects) {
            next unless exists $func_map->{$fe};
            my $n = $func_map->{$fe} || 0;
            my $w = $n / $max_total * 100;
            print qq{        <span class="bar-seg $fe" style="width: } . sprintf('%.1f', $w) . qq{%"></span>\n};
        }

        for my $fe (keys %{$func_map}) {
            next if grep { $_ eq $fe } @known_effects;
            $other_bucket += $func_map->{$fe} || 0;
        }
        if ($other_bucket > 0) {
            my $w = $other_bucket / $max_total * 100;
            print qq{        <span class="bar-seg other" style="width: } . sprintf('%.1f', $w) . qq{%"></span>\n};
        }

        print qq{      </div>\n};
        print qq{      <div class="bar-total">} . $total . qq{</div>\n};
        print qq{    </div>\n};
    }

    print qq{    <div class="legend">\n};
    my @legend_defs = (
        ['missense',   'missense'],
        ['nonsense',   'nonsense'],
        ['frameshift', 'frameshift'],
        ['splice',     'splice'],
        ['other',      'other'],
    );
    for my $ld (@legend_defs) {
        my ($class, $label) = @$ld;
        print qq{      <div class="legend-item"><span class="legend-color bar-seg $class"></span><span>$label</span></div>\n};
    }
    print qq{    </div>\n};

    #print qq{    <div class="small-note">Bar length is proportional to the total variant count of the gene within the scope (max gene = 100%). Colour shows functional_effect breakdown.</div>\n};
} else {
    print qq{    <div class="small-note">No variants in this scope.</div>\n};
}

print qq{  </div>\n}; # section

# -----------------------------
# JS: select change -> reload same CGI with params
# -----------------------------
print <<'HTML';
  <script>
    (function(){
      const srcSel = document.getElementById('srcSel');
      const organSel = document.getElementById('organSel');

      function reload() {
        const src = srcSel.value;
        const organ = organSel.value;
        const url = new URL(window.location.href);
        url.searchParams.set('src', src);
        if (organ) url.searchParams.set('organ', organ);
        else url.searchParams.delete('organ');
        // Reload within the same iframe
        window.location.href = url.toString();
      }

      srcSel.addEventListener('change', () => {
        // When src changes, reset organ and reload so a fresh candidate list is fetched
        const url = new URL(window.location.href);
        url.searchParams.set('src', srcSel.value);
        url.searchParams.delete('organ');
        window.location.href = url.toString();
      });

      organSel.addEventListener('change', reload);
    })();
  </script>
</body>
</html>
HTML

