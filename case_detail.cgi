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
my $case_id = $q->param('case_id');

print $q->header(-type => 'text/html', -charset => 'utf-8');

unless (defined $case_id && $case_id =~ /^\d+$/) {
    print "<html><body><p>Invalid case_id.</p></body></html>";
    exit;
}

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

# ---- case information ----
my $case_sth = $dbh->prepare(<<'SQL');
SELECT
    case_id, panel_name, panel_type, vendor, report_id, patient_id,
    sex, age, disease, disease_ontology, tissue_of_origin,
    pathology_diagnosis, specimen_id, test_type,
    date,                         -- added
    percent_tumor_nuclei, purity, msi_status,
    tmb_score, tmb_status, tmb_unit,
    non_human_content, other_info
FROM cases
WHERE case_id = ?
SQL
$case_sth->execute($case_id);
my $case = $case_sth->fetchrow_hashref;
$case_sth->finish;

unless ($case) {
    print "<html><body><p>Case not found (case_id = $case_id).</p></body></html>";
    $dbh->disconnect;
    exit;
}

# ---- Non-human contents ----
my $nh_sth = $dbh->prepare(<<'SQL');
SELECT organism, reads_per_million, status, sample
FROM non_human_contents
WHERE case_id = ?
ORDER BY reads_per_million DESC, organism
SQL
$nh_sth->execute($case_id);
my @nh_rows;
while (my $r = $nh_sth->fetchrow_hashref) {
    push @nh_rows, $r;
}
$nh_sth->finish;

# ---- Variants ----
my $var_sth = $dbh->prepare(<<'SQL');
SELECT
    gene, variant_type, variant_subtype,
    chrom, pos, pos2,
    cds_effect, protein_effect,
    strand, transcript,
    functional_effect, status, origin, classification,
    allele_fraction, depth, copy_number, cnv_ratio,
    other_gene, in_frame, supporting_read_pairs,
    tpm
FROM variants
WHERE case_id = ?
ORDER BY variant_type, gene, pos
SQL
$var_sth->execute($case_id);
my @variants;
while (my $r = $var_sth->fetchrow_hashref) {
    push @variants, $r;
}
$var_sth->finish;

$dbh->disconnect;

# ---- HTML output ----

print <<'HTML';
<!DOCTYPE html>
<html lang="ja">
<head>
  <meta charset="UTF-8">
  <title>Case detail</title>
  <style>
    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 1.5rem;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: radial-gradient(circle at top, #f5f3ff 0, #eff6ff 40%, #f3f4f6 100%);
      color: #111827;
    }
    a { color: #2563eb; text-decoration: none; }
    a:hover { text-decoration: underline; }

    .shell { max-width: 1200px; margin: 0 auto; }

    .card {
      background: #ffffffee;
      border-radius: 1rem;
      padding: 1.4rem 1.6rem 1.4rem;
      box-shadow:
        0 18px 45px rgba(15, 23, 42, 0.16),
        0 0 0 1px rgba(255, 255, 255, 0.85);
      backdrop-filter: blur(4px);
      margin-bottom: 1rem;
    }

    h1 {
      margin: 0 0 0.4rem;
      font-size: 1.3rem;
      font-weight: 650;
      letter-spacing: 0.03em;
    }
    h2 {
      margin: 1.2rem 0 0.4rem;
      font-size: 1rem;
      font-weight: 600;
      color: #1f2937;
    }

    .subtitle {
      font-size: 0.8rem;
      color: #6b7280;
      margin-bottom: 0.6rem;
    }

    .back-link {
      font-size: 0.8rem;
      margin-bottom: 0.8rem;
    }

    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 0.8rem;
      margin-top: 0.2rem;
    }
    th, td {
      border: 1px solid #e5e7eb;
      padding: 0.3rem 0.4rem;
      vertical-align: top;
    }
    th {
      background: #f9fafb;
      font-weight: 600;
      color: #374151;
      white-space: nowrap;
    }
    td.label {
      width: 20%;
      background: #f9fafb;
      font-weight: 600;
      color: #374151;
      white-space: nowrap;
    }
    .num { text-align: right; white-space: nowrap; }
    .mono {
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
    }
    .tag-panel {
      display: inline-flex;
      align-items: center;
      padding: 0.1rem 0.5rem;
      border-radius: 999px;
      background: #eff6ff;
      color: #1d4ed8;
      border: 1px solid #bfdbfe;
      font-size: 0.75rem;
    }
    .small-note {
      font-size: 0.75rem;
      color: #6b7280;
      margin-top: 0.3rem;
    }
    .scroll-box {
      max-height: 55vh;
      overflow: auto;
      border-radius: 0.5rem;
      border: 1px solid #e5e7eb;
      background: #f9fafb;
      padding: 0.3rem;
    }

    @media (max-width: 768px) {
      body { padding: 0.75rem; }
      .card { padding: 1rem 1rem 1rem; }
      td.label { width: 30%; }
    }
  </style>
</head>
<body>
  <div class="shell">
HTML

my $report_id  = _safe($case->{report_id});
my $panel_name = _safe($case->{panel_name});
my $date       = _safe($case->{date});

print qq{
    <div class="card">
      <div class="back-link">
        <a href="/panel/search.html">&laquo; Back to search</a>
      </div>
      <h1>Case detail</h1>
      <div class="subtitle">
        case_id = } . _safe($case->{case_id}) . qq{ / report_id = } . $report_id . qq{ /
        <span class="tag-panel">} . $panel_name . qq{</span>
      </div>

      <h2>Basic information</h2>
      <table>
        <tr>
          <td class="label">Report ID</td>
          <td>} . $report_id . qq{</td>
          <td class="label">Panel type</td>
          <td>} . _safe($case->{panel_type}) . qq{</td>
        </tr>
        <tr>
          <td class="label">Vendor</td>
          <td>} . _safe($case->{vendor}) . qq{</td>
          <td class="label">Patient ID</td>
          <td>} . _safe($case->{patient_id}) . qq{</td>
        </tr>
        <tr>
          <td class="label">Test date (date)</td>
          <td>} . $date . qq{</td>
          <td class="label">Specimen ID</td>
          <td>} . _safe($case->{specimen_id}) . qq{</td>
        </tr>
        <tr>
          <td class="label">Sex / Age</td>
          <td>} . _safe($case->{sex}) . qq{ / } . _safe($case->{age}) . qq{</td>
          <td class="label">Test type</td>
          <td>} . _safe($case->{test_type}) . qq{</td>
        </tr>
        <tr>
          <td class="label">Disease</td>
          <td>} . _safe($case->{disease}) . qq{</td>
          <td class="label">Tissue of origin</td>
          <td>} . _safe($case->{tissue_of_origin}) . qq{</td>
        </tr>
        <tr>
          <td class="label">Pathology diagnosis</td>
          <td colspan="3">} . _safe($case->{pathology_diagnosis}) . qq{</td>
        </tr>
        <tr>
          <td class="label">Tumor nuclei %</td>
          <td>} . _safe($case->{percent_tumor_nuclei}) . qq{</td>
          <td class="label">Purity</td>
          <td>} . _safe($case->{purity}) . qq{</td>
        </tr>
        <tr>
          <td class="label">MSI status</td>
          <td>} . _safe($case->{msi_status}) . qq{</td>
          <td class="label">TMB</td>
          <td>} . _safe($case->{tmb_score}) . qq{ } . _safe($case->{tmb_unit}) . qq{ (} . _safe($case->{tmb_status}) . qq{)</td>
        </tr>
        <tr>
          <td class="label">Non-human content (attr)</td>
          <td>} . _safe($case->{non_human_content}) . qq{</td>
          <td class="label">Other info</td>
          <td>} . _safe($case->{other_info}) . qq{</td>
        </tr>
      </table>
};

# --- non_human_contents table ---
print qq{
      <h2>Non-human contents</h2>
};

if (@nh_rows) {
    print qq{
      <div class="scroll-box">
        <table>
          <thead>
            <tr>
              <th>Organism</th>
              <th>Reads per million</th>
              <th>Status</th>
              <th>Sample</th>
            </tr>
          </thead>
          <tbody>
    };
    for my $nh (@nh_rows) {
        my $org  = _safe($nh->{organism});
        my $rpm  = defined $nh->{reads_per_million} ? sprintf('%.0f', $nh->{reads_per_million}) : '';
        my $stat = _safe($nh->{status});
        my $smp  = _safe($nh->{sample});
        print "<tr>";
        print "<td class=\"mono\">$org</td>";
        print "<td class=\"num\">$rpm</td>";
        print "<td>$stat</td>";
        print "<td class=\"mono\">$smp</td>";
        print "</tr>\n";
    }
    print qq{
          </tbody>
        </table>
      </div>
      <div class="small-note">Non-human reads detected by, e.g., FoundationOne(Liquid).</div>
    };
} else {
    print qq{
      <div class="small-note">No non-human-content records are registered for this case.</div>
    };
}

# --- variants table ---
my $var_count = scalar @variants;

print qq{
      <h2>Variants ($var_count rows)</h2>
      <div class="scroll-box">
        <table>
          <thead>
            <tr>
              <th>Gene</th>
              <th>Type</th>
              <th>Subtype</th>
              <th>Chr:Pos</th>
              <th>Pos2</th>
              <th>cDNA</th>
              <th>Protein</th>
              <th>Strand</th>
              <th>Transcript</th>
              <th>Other gene</th>
              <th>In-frame</th>
              <th>Supp pairs</th>
              <th>Func effect</th>
              <th>Status</th>
              <th>Origin</th>
              <th>Class</th>
              <th>AF</th>
              <th>Depth</th>
              <th>Copy#</th>
              <th>CNV ratio</th>
              <th>TPM</th>
            </tr>
          </thead>
          <tbody>
};

for my $v (@variants) {
    my $af    = defined $v->{allele_fraction} ? sprintf('%.4f', $v->{allele_fraction}) : '';
    my $depth = defined $v->{depth}           ? $v->{depth} : '';
    my $cn    = defined $v->{copy_number}     ? sprintf('%.2f', $v->{copy_number}) : '';
    my $cnr   = defined $v->{cnv_ratio}       ? sprintf('%.2f', $v->{cnv_ratio}) : '';
    my $tpm   = defined $v->{tpm}             ? sprintf('%.2f', $v->{tpm}) : '';

    my $chrpos = '';
    if (defined $v->{chrom} && length $v->{chrom}) {
        $chrpos = $v->{chrom};
        $chrpos .= ':' . $v->{pos} if defined $v->{pos};
    }
    my $pos2 = defined $v->{pos2} ? $v->{pos2} : '';

    my $supp_pairs = defined $v->{supporting_read_pairs} ? $v->{supporting_read_pairs} : '';

    print "<tr>";
    print "<td class=\"mono\">" . _safe($v->{gene}) . "</td>";
    print "<td>" . _safe($v->{variant_type}) . "</td>";
    print "<td>" . _safe($v->{variant_subtype}) . "</td>";
    print "<td class=\"mono\">" . _safe($chrpos) . "</td>";
    print "<td class=\"mono\">" . _safe($pos2) . "</td>";
    print "<td class=\"mono\">" . _safe($v->{cds_effect}) . "</td>";
    print "<td class=\"mono\">" . _safe($v->{protein_effect}) . "</td>";
    print "<td>" . _safe($v->{strand}) . "</td>";
    print "<td class=\"mono\">" . _safe($v->{transcript}) . "</td>";
    print "<td>" . _safe($v->{other_gene}) . "</td>";
    print "<td>" . _safe($v->{in_frame}) . "</td>";
    print "<td class=\"num\">" . _safe($supp_pairs) . "</td>";
    print "<td>" . _safe($v->{functional_effect}) . "</td>";
    print "<td>" . _safe($v->{status}) . "</td>";
    print "<td>" . _safe($v->{origin}) . "</td>";
    print "<td>" . _safe($v->{classification}) . "</td>";
    print "<td class=\"num\">" . _safe($af) . "</td>";
    print "<td class=\"num\">" . _safe($depth) . "</td>";
    print "<td class=\"num\">" . _safe($cn) . "</td>";
    print "<td class=\"num\">" . _safe($cnr) . "</td>";
    print "<td class=\"num\">" . _safe($tpm) . "</td>";
    print "</tr>\n";
}

print qq{
          </tbody>
        </table>
      </div>
    </div> <!-- card -->
  </div> <!-- shell -->
</body>
</html>
};


