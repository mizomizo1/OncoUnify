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
my $type  = lc($q->param('type') // '');
my $query = $q->param('q') // '';
my $limit = $q->param('limit') // 30;

$query =~ s/^\s+|\s+$//g;
$limit = 30 if $limit !~ /^\d+$/;
$limit = 100 if $limit > 100;

# JSON output
print $q->header(-type => 'application/json', -charset => 'utf-8');

sub json_escape {
  my ($s) = @_;
  $s //= '';
  $s =~ s/\\/\\\\/g;
  $s =~ s/"/\\"/g;
  $s =~ s/\r/\\r/g;
  $s =~ s/\n/\\n/g;
  $s =~ s/\t/\\t/g;
  return $s;
}

sub emit_json {
  my (@items) = @_;
  my $arr = join(",", map { '"' . json_escape($_) . '"' } @items);
  print qq|{"items":[${arr}]}|;
}

if ($query eq '') { emit_json(); exit; }

my $dbh = DBI->connect(
  "dbi:SQLite:dbname=$DB_FILE", "", "",
  {
    RaiseError => 1,
    sqlite_unicode => 1,
    AutoCommit => 1,
    ReadOnly => 1,
    sqlite_open_flags => DBD::SQLite::OPEN_READONLY(),
  }
) or do { emit_json(); exit; };

my @items;

if ($type eq 'gene') {
  # gene: prefix match (case-sensitivity depends on SQLite LIKE; ASCII names are fine)
  my $sth = $dbh->prepare(<<'SQL');
SELECT gene
FROM variants
WHERE gene IS NOT NULL AND TRIM(gene) != ''
  AND gene LIKE ?
GROUP BY gene
ORDER BY gene
LIMIT ?
SQL
  $sth->execute($query . "%", $limit);
  while (my ($v) = $sth->fetchrow_array) { push @items, $v; }
  $sth->finish;

} elsif ($type eq 'protein') {
  # protein: tolerate optional "p." prefix and collect candidates broadly
  my $norm = $query;
  $norm =~ s/^\s*p\.?\s*//i;  # strip leading 'p.' / 'p' prefix
  $norm =~ s/\s+//g;

  my @patterns;
  push @patterns, $query . "%";
  push @patterns, $norm . "%";
  push @patterns, "p." . $norm . "%";
  push @patterns, "p" . $norm . "%";

  my %seen; @patterns = grep { !$seen{$_}++ } @patterns;

  my $or = join(" OR ", map { "protein_effect LIKE ?" } @patterns);

  my $sql = "SELECT protein_effect
             FROM variants
             WHERE protein_effect IS NOT NULL AND TRIM(protein_effect) != ''
               AND ($or)
             GROUP BY protein_effect
             ORDER BY protein_effect
             LIMIT ?";

  my $sth = $dbh->prepare($sql);
  $sth->execute(@patterns, $limit);
  while (my ($v) = $sth->fetchrow_array) { push @items, $v; }
  $sth->finish;

} elsif ($type eq 'disease') {
  # disease: collect candidates from three case columns (prefix first, also substring)
  my $pat_prefix = $query . "%";
  my $pat_like   = "%" . $query . "%";

  my $sql = <<'SQL';
SELECT val FROM (
  SELECT disease AS val, 1 AS pri
    FROM cases
   WHERE disease IS NOT NULL AND TRIM(disease) != '' AND disease LIKE ?
  UNION
  SELECT tissue_of_origin AS val, 1 AS pri
    FROM cases
   WHERE tissue_of_origin IS NOT NULL AND TRIM(tissue_of_origin) != '' AND tissue_of_origin LIKE ?
  UNION
  SELECT pathology_diagnosis AS val, 1 AS pri
    FROM cases
   WHERE pathology_diagnosis IS NOT NULL AND TRIM(pathology_diagnosis) != '' AND pathology_diagnosis LIKE ?
  UNION
  SELECT disease AS val, 2 AS pri
    FROM cases
   WHERE disease IS NOT NULL AND TRIM(disease) != '' AND disease LIKE ?
  UNION
  SELECT tissue_of_origin AS val, 2 AS pri
    FROM cases
   WHERE tissue_of_origin IS NOT NULL AND TRIM(tissue_of_origin) != '' AND tissue_of_origin LIKE ?
  UNION
  SELECT pathology_diagnosis AS val, 2 AS pri
    FROM cases
   WHERE pathology_diagnosis IS NOT NULL AND TRIM(pathology_diagnosis) != '' AND pathology_diagnosis LIKE ?
)
GROUP BY val
ORDER BY pri, val
LIMIT ?
SQL

  my $sth = $dbh->prepare($sql);
  $sth->execute($pat_prefix, $pat_prefix, $pat_prefix,
                $pat_like,   $pat_like,   $pat_like,
                $limit);
  while (my ($v) = $sth->fetchrow_array) { push @items, $v; }
  $sth->finish;

} else {
  # unknown type
  $dbh->disconnect;
  emit_json();
  exit;
}

$dbh->disconnect;
emit_json(@items);

