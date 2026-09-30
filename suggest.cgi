#!/usr/bin/perl
# suggest.cgi — autocomplete endpoint returning {"items": [...]} (read-only).
#   type=gene     prefix match on variants.gene and variants.other_gene
#   type=protein  prefix match on the canonical protein key (variants.hgvs_p)
#   type=disease  prefix-then-substring match on the case-level disease columns
#   type=panel    list of panel names present in the database (q ignored)
# Configuration: ONCOUNIFY_DB (see panel_search.cgi).
use strict;
use warnings;
use utf8;
use CGI qw(:standard);
use DBI;
use DBD::SQLite ();
use JSON::PP ();

binmode(STDOUT, ':encoding(UTF-8)');

my $DB_FILE = $ENV{ONCOUNIFY_DB} || '/var/www/data/panels.db';

my $q = CGI->new;
my $type  = lc($q->param('type') // '');
my $query = $q->param('q') // '';
my $limit = $q->param('limit') // 30;
$query =~ s/^\s+|\s+$//g;
$limit = 30 if $limit !~ /^\d+$/ || $limit < 1;
$limit = 100 if $limit > 100;

print $q->header(-type => 'application/json', -charset => 'utf-8');
my $JSON = JSON::PP->new->utf8(0)->canonical(1);
sub emit { print $JSON->encode({ items => [@_] }); exit; }

sub _like { my $s = shift; $s =~ s/([\\%_])/\\$1/g; return $s; }

my %AA3 = (Ala=>'A',Arg=>'R',Asn=>'N',Asp=>'D',Cys=>'C',Gln=>'Q',Glu=>'E',Gly=>'G',His=>'H',Ile=>'I',
           Leu=>'L',Lys=>'K',Met=>'M',Phe=>'F',Pro=>'P',Ser=>'S',Thr=>'T',Trp=>'W',Tyr=>'Y',Val=>'V',Ter=>'*');
my $AA3_RE = join('|', sort { length($b) <=> length($a) } keys %AA3);

emit() if $query eq '' && $type ne 'panel';
emit() unless -r $DB_FILE;
my $dbh = eval {
    DBI->connect("dbi:SQLite:dbname=$DB_FILE", '', '', {
        RaiseError => 1, PrintError => 0, sqlite_unicode => 1, AutoCommit => 1, ReadOnly => 1,
        sqlite_open_flags => DBD::SQLite::OPEN_READONLY(),
    });
} or emit();

my $items = [];
if ($type eq 'gene') {
    my $p = _like(uc $query) . '%';
    $items = $dbh->selectcol_arrayref(<<'SQL', undef, $p, $p, $limit);
SELECT g FROM (
  SELECT gene AS g FROM variants WHERE gene LIKE ? ESCAPE '\'
  UNION
  SELECT other_gene AS g FROM variants WHERE other_gene LIKE ? ESCAPE '\'
) WHERE g IS NOT NULL AND g != '' ORDER BY g LIMIT ?
SQL
} elsif ($type eq 'protein') {
    my $s = $query;
    $s =~ s/\s+//g;
    $s =~ s/^\(?p\.//;
    $s =~ s/^p(?=[A-Z*(])//;
    $s =~ s/($AA3_RE)/$AA3{$1}/g;
    $items = $dbh->selectcol_arrayref(
        "SELECT DISTINCT hgvs_p FROM variants WHERE hgvs_p LIKE ? ESCAPE '\\' ORDER BY hgvs_p LIMIT ?",
        undef, 'p.' . _like($s) . '%', $limit);
} elsif ($type eq 'disease') {
    my ($pre, $sub) = (_like($query) . '%', '%' . _like($query) . '%');
    my @cols = qw(disease oncotree_code disease_ontology tissue_of_origin pathology_diagnosis);
    my $union = join("\n  UNION ALL\n", map {
        "SELECT $_ AS val, CASE WHEN $_ LIKE ? ESCAPE '\\' THEN 1 ELSE 2 END AS pri FROM cases WHERE $_ LIKE ? ESCAPE '\\'"
    } @cols);
    $items = $dbh->selectcol_arrayref(
        "SELECT val FROM ($union) WHERE val IS NOT NULL AND TRIM(val) != '' GROUP BY val ORDER BY MIN(pri), val LIMIT ?",
        undef, (map { ($pre, $sub) } @cols), $limit);
} elsif ($type eq 'panel') {
    $items = $dbh->selectcol_arrayref('SELECT DISTINCT panel_name FROM cases ORDER BY panel_name');
}
$dbh->disconnect;
emit(@$items);
