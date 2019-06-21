package Math::DifferenceSet::Planar::Data;

use strict;
use warnings;
use File::Share qw(dist_file);
use DBD::SQLite::Constants qw(SQLITE_OPEN_READONLY);
use Math::DifferenceSet::Planar::Schema;

our $VERSION = '0.002';

my $DATABASE = _database();
my $schema   = undef;

sub _database { dist_file('Math-DifferenceSet-Planar', 'pds.sqlite') }

sub _schema {
    if (!defined $schema) {
        $schema = Math::DifferenceSet::Planar::Schema->connect(
            "dbi:SQLite:$DATABASE", q[], q[], 
            { sqlite_open_flags => SQLITE_OPEN_READONLY },
        );
    }
    return $schema;
}

sub get {
    my ($class, $order, @columns) = @_;
    return _schema->resultset('DifferenceSet')->search(
        { order_ => $order },
        @columns ? { columns => \@columns } : ()
    )->single;
}

sub iterate {
    my ($class, $min, $max) = @_;
    my @sel = ();
    push @sel, '>=' => $min if defined $min;
    push @sel, '<=' => $max if defined $max;
    my $results = _schema->resultset('DifferenceSet')->search(
        @sel? { order_ => { @sel } }: undef
    );
    return sub { $results->next };
}

sub max_order {
    return _schema->resultset('DifferenceSet')->get_column('order_')->max;
}

1;
