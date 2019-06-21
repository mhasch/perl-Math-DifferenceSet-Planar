# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl 01_basics.t'

#########################

use strict;
use warnings;

use Test::More tests => 78;
BEGIN { use_ok('Math::DifferenceSet::Planar') };

#########################

sub pds_ok {
    my ($ds, $res) = @_;
    isa_ok($ds, 'Math::DifferenceSet::Planar');
    SKIP: {
        if (eval { $ds->can('residues') }) {
            my @r = $ds->residues;
            is("@r", $res);
        }
        else {
            skip 'no residues from this object', 1;
        }
    }
};

ok(Math::DifferenceSet::Planar->available(9));
ok(Math::DifferenceSet::Planar->available(3, 2));
ok(!Math::DifferenceSet::Planar->available(6));
ok(!Math::DifferenceSet::Planar->available(4, 1));
ok(!Math::DifferenceSet::Planar->available(0, 2));

my $ds = Math::DifferenceSet::Planar->new(9);
pds_ok($ds, '0 1 3 9 27 49 56 61 77 81');

$ds = Math::DifferenceSet::Planar->new(3, 2);
pds_ok($ds, '0 1 3 9 27 49 56 61 77 81');

is($ds->order, 9);
is($ds->order_base, 3);
is($ds->order_exponent, 2);
is($ds->modulus, 91);
is($ds->residue(0), 0);
is($ds->residue(9), 81);
is($ds->residue(10), undef);

my $ds1 = $ds->translate(10);
pds_ok($ds1, '10 11 13 19 37 59 66 71 87 0');

my $ds2 = $ds1->normalize;
pds_ok($ds2, '0 1 3 9 27 49 56 61 77 81');
$ds2 = $ds->normalize;
pds_ok($ds2, '0 1 3 9 27 49 56 61 77 81');

SequentialIteratorTest: {
    local $Math::DifferenceSet::Planar::_MAX_ENUM_ORDER = 0;
    my $it = $ds->iterate_rotators;
    is(ref($it), 'CODE');
    my @p = ();
    while (my $ro = $it->()) {
      push @p, $ro;
    }
    is("@p", '1 2 4 5 8 10 16 19 20 23 29 46');
}

is($ds->n_planes, 12);
my $it = $ds->iterate_rotators;
my $it2 = $ds->iterate_rotators;
is(ref($it), 'CODE');
my @p = ();
while (my $ro = $it->()) {
  push @p, $ro;
}
is("@p", '1 2 4 5 8 10 16 19 20 23 29 46');
undef $it;
is($it2->(), 1);

my $it3 = $ds->iterate_planes;
pds_ok($it3->(), '0 1 3 9 27 49 56 61 77 81');
pds_ok($it3->(), '0 1 12 15 25 48 57 65 85 87');
my $c = 2;
while ($it3->()) {
    ++$c;
}
is($c, 12);

my @m = $ds->multipliers;
is("@m", '1 3 9 27 81 61');

my $ds3 = $ds->multiply(1);
pds_ok($ds3, '0 1 3 9 27 49 56 61 77 81');
$ds3 = $ds->multiply(10);
pds_ok($ds3, '90 0 10 14 30 35 42 64 82 88');
$ds3 = $ds->multiply(3);
pds_ok($ds3, '0 1 3 9 27 49 56 61 77 81');
$ds3 = $ds->multiply(16);
pds_ok($ds3, '48 49 53 56 66 68 77 0 16 22');

$ds3 = eval { $ds->multiply(35) };
is($ds3, undef);
like($@, qr/^35: factor is not coprime to modulus/);

$ds3 = eval { Math::DifferenceSet::Planar->new(6) };
is($ds3, undef);
like($@, qr/^PDS\(6\) not available/);

$ds3 = eval { Math::DifferenceSet::Planar->new(4, 1) };
is($ds3, undef);
like($@, qr/^PDS\(4, 1\) not available/);

$ds3 = eval { Math::DifferenceSet::Planar->new(0, 2) };
is($ds3, undef);
like($@, qr/^PDS\(0, 2\) not available/);

$ds3 = eval { Math::DifferenceSet::Planar->from_residues(
    1, 2, 4, 10, 28, 50, 57, 62, 78, 82,
)};
pds_ok($ds3, '1 2 4 10 28 50 57 62 78 82');

$ds3 = eval { Math::DifferenceSet::Planar->from_residues($ds3->residues) };
pds_ok($ds3, '1 2 4 10 28 50 57 62 78 82');

$ds3 = eval { Math::DifferenceSet::Planar->from_residues(
    -9, 2, 4, 10, 28, 50, 57, 62, 78
)};
is($ds3, undef);
like($@, qr/^residue values inside range 0\.\.72 expected/);

$ds3 = eval { Math::DifferenceSet::Planar->from_residues(0) };
is($ds3, undef);
like($@, qr/^this implementation cannot handle order 0/);

$ds3 = eval { Math::DifferenceSet::Planar->from_residues(
    1, 2, 4, 8, 16, 32, 21
)};
is($ds3, undef);
like($@, qr/^this implementation cannot handle order 6/);

my $rc = Math::DifferenceSet::Planar->verify_residues(1, 2, 4);
is($rc, 1);
$rc = Math::DifferenceSet::Planar->verify_residues(1, 2, 4, 1);
is($rc, q[]);
$rc = Math::DifferenceSet::Planar->verify_residues(1, 2, 3);
is($rc, q[]);
$rc = Math::DifferenceSet::Planar->verify_residues(-1, 2, 4);
is($rc, undef);
$rc = Math::DifferenceSet::Planar->verify_residues(1, 2.5, 4);
is($rc, undef);
$rc = Math::DifferenceSet::Planar->verify_residues(1, 2, 7);
is($rc, undef);
$rc = Math::DifferenceSet::Planar->verify_residues(0, 1);
is($rc, undef);
$rc = Math::DifferenceSet::Planar->verify_residues(4, 2, 1);
is($rc, 1);

SKIP: {
    if (Math::DifferenceSet::Planar->available(2, 10)) {
        my $ds4 = Math::DifferenceSet::Planar->new(2, 10);
        my $nm  = eval { $ds4->n_planes };
        diag("PDS(2, 10) has $nm planes\n") if defined $nm;
        ok(!defined($nm) || $nm == 19800);
        ok( defined($nm) || $@ =~ /^planes not supported for orders > /);
    }
    else {
        skip 'medium size sets not supported', 2;
    }
};

SKIP: {
    if (Math::DifferenceSet::Planar->available(3, 9)) {
        my $ds4 = Math::DifferenceSet::Planar->new(3, 9);
        my $nm  = eval { $ds4->n_planes };
        diag("PDS(3, 9) has $nm planes\n") if defined $nm;
        ok(!defined($nm) || $nm == 14183424);
        ok( defined($nm) || $@ =~ /^planes not supported for orders > /);
    }
    else {
        skip 'large sets not supported', 2;
    }
};

my $sit = Math::DifferenceSet::Planar->iterate_available_sets(10, 20);
my $isa = 1;
my @ords = ();
while (my $ds = $sit->()) {
    $isa &&= $ds->isa('Math::DifferenceSet::Planar');
    push @ords, $ds->order;
}
ok($isa);
is("@ords", '11 13 16 17 19');

$sit = Math::DifferenceSet::Planar->iterate_available_sets;
$isa = 1;
@ords = ();
while (my $ds = $sit->()) {
    $isa &&= $ds->isa('Math::DifferenceSet::Planar');
    push @ords, $ds->order;
    last if @ords >= 5;
}
ok($isa);
is("@ords", '2 3 4 5 7');
undef $sit;

my $max = Math::DifferenceSet::Planar->available_max_order;
diag("max available order is $max");
like($max, qr/^[1-9][0-9]*\z/);

