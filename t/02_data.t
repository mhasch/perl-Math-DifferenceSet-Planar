# Copyright (c) 2019 Martin Becker, Blaubeuren.  All rights reserved.
# This package is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.

# Before 'make install' is performed this script should be runnable with
# 'make test'. After 'make install' it should work as 'perl 02_data.t'

#########################

use strict;
use warnings;
use File::Spec;

use Test::More tests => 28;
BEGIN { use_ok('Math::DifferenceSet::Planar::Data') };

#########################

my $data = Math::DifferenceSet::Planar::Data->new;
isa_ok($data, 'Math::DifferenceSet::Planar::Data');

my @dbs = Math::DifferenceSet::Planar::Data->list_databases;
ok(@dbs >= 1);
ok(1 == grep { $_ eq 'pds.db' } @dbs);

my $rec = $data->get(9);
isa_ok($rec, 'Math::DifferenceSet::Planar::Schema::Result::DifferenceSet');
is($rec->order, 9);
is($rec->base, 3);
is($rec->exponent, 2);
is($rec->modulus, 91);
is($rec->n_planes, 12);
like($rec->deltas, qr/^[0-9]+(?: [0-9]+){7}\z/);

$rec = $data->get(9, qw(modulus n_planes));
isa_ok($rec, 'Math::DifferenceSet::Planar::Schema::Result::DifferenceSet');
ok(!defined $rec->order);
ok(!defined $rec->base);
ok(!defined $rec->exponent);
is($rec->modulus, 91);
is($rec->n_planes, 12);
ok(!defined $rec->deltas);

my $it = $data->iterate(6, 8);
$rec = $it->();
is($rec->order, 7);
like($rec->deltas, qr/^[0-9]+(?: [0-9]+){5}\z/);
$rec = $it->();
is($rec->order, 8);
like($rec->deltas, qr/^[0-9]+(?: [0-9]+){6}\z/);
$rec = $it->();
ok(!$rec);

$it = $data->iterate_properties(8, 6);
$rec = $it->();
is($rec->order, 8);
ok(!defined $rec->deltas);
$rec = $it->();
is($rec->order, 7);
ok(!defined $rec->deltas);
$rec = $it->();
ok(!$rec);

