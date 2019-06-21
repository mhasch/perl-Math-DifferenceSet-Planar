package Math::DifferenceSet::Planar::Schema::Result::DifferenceSet;

=head1 NAME

Math::DifferenceSet::Planar::Schema::Result::DifferenceSet

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 VERSION

This documentation refers to version 0.002 of
Math::DifferenceSet::Planar::Schema::Result::DifferenceSet.

=cut

our $VERSION = '0.002';

=head1 TABLE: C<difference_set>

=cut

__PACKAGE__->table("difference_set");

=head1 ACCESSORS

=head2 order

  data_type: 'integer'
  is_nullable: 0

=head2 base

  data_type: 'integer'
  is_nullable: 0

=head2 exponent

  data_type: 'integer'
  is_nullable: 0

=head2 modulus

  data_type: 'integer'
  is_nullable: 0

=head2 n_planes

  data_type: 'integer'
  is_nullable: 0

=head2 deltas

  data_type: 'text'
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "order_",
  { accessor => "order", data_type => "integer", is_nullable => 0 },
  "base",
  { data_type => "integer", is_nullable => 0 },
  "exponent",
  { data_type => "integer", is_nullable => 0 },
  "modulus",
  { data_type => "integer", is_nullable => 0 },
  "n_planes",
  { data_type => "integer", is_nullable => 0 },
  "deltas",
  { data_type => "text", is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</order>

=back

=cut

__PACKAGE__->set_primary_key("order_");

1;
