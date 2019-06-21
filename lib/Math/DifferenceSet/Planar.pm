package Math::DifferenceSet::Planar;

use strict;
use warnings;
use Carp qw(croak);
use Math::DifferenceSet::Planar::Data;
use Math::Prime::Util qw(is_prime_power euler_phi);

use constant F_ORDER    => 0;
use constant F_BASE     => 1;
use constant F_EXPONENT => 2;
use constant F_MODULUS  => 3;
use constant F_N_PLANES => 4;
use constant F_RESIDUES => 5;
use constant F_ROTATORS => 6;

our $VERSION = '0.002';

our $_MAX_ENUM_ORDER = 1024;
our $_LOG_MAX_ORDER = 21 * log(2);

my %memo_n_planes = ();

sub _coprime {
    my ($a, $b) = @_;
    while (0 < $b) {
        ($a, $b) = ($b, $a % $b);
    }
    return $a == 1;
}

sub _multipliers {
    my ($base, $exponent, $modulus) = @_;
    my @mult = (1);
    for (my $x = $base; $x != 1; $x = $x * $base % $modulus) {
        push @mult, $x;
    }
    die "assertion failed: not 3*$exponent elements: @mult"
        if @mult != 3 * $exponent;
    return @mult;
}

sub _rotators {
    my ($this) = @_;
    my (          $base,  $exponent,  $modulus,  $rotators) =
        @{$this}[F_BASE, F_EXPONENT, F_MODULUS, F_ROTATORS];
    return $rotators if @{$rotators};
    return undef if $this->[F_ORDER] > $_MAX_ENUM_ORDER;
    my @mult = _multipliers($base, $exponent, $modulus);
    my @sieve = (1) x $modulus;
    @sieve[@mult] = ();
    @{$rotators} = (1);
    for (my $x = 2; $x < $modulus ; ++$x) {
        if ($sieve[$x]) {
            if (0 == $modulus % $x) {
                for (my $i = $x; $i < $modulus; $i += $x) {
                    undef $sieve[$i];
                }
                next;
            }
            @sieve[ map { $x * $_ % $modulus } @mult ] = ();
            push @{$rotators}, $x;
        }
    }
    return $rotators;
}

sub _sequential_rotators {
    my ($this) = @_;
    my (          $base,  $exponent,  $modulus,  $n_planes) =
        @{$this}[F_BASE, F_EXPONENT, F_MODULUS, F_N_PLANES];
    my @mult = _multipliers($base, $exponent, $modulus);
    shift @mult;
    my $mx = 0;
    my $x  = 0;
    return sub {
        return 0 if $mx >= $n_planes;
        RESIDUE:
        while (1) {
            ++$x;
            next if !_coprime($modulus, $x);
            foreach my $e (@mult) {
                next RESIDUE if $x * $e % $modulus < $x;
            }
            ++$mx;
            return $x;
        }
    };
}

# small integer exponentiation
sub _pow {
    my ($base, $exponent) = @_;
    return 0 if !$base || log(0+$base) * $exponent >= $_LOG_MAX_ORDER;
    my $power = 1;
    while ($exponent) {
        $power *= $base if 1 & $exponent;
        $exponent >>= 1 and $base *= $base;
    }
    return $power;
}

# re-arrange elements of a planar difference set
sub _sort_residues {
    my $modulus  = shift;
    my @residues = sort { $a <=> $b } @_;
    my $lo = $residues[-1] - $modulus;
    my $hi = $residues[0];
    my $mx = 0;
    while ($hi - $lo != 1 && ++$mx < @residues) {
        ($lo, $hi) = ($hi, $residues[$mx]);
    }
    return undef if $mx >= @residues;
    $mx += @residues if !$mx--;
    push @residues, splice @residues, 0, $mx if $mx;
    return \@residues;
}

sub _residues_from_deltas {
    my $sum = 0;
    my @residues = map { $sum += $_ } 0, 1, split / /, $_[0];
    return \@residues;
}

# print "ok" if Math::DifferenceSet::Planar->available(9);
# print "ok" if Math::DifferenceSet::Planar->available(3, 2);
sub available {
    my ($class, $base, $exponent) = @_;
    my $order = defined($exponent)? _pow($base, $exponent): $base;
    my $pds = $order && Math::DifferenceSet::Planar::Data->get(
        $order, 'base'
    );
    return !!$pds && (!defined($exponent) || $base == $pds->base);
}

# $ds = Math::DifferenceSet::Planar->new(9);
# $ds = Math::DifferenceSet::Planar->new(3, 2);
sub new {
    my ($class, $base, $exponent) = @_;
    my $order = defined($exponent)? _pow($base, $exponent): $base;
    my $pds = $order && Math::DifferenceSet::Planar::Data->get($order);
    if (!$pds || defined($exponent) && $base != $pds->base) {
        my $key = defined($exponent)? "$base, $exponent": $order;
        croak "PDS($key) not available";
    }
    my $residues = _residues_from_deltas($pds->deltas);
    return bless [
        $pds->order,
        $pds->base,
        $pds->exponent,
        $pds->modulus,
        $pds->n_planes,
        $residues,
        [],
    ], $class;
}

# $ds = Math::DifferenceSet::Planar->from_residues(
#   0, 1, 3, 9, 27, 49, 56, 61, 77, 81
# );
sub from_residues {
    my $class    = shift;
    my $order    = $#_;
    my ($base, $exponent);
    $exponent    = is_prime_power($order, \$base)
        or croak "this implementation cannot handle order $order";
    my $modulus  = $order * ($order + 1) + 1;
    if (grep { $_ < 0 || $modulus <= $_ } @_) {
        my $max = $modulus - 1;
        croak "residue values inside range 0..$max expected";
    }
    my $residues = _sort_residues($modulus, @_);
    croak "residues of PDS expected" if !$residues;
    my $n_planes =
        $memo_n_planes{$order} ||= euler_phi($modulus) / (3 * $exponent);
    return bless [
        $order,
        $base,
        $exponent,
        $modulus,
        $n_planes,
        $residues,
        [],
    ], $class;
}

# $ds = Math::DifferenceSet::Planar->verify_residues(
#   0, 1, 3, 9, 27, 49, 56, 61, 77, 81
# );
sub verify_residues {
    my ($class, @residues) = @_;
    my $order   = $#residues;
    return undef if $order <= 1;
    my $modulus = $order * ($order + 1) + 1;
    my $median  = ($modulus - 1) / 2;
    my $seen    = '0' x $median;
    foreach my $r1 (@residues) {
        return undef if $r1 < 0 || $modulus <= $r1 || $r1 != int $r1;
        foreach my $r2 (@residues) {
            last if $r1 == $r2;
            my $d = $r1 < $r2? $r2 - $r1: $modulus + $r2 - $r1;
            $d = $modulus - $d if $d > $median;
            return q[] if substr($seen, $d-1, 1)++;
        }
    }
    return $median == $seen =~ tr/1//;
}

# $it1 = Math::DifferenceSet::Planar->iterate_available_sets;
# $it2 = Math::DifferenceSet::Planar->iterate_available_sets(10, 20);
# while (my $ds = $it2->()) {
#   ...
# }
sub iterate_available_sets {
    my ($class, @minmax) = @_;
    my $dit = Math::DifferenceSet::Planar::Data->iterate(@minmax);
    return sub {
        my $pds = $dit->();
        return undef if !$pds;
        my $residues = _residues_from_deltas($pds->deltas);
        return bless [
            $pds->order,
            $pds->base,
            $pds->exponent,
            $pds->modulus,
            $pds->n_planes,
            $residues,
            [],
        ], $class;
    };
}

# $om = Math::DifferenceSet::Planar->available_max_order;
sub available_max_order { Math::DifferenceSet::Planar::Data->max_order }

# $o  = $ds->order;
# $p  = $ds->order_base;
# $n  = $ds->order_exponent;
# $m  = $ds->modulus;
# $np = $ds->n_planes;
# @r  = $ds->residues;
# $r0 = $ds->residue(0);
sub order          {    $_[0]->[F_ORDER   ]            }
sub order_base     {    $_[0]->[F_BASE    ]            }
sub order_exponent {    $_[0]->[F_EXPONENT]            }
sub modulus        {    $_[0]->[F_MODULUS ]            }
sub n_planes       {    $_[0]->[F_N_PLANES]            }
sub residues       { @{ $_[0]->[F_RESIDUES]          } }
sub residue        {    $_[0]->[F_RESIDUES]->[$_[1]]   }

# $ds1 = $ds->translate(1);
sub translate {
    my ($this, $delta) = @_;
    my $modulus = $this->[F_MODULUS];
    $delta %= $modulus;
    return $this if !$delta;
    my @residues = map { ($_ + $delta) % $modulus } @{$this->[F_RESIDUES]};
    my $that = bless [@{$this}], ref $this;
    $that->[F_RESIDUES] = \@residues;
    return $that;
}

# $ds2 = $ds->normalize;
sub normalize { $_[0]->translate(- $_[0]->[F_RESIDUES]->[0]) }

# $it  = $ds->iterate_rotators;
# while (my $m = $it->()) {
#   ...
# }
sub iterate_rotators {
    my ($this) = @_;
    my $rotators = $this->_rotators;
    return $this->_sequential_rotators if !$rotators;
    my $mx = 0;
    return sub { $mx < @{$rotators}? $rotators->[$mx++]: 0 };
}

# $it = $ds->iterate_planes;
# while (my $ds = $it->()) {
#   ...
# }
sub iterate_planes {
    my ($this) = @_;
    my $r_it = $this->iterate_rotators;
    return sub {
        my $r = $r_it->();
        return $r? $this->multiply($r)->normalize: undef;
    };
}

# @pm = $ds->multipliers;
sub multipliers {
    my ($this) = @_;
    my (          $base,  $exponent,  $modulus) =
        @{$this}[F_BASE, F_EXPONENT, F_MODULUS];
    return _multipliers($base, $exponent, $modulus);
}

# $ds3 = $ds->multiply($m);
sub multiply {
    my ($this, $factor) = @_;
    my $modulus = $this->[F_MODULUS];
    $factor %= $modulus;
    croak "$_[1]: factor is not coprime to modulus"
        if !_coprime($modulus, $factor);
    return $this if 1 == $factor;
    my $residues = _sort_residues(
        $modulus,
        map { $_ * $factor % $modulus } @{$this->[F_RESIDUES]}
    );
    die "assertion failed: multiplied set is no PDS\n" if !$residues;
    my $that = bless [@{$this}], ref $this;
    $that->[F_RESIDUES] = $residues;
    return $that;
}

1;
__END__

=encoding utf8

=head1 NAME

Math::DifferenceSet::Planar - object class for planar difference sets

=head1 SYNOPSIS

  use Math::DifferenceSet::Planar;

  $ds = Math::DifferenceSet::Planar->new(9);
  $ds = Math::DifferenceSet::Planar->new(3, 2);
  $ds = Math::DifferenceSet::Planar->from_residues(
    0, 1, 3, 9, 27, 49, 56, 61, 77, 81
  );
  print "ok" if Math::DifferenceSet::Planar->verify_residues(
    0, 1, 3, 9, 27, 49, 56, 61, 77, 81
  );
  $o  = $ds->order;
  $m  = $ds->modulus;
  @r  = $ds->residues;
  $r0 = $ds->residue(0);
  $np = $ds->n_planes;
  $p  = $ds->order_base;
  $n  = $ds->order_exponent;

  $ds1 = $ds->translate(1);
  $ds2 = $ds->normalize;
  $ds2 = $ds->translate(- $ds->residue(0));
  @pm  = $ds->multipliers;
  $it  = $ds->iterate_rotators;
  while (my $m = $it->()) {
    $ds3 = $ds->multiply($m)->normalize;
  }
  $it = $ds->iterate_planes;
  while (my $ds3 = $it->()) {
    # as above
  }

  print "ok" if Math::DifferenceSet::Planar->available(9);
  print "ok" if Math::DifferenceSet::Planar->available(3, 2);
  $it1 = Math::DifferenceSet::Planar->iterate_available_sets;
  $it2 = Math::DifferenceSet::Planar->iterate_available_sets(10, 20);
  while (my $ds = $it2->()) {
    $o = $ds->order;
    $m = $ds->modulus;
    print "$o\t$m\n";
  }
  $om = Math::DifferenceSet::Planar->available_max_order;

=head1 DESCRIPTION

A planar difference set in a modular integer ring E<8484>_n, or cyclic
planar difference set, is a subset D = {d_1, d_2, ..., d_k} of E<8484>_n
such that each nonzero element of E<8484>_n can be represented as a
difference (d_i - d_j) in exactly one way.

Necessarily, for such a set to exist, the modulus n has to be equal to
(k - 1) E<183> k + 1.  If (k - 1) is a prime power, planar difference
sets can be constructed from a finite field of order (k - 1).  It is
conjectured that no other planar difference sets exist.  If other families
of planar difference sets should be discovered, this library would be
due to be extended accordingly.

If D = {d_1, d_2, ..., d_k} E<8834> E<8484>_n is a difference set and
a is an element of E<8484>_n, D + a = {d_1 + a, d_2 + a, ..., d_k + a}
is also a difference set.  D + a is called a translate of D.  The set
of all translates of a planar difference set as lines and the elements
of E<8484>_n as points make up a finite projective plane (hence the name).

If t is an element of E<8484>_n coprime to n, D E<183> t = {d_1 E<183> t,
d_2 E<183> t, ..., d_k E<183> t} is also a difference set.  If D E<183> t
is a translate of D, t is called a multiplicator of D.  If t is coprime
to n but either identical to 1 (mod n) or not a multiplicator, it is
called a rotator.  Rotators of planar difference sets are also rotators
of planes as translates of a difference set are mapped to translates of
the rotated set.  A minimal set of rotators spanning all plane rotations
is called a rotator base.

Math::DifferenceSet::Planar provides examples of small cyclic planar
difference sets constructed from finite fields.  It is primarily intended
as a helper module for algorithms employing such sets.  It also allows
to iterate over all sets of a given size via translations and rotations,
and to verify whether an arbitrary set of modular integers is a cyclic
planar difference set.

Currently, only sets with k E<8804> 10910, or moduli E<8804> 119017191,
are supported.  These limits can be extended by installing a database
with more samples.

=head1 CLASS VARIABLE

=over 4

=item I<$VERSION>

C<$VERSION> is the version number of the module.

=back

=head1 CLASS METHODS

=head2 Constructors

=over 4

=item I<new>

If C<$k> is a prime power, C<Math::DifferenceSet::Planar-E<gt>new($k)>
returns a sample planar difference set object with C<$k + 1> elements,
unless C<$k> exceeds some implementation limitation (see below).

If C<$p> is a prime number and C<$n> is an integer E<gt> 0,
C<Math::DifferenceSet::Planar-E<gt>new($p, $n)> returns a sample planar
difference set object with C<$p ** $n + 1> elements, unless C<$p ** $n>
exceeds some implementation limitation.

If C<$k> is not a prime power, or C<$p> is not a prime, or the number
of elements would exceed the limitation, an exception is raised.

=item I<from_residues>

If C<@d> is a cyclic planar difference set, represented as
distinct non-negative integer numbers less than some modulus,
C<Math::DifferenceSet::Planar-E<gt>from_residues(@d)> returns a planar
difference set object with precisely these elements.

Note that arguments not verified to define a planar difference set may
yield a broken object with undefined behaviour.  Note also that this
method expects residues to be normalized, i.e. integer values from zero
to the modulus minus one.

The modulus itself is not a parameter, as it can be computed from the
number I<k> of arguments as I<m> = I<k>E<178> - I<k> + 1.

=back

=head2 Other class methods

=over 4

=item I<verify_residues>

If C<@d> is an array of integer values,
C<Math::DifferenceSet::Planar-E<gt>verify_residues(@d)> returns a true
value if those values define a cyclic planar difference set and are
normalized, i.e. non-negative and less than the modulus, otherwise a
false value.  Note that this check is somewhat expensive, but should
work regardless of the method the set was constructed by.  It may thus
be used to verify cyclic planar difference sets this module would not
be capable of generating itself.

=item I<available>

The class method C<Math::DifferenceSet::Planar-E<gt>available(@params)>
checks whether I<new> can be called with the same parameters, i.e. either
an order C<$k> or a prime C<$p> and an exponent C<$n> indicating a
prime power order that are available from the database of PDS samples.
It returns a true value if sets with the given parameters are found,
otherwise false.

=item I<iterate_available_sets>

The class method
C<Math::DifferenceSet::Planar-E<gt>iterate_available_sets> returns a code
reference that, repeatedly called, returns all sample planar difference
sets known to the module, one by one.  The iterator returns a false
value when it is exhausted.

C<Math::DifferenceSet::Planar-E<gt>iterate_available_sets($lo, $hi)>
returns an iterator over all samples with orders between C<$lo>
and C<$hi>.

=item I<available_max_order>

The class method C<Math::DifferenceSet::Planar-E<gt>available_max_order>
returns the order of the largest sample planar difference set known to
the module.

=back

=head1 OBJECT METHODS

=head2 Constructors

=over 4

=item I<translate>

  $ds = $ds->translate(1);
  $ds = $ds->normalize;
  $ds = $ds->multiply($m);
  $it = $ds->iterate_planes;

=back

=head2 Property Accessors

=over 4

=item I<order>

  $o  = $ds->order;
  $m  = $ds->modulus;
  @r  = $ds->residues;
  $r0 = $ds->residue(0);
  $np = $ds->n_planes;
  $p  = $ds->order_base;
  $n  = $ds->order_exponent;

  @pm = $ds->multipliers;
  $it = $ds->iterate_rotators;

=back

=head1 DIAGNOSTICS

=over 4

=item PDS(%s) not available

The class method I<new> was called with parameters this implementation
does not cover.  The parameters are repeated in the message.  To avoid
this exception, verify the parameters using the I<available> method
before calling I<new>.

=item this implementation cannot handle order %d

The class method I<from_residues> was called with a number of residues
not equal to a prime power plus one.  The number of arguments minus one
is repeated in the message.  The given arguments may or may not define a
planar difference set, but if they were (i.e. I<verify_residues> called
with the same arguments returned true), the prime power conjecture would
be proven wrong.  Many mathematical journals would certainly be keen to
publish this counter-example.  Alternatively, you may report a bug in
this module's bug tracker.  Please include all arguments.

=item residue values inside range 0..%d expected

The class method I<from_residues> was called with residues that were
not normalized, i.e. integer values from zero to the modulus minus one,
or some values were too large for a difference set of the given size.
The modulus matching the number of arguments, minus one, is indicated
in the message.

=item residues of PDS expected

The class method I<from_residues> was called with values that obviously
define no planar difference set.  Note that not all cases of wrong values
will be detected this way.  Dubious values should always be verified
before they are turned into an object.

=item %d: factor is not coprime to modulus

The object method I<multiply> was called with an argument that was not an
integer coprime to the modulus.  The argument is repeated in the message.
Factors not coprime to the modulus would not yield a new difference set.

=back

=head1 BUGS AND LIMITATIONS

As this library depends on a database with sample sets, it will not
generate arbitrarily large sets.  The database packaged with the base
module is good for sets with at most 10910 elements.  An extension
by factor 6 is in preparation.  For much larger sets, the API should
presumably be changed to use PDL vectors rather than plain perl arrays,
to improve efficiency.

To handle difference sets on groups other than cyclic groups, some more
slight API changes would be required.  It should accept group elements as
well as small integers as arguments.  Although lacking practical examples,
this is intended to be dealt with in a future release.

Bug reports and suggestions are welcome.
Please submit them through the CPAN RT,
L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math-DifferenceSet-Planar>.

=head1 SEE ALSO

=over 4

=item *

Math::DifferenceSet::Planar::Data - planar difference set storage.

=item *

Math::ModInt - modular integer arithmetic.

=item *

PDL - the Perl Data Language.

=item *

Moore, Emily H., Pollatsek, Harriet S., "Difference Sets", American
Mathematical Society, Providence, 2013, ISBN 978-0-8218-9176-6.

=item *

Dinitz, J.H., Stinson, D.R., "Contemporary Design Theory: A collection
of surveys", John Wiley and Sons, New York, 1992, ISBN 0-471-53141-3.

=item *

Gordon, Daniel M., "La Jolla Difference Set Repository".
https://www.dmgordon.org/diffset/

=back

=head1 AUTHOR

Martin Becker, E<lt>becker-cpan-mp I<at> cozap.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2019 by Martin Becker, Blaubeuren.  All rights reserved.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

=head1 DISCLAIMER OF WARRANTY

This library is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of merchantability
or fitness for a particular purpose.

=cut
