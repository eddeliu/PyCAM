package Point;
use strict;
use warnings;
use Carp;
use POSIX;
use Data::Dumper;
use overload
	'+' => \&plus,
	'!=' => \&isnot,
	'==' => \&is,
	'-' => \&minus,
	'*' => \&times,
	'""' => \&signature,
	'<' => \&less_than,
	'<=' => \&less_than_equal,
	'cmp' => \&compare;

#my $unit = 1000;
#my $precision = 12;
my $unit = 1;
my $precision = 8;
my $delta = 5 * 10**(-($precision+1));

sub new {
	my $class = shift;
	my $self = [@_];
	bless($self, $class);
	return $self;
}

sub plus {
	my $a = shift;
	my $b = shift;
	confess 'not same type of points' unless $#{$a} == $#{$b};
	my @coordinates;
	push @coordinates, $a->[$_] + $b->[$_] for (0..$#{$a});
	return new Point(@coordinates);
}

sub times {
	my $a = shift;
	my $b = shift;
	my $p;
	my $s;
	if(ref $a eq 'Point') {
		$p = $a;
		$s = $b;
	} else {
		$p = $b;
		$s = $a;
	}
	my @coordinates;
	push @coordinates, $p->[$_] * $s for (0..$#{$p});
	return new Point(@coordinates);
}

#TODO: rewrite faster version ?
sub minus {
	my $a = shift;
	my $b = shift;
	return $a + -1*$b;
}

sub get_x { return $_[0]->[0] }

sub get_y { return $_[0]->[1] }

sub get_z { return $_[0]->[2] }

sub set_x {
	my $self = shift;
	my $x = shift;
	$self->[0] = $x;
}

sub set_y {
	my $self = shift;
	my $y = shift;
	$self->[1] = $y;
}

sub get_coordinates { return @{$_[0]} }

#Taken from http://www.perlmonks.org/?node_id=656450

sub rotate {
	my $self = shift;
	my $angle = shift;
	confess 'no angle in point rotation' unless defined $angle;
	confess 'rotate only works on 2d points' unless @{$self}==2;
	my $x = $self->[0];
	my $y = $self->[1];
	#TODO: change -angle to +angle and verify everywhere that it's ok
	my $c = cos(-$angle);
	my $s = sin(-$angle);
	return new Point($x*$c - $y*$s, $x*$s + $y*$c);
}

sub in_unit_square {
	my $self = shift;
	confess 'in_unit_square only works on 2d points' unless @{$self}==2;
	my $x = $self->get_x();
	my $y = $self->get_y();
	#TODO: this is a ugly hack to get things working despite
	#the awful hell that are floating points numbers
	return ((abs($x)<=$unit + 0.00000001) and (abs($y)<=$unit + 0.00000001));
}

sub isnot {
	my $a = shift;
	my $b = shift;
	return not (is($a,$b));
}

sub is {
	my $a = shift;
	my $b = shift;
	confess 'undefined point in point comparison' unless defined $b;
	confess 'it only works on 2d points' unless @{$a}==2 and @{$b}==2;
	for my $index (0..$#{$a}) {
		my $c1 = $a->[$index];
		my $c2 = $b->[$index];
		if ((round_coordinate($c1)!=round_coordinate($c2)) and (round_coordinate($c1+$delta)!=round_coordinate($c2)) and round_coordinate($c1) != round_coordinate($c2+$delta)) {
			return 0;
		}
	}
	return 1;
}

sub signature {
	my $self = shift;
	my @string;
	for my $coordinate (@{$self}) {
		my $rounded_coordinate = round_coordinate($coordinate);
		push @string, $rounded_coordinate;
	}
	return join(',', @string);
}

sub cross_product {
	my $a = shift;
	my $b = shift;
	my ($x1, $y1, $z1) = @{$a};
	my ($x2, $y2, $z2) = @{$b};
	confess 'works only on 3d points' unless defined $z1 and defined $z2;
	my @coordinates;
	$coordinates[0] = $y1*$z2 - $z1*$y2;
	$coordinates[1] = $z1*$x2 - $x1*$z2;
	$coordinates[2] = $x1*$y2 - $y1*$x2;
	return new Point(@coordinates);
}

sub check_order {
	my $self = shift;
	my $direction = shift;
	confess 'works only on 3d points' unless defined $self->[2];
	#my $normalized_point = $self * (1/$d);
	my $product = $self->scalar_product($direction);
	return ($product > 0);
}

sub scalar_product {
	my $a = shift;
	my $b = shift;
	confess 'not same type of points' unless $#{$a} == $#{$b};
	my $r = 0;
	$r += $a->[$_] * $b->[$_] for(0..$#{$a});
	return $r;
}

sub less_than {
	my $a = shift;
	my $b = shift;
	confess 'not same type of points' unless $#{$a} == $#{$b};
	for my $index (0..$#{$a}) {
		return 1 if $a->[$index] < $b->[$index];
		return 0 if $a->[$index] > $b->[$index];
	}
	return 0;
}

sub compare {
	my $a = shift;
	my $b = shift;
	confess 'not same type of points' unless $#{$a} == $#{$b};
	for my $index (0..$#{$a}) {
		return -1 if $a->[$index] < $b->[$index];
		return 1 if $a->[$index] > $b->[$index];
	}
	return 0;
}

sub y_symmetry {
	my $self = shift;
	return new Point(-$self->[0], $self->[1]);
}

sub on_unit_circle {
	my $angle = shift;
	return new Point(cos($angle)*$unit, sin($angle)*$unit);
}

sub unit {
	return $unit;
}

sub precision {
	return $precision;
}

sub delta {
	return $delta;
}

sub round_coordinate {
	my $c = shift;
	my $r = sprintf("%.${precision}f", $c);
	$r = 0 if $r == -0;
	return $r;
}

sub round {
	my $self = shift;
	confess 'wrong type in Point::round' unless (ref $self eq 'Point');
	for my $index (0..$#{$self}) {
		$self->[$index] = round_coordinate($self->[$index]);
	}
}

sub fullprint {
	my $self = shift;
	print STDERR join(' ', map {sprintf('%.20f', $_)} @{$self}) . "\n";
}

1;
