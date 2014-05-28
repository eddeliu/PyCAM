#!/usr/bin/env perl

use strict;
use warnings;
use Point;

die 'give points number' unless defined $ARGV[0];
my $points = $ARGV[0];

my %existing_points;

#create points
my $ok = 0;
while ($ok != $points) {
	my $point = new Point(int rand(10), int rand(10));
	next if exists $existing_points{$point};
	$existing_points{$point} = $point;
	$ok++;
}

#compute and print dot
my @points = values %existing_points;
my $size = @points;
print "$size\n";
for my $i1 (0..$#points) {
	for my $i2 (0..$#points) {
		my ($p1, $p2) = map {$points[$_]} ($i1, $i2);
		my ($x1, $x2) = map {$_->get_x()} ($p1, $p2);
		my ($y1, $y2) = map {$_->get_y()} ($p1, $p2);
		my $d = sqrt(($y2-$y1)*($y2-$y1)+($x2-$x1)*($x2-$x1));
		print "$d ";
	}
    print "\n";
}

#compute and print svg
print STDERR "<svg width=\"100\" height=\"100\">\n";
for my $i (0..$#points) {
	my $x = $points[$i]->get_x()*10;
	my $y = $points[$i]->get_y()*10;
	print STDERR "<text x=\"$x\" y=\"$y\">$i</text>\n";
}
print STDERR "</svg>\n";
