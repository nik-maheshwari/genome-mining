use strict;
use warnings;

use lib "./"; #change this path to location of Bridge.pm
use Bridge qw(&brid);
my ($bridge,@pos,%temp1)= &brid("MNKLNSNAVVSLNEVSDSELDTILGGNRWWQGVVPTVSYECRMNSWQHVFTCC");
print $bridge;
