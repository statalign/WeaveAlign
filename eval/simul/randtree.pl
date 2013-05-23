use Bio::TreeIO;
use Bio::Tree::RandomFactory;

use warnings;

@ARGV == 1 or print('Use: randtree.pl no_of_taxa\n'), exit;

# initialize a TreeIO writer to output the trees as we create them
my $str;
open $fh, ">", \$str;
my $out = Bio::TreeIO->new(-format => 'newick',
                           -fh   => $fh);
my $n = $ARGV[0];
my @taxa = map sprintf("tax%03d",$_), 0..$n-1;
my $factory = new Bio::Tree::RandomFactory(-taxa => \@taxa);
$out->write_tree($factory->next_tree);
$str =~ s/Node\d+//g;
print "$str\n";