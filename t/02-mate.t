#!perl
#                           -*- mode: CPerl -*-
use 5.006;
use strict;
use warnings; # FATAL => 'all';
use Test::More;
use aliased 'PerlGP::Lite::Individual';
use PerlGP::Lite::GPMisc qw/copies/;

plan tests => 6;

my $filestem = 'tmp-indiv';
# remove the individual storage files
unlink glob "$filestem.*";

my $functions = {
		 ROOT => [ 'q({BAR})' ],
		 BAR => [ '{BAR}{BARX}', '{BARX}' ],
		};
my $terminals = {
		 ROOT => [ 'q()' ],
		 BAR => [ 'a' .. 'e' ],
		 BARX => [ 'A' .. 'E' ],
		};

my %params = (
	      grammar_branching => $functions,
	      grammar_nonbranching => $terminals,
	      crossover_probability => 0.5,
	      mutation_probability => 0.5,
	      minimum_genome_size => 20,
	     );

my $parent1 = Individual->new(filestem => $filestem.'.p1', %params);
my $parent2 = Individual->new(filestem => $filestem.'.p2', %params);
my $child1 =  Individual->new(filestem => $filestem.'.c1', %params);
my $child2 =  Individual->new(filestem => $filestem.'.c2', %params);

my $p1code = $parent1->getCode();
diag("p1: $p1code");
my $p2code = $parent2->getCode();
diag("p2: $p2code");

$parent1->crossover($parent2, $child1, $child2);

my $c1code = $child1->getCode();
diag("c1: $c1code");
my $c2code = $child2->getCode();
diag("c2: $c2code");

isnt($p1code, $c1code, "p1 c1 different");
isnt($p1code, $c2code, "p1 c2 different");
isnt($p2code, $c1code, "p2 c1 different");
isnt($p2code, $c2code, "p2 c2 different");

$child1->mutate();
my $c1code2 = $child1->getCode();
diag("mutate c1: $c1code -> $c1code2");
isnt($c1code2, $c1code, "c1 mutated ok");

$child2->mutate();
my $c2code2 = $child2->getCode();
diag("mutate c2: $c2code -> $c2code2");
isnt($c2code2, $c2code, "c2 mutated ok");


diag("NOTE: TESTS CAN FAIL BY CHANCE - not necessarily a bug!");

# remove the individual storage files
unlink glob "$filestem.*";

