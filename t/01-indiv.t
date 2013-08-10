#!perl -T
#                           -*- mode: CPerl -*-
use 5.006;
use strict;
use warnings; # FATAL => 'all';
use Test::More;
use aliased 'PerlGP::Lite::Individual';
use PerlGP::Lite::GPMisc qw/copies/;

plan tests => 5;

my $filestem = 'tmp-indiv';
my $functions = {
		 ROOT => [ 'something({FOO}, {BAR})' ],
		 BAR => [ '{BAR}{BARX}', '{BARX}' ],
		};
my $terminals = {
		 ROOT => [ 'nothing()' ],
		 FOO => [ 'foo1', 'foo2' ],
		 BAR => [ 'bar' ],
		 BARX => [ 'eee', 'ooo', copies(2, 'aaa') ],
		};

my $individual = Individual->new
  (
   filestem => $filestem,
   grammar_branching => $functions,
   grammar_nonbranching => $terminals,
  );


isa_ok($individual, 'PerlGP::Lite::Individual');
is($individual->{MacroMutationDepthBias}, 2, "PerlGP default attribs raw access");
is($individual->PointMutationFrac, 0.8, "PerlGP default attribs accessors");
$individual->PointMutationFrac(0.9);
is($individual->PointMutationFrac, 0.9, "PerlGP default attribs setter");
my $code = $individual->getCode();
diag($code);
like($code, qr/thing/, "some code was generated");


# remove the individual storage files
unlink "$filestem.dir", "$filestem.pag";

