#                           -*- mode: CPerl -*-
package PerlGP::Lite::Individual;

use 5.006;
use strict;
use warnings; # FATAL => 'all';
use Mouse;
use SDBM_File;
# use DB_File;
use PerlGP::Lite::GPMisc qw/pickrandom poisson/;
use Fcntl;
use Digest::MD5 qw(md5_hex);
use Carp;
# for the alarm() calls
use POSIX ':signal_h';
sigaction SIGALRM, new POSIX::SigAction sub { die "alarmed (SigAction handler)" }
  or die "Error setting SIGALRM handler: $!\n";
# which replaces the next line which was unreliable and didn't work with 5.8.0
# actually the new SigAction stuff also causes panic: leave_scope inconsistency.
# but install the handler the old way for versions < 5.8.0
$SIG{ALRM} = sub { die "alarmed (\$SIG{ALRM} handler)" } if ($] < 5.008);
# to prevent deep recursion warnings during _tree_type_size
# and _fix_nodes and _crossover with deep trees
$SIG{__WARN__} = sub { print STDERR @_ unless "@_" =~ /Deep recursion/; };

=head1 NAME

PerlGP::Lite - The great new PerlGP::Lite!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

  TBC - see Lite.pm maybe

=head1 ATTRIBUTES

=head2 filestem

This is the full path (exluding filename extension) of the file
location of the genome-as-hash DBM file.

Can be relative (take care) or absolute!

=cut

has 'DBFileStem' => (
    init_arg => 'filestem',
    is => 'ro',
    isa => 'Str',
    required => 1,
    );


=head2 grammar_branching

Provide a hashref of %F as provided in old-style Grammar.pm
See http://perlgp.org/docs/manual/manual/node25.html#Grammar

=cut

has 'Functions' => (
    init_arg => 'grammar_branching',
    is => 'ro',
    isa => 'HashRef',
    required => 1,
    );

=head2 grammar_nonbranching

Provide a hashref of %T as provided in old-style Grammar.pm
See http://perlgp.org/docs/manual/manual/node25.html#Grammar

=cut

has 'Terminals' => (
    init_arg => 'grammar_nonbranching',
    is => 'ro',
    isa => 'HashRef',
    required => 1,
    );

=cut

=head2 quick_crossover_probability

QuickXoverProb

=cut

has 'QuickXoverProb' =>
  ( init_arg => 'quick_crossover_probability',
    is => 'rw',
    isa => 'Num',
    default => 0.9,
  );

=head2 crossover_probability

NodeXoverProb

=cut

has 'NodeXoverProb' =>
  ( init_arg => 'crossover_probability',
    is => 'rw',
    isa => 'Num',
    default => 0.01,
  );

=head2 crossover_depth

XoverDepthBias

greater than 1: more rooty
less than 1: more leafy
(must always be positive)

=cut

has 'XoverDepthBias' =>
  ( init_arg => 'crossover_depth',
    is => 'rw',
    isa => 'Num',
    default => 3,
  );

=head2 mutation_probability

NodeMutationProb

=cut

has 'NodeMutationProb' =>
  ( init_arg => 'mutation_probability',
    is => 'rw',
    isa => 'Num',
    default => 0.01,
  );

=head2 mutation_depth

MacroMutationDepthBias

>1 root bias
<1 leaf bias
(must always be positive)

=cut

has 'MacroMutationDepthBias' =>
  ( init_arg => 'macro_mutation_depth',
    is => 'rw',
    isa => 'Num',
    default => 2,
  );

=head2 minimum_genome_size

MinTreeNodes

=cut

has 'MinTreeNodes' =>
  ( init_arg => 'minimum_genome_size',
    is => 'rw',
    isa => 'Int',
    default => 0,
  );

=head2 numeric_mutation_probability

NumericMutationFrac

=cut

has 'NumericMutationFrac' =>
  ( init_arg => 'numeric_mutation_probability',
    is => 'rw',
    isa => 'Num',
    default => 0,
  );

=head2 numeric_mutation_types

NumericAllowNTypes

e.g. { CONST => 1.0 }

=cut

has 'NumericAllowNTypes' =>
  (
    init_arg => 'numeric_mutation_types',
    is => 'rw',
    isa => 'HashRef',
    default => sub { return {} },
  );

sub BUILD {
  my $self = shift;
  $self->_init();
}


# End of moose/mouse magic
__PACKAGE__->meta->make_immutable();
#############


# my $DBTYPE = 'SDBM_File';
my $DBTYPE = 'SDBM_File';
my $too_many_tries = 5000;

sub _init {
  my ($self) = @_;

 #################################################################
 # Paramater defaults, do not change them here, change the ones  #
 # you need in classes inheriting from this (e.g. Individual.pm) #
 #################################################################

  my %defaults = (

 ###### MUTATION ######
 # Probabilty that each node is mutated
		  ### NodeMutationProb => 1/100,
 # Or you can override this with a fixed number of mutations (beware bloat!)
 # and probability of this happening
                  FixedMutations => 0,
                  FixedMutationProb => 0,

 # Fraction of point mutations (the rest are macro mutations)
		  PointMutationFrac => 0.8,

 # Fraction of point mutations which involve constants that are randomly
 # adjusted rather than discretely
		  # NumericMutationFrac => 0.0,
 # Which types are allowed or ignored when doing numeric mutation
 # (for example you might not want some integers becoming floating point)
 # specify the following like this: { VARX=>1, CONSTY=>1 }
		  NumericIgnoreNTypes => {},
 # and/or the following like this: { VARX=>0.5, CONSTY=>0.2 }
 # where the value is the maximum change amount (*= or /= 1.5 and 1.2)
 # default is 0.1 (num = num*1.1 or num = num/1.1)
		  # NumericAllowNTypes => {},
 # the regular expression that defines what a number for numeric mutation is
 # for example you could set this to integers only: qr/^\d+$/ so that you
 # get 'one shot' numeric mutation
		  NumericMutationRegex =>
		  qr/^[+-]?\d+(?:\.\d+)?([eE][+-]?\d+)?$/,

 # When set, require that all mutations make some visible difference to code
		  NoNeutralMutations => 0,

 ### Depth bias: 1 means no bias, just pick a random node
 ### higher numbers mean that the nodes closer to the root are favoured
 ### smaller positive nonzero numbers mean nodes closer to the leaves
 # Depth bias for point mutations
		  PointMutationDepthBias => 0.5,
 # Depth bias for macro mutations
 #		  MacroMutationDepthBias => 2,

 # What types of macro mutation are used
 # (you can bias certain types by specifying them more than once)
 # 'encapsulate_subtree' and 'point_mutate_deep' are not used by default
		  MacroMutationTypes => [ qw(swap_subtrees copy_subtree
                                             replace_subtree insert_internal
                                             delete_internal) ],

 # If using 'encapsulate_subtree', don't encapsulate subtrees with these
 # node types (specify node types as keys, e.g. { NUM=>1, SUB=>1 })
		  EncapsulateIgnoreNTypes => {},
 # Maximum number of nodes allowed to encapsulate (fraction of total nodes)
		  EncapsulateFracMax => 0.25,
 # Probability that already encapsulated subtrees are used during
 # tree generation as functions or terminals (can make too-small trees)
		  UseEncapsTerminalsFrac => 0.0,

 # Logging mutation data to a file
		  MutationLogFile => 'results/mutation.log',
		  MutationLogProb => 1/50,

 ###### CROSSOVER ######
 # Per-node crossover-point selection probability
		  ### NodeXoverProb => 1/50,
 # Or you can override this with a fixed number of crossovers (beware bloat!)
 # and probability of this happening
		  FixedXovers => 0,
		  FixedXoverProb => 0,

 ### Crossover bias.  small (say 0.01) -> little bias
 ###                  large (say 10)   -> a lot of bias
 # Bias towards subtrees of the same size
		  XoverSizeBias => 1,
 # Bias towards subtrees with similar contents (by crude identity measure)
		  XoverHomologyBias => 1,
 # Quick Homologous Crossover (new in version 1.1)
 # 0 = off, 1 = always on
                #  QuickXoverProb => 0.9,


 # Depth bias for crossover point selection (see MacroMutationDepthBias)
 #		  XoverDepthBias => 3,

 # Only do asexual reproduction (simple copy of genomes)
		  AsexualOnly => 0,

 # Logging crossover data to a file
		  XoverLogFile => 'results/crossover.log',
		  XoverLogProb => 1/50,


 ###### Tree/subtree size/shape parameters
 # Maximum number of nodes allowed in tree as a whole
 # (random terminal nodes are used after this limit is reached)
		  MaxTreeNodes => 10000,
 # Minimum number of nodes for FRESHLY generated trees (_init_tree())
 # If the tree is too small, it will try again, until the tree is big enough.
 # So this could TAKE TIME...
                  ### MinTreeNodes => 0,

 ### The following five tree depth options must be changed if
 ### non-naturally terminating grammars are used
 # Maximum allowed depth of *new* trees or subtrees
		  TreeDepthMax => 20,
 # Probability that a terminal node is added during tree or subtree generation
		  TerminateTreeProb => 0,
 # Minimum depth before "TerminateTreeProb" takes effect
		  TreeDepthMin => 1,
 # Mean and maximum (cap) of Poisson distribution of new subtree depths
		  NewSubtreeDepthMean => 20,
		  NewSubtreeDepthMax => 20,

 # During new tree/subtree generation, force a certain fraction
 # of terminals to come only from existing terminals in the tree
 # (possibly useful if you're using numeric mutation)
		  UseExistingTerminalsFrac => 0.0,

 ###### Other things ######
 # Tells the getSize() method to ignore all nodes below nodes of types
 # specified (in the form { NTYPEX=>1, NTYPEY=>1 })
 # can be useful if your fitness function uses the size of the tree
		  GetSizeIgnoreNTypes => {},

 # if there's a syntax or other error in your evolved
 # subroutines, then PerlGP will sleep for a while before
 # reinitialising that individual.
 # set this to zero if you know what you're doing!
		  SleepAfterSubEvalError => 15,

		 );

  my $meta = __PACKAGE__->meta;
  foreach my $key (keys %defaults) {
    my $default = $defaults{$key};
    my $type = 'Str';
    if (my $reftype = ref($default)) {
      # $default = sub { $default };
      $type = ucfirst(lc($reftype)).'Ref';
    }
    $meta->add_attribute( $key => ( is => 'rw', isa => $type ) );
    $self->$key($default);
#old way:    $self->{$key} = $defaults{$key};
  }

}


# provide basic Get and Set routines using AUTOLOAD
#
# implicit Get if paramater exists:
# my $x = $obj->Xcoord;
# implicit Set if paramater exists, and you can chain them (returns $self)
# $obj->Xcoord(12.0)->Ycoord(3.5);
#
#sub AUTOLOAD {
#  my ($self, @args) = @_;
#
#  my ($mystery) = $AUTOLOAD =~ /(\w+)$/i;
#
#  if (exists $self->{$mystery}) {
#    if (@args > 1) {
#      # Set: multi-valued using array reference
#      $self->{$mystery} = [ @args ];
#    } elsif (@args) {
#      # Set: single valued
#      $self->{$mystery} = $args[0];
#    } else {
#      # Get: return value
#      return $self->{$mystery};
#    }
#  } else {
#    croak "unknown method/parameter $mystery accessed via AUTOLOAD for $self";
#  }
#  return $self;
#}


sub reInitialise {
  my $self = shift;

#  $self->evalEvolvedSubs();
#  $self->evolvedInit();
}

# you MUST override this method in your evolved code (in Grammar.pm)
# $input is the training or testing data structure
sub evaluateOutput {
  my ($self, $input) = @_;
  my $output = 'can be a scalar value or reference to data structure';
  return $output;
}

# you can override this method in your evolved code
# to give you evolvable parameters, for example
sub evolvedInit {
  my $self = shift;
  return;
}

# you can override this method to give info (like evolved params)
# which is used (at least) for the tournament log
sub extraLogInfo {
  my $self = shift;
  return $self->DBFileStem();
}


# tie the hash $self->{genome} to the file on disk
# return pointer to hash
# you must use an equal number of tieGenome and untieGenome !!
# nested tieGenomes are ignored using a depth counting system
sub tieGenome {
  my ($self, $debug) = @_;
  unless ($self->{tie_level}) {
    # warn "tie $debug\n" if ($debug);
    $self->{genome} = {};
    tie %{$self->{genome}},
      $DBTYPE, $self->{DBFileStem}, O_RDWR | O_CREAT, 0644;
    # completely initialise tree if the hash wasn't untied properly (see below)
    $self->_init_tree()
      if (defined $self->{genome}{'tied'} && $self->{genome}{'tied'});
    $self->{genome}{'tied'} = 1;
  }
  $self->{tie_level}++;
  return $self->{genome};
}

# rewrite the hash back to disk (optimises space usage)
sub retieGenome {
  my $self = shift;
  my $genome = $self->tieGenome('retie');
  my %oldgenome = %$genome;
  unlink glob("$self->{DBFileStem}*");
  tie %$genome, $DBTYPE, $self->{DBFileStem}, O_RDWR | O_CREAT, 0644;
  %$genome = %oldgenome;
  $self->untieGenome();
}

sub untieGenome {
  my $self = shift;
  $self->{tie_level}--;
  if ($self->{tie_level} == 0) {
    # this is untying properly (see above)
    $self->{genome}{'tied'} = 0;
    untie %{$self->{genome}};
  }
}

sub _tree_error {
  my ($self, $node, $msg) = @_;
  warn "node $node in genome not found in tree during $msg\n";
  # $self->_display_tree('root');
  $self->_init_tree();
  die "died after initialising genome\n";
}

sub initTree {
  my $self = shift;
  $self->tieGenome();
  $self->_init_tree();
  $self->untieGenome();
}

# assume genome is tied
sub _init_tree {
  my $self = shift;
#warn "DEBUG: initialising tree of ".$self->DBFileStem."\n";
  my $genome = $self->{genome};
  do {
    %{$genome} = ( 'tied'=>1, root=>'{nodeROOT0}', nodeROOT0=>'');
    # it's necessary to set nodeROOT0 twice, believe!
    $genome->{nodeROOT0} =
      $self->_grow_tree(depth=>0,
			TreeDepthMax=>$self->{TreeDepthMax},
			type=>'ROOT');
  } until (keys(%$genome)-2 >= $self->MinTreeNodes());
  $self->{sizeCalcRequired} = 1;
  # can't use $self->getSize because it does a tie() which can
  # create an endless recursive loop
}

sub initFitness {
  my $self = shift;
  my $genome = $self->tieGenome('initFitness');
  delete $self->{memory}{fitness};
  delete $genome->{fitness};
  $self->untieGenome();
}

sub eraseMemory {
  my $self = shift;
  $self->{memory} = { };
}

sub memory {
  my ($self, @args) = @_;
  if (@args == 1) {
    return $self->{memory}{$args[0]}; # was $arg[0] in old version BUG alert!
  } elsif (@args % 2 == 0) {
    while (@args) {
      my ($key, $value) = splice @args, 0, 2;
      $self->{memory}{$key} = $value;
    }
  } else {
    die "Individual::memory() called with odd number of elements in hash\n";
  }
}

sub getMemory {
  my ($self, $arg) = @_;
  return $self->memory($arg);
}

sub setMemory {
  my ($self, @args) = @_;
  $self->memory(@args);
}

sub getSize {
  my ($self) = @_;
  my $result;
  my $genome = $self->tieGenome('getSize');
  if ($self->{GetSizeIgnoreNTypes} &&
      keys %{$self->{GetSizeIgnoreNTypes}} && $genome->{root}) {
    $result = $self->_tree_type_size('root', undef, undef,
				     $self->{GetSizeIgnoreNTypes});
  } else {
    $result = scalar grep /^node/, keys %$genome;
  }
  $self->untieGenome();
  return $result;
}

# get and set Fitness routine
# fitness is stored in $self->{fitness} AND $self->{genome}{fitness}
# and some methods in this class access them directly
sub Fitness {
  my ($self, $setval) = @_;
  my $res;
  if (defined $setval) {
    my $genome = $self->tieGenome('setFitness');
    $res = $self->{memory}{fitness} = $genome->{fitness} = $setval;
    $self->untieGenome();
  } else {
    if (defined $self->{memory}{fitness}) { # get quick memory version
      $res = $self->{memory}{fitness};
    } else { # slower disk retrieval
      my $genome = $self->tieGenome('getFitness');
      $res = $genome->{fitness};
      $self->{memory}{fitness} = $res; # set the memory version
      $self->untieGenome();
    }
  }
  return $res;
}

# get or increment Age (age is only stored in memory)
sub Age {
  my ($self, $incr) = @_;
  my $res;
  $self->{memory}{age} = 0 unless (defined $self->{memory}{age});
  if (defined $incr) {
    $res = $self->{memory}{age} += $incr;
  } else {
    $res = $self->{memory}{age};
  }
  return $res;
}

sub getCode {
  my $self = shift;
  my $code = $self->_expand_tree();

  if (0) {
    # put the foreach loop before the code or in place of '__loophead__'
    my $loophead = 'foreach my $input (@$inputs) {
  undef $output;';
    my $looptail = '  push @results, $output;
}';
    $code =~ s{__loophead__}{$loophead};
    $code =~ s{__looptail__}{$looptail};
  }

  $self->{last_code} = $code;
  return $code;
}

sub evalEvolvedSubs {
  my $self = shift;
  my $code = $self->getCode();
  local $SIG{__WARN__} =
    sub { die @_ unless ($_[0] =~ /Subroutine.+redefined/ ||
			 $_[0] =~ /Attempt to free unreferenced scalar/) };
  eval $code;
  if ($@) {
    print STDERR "$@ during eval of gp subroutines - will reinit after sleep - code follows:\n\n$code";
    sleep $self->SleepAfterSubEvalError(); # because of the recursion!
    $self->initTree(); # so that next time at least the error won't be there
    $self->evalEvolvedSubs();
  }
}

sub _random_terminal {
  my ($self, $type) = @_;
  die "can't find terminal of type $type"
    unless (defined $self->{Terminals}{$type});

  return $self->_random_existing_terminal($type) ||
    pickrandom($self->{Terminals}{$type});
}

sub _random_existing_terminal {
  # assume tied.
  my ($self, $type, $encapsulated) = @_;
  my $prob = $self->{UseExistingTerminalsFrac} || 0;
  $prob = $self->{UseEncapsTerminalsFrac} || 0 if ($encapsulated);
  if ($prob && rand(1) < $prob) {
    my %seenalready;
    my @termnodes = grep {
      /^node$type\d/ &&
	$self->{genome}{$_} &&	# not a null string (or zero)
	  # starts with ;; if encapsulated node requested
	  (!$encapsulated || $self->{genome}{$_} =~ /^;;/) &&
	  $self->{genome}{$_} !~ /\{node[A-Z]+\d+\}/ && # and has no subnodes
	    !$seenalready{$self->{genome}{$_}}++ # just one of each!
                                                 # or nodes can saturate
	  } keys %{$self->{genome}};
    if (@termnodes) {
      my $terminalcopy = pickrandom(\@termnodes);
      $terminalcopy = $self->{genome}{$terminalcopy};
      return $terminalcopy;
    }
  }
  # returns '' if nothing found
  return '';
}

sub _random_function {
  my ($self, $type) = @_;
  if (!defined $self->{Functions}{$type}) {
    return $self->_random_terminal($type);
  } elsif ($self->{UseEncapsTerminalsFrac} &&
	   rand(1) < $self->{UseEncapsTerminalsFrac}) {
    return $self->_random_existing_terminal($type, 'encaps_only') ||
      pickrandom($self->{Functions}{$type});
  } else {
    return pickrandom($self->{Functions}{$type});
  }
}

sub _grow_tree {
  my $self = shift;
  my %p = @_;

  return $self->_random_terminal($p{type})
    if ($p{depth} >= $p{TreeDepthMax} ||
	($p{depth} >= $self->{TreeDepthMin} &&
	 rand()<$self->{TerminateTreeProb}));

  my $node = $self->_random_function($p{type});
  while ($node =~ m/{([A-Z]+)}/g) {
    my $ntype = $1;
    # use random node numbers so that we can do
    # quick homologous crossover
    my $nodeid = $self->nid();
    $nodeid = $self->nid() while (defined $self->{genome}{"node$ntype$nodeid"});
    my $newnode = "node$ntype$nodeid";
    $node =~ s/{$ntype}/{$newnode}/;
    $self->{genome}{$newnode} = '';
    $self->{genome}{$newnode} =
      $self->_grow_tree(depth=>$p{depth}+1,
			TreeDepthMax=>$p{TreeDepthMax}, type=>$ntype);
  }
  $self->{sizeCalcRequired} = 1;
  return $node;
}

sub nid {
  my $self = shift;
  return int(10*rand($self->{MaxTreeNodes}));
}

sub _expand_tree {
  my ($self, $node) = @_;
  my $genome = $self->tieGenome('expand');
  $node = 'root' unless ($node);

  $self->_init_tree() unless (defined $genome->{$node});

  my $x=0;
  my $maxn = $self->{MaxTreeNodes};
  my $code = $genome->{$node};
  while ($x++<$maxn && $code =~ s/{(node[A-Z]+\d+)}/$genome->{$1}/) {
  }

  if ($x>=$maxn) {
    $code =~ s/{node([A-Z]+)\d+}/$self->_random_terminal($1)/ge;
  }

  # delete encapsulated subtree size prefixes
  $code =~ s/;;\d+;;//g;

  $self->untieGenome();
  return $code;
}

sub crossover {
  my ($self, $mate, $recip1, $recip2) = @_;
  $self->reInitialise(); # if implemented, redefined xover params
  my $mygenome = $self->tieGenome('crossme');
  my $mategenome = $mate->tieGenome('crossmate');

  # how many crossovers will we do?
  my $numtodo = 0;
  unless ($self->FixedXovers()) {
    my $imax = keys %$mygenome;
    for (my $i=0; $i<$imax; $i++) {
      $numtodo++ if (rand() < $self->{NodeXoverProb});
    }
  } else {
    $numtodo = $self->FixedXovers()
      if (rand() < $self->FixedXoverProb());
  }

  # recursively determine all subtree sizes and node types for self and mate
  my @mynodes = grep !/ROOT/, grep /^node/, keys %$mygenome;
  my %mysizes; my %mytypes;
  $self->_tree_type_size('root', \%mysizes, \%mytypes);
  # sort nodes rootiest first
  @mynodes = sort { $mysizes{$b} <=> $mysizes{$a} } @mynodes;

  my @matenodes = grep !/ROOT/, grep /^node/, keys %$mategenome;
  my %matesizes; my %matetypes;
  $mate->_tree_type_size('root', \%matesizes, \%matetypes);
  @matenodes = sort { $matesizes{$b} <=> $matesizes{$a} } @matenodes;

  my ($samples, $maxsamples, $xovercount) = (0, $numtodo*100, 0);

  my (%myxsubnodes, %matexsubnodes); # these are the children of xover nodes
  my (%myxpair, %matexpair); # these are the xover nodes themselves (value=partners)
  my ($asexual, %pcid); # store the percent ids for logging

  if ($self->{AsexualOnly} || $numtodo == 0) {
    ($xovercount, $numtodo) = (1, 1);
    my ($myxover) = grep /ROOT/, grep /^node/, keys %$mygenome;
    my ($matexover) = grep /ROOT/, grep /^node/, keys %$mategenome;
    die "problem with asexual 'crossover' - no ROOT nodes"
      unless ($myxover && $matexover);
    $myxpair{$myxover} = $matexover;
    $matexpair{$matexover} = $myxover;
    $asexual = 1;
  }

  while ($xovercount < $numtodo) {
    # select one of my nodes - with optional bias towards root (XoverDepthBias > 1) or leaves (XoverDepthBias < 1)
    my $mynode = $mynodes[int(@mynodes * rand(1)**$self->{XoverDepthBias})];


    # select one of the mate's nodes
    my $matenode;

    # if a node with the same name exists in the mate,
    # use it (if QuickXoverProb is set accordingly)
    if (exists $matesizes{$mynode} &&
	$self->{QuickXoverProb} && rand() < $self->{QuickXoverProb}) {
      $matenode = $mynode;
    } else {
      # otherwise pick a random node of the same type
      # first restrict to same type
      my $mytype = $mytypes{$mynode};
      my @sametypematenodes = grep { $matetypes{$_} eq $mytype } @matenodes;
      # now pick one
      $matenode = $sametypematenodes[int(@sametypematenodes * rand(1)**$mate->{XoverDepthBias})];

      # do reverse quick homol check too
      if ($matenode && exists $mysizes{$matenode} &&
	  $mate->{QuickXoverProb} && rand() < $mate->{QuickXoverProb}) {
	$mynode = $matenode;
      }
    }

    if ($samples++>$maxsamples) {
      warn "xover too many tries ($xovercount xover pairs found out of $numtodo)\n";
      last;
    }

    # nodes have the same type structure (should already be) and are not in subtrees of previously
    # picked xover nodes - or have been used before
    if (defined $matenode &&
	$mytypes{$mynode} eq $matetypes{$matenode} &&
	!exists $myxsubnodes{$mynode} && !exists $matexsubnodes{$matenode} &&
	!exists $myxpair{$mynode} && !exists $matexpair{$matenode}) {

      my $smaller = $mysizes{$mynode} < $matesizes{$matenode} ?
	$mysizes{$mynode} : $matesizes{$matenode};  # size of the smallest
      my $bigger = $mysizes{$mynode} > $matesizes{$matenode} ?
	$mysizes{$mynode} : $matesizes{$matenode};  # size of the biggest

      my $id;
      # if we accept the two subtrees as similar size:
      if (rand() > (abs($mysizes{$mynode} - $matesizes{$matenode})/$bigger)
	            ** (1/$self->{XoverSizeBias}) &&
	  # and content (crude homology):
	  rand() < ((($id = $self->_tree_id($mate, $mynode, $matenode)) + 0.1)
		    /$smaller)**$self->{XoverHomologyBias}) {

	my $pcid = int(100*$id/$smaller);

	# check to see if these nodes have other xover points in their subtrees
	my @myxsubnodes = $self->_get_subnodes($mynode);
	my @matexsubnodes = $mate->_get_subnodes($matenode);
	my $problem = 0;
	grep { $problem++ if (exists $myxpair{$_}); } @myxsubnodes;
	grep { $problem++ if (exists $matexpair{$_}); } @matexsubnodes;
	if ($problem == 0) {
	  # remember all the subnodes for the future
	  grep { $myxsubnodes{$_} = 1 } @myxsubnodes;
	  grep { $matexsubnodes{$_} = 1 } @matexsubnodes;

	  # and store the crossover point relationships
          $myxpair{$mynode} = $matenode;
	  $matexpair{$matenode} = $mynode;
	  $pcid{$mynode} = $pcid;
          $xovercount++;
	  # warn $mynode eq $matenode ? "quick pair\n" : "random crossover pair: $mynode ($mysizes{$mynode}) $matenode ($matesizes{$matenode})\n";
	}
      }
    }
  }

  # copy and paste alert (from a few blocks above)
  if ($xovercount == 0) {
    warn "no crossover points (of $numtodo) found - doing asexual\n";
    my ($myxover) = grep /ROOT/, grep /^node/, keys %$mygenome;
    my ($matexover) = grep /ROOT/, grep /^node/, keys %$mategenome;
    die "problem with asexual 'crossover' - no ROOT nodes"
      unless ($myxover && $matexover);
    $myxpair{$myxover} = $matexover;
    $matexpair{$matexover} = $myxover;
    $asexual = 1;
    $xovercount++;
  }

  # warn "going to do $xovercount of planned $numtodo crossovers\n";
  # now we actually do the cross over(s)
  if ($xovercount > 0) {

    # do some logging
    if ($self->XoverLogProb() &&
	rand() < $self->XoverLogProb() &&
	$self->XoverLogFile()) {
      if (open(FILE, ">>$self->{XoverLogFile}")) {
	if ($asexual) {
	  print FILE "doing 0 asexual\n";
	} else {
	  printf FILE "doing %d xovers\n", $xovercount;
	  foreach my $mynode (keys %myxpair) {
	    printf FILE "nodes %s %s sizes %d %d identity %d\n",
	      $mynode, $myxpair{$mynode},
		$mysizes{$mynode}, $matesizes{$myxpair{$mynode}},
		  $pcid{$mynode};
	  }
	}
	close(FILE);
      }
    }

    # do the crossovers
    $self->_start_crossover($mate, $recip1, \%myxpair);
    $mate->_start_crossover($self, $recip2, \%matexpair);
  }
  $self->untieGenome(); $mate->untieGenome();
  $recip1->{sizeCalcRequired} = 1;
  $recip2->{sizeCalcRequired} = 1;
}

sub _start_crossover {
  my ($self, $mate, $recipient, $selfxnode, $matexnode) = @_;
  my $rgenome = $recipient->tieGenome('crossrecip');
  %$rgenome = (); # wipe recipient's genome
  $recipient->retieGenome(); # actually wipe the DBM file
  $recipient->eraseMemory(); # erase the quick memory copies (fitness,age,...)
  $rgenome->{'tied'} = 1;
  $self->_crossover($mate, $recipient, 'root', '', $selfxnode);
  $recipient->_fix_nodes('root');
  $recipient->untieGenome();
}

sub _has_children {
  my ($self, $node) = @_;

  $self->{sizeCalcRequired} = 1 if (! exists $self->{sizeCalcRequired});
  $self->{nodesLookup} = {} if (! exists $self->{nodesLookup});

  my $children;
  if (($self->{sizeCalcRequired} == 1) || ! (exists $self->{nodesLookup}{$node})) {
    $children = $self->{genome}{$node} =~ qr/{node[A-Z]+\d+}/;
  }
  else {
    $children = @{$self->{nodesLookup}{$node}{subnodes}};
  }
  return $children;
}

sub _children {
  my ($self, $node) = @_;

  $self->{sizeCalcRequired} = 1 if (! exists $self->{sizeCalcRequired});
  $self->{nodesLookup} = {} if (! exists $self->{nodesLookup});

  my $children;
  if (($self->{sizeCalcRequired} == 1) || ! (exists $self->{nodesLookup}{$node})) {
    $children = $self->{genome}{$node} =~ qr/{node[A-Z]+\d+}/;
  }
  else {
    $children = $self->{nodesLookup}{$node}{subnodes};
  }
  return $children;
}

sub _tree_id_lookup_no_rec {
  my ($self, $mate, $mystartnode, $matestartnode) = @_;
  my @nodestack = [];
  my $sum = 0;
  my $mynode;
  my $matenode;
  push @nodestack, $mystartnode;
  push @nodestack, $matestartnode;
  while ($matenode = pop(@nodestack)) {
    print "Mate:".$matenode.":";
    $mynode = pop(@nodestack);
    print ", My:".$mynode."\n";
    if (@{$self->{nodesLookup}{$mynode}{subnodes}}) {
      if ($mate->_children($matenode)) {
        my $mycopy = $self->{genome}{$mynode};
        my $matecopy = $mate->{genome}{$matenode};
        $mycopy =~ s/{node([A-Z]+)\d+}/{$1}/g;
        $matecopy =~ s/{node([A-Z]+)\d+}/{$1}/g;
        if ($mycopy eq $matecopy) {
	  $sum += 1;
          my @mysubnodes = $self->_children($mynode);
          my @matesubnodes = $mate->_children($matenode);
          my $i;
          for ($i=0; $i<@mysubnodes; $i++) {
            push @nodestack, $mysubnodes[$i];
            push @nodestack, $matesubnodes[$i];
          }
        }
      } 
    }
    else { # they could both be terminal nodes
      $sum += 1 if ($self->{genome}{$mynode} eq $mate->{genome}{$matenode});
    }
  }
  return $sum;
}

sub _tree_id_lookup {
  my ($self, $mate, $mynode, $matenode) = @_;
  if (@{$self->{nodesLookup}{$mynode}{subnodes}}) {
    if ($mate->_has_children($matenode)) {
      my $mycopy = $self->{genome}{$mynode};
      my $matecopy = $mate->{genome}{$matenode};
      $mycopy =~ s/{node([A-Z]+)\d+}/{$1}/g;
      $matecopy =~ s/{node([A-Z]+)\d+}/{$1}/g;
      if ($mycopy eq $matecopy) {
	my $sum = 1;
	my @mysubnodes = $self->{genome}{$mynode} =~ /{(node[A-Z]+\d+)}/g;
	my @matesubnodes = $mate->{genome}{$matenode} =~ /{(node[A-Z]+\d+)}/g;
	my $i;
	for (my $i=0; $i<@mysubnodes; $i++) {
	  $sum += $self->_tree_id_lookup($mate, $mysubnodes[$i], $matesubnodes[$i]);
	}
	return $sum;
      }
      else {
	return 0;
      }
    }
    else {
      return 0;
    }
  }
  else { # they could both be terminal nodes
    return 1 if ($self->{genome}{$mynode} eq $mate->{genome}{$matenode});
  }
  return 0;
}

# assume tied
sub _tree_id_search {
  my ($self, $mate, $mynode, $matenode) = @_;
  if ($self->{genome}{$mynode} =~ qr/{node[A-Z]+\d+}/) {
    if ($mate->_has_children($matenode)) {
      my $mycopy = $self->{genome}{$mynode};
      my $matecopy = $mate->{genome}{$matenode};
      $mycopy =~ s/{node([A-Z]+)\d+}/{$1}/g;
      $matecopy =~ s/{node([A-Z]+)\d+}/{$1}/g;
      if ($mycopy eq $matecopy) {
	my $sum = 1;
	my @mysubnodes = $self->{genome}{$mynode} =~ /{(node[A-Z]+\d+)}/g;
	my @matesubnodes = $mate->{genome}{$matenode} =~ /{(node[A-Z]+\d+)}/g;
	my $i;
	for (my $i=0; $i<@mysubnodes; $i++) {
	  $sum += $self->_tree_id_search($mate, $mysubnodes[$i], $matesubnodes[$i]);
	}
	return $sum;
      }
      else {
	return 0;
      }
    }
    else {
      return 0;
    }
  }
  else { # they could both be terminal nodes
    return 1 if ($self->{genome}{$mynode} eq $mate->{genome}{$matenode});
  }
  return 0;
}

# in _tree_id assume genomes are tied for speed
sub _tree_id {
  my ($self, $mate, $mynode, $matenode) = @_;

  $self->{sizeCalcRequired} = 1 if (! exists $self->{sizeCalcRequired});
  $self->{nodesLookup} = {} if (! exists $self->{nodesLookup});

  my $sum;
  if (($self->{sizeCalcRequired} == 1) || ! (exists $self->{nodesLookup}{$mynode})) {
    $sum = $self->_tree_id_search($mate, $mynode, $matenode);
  }
  else {
    $sum = $self->_tree_id_lookup($mate, $mynode, $matenode);
  }
}

# assume genomes are tied for speed
sub _crossover {
  my ($self, $mate, $recip,
      $selfnode, $matenode, $myxpoint) = @_;
  my $newnode;

  if ($selfnode) {
    $newnode = ($recip->{genome}{$selfnode} = $self->{genome}{$selfnode});
  } else {
    $newnode = ($recip->{genome}{$matenode.'x'} = $mate->{genome}{$matenode});
  }

  my @subnodes = $newnode =~ /{(node[A-Z]+\d+)}/g;
  if ($matenode) {
    $recip->{genome}{$matenode.'x'} =~ s/{(node[A-Z]+\d+)}/{$1x}/g;
  }

  foreach my $subnode (@subnodes) {
    my ($mynext, $matenext) = ($subnode, '');
    if ($selfnode && exists $myxpoint->{$subnode}) {
      ($mynext, $matenext) = ('', $myxpoint->{$subnode});
      $recip->{genome}{$selfnode} =~ s/{$subnode}/{$myxpoint->{$subnode}x}/;
    } elsif ($matenode) {
      ($mynext, $matenext) = ('', $subnode);
    }
    $self->_crossover($mate, $recip, $mynext, $matenext, $myxpoint);
  }
  $recip->{sizeCalcRequired} = 1;
}

# assume genome is tied
sub _fix_nodes {
  my ($self, $node) = @_;

  my @subnodes = $self->{genome}{$node} =~ /{(node[A-Z]+\d+x?)}/g;
  foreach my $subnode (@subnodes) {
    if ($subnode =~ /x$/) {
      my $fixed = $subnode;
      $fixed =~ s/x$//;
      $fixed =~ s/(\d+)$/$self->nid()/e while (defined $self->{genome}{$fixed});
      $self->{genome}{$node} =~ s/$subnode/${fixed}/;
      # ${fixed} is like this because emacs perl mode was playing up...!
      $self->{genome}{$fixed} = $self->{genome}{$subnode};
      delete $self->{genome}{$subnode};
      $subnode = $fixed;
    }
    $self->_fix_nodes($subnode);
  }
  $self->{sizeCalcRequired} = 1;
}

sub _random_node {
  my ($self, %p) = @_;

  my $depth_bias = $p{depth_bias} || 1;

  my %mysizes;
  my %mytypes;
  my $start_node = $p{start_node} || 'root';
  $self->_tree_type_size($start_node, \%mysizes, \%mytypes);
  my @mynodes = grep /^node(?!ROOT)/, keys %mysizes;

  if (defined $p{node_type}) {
    @mynodes = grep { $mytypes{$_} eq $p{node_type} } @mynodes;
  }
  if (defined $p{not_this_node}) {
    @mynodes = grep { $_ ne $p{not_this_node} } @mynodes;
  }
  if (defined $p{not_this_subtree}) {
    my %subnodes;
    grep { $subnodes{$_} = 1 } $self->_get_subnodes($p{not_this_subtree});
    $subnodes{$p{not_this_subtree}} = 1;
    @mynodes = grep { !exists $subnodes{$_} } @mynodes;
  }

  # sort rootiest first
  @mynodes = sort { $mysizes{$b} <=> $mysizes{$a} } @mynodes;

  # return one biased to beginning (bias>1) or end (bias<1) of list
  return $mynodes[int(@mynodes * rand(1)**$depth_bias)];
}


sub mutate {
  my ($self, %p) = @_;
  $self->reInitialise();
  my $genome = $self->tieGenome('mutate');

  my $i;
  $self->{mutednodes} = {};

  my $mutes_to_do = 0;
  unless ($self->FixedMutations()) {
    my $imax = keys %$genome;
    for (my $i=0; $i<$imax; $i++) {
      $mutes_to_do++ if (rand() < $self->{NodeMutationProb});
    }
  } else {
    $mutes_to_do = $self->FixedMutations()
      if (rand() < $self->FixedMutationProb());

    my $imax = keys %$genome;
    for (my $i=0; $i<$imax; $i++) {
      $i++ if (rand() < $self->{NodeMutationProb});
    }
  }

  for (my $i=0; $i<$mutes_to_do; $i++) {
    my $unmutedcode = $self->{NoNeutralMutations} ? $self->_expand_tree() : '';

    if (rand() < $self->{PointMutationFrac}) {
      $self->point_mutate_shallow();
    } else {
      $self->macro_mutate();
    }

    # if necessary, check that the mutation actually changed the program
    if ($self->{NoNeutralMutations}) {
      my $mutedcode = $self->_expand_tree();
      if ($mutedcode eq $unmutedcode) {
	# allow another pass through mutation loop
	$mutes_to_do++ if ($mutes_to_do<$too_many_tries);
      }
    }
  }

  # reset fitness if we did any mutations
  if ($mutes_to_do) {
    $self->initFitness();
  }

  # do some logging
  if ($mutes_to_do && $self->MutationLogProb() &&
      rand() < $self->MutationLogProb() &&
      $self->MutationLogFile()) {
    if (open(FILE, ">>$self->{MutationLogFile}")) {
      printf FILE "did %d of %d\n",
	scalar(keys %{$self->{mutednodes}}), $mutes_to_do;
      foreach my $node (keys %{$self->{mutednodes}}) {
	my $size = 0;
	# node might not exist any more
	if (defined $genome->{$node}) {
	  $size = $self->_tree_type_size($node);
	}
	printf FILE "node $node size %d %s\n", $size,
	  $self->{mutednodes}{$node};
      }
      close(FILE);
    }
  }

  $self->untieGenome();
}

sub point_mutate_shallow {
  my $self = shift;
  $self->point_mutate($self->{PointMutationDepthBias});
}

sub point_mutate_deep {
  my $self = shift;
  $self->point_mutate($self->{MacroMutationDepthBias});
}

sub point_mutate {
  my ($self, $depth_bias) = @_;
  my $genome = $self->{genome};
  $depth_bias = $self->{PointMutationDepthBias}
    unless defined ($depth_bias);
  my $mutnode =
    $self->_random_node(depth_bias=>$depth_bias);
  return unless ($mutnode);

  my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

  # if it's an internal node
  if ($genome->{$mutnode} =~ qr/{node[A-Z]+\d+}/) {
    my @subnodes = $genome->{$mutnode} =~ /{(node[A-Z]+\d+)}/g;
    my @subtypes = $genome->{$mutnode} =~ /{node([A-Z]+)\d+}/g;
    my $z = 0; my @newtypes;
    while ($z++<$too_many_tries) {
      my $newnode = $self->_random_function($ntype);
      @newtypes = $newnode =~ /{([A-Z]+)}/g;
      if ("@subtypes" eq "@newtypes") {
	$newnode =~ s/{[A-Z]+}/'{'.(shift @subnodes).'}'/ge;
	$genome->{$mutnode} = $newnode;
	$self->{mutednodes}{$mutnode} = 'point_internal';
	last;
      }
    }
    # Changed the structure, so size may need recalculating
    $self->{sizeCalcRequired} = 1;
  } else { # it's a terminal node
    # is it a number? (doesn't catch ".123")
    # and are we allowed to mutate them by multiplication?
    if (rand() < $self->{NumericMutationFrac} &&
	!$self->{NumericIgnoreNTypes}{$ntype} &&
	($self->{NumericAllowNTypes}{$ntype} ||
	 keys %{$self->{NumericAllowNTypes}} == 0) &&
	$genome->{$mutnode} =~ $self->{NumericMutationRegex}) {
      my $amount = $self->{NumericAllowNTypes}{$ntype} || 0.1;
      if (rand() < 0.5) {
	$genome->{$mutnode} *= 1+rand($amount); # small random change
      } else {
	$genome->{$mutnode} /= 1+rand($amount); # and the other way
      }
      # now keep it with a sensible number of significant figures
      $genome->{$mutnode} = sprintf "%.4g", $genome->{$mutnode};
      $self->{mutednodes}{$mutnode} = 'point_numeric';
    } else { # normal mutation by replacement
      $genome->{$mutnode} = $self->_random_terminal($ntype);
      $self->{mutednodes}{$mutnode} = 'point_terminal';
    }
  }
}

sub macro_mutate {
  my $self = shift;
  my $macro_type = pickrandom($self->{MacroMutationTypes});
  $self->$macro_type();
}

sub replace_subtree {
  my ($self) = @_;
  my $genome = $self->{genome};
  my $mutnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias});
  return unless ($mutnode);
  my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

  # replace subtree with random subtree tree
  $self->_del_subtree($mutnode);
  my $newdepth = poisson($self->{NewSubtreeDepthMean},
				 $self->{NewSubtreeDepthMax});
  $genome->{$mutnode} = '';
  $genome->{$mutnode} =
    $self->_grow_tree(depth=>0, TreeDepthMax=>$newdepth, type=>$ntype);
  $self->{mutednodes}{$mutnode} = 'replace_subtree';
  $self->{sizeCalcRequired} = 1;
}

sub insert_internal {
  my ($self) = @_;
  my $genome = $self->{genome};
  my $mutnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias});
  return unless ($mutnode);
  my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

  my $nodeid = $self->nid();
  $nodeid = $self->nid() while (defined $genome->{"node$ntype$nodeid"});
  my $newnode = "node$ntype$nodeid";
  my $copy = $genome->{$mutnode};
  my $insert;
  my $z = 0;
  do {
    $insert = $self->_random_function($ntype);
  } until ($z++>$too_many_tries || $insert =~ /{$ntype}/);
  if ($z<$too_many_tries) { # otherwise don't bother
    # copy mutnode to newnode so we can insert the new node at mutnode
    $genome->{$newnode} = $genome->{$mutnode};
    $genome->{$mutnode} = $insert;

    $self->{mutednodes}{$mutnode} = 'insert_internal';
    $self->{mutednodes}{$newnode} = 'insert_internal';

    my @posslinknodes = $genome->{$mutnode} =~ /{($ntype)}/g;
    my $linkupnode = int(rand(scalar @posslinknodes));
    $z = 0;
    while ($genome->{$mutnode} =~ /{([A-Z]+)}/) {
      my $subtype = $1;
      if ($subtype eq $ntype && $z++ == $linkupnode) { 
	# reconnect to original (copied) node
	$genome->{$mutnode} =~ s/{$ntype}/{$newnode}/;
      } else {
	# add new subtree of depth 0,1 or 2
	my $newdepth = poisson($self->{NewSubtreeDepthMean}, $self->{NewSubtreeDepthMax});
	$nodeid = $self->nid();
	$nodeid = $self->nid() while (defined $genome->{"node$subtype$nodeid"});
	my $subnode = "node$subtype$nodeid";
	$genome->{$mutnode} =~ s/{$subtype}/{$subnode}/;
	$genome->{$subnode} = '';
	$genome->{$subnode} =
	  $self->_grow_tree(depth=>0, TreeDepthMax=>$newdepth, type=>$subtype);
      }
    }
  }
  $self->{sizeCalcRequired} = 1;
}

sub delete_internal {
  my ($self) = @_;
  my $genome = $self->{genome};
  my $mutnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias});
  return unless ($mutnode);
  my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

  my $secondnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias},
			start_node=>$mutnode,
			node_type=>$ntype,
			not_this_node=>$mutnode);
  if ($secondnode) {
    my @ripnodes = $genome->{$mutnode} =~ /{(node[A-Z]+\d+)}/g;
    $genome->{$mutnode} = $genome->{$secondnode};
    $self->_xcopy_subtree($mutnode);
    grep { $self->_del_subtree($_) } @ripnodes;
    grep { delete $genome->{$_} } @ripnodes;
    $self->_fix_nodes($mutnode);
    $self->{mutednodes}{$mutnode} = 'delete_internal';
    $self->{mutednodes}{$secondnode} = 'delete_internal';
  }
  $self->{sizeCalcRequired} = 1;
}

sub copy_subtree {
  my ($self) = @_;
  my $genome = $self->{genome};
  my $mutnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias});
  return unless ($mutnode);
  my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

  # get another node that isn't in the subtree of mutnode
  my $secondnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias},
			node_type=>$ntype,
			not_this_subtree=>$mutnode);

  if ($secondnode) {
    # now check that mutnode isn't in the subtree of secondnode
    my %secondsubnodes;
    grep { $secondsubnodes{$_} = 1 } $self->_get_subnodes($secondnode);
    if (not exists $secondsubnodes{$mutnode}) {
      $self->_del_subtree($secondnode);
      $genome->{$secondnode} = $genome->{$mutnode};
      $self->_xcopy_subtree($secondnode);
      $self->_fix_nodes('root');
      $self->{mutednodes}{$mutnode} = 'copy_subtree';
      $self->{mutednodes}{$secondnode} = 'copy_subtree';
    }
  }
  $self->{sizeCalcRequired} = 1;
}

sub swap_subtrees {
  my ($self) = @_;
  my $genome = $self->{genome};
  my $mutnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias});
  return unless ($mutnode);
  my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

  # get another node that isn't in the subtree of mutnode
  my $secondnode =
    $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias},
			node_type=>$ntype,
			not_this_subtree=>$mutnode);

  if ($secondnode) {
    # now check that mutnode isn't in the subtree of secondnode
    my %secondsubnodes;
    grep { $secondsubnodes{$_} = 1 } $self->_get_subnodes($secondnode);
    if (not exists $secondsubnodes{$mutnode}) {
      # simply swap the contents of the two nodes and
      # everything will still be connected to something
      ($genome->{$mutnode}, $genome->{$secondnode}) =
	($genome->{$secondnode}, $genome->{$mutnode});

      $self->{mutednodes}{$mutnode} = 'swap_subtrees';
      $self->{mutednodes}{$secondnode} = 'swap_subtrees';
    }
  }
  $self->{sizeCalcRequired} = 1;
}

sub encapsulate_subtree {
  my ($self) = @_;
  my $genome = $self->{genome};
  my $genomesize = scalar keys %$genome;

  my $z = 0;
  while ($z++ < $too_many_tries) {
    my $mutnode =
      $self->_random_node(depth_bias=>$self->{MacroMutationDepthBias});
    return unless ($mutnode);
    my ($ntype) = $mutnode =~ qr/node([A-Z]+)\d+/;

    # some individuals don't allow certain subtrees to be frozen
    next if ($self->{EncapsulateIgnoreNTypes}{$ntype});

    my $subsize = $self->_tree_type_size($mutnode);
    # don't allow large fractions of the tree to be encapsulated
    next if ($subsize/$genomesize > $self->{EncapsulateFracMax});

    my $code = $self->_expand_tree($mutnode);
    $code = $self->simplify($code);

    # the length limit is a 'feature' of SDBM_File
    if (length($code) < 1000) {
      # replace subtree with the code that was the subtree
      $self->_del_subtree($mutnode);
      $genome->{$mutnode} = ";;$subsize;;$code";
      $self->{mutednodes}{$mutnode} = 'encapsulate_subtree';
      last; # we're done!
    }
  }
  $self->{sizeCalcRequired} = 1;
}

sub simplify {
  my ($self, $code) = @_;
  return $code;
}

# assume tied
sub _xcopy_subtree {
  my ($self, $node) = @_;
  my $genome = $self->{genome};

  my @subnodes = $genome->{$node} =~ /{(node[A-Z]+\d+)}/g;
  $genome->{$node} =~ s/{(node[A-Z]+\d+)}/{$1x}/g;
  foreach my $subnode (@subnodes) {
    $genome->{$subnode.'x'} = $genome->{$subnode};
    $self->_xcopy_subtree($subnode.'x');
  }
  $self->{sizeCalcRequired} = 1;
}

sub _get_subnodes_lookup {
  my ($self, $node) = @_;

  my @subnodes = @{$self->{nodesLookup}{$node}{subnodes}};
  my @retnodes = ();
  foreach my $subnode (@subnodes) {
      push @retnodes, $self->_get_subnodes_lookup($subnode);
  }
  return (@subnodes, @retnodes);
}

# assume tied
sub _get_subnodes_search {
  my ($self, $node) = @_;

  $self->{sizeCalcRequired} = 1 if (! exists $self->{sizeCalcRequired});
  $self->{nodesLookup} = {} if (! exists $self->{nodesLookup});

  my @subnodes = $self->{genome}{$node} =~ /{(node[A-Z]+\d+)}/g;
  my @retnodes = ();
  foreach my $subnode (@subnodes) {
      push @retnodes, $self->_get_subnodes_search($subnode);
  }
  return (@subnodes, @retnodes);
}

# assume tied
sub _get_subnodes {
  my ($self, $node) = @_;

  $self->{sizeCalcRequired} = 1 if (! exists $self->{sizeCalcRequired});
  $self->{nodesLookup} = {} if (! exists $self->{nodesLookup});

  my @subnodes;
  my @retnodes;
  if (($self->{sizeCalcRequired} == 1) || ! (exists $self->{nodesLookup}{$node})) {
    (@subnodes, @retnodes) = $self->_get_subnodes_search($node);
  }
  else {
    (@subnodes, @retnodes) = $self->_get_subnodes_lookup($node);
  }
  return (@subnodes, @retnodes);
}


# assume tied
sub _del_subtree {
  my ($self, $node) = @_;
  my $genome = $self->{genome};

  my @subnodes = $genome->{$node} =~ /{(node[A-Z]+\d+)}/g;
  if (@subnodes == 0) {
    delete $genome->{$node};
  }
  foreach my $subnode (@subnodes) {
    $self->_del_subtree($subnode);
    delete $genome->{$subnode};
  }
  $self->{sizeCalcRequired} = 1;
}

# assume tied
sub _tree_type_size { 
  my ($self, $node, $sizes, $types, $ignore, $stack) = @_;

  $self->{sizeCalcRequired} = 1 if (! exists $self->{sizeCalcRequired});
  $self->{nodesLookup} = {} if (! exists $self->{nodesLookup});

  my $nodetype;
  my $sum;
  if (($self->{sizeCalcRequired} == 1) || ! (exists $self->{nodesLookup}{$node})) {
    my @subnodes;
    my $genome = $self->{genome};

    $self->_tree_error($node, '_tree_type_size')
      if (!defined $genome->{$node});

    $nodetype = $node;
    $nodetype =~ s/node|\d+//g;
    @subnodes = $genome->{$node} =~ /{(node[A-Z]+\d+)}/g;
    if (@subnodes == 0) { # this is a leaf
      $sizes->{$node} = 1 if (defined $sizes);
      $types->{$node} = $nodetype if (defined $types);
      $self->{nodesLookup}{$node} = { sum=>1, nodetype=>$nodetype, subnodes=>[] };
      # is this an encapsulated subtree?  if so return the size of
      # the original - else 1
      return $genome->{$node}  =~ /^;;(\d+);;/ && $1 || 1;
    }
    # otherwise sum up the sizes of the subtrees
    # however if the $ignore hashref is defined then
    # don't recurse into subtree if the nodetype is in that hash
    $sum = 1;
    unless (defined $ignore && $ignore->{$nodetype}) {
      foreach my $subnode (@subnodes) {
        $sum += $self->_tree_type_size($subnode, $sizes, $types, $ignore);
      }
    }
    $self->{nodesLookup}{$node} = { sum=>$sum, nodetype=> $nodetype, subnodes=> [@subnodes] };
  }
  else
  {
    $sum = $self->{nodesLookup}{$node}{sum};
    $nodetype = $self->{nodesLookup}{$node}{nodetype};
  }
  $sizes->{$node} = $sum if (defined $sizes);
  $types->{$node} = $nodetype if (defined $types);
  $self->{sizeCalcRequired} = 0 if ($node eq 'root');
  return $sum;
}

sub _display_tree {
  my ($self, $node) = @_;

  printf STDERR "%-12s = '%s'\n", $node, $self->{genome}{$node};
  my @subnodes = $self->{genome}{$node} =~ /{(node[A-Z]+\d+x?)}/g;
  grep $self->_display_tree($_), @subnodes;
}


sub saveCode {
  my ($self, %p) = @_;
  die unless ($p{Filename});

  if (open(FILE, ">$p{Filename}")) {
    my $fitness = defined $self->Fitness() ? $self->Fitness() : 'unevaluated';
    my $age = $self->Age();
    my $code = $self->getCode(); # get this first, in case there is no code
    my $size = $self->getSize();

    printf FILE "# Experiment: $self->{ExperimentId}\n";
    if (defined $p{Tournament}) {
      printf FILE "# Tournament: %d\n", $p{Tournament};
    }
    print FILE "# Fitness:  $fitness\n";
    print FILE "# Age:      $age\n";
    print FILE "# CodeSize: $size\n";
    print FILE "# Code follows:\n$code\n";
    close(FILE);
  }
}

sub save {
  my ($self, %p) = @_;
  die unless ($p{FileStem});

  my %copygenome = ();
  tie %copygenome, $DBTYPE, $p{FileStem}, O_RDWR | O_CREAT, 0644;
  my $genome = $self->tieGenome('save');
  foreach my $key (keys %$genome) {
    $copygenome{$key} = $genome->{$key};
  }
  $self->untieGenome();
  $copygenome{'tied'} = 0;
  untie %copygenome;
}

sub load {
  my ($self, %p) = @_;
  die unless ($p{FileStem});

  my %copygenome;
  tie %copygenome, $DBTYPE, $p{FileStem}, O_RDONLY, 0644;
  my $genome = $self->tieGenome('load');
  %$genome = ();
  foreach my $key (keys %copygenome) {
    $genome->{$key} = $copygenome{$key};
  }
  $self->untieGenome();
  untie %copygenome;
}

1;

=head1 AUTHOR

Bob MacCallum, C<< <uncoolbob at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-perlgp-lite at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=PerlGP-Lite>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc PerlGP::Lite


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=PerlGP-Lite>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/PerlGP-Lite>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/PerlGP-Lite>

=item * Search CPAN

L<http://search.cpan.org/dist/PerlGP-Lite/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2013 Bob MacCallum.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

1; # End of PerlGP::Lite
