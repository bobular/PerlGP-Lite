package GPMisc;
#GPL#

sub pickrandom {
    my $ref = shift;
    return $ref->[int rand @$ref];
}

sub getTime {
    my @a = times();
    return $a[0];
}

sub poisson
{
  my ($mean, $cap) = @_;
  my $prob = 1/$mean;
  my $i=0;
  while (rand()>$prob)
  {
    $i++;
  }
  return $cap if ($cap && $i>$cap);
  return $i;
}

sub gaussian {
  my ($x, $m, $s) = @_;
  $m = 0 unless (defined $m);
  $s = 1 unless (defined $s);

  return exp( -1 * ($x-$m)*($x-$m)/(2*$s*$s) ) / ( $s * 2.506628274631 );
  # 2.506628274631 = sqrt(2*pi)
}

# returns random numbers between -$r and +$r
# from a gaussian distribution with stdev $s
sub pseudo_gaussian_rand {
  my ($r, $s) = @_;
  $s = 1 unless ($s);
  my $x;
  while (1) {
    $x = rand($r*2)-$r;
    if (rand() < gaussian($x, 0, $s)) {
      return $x;
    }
  }
}

sub pid_filename {
  my ($scriptname) = $0 =~ m|([^/]+)$|;
  $scriptname =~ s/\.\w+$//;
  return "/tmp/${scriptname}_$ENV{USER}.pid";
}

sub not_running {

  my $pidfile = pid_filename();

  if (-e $pidfile) {
    open(PIDLOCK, $pidfile) || die "can't open pidfile $pidfile";
    my $pid = <PIDLOCK>;
    close(PIDLOCK);
    chomp($pid);
    return 0 if (kill 0, $pid);
  }
  open(PIDLOCK, ">$pidfile") || die "can't write $pidfile";
  print PIDLOCK "$$\n";
  close PIDLOCK;

  return 1;
}

sub unique
{
  my %got = ();
  my @got = ();
  my $thing;
  foreach $thing (@_)
  {
    push @got, $thing if (not defined $got{$thing});
    $got{$thing} = 1;
  }
  @got;
}

return 1;
