use 5.006;
use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'PerlGP::Lite',
    AUTHOR           => q{Bob MacCallum <uncoolbob@gmail.com>},
    VERSION_FROM     => 'lib/PerlGP/Lite.pm',
    ABSTRACT_FROM    => 'lib/PerlGP/Lite.pm',
    LICENSE          => 'Artistic_2_0',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
	'Mouse' => 0,
	'aliased' => 0,
	'SDBM_File' => 0,
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'PerlGP-Lite-*' },
);
