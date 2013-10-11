my $r=new Ring(qw(a b c d e f));
my ($a,$b,$c,$d,$e,$f)=$r->variables;
my $p1 = $f^3;
my $p2 = $e*($f)^2;
my $p3 = ($e^2)*$f;
my $p4 = $b*$c*$f - $a*$d*$f;
my $p5 = $d*$e + $c*$f;
my $p6 = $b*$e + $a*$f;
my $p7 = $e^3;
my $i = new Ideal(GENERATORS=>[$p1,$p2,$p3,$p4,$p5,$p6,$p7]);
my @pd = $i->PRIMARY_DECOMPOSITION;
compare_object('1', $pd[0]) && compare_object('2', $pd[1]);


