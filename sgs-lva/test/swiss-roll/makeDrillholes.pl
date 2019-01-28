# usage: perl makeDrillholes.pl label_file.csv
# with label_file.csv as:
#
# x,y,z,val
# 0.1,0.2,0.3,30.0
# ...
#

use Math::Vector::Real::kdTree;
use Math::Vector::Real;


$xmn=0.0;
$ymn=0.0;
$zmn=0.0;

$nx=40.0;
$ny=40.0;
$nz=40.0;

$xsiz=1.0/$nx;
$ysiz=1.0/$ny;
$zsiz=1.0/$nz;

$drillspacex=8;
$drillspacey=8;
$drillspacez=4;


open(IN,"<$ARGV[0]");
$header=<IN>;
print "x,y,z,val\n";
$const=3.14159265359/180;

%values=();
@v=();

while(<IN>){

	#$_=~/([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),(0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+)\n/;

	$_=~/([^,]+),([^,]+),([^,]+),([^,]+)/;

	$x=$1; $y=$2; $z=$3;
	$val=$4;

	$values{V($x,$y,$z)}="$val";
	
	push(@v,V($x,$y,$z));
}

my $tree = Math::Vector::Real::kdTree->new(@v);

#my $ix = $tree->find_nearest_vector(V(0.1020039273176382,-0.1962608937928392,0.4719450003066917));
#my $ix = $tree->find_nearest_vector(V(0.412527332681631,0.29568461144239,0.4928001779865698));
 
#print "nearest vector is $ix, $v[$ix], with values $values{$v[$ix]}\n";

for($i=0;$i<$nx;$i=$i+$drillspacex){
	$x=$xmn+$i*$xsiz;
for($j=0;$j<$ny;$j=$j+$drillspacey){
	$y=$ymn+$j*$ysiz;
for($k=0;$k<$nz;$k=$k+$drillspacez){
	$z=$zmn+$k*$zsiz;
	
	$ix = $tree->find_nearest_vector(V($x,$y,$z));
	print "$x,$y,$z,$values{$v[$ix]}";
	#print $v[$ix],",",$values{$v[$ix]},"\n";

}
}
}

