# usage: perl script.pl lva_file.csv
# with lva_file.csv as:
#
# x,y,z,a1,a2,a3,r1,r2
# 0.1,0.2,0.3,30.0,40.0,0.0,10.0,1.0
# ...
#
open(IN,"<$ARGV[0]");
$header=<IN>;
print "x,y,z,x1,x2,x3,y1,y2,y3,z1,z2,z3\n";
$const=3.14159265359/180;
while(<IN>){

	#$_=~/([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),(0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+)\n/;

	$_=~/([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),([^,]+),([^,]+)/;

	$x=$1; $y=$2; $z=$3;
	$a1=$4; $a2=$5; $a3=$6;
	$r1=$7; $r2=$8;

	$a2=0.0;
	$a3=0.0;

	#print $_;
	#print "$x,$y,$z,$a1,$a2,$a3,$r1,$r2\n";

	if($a1 >= 0.0 && $a1 < 270.0){
		$a1 = ( 90.0 - $a1) ;
	}
	else{
		$a1 = (450.0 - $a1) ;
	}
	$a2 = -1.0 * $a2;
	$a3 = 1.0 * $a3;
	$a=$a1*$const;
	$b=$a2*$const;
	$t=$a3*$const;

	$r1 = 1.0 / $r1;
	$r2 = 1.0 / $r2;

	$v100x=cos($b) * cos($a);
	$v010x=cos($b) * sin($a);
	$v001x=-sin($b);

	$v100y=$r1*(-cos($t)*sin($a)+sin($t)*sin($b)*cos($a));
	$v010y=$r1*( cos($t)*cos($a)+sin($t)*sin($b)*sin($a));
	$v001y=$r1*( sin($t)*cos($b));

	$v100z=$r2*( sin($t)*sin($a)+cos($t)*sin($b)*cos($a));
	$v010z=$r2*(-sin($t)*cos($a)+cos($t)*sin($b)*sin($a));
	$v001z=$r2*( cos($t)*cos($b));

	print "$x,$y,$z,$v100x,$v100y,$v100z,$v010x,$v010y,$v010z,$v001x,$v001y,$v001z\n";
}
