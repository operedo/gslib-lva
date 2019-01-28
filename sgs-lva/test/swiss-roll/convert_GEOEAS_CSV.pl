$xmn=0.0;
$ymn=0.0;
$zmn=0.0;

$nx=40.0;
$ny=40.0;
$nz=40.0;

$xsiz=1.0/$nx;
$ysiz=1.0/$ny;
$zsiz=1.0/$nz;

$loc=1;
while(<>){
	$_=~/(.*)/;
	$val=$1;
	chomp($val);
	$iz = 1 + int( ($loc-1)/($nx*$ny) );
	$iy = 1 + int( ($loc-($iz-1)*$nx*$ny)/$nx );
	$ix = $loc - ($iz-1)*$nx*$ny - ($iy-1)*$nx;
	
	$x=$xmn+$ix*$xsiz;
	$y=$ymn+$iy*$ysiz;
	$z=$zmn+$ik*$zsiz;

	print "$x,$y,$z,$val\n";

	$loc=$loc+1;
}
