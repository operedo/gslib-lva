use POSIX;

$xmn=0.0;
$ymn=0.0;
$zmn=0.0;

$nx=40.0;
$ny=40.0;
$nz=40.0;

$xsiz=1.0/$nx;
$ysiz=1.0/$ny;
$zsiz=1.0/$nz;

$loc=1.0;

for($iz=0;$iz<$nz;$iz++){
for($iy=0;$iy<$ny;$iy++){
for($ix=0;$ix<$nx;$ix++){
	$x=$xmn+($ix)*$xsiz;
	$y=$ymn+($iy)*$ysiz;
	$z=$zmn+($iz)*$zsiz;
	$val=<>;
	chomp($val);
	print "$x,$y,$z,$val\n";
		
}
}
}


#
#while(<>){
#	$_=~/(.*)/;
#	$val=$1;
#	chomp($val);
#	$iz = 1.0 + floor( ($loc-1.0)/($nx*$ny) );
#	$iy = 1.0 + floor( (($loc-($iz-1.0)*$nx*$ny))/($nx) );
#	$ix = $loc - ($iz-1.0)*$nx*$ny - ($iy-1.0)*$nx -0.0;
#	
#	$x=$xmn+($ix)*$xsiz;
#	$y=$ymn+($iy)*$ysiz;
#	$z=$zmn+($iz)*$zsiz;
#
#	print "$ix,$iy,$iz : $x,$y,$z,$val\n";
#
#	$loc=$loc+1;
#}
