set sel [atomselect top "all"]
set com [measure center $sel weight mass]
set matrix [trans x 90]
$sel moveby [vecscale -1.0 $com]
$sel move $matrix
$sel moveby $com 