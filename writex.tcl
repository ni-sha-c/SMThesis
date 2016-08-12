mol new poiseuilleD40.lammpstrj
pbc wrap -all
for {set i 10} {$i <= 40} {incr i} {
set liquid [atomselect top "(type==5) || (type==6) && (z < 30.0)" frame $i]
set xliquid [$liquid get x]
set filename "xcoord$i.dat"
set fileid [open $filename "w"]
puts -nonewline $fileid $xliquid
close $fileid
$liquid delete
}

