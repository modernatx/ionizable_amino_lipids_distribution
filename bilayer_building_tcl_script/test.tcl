
#  Slurp up the data file
set fp [open "../output/gsmd.0.xsc" r]
set file_data [read $fp]
set data [split $file_data "\n"]
set r1 [lindex $data 2 1]
set r2 [lindex $data 2 5]
set r3 [lindex $data 2 9]
#foreach line $data {
#     puts $line
#}
close $fp
