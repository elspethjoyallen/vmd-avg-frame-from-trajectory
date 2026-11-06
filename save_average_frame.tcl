# Save the trajectory frame closest to the aligned average structure.
#
# Expected input files in the current working directory:
#   <system>_md.prmtop
#   <system>_step<step><replicate>.nc
#
# For example, system beta2_tyrosine, steps 7-10, and replicates a-d use:
#   beta2_tyrosine_md.prmtop
#   beta2_tyrosine_step7a.nc ... beta2_tyrosine_step10d.nc

proc prompt {label} {
    puts -nonewline "$label: "
    flush stdout
    if {[gets stdin value] < 0} {
        error "Could not read input for '$label'. Run VMD from a terminal so the prompts can read your responses."
    }
    return [string trim $value]
}

proc parse_step_range {value} {
    if {[regexp {^([0-9]+)$} $value -> first]} {
        return [list $first $first]
    }
    if {![regexp {^([0-9]+)-([0-9]+)$} $value -> first last]} {
        error "Enter steps as one number or an inclusive range, for example 7-10."
    }
    if {$first > $last} {
        error "The first step must be less than or equal to the last step."
    }
    return [list $first $last]
}

proc parse_replicate_range {value} {
    if {![regexp {^([A-Za-z])(?:-([A-Za-z]))?$} $value -> first last]} {
        error "Enter replicates as one letter or an inclusive letter range, for example a-d."
    }
    set first [string tolower $first]
    if {$last eq ""} {
        set last $first
    } else {
        set last [string tolower $last]
    }
    scan $first %c first_code
    scan $last %c last_code
    if {$first_code < 97 || $last_code > 122 || $first_code > $last_code} {
        error "Replicates must be an ascending range from a to z, such as a-d."
    }
    set replicates {}
    for {set code $first_code} {$code <= $last_code} {incr code} {
        lappend replicates [format %c $code]
    }
    return $replicates
}

proc parse_residue_ranges {value} {
    set tokens [regexp -all -inline {\S+} $value]
    if {[llength $tokens] == 0} {
        error "Enter at least one residue number or range."
    }

    set normalized {}
    set i 0
    set token_count [llength $tokens]
    while {$i < $token_count} {
        set first [lindex $tokens $i]
        if {![regexp {^-?[0-9]+$} $first]} {
            error "Residue selections must use numbers and the word 'to', for example: 37 to 59 66 to 90."
        }
        lappend normalized $first
        incr i

        if {$i < $token_count && [string tolower [lindex $tokens $i]] eq "to"} {
            incr i
            if {$i >= $token_count} {
                error "A residue range ending in 'to' needs an ending residue number."
            }
            set last [lindex $tokens $i]
            if {![regexp {^-?[0-9]+$} $last]} {
                error "Residue selections must use numbers and the word 'to', for example: 37 to 59 66 to 90."
            }
            if {$first > $last} {
                error "Each residue range must be in ascending order."
            }
            lappend normalized to $last
            incr i
        }
    }
    return [join $normalized " "]
}

puts "Save the trajectory frame closest to the aligned average structure"
puts "Input files must be in the current working directory."

set system [prompt "System name (filename prefix)"]
if {$system eq "" || ![regexp {^[A-Za-z0-9_.-]+$} $system] || $system eq "." || $system eq ".."} {
    error "Use a non-empty system name containing only letters, numbers, periods, underscores, or hyphens."
}

set step_input [prompt "Steps to include (for example, 7-10)"]
lassign [parse_step_range $step_input] first_step last_step

set replicate_input [prompt "Replicates to include (for example, a-d)"]
set replicates [parse_replicate_range $replicate_input]

set stride_input [prompt "Trajectory stride (positive integer; for example, 100)"]
if {![regexp {^[0-9]+$} $stride_input] || $stride_input < 1} {
    error "Stride must be a positive integer."
}
set stride $stride_input

puts "For alignment, prefer the most stable protein regions (for GPCRs, the transmembrane domains are a common choice). Avoid flexible loops or tails where possible."
set residue_input [prompt "Residue ranges (for example, 37 to 59 66 to 90)"]
set residue_ranges [parse_residue_ranges $residue_input]
set fit_sel_text "backbone and (resid $residue_ranges)"
set rms_sel_text $fit_sel_text

set topology_file "${system}_md.prmtop"
if {![file isfile $topology_file]} {
    error "Topology file not found: [file normalize $topology_file]"
}

set trajectory_files {}
for {set step $first_step} {$step <= $last_step} {incr step} {
    foreach replicate $replicates {
        set trajectory_file [format "%s_step%d%s.nc" $system $step $replicate]
        if {![file isfile $trajectory_file]} {
            error "Trajectory file not found: [file normalize $trajectory_file]"
        }
        lappend trajectory_files $trajectory_file
    }
}

set output_file "${system}_avg_frame.pdb"
if {[file exists $output_file]} {
    set overwrite [string tolower [prompt "Output $output_file already exists. Overwrite? (y/n)"]]
    if {$overwrite ni {y yes}} {
        error "Stopped without overwriting the existing output file."
    }
}

puts "Loading topology: $topology_file"
set molid [mol new $topology_file type parm7 waitfor all]
foreach trajectory_file $trajectory_files {
    puts "Loading trajectory: $trajectory_file (stride $stride)"
    if {[catch {
        mol addfile $trajectory_file type netcdf first 0 last -1 step $stride waitfor all molid $molid
    } message]} {
        mol delete $molid
        error "Could not load $trajectory_file: $message"
    }
}

set nframes [molinfo $molid get numframes]
puts "Total frames loaded: $nframes"
if {$nframes < 1} {
    mol delete $molid
    error "No trajectory frames loaded."
}

set ref_fit [atomselect $molid $fit_sel_text frame 0]
set nfit [$ref_fit num]
if {$nfit < 1} {
    $ref_fit delete
    mol delete $molid
    error "The alignment selection matched no atoms: $fit_sel_text"
}

puts "Aligning all frames to frame 0 using: $fit_sel_text"
for {set i 0} {$i < $nframes} {incr i} {
    set mob_fit [atomselect $molid $fit_sel_text frame $i]
    set all_atoms [atomselect $molid "all" frame $i]
    if {[$mob_fit num] != $nfit} {
        $mob_fit delete
        $all_atoms delete
        $ref_fit delete
        mol delete $molid
        error "The alignment selection has inconsistent atom counts at frame $i."
    }
    set trans_mat [measure fit $mob_fit $ref_fit]
    $all_atoms move $trans_mat
    $mob_fit delete
    $all_atoms delete
}
$ref_fit delete

puts "Computing the average coordinates..."
set avg_sel [atomselect $molid $rms_sel_text frame 0]
set natoms [$avg_sel num]
set sum_coords {}
for {set j 0} {$j < $natoms} {incr j} {
    lappend sum_coords [list 0.0 0.0 0.0]
}

for {set i 0} {$i < $nframes} {incr i} {
    set sel [atomselect $molid $rms_sel_text frame $i]
    set coords [$sel get {x y z}]
    $sel delete
    for {set j 0} {$j < $natoms} {incr j} {
        set oldsum [lindex $sum_coords $j]
        set xyz [lindex $coords $j]
        lset sum_coords $j [list \
            [expr {[lindex $oldsum 0] + [lindex $xyz 0]}] \
            [expr {[lindex $oldsum 1] + [lindex $xyz 1]}] \
            [expr {[lindex $oldsum 2] + [lindex $xyz 2]}]]
    }
}

set avg_coords {}
for {set j 0} {$j < $natoms} {incr j} {
    set sum [lindex $sum_coords $j]
    lappend avg_coords [list \
        [expr {[lindex $sum 0] / double($nframes)}] \
        [expr {[lindex $sum 1] / double($nframes)}] \
            [expr {[lindex $sum 2] / double($nframes)}]]
}
$avg_sel delete

puts "Finding the loaded frame closest to the average..."
set best_frame -1
set best_rmsd 1.0e300
for {set i 0} {$i < $nframes} {incr i} {
    set sel [atomselect $molid $rms_sel_text frame $i]
    set coords [$sel get {x y z}]
    $sel delete
    set sum_squared_displacement 0.0
    for {set j 0} {$j < $natoms} {incr j} {
        set xyz [lindex $coords $j]
        set mean_xyz [lindex $avg_coords $j]
        set dx [expr {[lindex $xyz 0] - [lindex $mean_xyz 0]}]
        set dy [expr {[lindex $xyz 1] - [lindex $mean_xyz 1]}]
        set dz [expr {[lindex $xyz 2] - [lindex $mean_xyz 2]}]
        set sum_squared_displacement [expr {$sum_squared_displacement + $dx*$dx + $dy*$dy + $dz*$dz}]
    }
    set rmsd [expr {sqrt($sum_squared_displacement / double($natoms))}]
    if {$rmsd < $best_rmsd} {
        set best_rmsd $rmsd
        set best_frame $i
    }
}

set output_sel [atomselect $molid all frame $best_frame]
$output_sel writepdb $output_file
$output_sel delete
mol delete $molid

puts "Saved frame $best_frame (0-based VMD frame index; RMSD $best_rmsd angstrom) to $output_file"
