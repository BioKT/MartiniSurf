# Render first 5 frames from DNA_surface.pdb with VMD + TachyonInternal.
# Usage:
#   vmd -dispdev text -e render_dna_surface_frames_vmd.tcl -args <pdb_path> <out_dir>

proc has_atoms {mol selection} {
    set sel [atomselect $mol $selection]
    set n [$sel num]
    $sel delete
    return $n
}

proc render_one {out_base} {
    set tga_out "${out_base}.tga"
    set png_out "${out_base}.png"
    if {[catch {render TachyonInternal $tga_out} err]} {
        puts "render failed: $err"
        return ""
    }

    if {[catch {exec magick convert $tga_out $png_out} _] != 0} {
        if {[catch {exec convert $tga_out $png_out} _] != 0} {
            puts "rendered $tga_out (PNG conversion tool not found)"
            return $tga_out
        } else {
            file delete -force $tga_out
            puts "rendered $png_out"
            return $png_out
        }
    } else {
        file delete -force $tga_out
        puts "rendered $png_out"
        return $png_out
    }
}

if {[llength $argv] < 1} {
    puts "Usage: vmd -dispdev text -e render_dna_surface_frames_vmd.tcl -args <pdb_path> <out_dir>"
    quit
}

set pdb_path [lindex $argv 0]
set out_dir [file join [pwd] "dna_surface_frames"]
if {[llength $argv] >= 2} {
    set out_dir [lindex $argv 1]
}
file mkdir $out_dir

if {![file exists $pdb_path]} {
    puts "PDB not found: $pdb_path"
    quit
}

mol new $pdb_path type pdb waitfor all
set m [molinfo top]
set nframes [molinfo $m get numframes]
set nrender [expr {$nframes < 5 ? $nframes : 5}]
puts "Loaded $pdb_path with $nframes frames; rendering $nrender"

# Display settings (same style as example renders)
display projection Orthographic
display shadows on
display ambientocclusion on
display aoambient 0.30
display aodirect 0.45
display depthcue off
display resize 2800 2200
color Display Background white
axes location Off

mol delrep 0 $m

set exclude_solvent_ions "resname W WF SOL NA CL K CA MG ZN LI RB CS BA SR F BR I"

# Positive ions in blue
set ion_pos_sel "(resname NA or (resname ION and name NA+))"
if {[has_atoms $m $ion_pos_sel] > 0} {
    mol representation VDW 1.0 18.0
    mol color ColorID 0
    mol selection $ion_pos_sel
    mol material AOEdgy
    mol addrep $m
}

# Negative ions in red
set ion_neg_sel "(resname CL or (resname ION and name CL-))"
if {[has_atoms $m $ion_neg_sel] > 0} {
    mol representation VDW 1.0 18.0
    mol color ColorID 1
    mol selection $ion_neg_sel
    mol material AOEdgy
    mol addrep $m
}

# Surface C/C1 in gray
set surf_c_sel "(resname SRF GRA) and name C C1 and not ($exclude_solvent_ions)"
if {[has_atoms $m $surf_c_sel] > 0} {
    mol representation VDW 1.0 20.0
    mol color ColorID 2
    mol selection $surf_c_sel
    mol material AOEdgy
    mol addrep $m
}

# Surface P4 in orange
set surf_p4_sel "(resname SRF GRA) and name P4 and not ($exclude_solvent_ions)"
if {[has_atoms $m $surf_p4_sel] > 0} {
    mol representation VDW 1.4 20.0
    mol color ColorID 3
    mol selection $surf_p4_sel
    mol material AOEdgy
    mol addrep $m
}

# Linker in yellow
set linker_sel "resname ALK EPOX MOL1 and not ($exclude_solvent_ions)"
if {[has_atoms $m $linker_sel] > 0} {
    mol representation VDW 1.2 24.0
    mol color ColorID 0
    mol selection $linker_sel
    mol material Opaque
    mol addrep $m
}

# DNA in Name coloring
set dna_sel "resname DA DC DG DT and not ($exclude_solvent_ions)"
if {[has_atoms $m $dna_sel] > 0} {
    mol representation VDW 1.2 20.0
    mol color Name
    mol selection $dna_sel
    mol material AOEdgy
    mol addrep $m
}

for {set i 0} {$i < $nrender} {incr i} {
    animate goto $i

    # Top view
    display resetview
    rotate z by -90
    scale by 1.8
    set top_base [file join $out_dir [format "DNA_surface_frame%02d_top" [expr {$i+1}]]]
    render_one $top_base

    # Side view
    display resetview
    rotate y by 90
    scale by 1.8
    set side_base [file join $out_dir [format "DNA_surface_frame%02d_side" [expr {$i+1}]]]
    render_one $side_base

    puts "done frame [expr {$i+1}]"
}

mol delete $m
puts "All frame renders completed in $out_dir"
quit
