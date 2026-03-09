# Render MartiniSurf examples with VMD + TachyonInternal.
# Usage:
#   vmd -dispdev text -e render_examples_vmd.tcl -args <examples_dir> <output_dir>

proc pick_system_gro {example_dir} {
    set candidates [list \
        [file join $example_dir Simulation_Files 2_system final_system.gro] \
        [file join $example_dir Simulation_Files 2_system system_final.gro] \
        [file join $example_dir Simulation_Files 2_system system.gro] \
        [file join $example_dir Simulation_Files 2_system immobilized_system.gro] \
    ]
    foreach c $candidates {
        if {[file exists $c]} {
            return $c
        }
    }
    return ""
}

proc has_atoms {mol selection} {
    set sel [atomselect $mol $selection]
    set n [$sel num]
    $sel delete
    return $n
}

proc render_view {ex_name out_dir view_tag view_setup_cmd} {
    # Apply camera transform for this view.
    display resetview
    eval $view_setup_cmd

    set tga_out [file join $out_dir "${ex_name}_${view_tag}.tga"]
    set png_out [file join $out_dir "${ex_name}_${view_tag}.png"]

    render TachyonInternal $tga_out

    # Convert to PNG when available.
    if {[catch {exec magick convert $tga_out $png_out} _] != 0} {
        if {[catch {exec convert $tga_out $png_out} _] != 0} {
            puts "$ex_name ($view_tag): rendered $tga_out (PNG conversion tool not found)"
            return ""
        } else {
            file delete -force $tga_out
            puts "$ex_name ($view_tag): rendered $png_out"
            return $png_out
        }
    } else {
        file delete -force $tga_out
        puts "$ex_name ($view_tag): rendered $png_out"
        return $png_out
    }
}

proc anchor_index_selection {example_dir} {
    set ndx_path [file join $example_dir Simulation_Files 0_topology index.ndx]
    if {![file exists $ndx_path]} {
        return ""
    }

    set fh [open $ndx_path r]
    set in_anchor_group 0
    set atom_ids {}

    while {[gets $fh line] >= 0} {
        set s [string trim $line]
        if {$s eq ""} {
            continue
        }
        if {[regexp {^\[\s*([^\]]+)\s*\]$} $s -> group_name]} {
            if {[string match "Anchor_*" $group_name]} {
                set in_anchor_group 1
            } else {
                set in_anchor_group 0
            }
            continue
        }
        if {$in_anchor_group} {
            foreach tok [split $s] {
                if {[string is integer -strict $tok]} {
                    # GROMACS index atoms are 1-based; VMD "index" is 0-based.
                    lappend atom_ids [expr {$tok - 1}]
                }
            }
        }
    }
    close $fh

    if {[llength $atom_ids] == 0} {
        return ""
    }
    return "index [join [lsort -integer -unique $atom_ids] { }]"
}

proc render_example {example_dir out_dir} {
    set ex_name [file tail $example_dir]
    set gro [pick_system_gro $example_dir]
    if {$gro eq ""} {
        puts "$ex_name: skipped (no system GRO found)"
        return
    }

    mol new $gro type gro waitfor all
    set m [molinfo top]

    # Keep complex components contiguous across periodic boundaries.
    if {[catch {package require pbctools} _] == 0} {
        set center_sel_text "resname ALK EPOX MOL1 DA DC DG DT NAD ETO or (name BB and resname ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL)"
        if {[has_atoms $m $center_sel_text] == 0} {
            set center_sel_text "all"
        }
        pbc unwrap -all -molid $m
        pbc wrap -all -molid $m -compound residue -center com -centersel "$center_sel_text"
    }

    # Global display/ray settings.
    display projection Orthographic
    display shadows on
    display ambientocclusion on
    display aoambient 0.30
    display aodirect 0.45
    display depthcue off
    display resize 2800 2200
    color Display Background white
    color change rgb 31 1.00 0.12 0.12
    color change rgb 29 1.00 0.90 0.00
    color change rgb 22 0.28 0.56 0.86
    color change rgb 23 0.36 0.85 0.93
    color change rgb 24 0.95 0.62 0.80
    color change rgb 25 0.28 0.70 0.72
    color change rgb 26 1.00 0.45 0.00
    color change rgb 27 0.05 0.20 0.60
    color change rgb 28 1.00 1.00 1.00
    axes location Off

    # Reset default reps.
    mol delrep 0 $m

    # Exclusions shared by all selections.
    set exclude_solvent_ions "resname W WF SOL NA CL K CA MG ZN LI RB CS BA SR F BR I"
    set protein_res "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"
    # Anchor residues from example 05 to highlight across protein renders.
    set protein_anchor_resids "8 10 11 1025 1027 1028"
    set anchor_idx_sel [anchor_index_selection $example_dir]
    if {$anchor_idx_sel ne ""} {
        set anchor_core_sel "(same residue as ($anchor_idx_sel)) and (resname $protein_res or resname DA DC DG DT)"
    } else {
        set anchor_core_sel "resname $protein_res and resid $protein_anchor_resids"
    }

    # Rep 1: protein backbone beads only (BB), uniform white.
    set protein_bb_sel "name BB and resname $protein_res and not ($exclude_solvent_ions) and not ($anchor_core_sel)"
    if {[has_atoms $m $protein_bb_sel] > 0} {
        mol representation VDW 1.4 24.0
        mol color ColorID 28
        mol selection $protein_bb_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Rep 1b: highlight anchor protein residues in red, BB bead only.
    set protein_anchor_sel "($anchor_core_sel) and resname $protein_res and name BB and not ($exclude_solvent_ions)"
    if {[has_atoms $m $protein_anchor_sel] > 0} {
        mol representation VDW 1.4 24.0
        mol color ColorID 1
        mol selection $protein_anchor_sel
        mol material Glossy
        mol addrep $m
    }

    # Rep 2: surface C-type beads in orange (supports C and C1 naming).
    set surf_c1_sel "(resname SRF GRA) and name C C1 and not ($exclude_solvent_ions)"
    if {[has_atoms $m $surf_c1_sel] > 0} {
        mol representation VDW 1.0 20.0
        mol color ColorID 3
        mol selection $surf_c1_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Rep 3: surface P4 beads in orange.
    set surf_p4_sel "(resname SRF GRA) and name P4 and not ($exclude_solvent_ions)"
    if {[has_atoms $m $surf_p4_sel] > 0} {
        mol representation VDW 1.4 20.0
        mol color ColorID 3
        mol selection $surf_p4_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Rep 4: linkers in the same orange as the surface.
    set linker_sel "resname ALK EPOX MOL1 and not ($exclude_solvent_ions)"
    if {[has_atoms $m $linker_sel] > 0} {
        mol representation VDW 1.0 24.0
        mol color ColorID 3
        mol selection $linker_sel
        mol material Glossy
        mol addrep $m
    }

    # Rep 5: DNA backbone beads in blue.
    set dna_bb_sel "resname DA DC DG DT and name BB1 BB2 BB3 and not ($exclude_solvent_ions)"
    if {[has_atoms $m $dna_bb_sel] > 0} {
        mol representation VDW 1.0 20.0
        mol color ColorID 22
        mol selection $dna_bb_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Rep 5b: DNA side beads in turquoise.
    set dna_sc_sel "resname DA DC DG DT and name SC1 SC2 SC3 SC4 and not ($exclude_solvent_ions)"
    if {[has_atoms $m $dna_sc_sel] > 0} {
        mol representation VDW 1.0 20.0
        mol color ColorID 23
        mol selection $dna_sc_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Rep 6: ETO substrate in blue.
    set eto_sel "resname ETO and not ($exclude_solvent_ions)"
    if {[has_atoms $m $eto_sel] > 0} {
        mol representation VDW 1.4 20.0
        mol color ColorID 0
        mol selection $eto_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Rep 7: cofactors/substrates/other molecules (excluding surface/linker/protein/DNA/ETO/ions/water), color by name.
    set cof_sub_sel "(not resname $protein_res) and not resname DA DC DG DT SRF GRA ALK EPOX MOL1 ETO W WF SOL NA CL K CA MG ZN LI RB CS BA SR F BR I"
    if {[has_atoms $m $cof_sub_sel] > 0} {
        mol representation VDW 1.4 20.0
        mol color Name
        mol selection $cof_sub_sel
        mol material AOEdgy
        mol addrep $m
    }

    # Fallback: if no specific reps matched, render all non-water/non-ion atoms by Name.
    if {[molinfo $m get numreps] == 0} {
        set fallback_sel "not ($exclude_solvent_ions)"
        if {[has_atoms $m $fallback_sel] > 0} {
            mol representation VDW 1.0 20.0
            mol color Name
            mol selection $fallback_sel
            mol material AOEdgy
            mol addrep $m
        }
    }

    # Top view: keep Z up and Y to the right.
    set top_png [render_view $ex_name $out_dir "top" {rotate z by -90; scale by 1.8}]
    # Side view: Z positive up and X positive to the right.
    if {$ex_name eq "07_protein_nad"} {
        # For figure 07, flip in-plane X direction while keeping +Z up.
        set side_png [render_view $ex_name $out_dir "side" {rotate x by -90; rotate z by 180; scale by 1.8}]
    } else {
        set side_png [render_view $ex_name $out_dir "side" {rotate x by -90; scale by 1.8}]
    }

    mol delete $m
}

# Args
set examples_dir [pwd]
set out_dir [file join $examples_dir images]
if {[llength $argv] >= 1} {
    set examples_dir [lindex $argv 0]
}
if {[llength $argv] >= 2} {
    set out_dir [lindex $argv 1]
}

file mkdir $out_dir

set dirs [lsort [glob -nocomplain -directory $examples_dir -types d {[0-9][0-9]_*}]]

if {[llength $dirs] == 0} {
    puts "No numbered example directories found under $examples_dir"
    quit
}

foreach d $dirs {
    render_example $d $out_dir
}

quit
