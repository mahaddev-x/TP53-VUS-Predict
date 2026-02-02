# TP53 Variant Structural Visualization
# PDB: 1TUP | Top 5 pathogenic variants (ESM-2 + AlphaMissense)
# Open in PyMOL: File > Run Script > tp53_variant_views.pml
# Or: pymol tp53_variant_views.pml

# ?? Load structure ????????????????????????????????????????????
fetch 1TUP, p53_dna, type=pdb
remove resn HOH
remove solvent

# ?? Base styling ??????????????????????????????????????????????
hide everything, all

# Protein: light gray cartoon
select protein_all, chain A+B+C and polymer.protein
show cartoon, protein_all
color gray80, protein_all
color gray90, chain B and polymer.protein

# DNA: orange
select dna_all, chain E+F and polymer.nucleic
show cartoon, dna_all
color orange, dna_all
set cartoon_ring_mode, 3
set cartoon_ring_color, orange, dna_all

# Zinc ions
select zinc_ions, resn ZN
show spheres, zinc_ions
color slate, zinc_ions
set sphere_scale, 0.5, zinc_ions

# Global render settings
set ray_opaque_background, 1
bg_color white
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 1
set spec_reflect, 0.3
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set label_font_id, 7
set label_size, 14
set label_color, black
set float_labels, 1

# ?? Variant 1: p.Leu257Arg (L257R) ??????????????????
# Hydrophobic core disruption
# ESM-2 LLR: -13.88 | AlphaMissense: 0.9938

hide sticks, all
hide labels, all
show cartoon, protein_all
show cartoon, dna_all
show spheres, zinc_ions

select wt_residue, chain B and resi 257
show sticks, wt_residue
color green, wt_residue
set stick_radius, 0.2, wt_residue

select nearby_res, (byres (all within 4.0 of wt_residue)) and polymer.protein and not wt_residue
show sticks, nearby_res
color lightblue, nearby_res
set stick_radius, 0.12, nearby_res

label wt_residue and name CA, "L257R"
label nearby_res and name CA, resn + resi

select focus, wt_residue or nearby_res
orient focus
zoom focus, 6
turn y, 15

scene Variant_L257R, store

ray 2400, 1800
png F:/SideProjects/MutationResearch/pymol_views/tp53_L257R_variant_view_pymol.png, dpi=300
print('Exported: tp53_L257R_variant_view_pymol.png')

delete dna_hbonds
delete zn_bonds
delete wt_residue
delete nearby_res
delete dna_contacts
delete focus

# ?? Variant 2: p.Val157Asp (V157D) ??????????????????
# Hydrophobic core disruption
# ESM-2 LLR: -13.69 | AlphaMissense: 0.9992

hide sticks, all
hide labels, all
show cartoon, protein_all
show cartoon, dna_all
show spheres, zinc_ions

select wt_residue, chain B and resi 157
show sticks, wt_residue
color green, wt_residue
set stick_radius, 0.2, wt_residue

select nearby_res, (byres (all within 4.0 of wt_residue)) and polymer.protein and not wt_residue
show sticks, nearby_res
color lightblue, nearby_res
set stick_radius, 0.12, nearby_res

label wt_residue and name CA, "V157D"
label nearby_res and name CA, resn + resi

select focus, wt_residue or nearby_res
orient focus
zoom focus, 6
turn y, 15

scene Variant_V157D, store

ray 2400, 1800
png F:/SideProjects/MutationResearch/pymol_views/tp53_V157D_variant_view_pymol.png, dpi=300
print('Exported: tp53_V157D_variant_view_pymol.png')

delete dna_hbonds
delete zn_bonds
delete wt_residue
delete nearby_res
delete dna_contacts
delete focus

# ?? Variant 3: p.Arg248Pro (R248P) ??????????????????
# DNA minor groove contact lost
# ESM-2 LLR: -12.63 | AlphaMissense: 0.9994

hide sticks, all
hide labels, all
show cartoon, protein_all
show cartoon, dna_all
show spheres, zinc_ions

select wt_residue, chain B and resi 248
show sticks, wt_residue
color green, wt_residue
set stick_radius, 0.2, wt_residue

select nearby_res, (byres (all within 4.0 of wt_residue)) and polymer.protein and not wt_residue
show sticks, nearby_res
color lightblue, nearby_res
set stick_radius, 0.12, nearby_res

select dna_contacts, (chain E+F and polymer.nucleic) within 4.0 of wt_residue
show sticks, dna_contacts
color brightorange, dna_contacts
set stick_radius, 0.15, dna_contacts
distance dna_hbonds, wt_residue, dna_contacts, 4.0, 2
set dash_color, red, dna_hbonds
set dash_width, 2.5, dna_hbonds
hide labels, dna_hbonds

label wt_residue and name CA, "R248P [DNA CONTACT]"
label nearby_res and name CA, resn + resi

select focus, wt_residue or nearby_res or dna_contacts
orient focus
zoom focus, 6
turn y, 15

scene Variant_R248P, store

ray 2400, 1800
png F:/SideProjects/MutationResearch/pymol_views/tp53_R248P_variant_view_pymol.png, dpi=300
print('Exported: tp53_R248P_variant_view_pymol.png')

delete dna_hbonds
delete zn_bonds
delete wt_residue
delete nearby_res
delete dna_contacts
delete focus

# ?? Variant 4: p.Cys176Arg (C176R) ??????????????????
# Zinc coordination abolished
# ESM-2 LLR: -12.53 | AlphaMissense: 0.9999

hide sticks, all
hide labels, all
show cartoon, protein_all
show cartoon, dna_all
show spheres, zinc_ions

select wt_residue, chain B and resi 176
show sticks, wt_residue
color green, wt_residue
set stick_radius, 0.2, wt_residue

select nearby_res, (byres (all within 4.0 of wt_residue)) and polymer.protein and not wt_residue
show sticks, nearby_res
color lightblue, nearby_res
set stick_radius, 0.12, nearby_res

distance zn_bonds, wt_residue, zinc_ions, 3.0, 2
set dash_color, slate, zn_bonds
set dash_width, 3.0, zn_bonds
hide labels, zn_bonds

label wt_residue and name CA, "C176R [Zn2+ SITE]"
label nearby_res and name CA, resn + resi

select focus, wt_residue or nearby_res
orient focus
zoom focus, 6
turn y, 15

scene Variant_C176R, store

ray 2400, 1800
png F:/SideProjects/MutationResearch/pymol_views/tp53_C176R_variant_view_pymol.png, dpi=300
print('Exported: tp53_C176R_variant_view_pymol.png')

delete dna_hbonds
delete zn_bonds
delete wt_residue
delete nearby_res
delete dna_contacts
delete focus

# ?? Variant 5: p.Arg280Ile (R280I) ??????????????????
# DNA major groove contact lost
# ESM-2 LLR: -12.47 | AlphaMissense: 0.9996

hide sticks, all
hide labels, all
show cartoon, protein_all
show cartoon, dna_all
show spheres, zinc_ions

select wt_residue, chain B and resi 280
show sticks, wt_residue
color green, wt_residue
set stick_radius, 0.2, wt_residue

select nearby_res, (byres (all within 4.0 of wt_residue)) and polymer.protein and not wt_residue
show sticks, nearby_res
color lightblue, nearby_res
set stick_radius, 0.12, nearby_res

select dna_contacts, (chain E+F and polymer.nucleic) within 4.0 of wt_residue
show sticks, dna_contacts
color brightorange, dna_contacts
set stick_radius, 0.15, dna_contacts
distance dna_hbonds, wt_residue, dna_contacts, 4.0, 2
set dash_color, red, dna_hbonds
set dash_width, 2.5, dna_hbonds
hide labels, dna_hbonds

label wt_residue and name CA, "R280I [DNA CONTACT]"
label nearby_res and name CA, resn + resi

select focus, wt_residue or nearby_res or dna_contacts
orient focus
zoom focus, 6
turn y, 15

scene Variant_R280I, store

ray 2400, 1800
png F:/SideProjects/MutationResearch/pymol_views/tp53_R280I_variant_view_pymol.png, dpi=300
print('Exported: tp53_R280I_variant_view_pymol.png')

delete dna_hbonds
delete zn_bonds
delete wt_residue
delete nearby_res
delete dna_contacts
delete focus

# ?? Overview: All 5 variant sites ?????????????????????????????
hide sticks, all
hide labels, all
show cartoon, protein_all
show cartoon, dna_all
show spheres, zinc_ions
select var_L257R, chain B and resi 257
show spheres, var_L257R
color red, var_L257R
set sphere_scale, 0.8, var_L257R
label var_L257R and name CA, "L257R"
select var_V157D, chain B and resi 157
show spheres, var_V157D
color magenta, var_V157D
set sphere_scale, 0.8, var_V157D
label var_V157D and name CA, "V157D"
select var_R248P, chain B and resi 248
show spheres, var_R248P
color yellow, var_R248P
set sphere_scale, 0.8, var_R248P
label var_R248P and name CA, "R248P"
select var_C176R, chain B and resi 176
show spheres, var_C176R
color cyan, var_C176R
set sphere_scale, 0.8, var_C176R
label var_C176R and name CA, "C176R"
select var_R280I, chain B and resi 280
show spheres, var_R280I
color orange, var_R280I
set sphere_scale, 0.8, var_R280I
label var_R280I and name CA, "R280I"
orient chain B or chain E+F
zoom chain B or chain E+F, 5
turn y, 30
scene Overview_All_Variants, store
ray 2400, 1800
png F:/SideProjects/MutationResearch/pymol_views/tp53_all_variants_overview_pymol.png, dpi=300

save F:/SideProjects/MutationResearch/pymol_views/tp53_variants.pse
print('Session saved: tp53_variants.pse')
print('Use Scene menu to switch between variant views')
