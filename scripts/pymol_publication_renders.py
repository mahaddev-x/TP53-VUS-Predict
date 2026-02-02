#!/usr/bin/env python3
"""
TP53 Publication-Quality Structural Renders (PyMOL cmd API)
===========================================================
Generates high-DPI publication renders for the top 5 TP53 mutations
using PyMOL's headless rendering mode.

Run with:
    micromamba run -n pymol pymol -cqr scripts/pymol_publication_renders.py

PDB 1TUP: Cho et al. (1994) Science 265:346-355
  - Chain B: p53 core domain (DNA-bound monomer)
  - Chains E, F: 21-bp DNA duplex
  - Zn2+: structural zinc ion
"""

import os
import sys
from pathlib import Path

from pymol import cmd

# ═════════════════════════════════════════════════════════════════════════════
# Configuration
# ═════════════════════════════════════════════════════════════════════════════

# PyMOL -cqr may resolve __file__ relative to its own install directory,
# so we derive the project root from the current working directory instead.
# The run command should always be invoked from the project root:
#   cd F:\SideProjects\MutationResearch && pymol -cqr scripts/...
PROJECT_DIR = Path(os.getcwd())
FIGURES_DIR = PROJECT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

PDB_ID = "1TUP"
PROTEIN_CHAIN = "B"
DNA_CHAINS = "E+F"
CONTACT_DISTANCE = 4.0

VARIANTS = [
    {
        "name": "L257R",
        "resi": 257,
        "wt_one": "L",
        "mut_one": "R",
        "wt_three": "LEU",
        "mut_three": "ARG",
        "mechanism": "Hydrophobic core disruption",
        "detail": "Buried Leu in beta-sandwich; Arg introduces charge into core",
    },
    {
        "name": "V157D",
        "resi": 157,
        "wt_one": "V",
        "mut_one": "D",
        "wt_three": "VAL",
        "mut_three": "ASP",
        "mechanism": "Hydrophobic core disruption",
        "detail": "Beta-strand in sandwich core; Val->Asp introduces charge",
    },
    {
        "name": "R248P",
        "resi": 248,
        "wt_one": "R",
        "mut_one": "P",
        "wt_three": "ARG",
        "mut_three": "PRO",
        "mechanism": "DNA minor groove contact lost",
        "detail": "R248 inserts into DNA minor groove; Pro eliminates contact + kinks backbone",
    },
    {
        "name": "C176R",
        "resi": 176,
        "wt_one": "C",
        "mut_one": "R",
        "wt_three": "CYS",
        "mut_three": "ARG",
        "mechanism": "Zinc coordination abolished",
        "detail": "C176 chelates structural Zn2+; Arg cannot coordinate zinc",
    },
    {
        "name": "R280I",
        "resi": 280,
        "wt_one": "R",
        "mut_one": "I",
        "wt_three": "ARG",
        "mut_three": "ILE",
        "mechanism": "DNA major groove contact lost",
        "detail": "R280 H-bonds DNA major groove base; Ile is hydrophobic",
    },
]

# Mapping from 1-letter to 3-letter codes for the mutagenesis wizard
ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


# ═════════════════════════════════════════════════════════════════════════════
# Scene setup
# ═════════════════════════════════════════════════════════════════════════════

def setup_scene():
    """Fetch 1TUP, remove water, apply base styling and render settings."""
    print("[1/3] Setting up scene...")

    cmd.reinitialize()
    cmd.fetch(PDB_ID, name="p53_dna", type="pdb")
    cmd.remove("resn HOH")
    cmd.remove("solvent")

    # Base representation
    cmd.hide("everything", "all")

    # Protein: gray cartoon
    cmd.select("protein_all", "chain A+B+C and polymer.protein")
    cmd.show("cartoon", "protein_all")
    cmd.color("gray80", "protein_all")
    cmd.color("gray90", f"chain {PROTEIN_CHAIN} and polymer.protein")

    # DNA: orange cartoon with ring fill
    cmd.select("dna_all", f"chain {DNA_CHAINS} and polymer.nucleic")
    cmd.show("cartoon", "dna_all")
    cmd.color("orange", "dna_all")
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_color", "orange", "dna_all")

    # Zinc ions: slate spheres
    cmd.select("zinc_ions", "resn ZN")
    cmd.show("spheres", "zinc_ions")
    cmd.color("slate", "zinc_ions")
    cmd.set("sphere_scale", 0.5, "zinc_ions")

    # Publication render settings
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 4)
    cmd.set("ray_trace_mode", 1)
    cmd.set("spec_reflect", 0.3)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("label_font_id", 7)
    cmd.set("label_size", 14)
    cmd.set("label_color", "black")
    cmd.set("float_labels", 1)

    print("  Scene ready: protein (gray cartoon), DNA (orange), Zn (slate spheres)")


# ═════════════════════════════════════════════════════════════════════════════
# Mutagenesis
# ═════════════════════════════════════════════════════════════════════════════

def mutate_residue(obj_name, chain, resi, mut_three):
    """
    Create a copy of obj_name, apply a point mutation using the mutagenesis
    wizard, and return the name of the mutant object.

    The mutagenesis wizard replaces the residue side chain in-place on the
    duplicate object so we can overlay wild-type and mutant.
    """
    mut_obj = f"{obj_name}_mut_{resi}"
    cmd.create(mut_obj, obj_name)

    # Use the mutagenesis wizard to mutate the residue
    cmd.wizard("mutagenesis")
    wizard = cmd.get_wizard()

    # Select the residue to mutate on the duplicate
    cmd.get_wizard().do_select(f"/{mut_obj}//{chain}/{resi}")
    cmd.get_wizard().set_mode(mut_three)

    # Apply the mutation (pick the most common rotamer — index 1)
    cmd.get_wizard().apply()
    cmd.set_wizard()  # dismiss wizard

    print(f"    Mutated {obj_name} resi {resi} -> {mut_three} in {mut_obj}")
    return mut_obj


# ═════════════════════════════════════════════════════════════════════════════
# Per-variant render
# ═════════════════════════════════════════════════════════════════════════════

def render_variant(var):
    """Render a single variant showing WT and mutant overlaid with context."""
    name = var["name"]
    resi = var["resi"]
    mut_three = var["mut_three"]
    wt_three = var["wt_three"]
    mechanism = var["mechanism"]

    print(f"\n  Rendering {name} ({mechanism})...")

    # ── Reset view to base scene ──────────────────────────────────────────
    cmd.hide("sticks", "all")
    cmd.hide("labels", "all")
    cmd.hide("dashes", "all")
    cmd.show("cartoon", "protein_all")
    cmd.show("cartoon", "dna_all")
    cmd.show("spheres", "zinc_ions")

    # Delete any leftover objects from previous iterations
    for obj in cmd.get_names("objects"):
        if obj.startswith("p53_dna_mut_") or obj in (
            "wt_res", "mut_res", "nearby", "dna_contacts",
            "dna_hbonds", "zn_bonds", "focus",
        ):
            cmd.delete(obj)

    # ── Wild-type residue: green sticks ───────────────────────────────────
    cmd.select("wt_res", f"chain {PROTEIN_CHAIN} and resi {resi} and p53_dna")
    cmd.show("sticks", "wt_res")
    cmd.color("green", "wt_res")
    cmd.set("stick_radius", 0.2, "wt_res")

    # ── Create mutant copy and show mutated residue ───────────────────────
    mut_obj = mutate_residue("p53_dna", PROTEIN_CHAIN, resi, mut_three)

    # Hide everything on the mutant object except the mutated residue
    cmd.hide("everything", mut_obj)
    cmd.select("mut_res", f"{mut_obj} and chain {PROTEIN_CHAIN} and resi {resi}")
    cmd.show("sticks", "mut_res")
    cmd.color("salmon", "mut_res")
    cmd.set("stick_radius", 0.2, "mut_res")

    # ── Neighboring residues: light-blue thin sticks ──────────────────────
    cmd.select(
        "nearby",
        f"(byres (all within {CONTACT_DISTANCE} of wt_res)) "
        f"and polymer.protein and chain {PROTEIN_CHAIN} "
        f"and not wt_res and p53_dna"
    )
    cmd.show("sticks", "nearby")
    cmd.color("lightblue", "nearby")
    cmd.set("stick_radius", 0.12, "nearby")

    # ── DNA polar contacts (red dashes) ───────────────────────────────────
    cmd.select(
        "dna_contacts",
        f"(chain {DNA_CHAINS} and polymer.nucleic) "
        f"within {CONTACT_DISTANCE} of wt_res"
    )
    dna_count = cmd.count_atoms("dna_contacts")
    if dna_count > 0:
        cmd.show("sticks", "dna_contacts")
        cmd.color("brightorange", "dna_contacts")
        cmd.set("stick_radius", 0.15, "dna_contacts")
        cmd.distance("dna_hbonds", "wt_res", "dna_contacts", CONTACT_DISTANCE, 2)
        cmd.set("dash_color", "red", "dna_hbonds")
        cmd.set("dash_width", 2.5, "dna_hbonds")
        cmd.hide("labels", "dna_hbonds")

    # ── Zinc coordination (slate dashes) ──────────────────────────────────
    zn_count = cmd.count_atoms(f"zinc_ions within 3.5 of wt_res")
    if zn_count > 0:
        cmd.distance("zn_bonds", "wt_res", "zinc_ions", 3.0, 2)
        cmd.set("dash_color", "slate", "zn_bonds")
        cmd.set("dash_width", 3.0, "zn_bonds")
        cmd.hide("labels", "zn_bonds")

    # ── Labels ────────────────────────────────────────────────────────────
    cmd.label(f"wt_res and name CA", f'"{wt_three} {resi} (WT)"')
    cmd.label(f"mut_res and name CA", f'"{mut_three} {resi} (MUT)"')

    # ── Camera: orient and zoom ───────────────────────────────────────────
    if dna_count > 0:
        cmd.select("focus", "wt_res or mut_res or nearby or dna_contacts")
    else:
        cmd.select("focus", "wt_res or mut_res or nearby")
    cmd.orient("focus")
    cmd.zoom("focus", 6)
    cmd.turn("y", 15)

    # ── Store scene, ray-trace, save ──────────────────────────────────────
    cmd.scene(f"Variant_{name}", "store")

    output_path = str(FIGURES_DIR / f"{name}_impact.png").replace("\\", "/")
    cmd.ray(2400, 2400)
    cmd.png(output_path, dpi=300)
    print(f"    Saved: {name}_impact.png")

    # ── Cleanup mutant object ─────────────────────────────────────────────
    cmd.delete(mut_obj)
    cmd.delete("wt_res")
    cmd.delete("mut_res")
    cmd.delete("nearby")
    cmd.delete("dna_contacts")
    cmd.delete("dna_hbonds")
    cmd.delete("zn_bonds")
    cmd.delete("focus")


# ═════════════════════════════════════════════════════════════════════════════
# Overview render — all 5 sites as colored spheres
# ═════════════════════════════════════════════════════════════════════════════

def render_overview():
    """Render an overview showing all 5 mutation sites on the full complex."""
    print("\n  Rendering overview (all 5 sites)...")

    # Reset to base
    cmd.hide("sticks", "all")
    cmd.hide("labels", "all")
    cmd.hide("dashes", "all")
    cmd.show("cartoon", "protein_all")
    cmd.show("cartoon", "dna_all")
    cmd.show("spheres", "zinc_ions")

    site_colors = ["red", "magenta", "yellow", "cyan", "orange"]

    for var, color in zip(VARIANTS, site_colors):
        resi = var["resi"]
        sel_name = f"site_{var['name']}"
        cmd.select(sel_name, f"chain {PROTEIN_CHAIN} and resi {resi} and p53_dna")
        cmd.show("spheres", sel_name)
        cmd.color(color, sel_name)
        cmd.set("sphere_scale", 0.8, sel_name)
        cmd.label(f"{sel_name} and name CA", f'"{var["name"]}"')

    cmd.orient(f"chain {PROTEIN_CHAIN} or chain {DNA_CHAINS}")
    cmd.zoom(f"chain {PROTEIN_CHAIN} or chain {DNA_CHAINS}", 5)
    cmd.turn("y", 30)

    cmd.scene("Overview_All_Sites", "store")

    output_path = str(FIGURES_DIR / "overview_impact.png").replace("\\", "/")
    cmd.ray(2400, 2400)
    cmd.png(output_path, dpi=300)
    print("    Saved: overview_impact.png")

    # Cleanup site selections
    for var in VARIANTS:
        cmd.delete(f"site_{var['name']}")


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print(" TP53 Publication-Quality Renders (PyMOL cmd API)")
    print(" Top 5 pathogenic variants | PDB: 1TUP")
    print("=" * 72)

    # Step 1: Setup
    setup_scene()

    # Step 2: Per-variant renders
    print("\n[2/3] Rendering individual variant views...")
    for var in VARIANTS:
        render_variant(var)

    # Step 3: Overview
    print("\n[3/3] Rendering overview...")
    render_overview()

    # Save session file
    session_path = str(FIGURES_DIR / "tp53_publication.pse").replace("\\", "/")
    cmd.save(session_path)
    print(f"\n  Session saved: tp53_publication.pse")

    # Summary
    print("\n" + "=" * 72)
    print(" Render complete. Output files:")
    print("=" * 72)
    for var in VARIANTS:
        fpath = FIGURES_DIR / f"{var['name']}_impact.png"
        if fpath.exists():
            size_kb = fpath.stat().st_size / 1024
            print(f"  {fpath.name:<30} ({size_kb:.0f} KB)")
        else:
            print(f"  {fpath.name:<30} (not found)")
    overview = FIGURES_DIR / "overview_impact.png"
    if overview.exists():
        size_kb = overview.stat().st_size / 1024
        print(f"  {overview.name:<30} ({size_kb:.0f} KB)")
    print(f"  {'tp53_publication.pse':<30} (session)")
    print("=" * 72)


# When run via `pymol -cqr`, __name__ is "__main__"
if __name__ == "__main__":
    main()
