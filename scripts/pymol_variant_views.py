#!/usr/bin/env python3
"""
TP53 Variant Structural Visualization (PDB: 1TUP)
==================================================
Two-mode script:
  1) Generates a .pml script for PyMOL GUI (always)
  2) Creates matplotlib structural context plots using BioPython (always)

If PyMOL Python bindings are available, also ray-traces PNGs directly.

PDB 1TUP: Cho et al. (1994) Science 265:346-355
  - Chains A, B, C: p53 core domain (residues ~94-312)
  - Chains E, F: 21-bp DNA duplex
  - Chain B is the monomer bound sequence-specifically to DNA
"""

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D

from Bio.PDB import PDBParser, PDBList, NeighborSearch, Selection

warnings.filterwarnings("ignore", category=UserWarning)

WORK_DIR = Path(__file__).parent
OUTPUT_DIR = WORK_DIR / "pymol_views"
OUTPUT_DIR.mkdir(exist_ok=True)

PDB_ID = "1TUP"
PROTEIN_CHAIN = "B"  # DNA-bound monomer
DNA_CHAINS = {"E", "F"}
ALL_PROTEIN_CHAINS = {"A", "B", "C"}
CONTACT_DISTANCE = 4.0

VARIANTS = [
    {"name": "L257R", "pchange": "p.Leu257Arg", "resi": 257,
     "wt": "LEU", "mut": "ARG", "esm2": -13.88, "am": 0.9938,
     "mechanism": "Hydrophobic core disruption",
     "detail": "Buried in beta-sandwich; Leu->Arg introduces charge into core"},
    {"name": "V157D", "pchange": "p.Val157Asp", "resi": 157,
     "wt": "VAL", "mut": "ASP", "esm2": -13.69, "am": 0.9992,
     "mechanism": "Hydrophobic core disruption",
     "detail": "Beta-strand in sandwich core; Val->Asp introduces charge"},
    {"name": "R248P", "pchange": "p.Arg248Pro", "resi": 248,
     "wt": "ARG", "mut": "PRO", "esm2": -12.63, "am": 0.9994,
     "mechanism": "DNA minor groove contact lost",
     "detail": "R248 inserts into DNA minor groove; Pro eliminates contact + kinks backbone"},
    {"name": "C176R", "pchange": "p.Cys176Arg", "resi": 176,
     "wt": "CYS", "mut": "ARG", "esm2": -12.53, "am": 0.9999,
     "mechanism": "Zinc coordination abolished",
     "detail": "C176 chelates structural Zn2+; Arg cannot coordinate zinc"},
    {"name": "R280I", "pchange": "p.Arg280Ile", "resi": 280,
     "wt": "ARG", "mut": "ILE", "esm2": -12.47, "am": 0.9996,
     "mechanism": "DNA major groove contact lost",
     "detail": "R280 H-bonds DNA major groove base; Ile is hydrophobic"},
]


# ═════════════════════════════════════════════════════════════════════════════
# PART 1: BioPython structural analysis + matplotlib visualization
# ═════════════════════════════════════════════════════════════════════════════

def download_pdb(pdb_id: str) -> str:
    """Download PDB file, return path."""
    pdb_dir = WORK_DIR / "pdb_files"
    pdb_dir.mkdir(exist_ok=True)
    expected = pdb_dir / f"pdb{pdb_id.lower()}.ent"
    if expected.exists():
        print(f"  [Cache hit] {expected.name}")
        return str(expected)
    pdbl = PDBList(verbose=False)
    path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(pdb_dir), file_format="pdb")
    print(f"  Downloaded {pdb_id} -> {Path(path).name}")
    return path


def analyze_structure(pdb_path: str):
    """Parse PDB and compute structural context for each variant."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, pdb_path)
    model = structure[0]

    # Collect all atoms for neighbor search
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)

    # Identify DNA atoms and zinc atoms
    dna_atoms = []
    zinc_atoms = []
    protein_chain_b_atoms = []

    for chain in model:
        cid = chain.id
        for residue in chain:
            hetflag = residue.id[0]
            if cid in DNA_CHAINS and hetflag == " ":
                dna_atoms.extend(residue.get_atoms())
            if cid == PROTEIN_CHAIN and hetflag == " ":
                protein_chain_b_atoms.extend(residue.get_atoms())
            # Zinc
            if residue.get_resname() == "ZN":
                zinc_atoms.extend(residue.get_atoms())

    dna_ns = NeighborSearch(dna_atoms) if dna_atoms else None
    zinc_ns = NeighborSearch(zinc_atoms) if zinc_atoms else None

    print(f"  Chain {PROTEIN_CHAIN}: {len(protein_chain_b_atoms)} atoms")
    print(f"  DNA (chains E+F): {len(dna_atoms)} atoms")
    print(f"  Zinc ions: {len(zinc_atoms)} atoms")

    results = []

    for var in VARIANTS:
        resi = var["resi"]
        name = var["name"]
        print(f"\n  Analyzing {var['pchange']} ({name})...")

        # Find the target residue in chain B
        target_residue = None
        chain = model[PROTEIN_CHAIN]
        for res in chain:
            if res.id[1] == resi and res.id[0] == " ":
                target_residue = res
                break

        if target_residue is None:
            print(f"    WARNING: Residue {resi} not found in chain {PROTEIN_CHAIN}")
            # Try other protein chains
            for alt_cid in ALL_PROTEIN_CHAINS:
                for res in model[alt_cid]:
                    if res.id[1] == resi and res.id[0] == " ":
                        target_residue = res
                        print(f"    Found in chain {alt_cid}")
                        break
                if target_residue:
                    break

        if target_residue is None:
            print(f"    SKIPPED: not in structure")
            results.append(None)
            continue

        resname = target_residue.get_resname()
        print(f"    Residue: {resname} {resi}")

        # Get CA atom position for center
        ca_atom = None
        for a in target_residue:
            if a.get_name() == "CA":
                ca_atom = a
                break
        center = ca_atom.get_vector().get_array() if ca_atom else list(target_residue.get_atoms())[0].get_vector().get_array()

        # Find nearby residues (all atoms within CONTACT_DISTANCE of any atom in target)
        nearby_residues = set()
        target_atoms = list(target_residue.get_atoms())
        for atom in target_atoms:
            neighbors = ns.search(atom.get_vector().get_array(), CONTACT_DISTANCE)
            for nb in neighbors:
                parent = nb.get_parent()
                if parent != target_residue and parent.id[0] == " ":
                    nearby_residues.add(parent)

        # DNA contact check
        min_dna_dist = float("inf")
        dna_contact_atoms = []
        dna_contact_residues = set()
        if dna_ns:
            for atom in target_atoms:
                dna_nearby = dna_ns.search(atom.get_vector().get_array(), CONTACT_DISTANCE)
                for da in dna_nearby:
                    dist = atom - da
                    if dist < min_dna_dist:
                        min_dna_dist = dist
                    dna_contact_atoms.append(da)
                    dna_contact_residues.add(da.get_parent())
            if not dna_contact_atoms:
                # Find nearest DNA distance
                for atom in target_atoms:
                    dna_far = dna_ns.search(atom.get_vector().get_array(), 25.0)
                    for da in dna_far:
                        dist = atom - da
                        if dist < min_dna_dist:
                            min_dna_dist = dist

        is_dna_contact = len(dna_contact_atoms) > 0

        # Zinc contact check
        min_zn_dist = float("inf")
        is_zinc_site = False
        if zinc_ns:
            for atom in target_atoms:
                zn_nearby = zinc_ns.search(atom.get_vector().get_array(), CONTACT_DISTANCE)
                if zn_nearby:
                    is_zinc_site = True
                    for za in zn_nearby:
                        dist = atom - za
                        if dist < min_zn_dist:
                            min_zn_dist = dist

        # Report
        tag = ""
        if is_dna_contact:
            tag = f"** DNA CONTACT ** ({len(dna_contact_atoms)} atoms, min {min_dna_dist:.1f}A)"
        elif is_zinc_site:
            tag = f"** ZINC SITE ** (min {min_zn_dist:.1f}A)"
        else:
            tag = f"Nearest DNA: {min_dna_dist:.1f}A"
        print(f"    Nearby residues (<{CONTACT_DISTANCE}A): {len(nearby_residues)}")
        print(f"    {tag}")

        # Collect atom coordinates for plotting
        result = {
            "variant": var,
            "target_residue": target_residue,
            "target_atoms": target_atoms,
            "center": center,
            "nearby_residues": nearby_residues,
            "is_dna_contact": is_dna_contact,
            "dna_contact_atoms": dna_contact_atoms,
            "dna_contact_residues": dna_contact_residues,
            "min_dna_dist": min_dna_dist,
            "is_zinc_site": is_zinc_site,
            "min_zn_dist": min_zn_dist,
            "tag": tag,
        }
        results.append(result)

    return results, model


def create_matplotlib_views(results, model):
    """Create per-variant 3D structural context plots with matplotlib."""
    print("\n  Generating matplotlib structural views...")

    for i, result in enumerate(results):
        if result is None:
            continue

        var = result["variant"]
        name = var["name"]
        center = result["center"]

        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection="3d")

        # Collect atoms in a 12A radius for context
        view_radius = 12.0

        # ── DNA backbone in view ──────────────────────────────────────────
        dna_coords = []
        for cid in DNA_CHAINS:
            chain = model[cid]
            for res in chain:
                if res.id[0] != " ":
                    continue
                for atom in res:
                    coord = atom.get_vector().get_array()
                    if np.linalg.norm(coord - center) <= view_radius:
                        dna_coords.append(coord)

        if dna_coords:
            dna_arr = np.array(dna_coords)
            ax.scatter(dna_arr[:, 0], dna_arr[:, 1], dna_arr[:, 2],
                       c="orange", s=8, alpha=0.3, label="DNA")

        # ── Protein backbone (CA trace) in view ───────────────────────────
        chain_b = model[PROTEIN_CHAIN]
        backbone_coords = []
        for res in chain_b:
            if res.id[0] != " ":
                continue
            for atom in res:
                if atom.get_name() == "CA":
                    coord = atom.get_vector().get_array()
                    if np.linalg.norm(coord - center) <= view_radius:
                        backbone_coords.append(coord)

        if backbone_coords:
            bb_arr = np.array(backbone_coords)
            ax.plot(bb_arr[:, 0], bb_arr[:, 1], bb_arr[:, 2],
                    c="lightgray", linewidth=1.5, alpha=0.5)
            ax.scatter(bb_arr[:, 0], bb_arr[:, 1], bb_arr[:, 2],
                       c="lightgray", s=10, alpha=0.3, label="Protein backbone")

        # ── Nearby residues: blue sticks ──────────────────────────────────
        for nearby_res in result["nearby_residues"]:
            res_coords = []
            for atom in nearby_res:
                coord = atom.get_vector().get_array()
                if np.linalg.norm(coord - center) <= view_radius:
                    res_coords.append(coord)
            if res_coords:
                nr_arr = np.array(res_coords)
                ax.scatter(nr_arr[:, 0], nr_arr[:, 1], nr_arr[:, 2],
                           c="cornflowerblue", s=25, alpha=0.6, edgecolors="navy",
                           linewidths=0.3)

        # ── Target residue: green (large) ─────────────────────────────────
        target_coords = []
        for atom in result["target_atoms"]:
            target_coords.append(atom.get_vector().get_array())
        tgt_arr = np.array(target_coords)
        ax.scatter(tgt_arr[:, 0], tgt_arr[:, 1], tgt_arr[:, 2],
                   c="limegreen", s=120, edgecolors="darkgreen",
                   linewidths=1.5, zorder=10, label=f"{name} (wild-type)")

        # ── DNA contact atoms: red highlight ──────────────────────────────
        if result["is_dna_contact"] and result["dna_contact_atoms"]:
            dc_coords = [a.get_vector().get_array() for a in result["dna_contact_atoms"]]
            dc_arr = np.array(dc_coords)
            ax.scatter(dc_arr[:, 0], dc_arr[:, 1], dc_arr[:, 2],
                       c="red", s=60, marker="D", edgecolors="darkred",
                       linewidths=1, zorder=9, label="DNA contact atoms")
            # Draw lines from variant to DNA contacts
            for dc in dc_coords:
                ax.plot([center[0], dc[0]], [center[1], dc[1]], [center[2], dc[2]],
                        "r--", linewidth=1.2, alpha=0.7)

        # ── Zinc ion ─────────────────────────────────────────────────────
        for chain in model:
            for res in chain:
                if res.get_resname() == "ZN":
                    for atom in res:
                        zn_coord = atom.get_vector().get_array()
                        if np.linalg.norm(zn_coord - center) <= view_radius:
                            ax.scatter(*zn_coord, c="slategray", s=200,
                                       marker="o", edgecolors="black",
                                       linewidths=2, zorder=11, label="Zn2+")
                            if result["is_zinc_site"]:
                                ax.plot([center[0], zn_coord[0]],
                                        [center[1], zn_coord[1]],
                                        [center[2], zn_coord[2]],
                                        "k-", linewidth=2, alpha=0.8)

        # ── Draw 4A sphere wireframe ──────────────────────────────────────
        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 10)
        sx = center[0] + CONTACT_DISTANCE * np.outer(np.cos(u), np.sin(v))
        sy = center[1] + CONTACT_DISTANCE * np.outer(np.sin(u), np.sin(v))
        sz = center[2] + CONTACT_DISTANCE * np.outer(np.ones_like(u), np.cos(v))
        ax.plot_wireframe(sx, sy, sz, color="gray", alpha=0.08, linewidth=0.3)

        # ── Labels and formatting ─────────────────────────────────────────
        contact_label = ""
        badge_color = "green"
        if result["is_dna_contact"]:
            contact_label = "  [DNA CONTACT SITE]"
            badge_color = "red"
        elif result["is_zinc_site"]:
            contact_label = "  [Zn2+ COORDINATION SITE]"
            badge_color = "slategray"

        ax.set_title(
            f"TP53 {var['pchange']} ({name}){contact_label}\n"
            f"{var['mechanism']}\n"
            f"ESM-2 LLR: {var['esm2']:.2f}  |  AlphaMissense: {var['am']:.4f}",
            fontsize=13, fontweight="bold", color=badge_color
        )

        ax.set_xlabel("X (A)")
        ax.set_ylabel("Y (A)")
        ax.set_zlabel("Z (A)")

        # Remove duplicate legend entries
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(),
                  loc="upper left", fontsize=8, framealpha=0.9)

        # Add text annotation box
        textstr = (
            f"Wild-type: {var['wt']} {var['resi']}\n"
            f"Mutant: {var['mut']}\n"
            f"Nearby residues: {len(result['nearby_residues'])}\n"
            f"Min DNA distance: {result['min_dna_dist']:.1f} A"
        )
        if result["is_zinc_site"]:
            textstr += f"\nMin Zn distance: {result['min_zn_dist']:.1f} A"

        props = dict(boxstyle="round", facecolor="wheat", alpha=0.8)
        fig.text(0.02, 0.02, textstr, fontsize=9, verticalalignment="bottom",
                 bbox=props, family="monospace")

        plt.tight_layout()
        output_path = OUTPUT_DIR / f"tp53_{name}_variant_view.png"
        fig.savefig(str(output_path), dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"    {name}: {output_path.name}")

    # ── Overview panel: all 5 positions ───────────────────────────────────
    print("\n  Generating overview panel...")
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection="3d")

    # Full DNA
    dna_all_coords = []
    for cid in DNA_CHAINS:
        for res in model[cid]:
            if res.id[0] != " ":
                continue
            for atom in res:
                if atom.get_name() in ("P", "O5'", "O3'", "C3'", "C5'"):
                    dna_all_coords.append(atom.get_vector().get_array())
    if dna_all_coords:
        dna_arr = np.array(dna_all_coords)
        ax.scatter(dna_arr[:, 0], dna_arr[:, 1], dna_arr[:, 2],
                   c="orange", s=12, alpha=0.4, label="DNA backbone")

    # Full protein chain B backbone
    bb_all = []
    for res in model[PROTEIN_CHAIN]:
        if res.id[0] != " ":
            continue
        for atom in res:
            if atom.get_name() == "CA":
                bb_all.append(atom.get_vector().get_array())
    if bb_all:
        bb_arr = np.array(bb_all)
        ax.plot(bb_arr[:, 0], bb_arr[:, 1], bb_arr[:, 2],
                c="lightgray", linewidth=1, alpha=0.4)

    # Zinc
    for chain in model:
        for res in chain:
            if res.get_resname() == "ZN":
                for atom in res:
                    zn = atom.get_vector().get_array()
                    ax.scatter(*zn, c="slategray", s=150, marker="o",
                               edgecolors="black", linewidths=2, zorder=10)

    # Plot each variant site
    colors = ["red", "magenta", "gold", "cyan", "darkorange"]
    for result, color in zip(results, colors):
        if result is None:
            continue
        var = result["variant"]
        c = result["center"]
        ax.scatter(*c, c=color, s=250, edgecolors="black",
                   linewidths=2, zorder=15, marker="*")
        tag = var["name"]
        if result["is_dna_contact"]:
            tag += "\n[DNA]"
        elif result["is_zinc_site"]:
            tag += "\n[Zn2+]"
        ax.text(c[0] + 1.5, c[1] + 1.5, c[2] + 1.5, tag,
                fontsize=10, fontweight="bold", color=color,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor=color, alpha=0.85))

    ax.set_title(
        "TP53 p53-DNA Complex (1TUP): All 5 Pathogenic Variant Sites\n"
        "Chain B (DNA-bound monomer) | DNA duplex (chains E+F)",
        fontsize=14, fontweight="bold"
    )
    ax.set_xlabel("X (A)")
    ax.set_ylabel("Y (A)")
    ax.set_zlabel("Z (A)")

    # Legend
    legend_patches = []
    for var, color in zip(VARIANTS, colors):
        legend_patches.append(mpatches.Patch(color=color, label=f"{var['name']}: {var['mechanism']}"))
    legend_patches.append(mpatches.Patch(color="orange", label="DNA"))
    legend_patches.append(mpatches.Patch(color="slategray", label="Zn2+ ion"))
    ax.legend(handles=legend_patches, loc="upper right", fontsize=8, framealpha=0.9)

    plt.tight_layout()
    overview_path = OUTPUT_DIR / "tp53_all_variants_overview.png"
    fig.savefig(str(overview_path), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    Overview: {overview_path.name}")

    return results


# ═════════════════════════════════════════════════════════════════════════════
# PART 2: Generate .pml script for PyMOL GUI
# ═════════════════════════════════════════════════════════════════════════════

def generate_pml_script(results):
    """Write a .pml file that can be opened in PyMOL to produce ray-traced views."""
    pml_path = OUTPUT_DIR / "tp53_variant_views.pml"

    lines = []
    lines.append("# TP53 Variant Structural Visualization")
    lines.append("# PDB: 1TUP | Top 5 pathogenic variants (ESM-2 + AlphaMissense)")
    lines.append("# Open in PyMOL: File > Run Script > tp53_variant_views.pml")
    lines.append("# Or: pymol tp53_variant_views.pml")
    lines.append("")
    lines.append("# ── Load structure ────────────────────────────────────────────")
    lines.append("fetch 1TUP, p53_dna, type=pdb")
    lines.append("remove resn HOH")
    lines.append("remove solvent")
    lines.append("")
    lines.append("# ── Base styling ──────────────────────────────────────────────")
    lines.append("hide everything, all")
    lines.append("")
    lines.append("# Protein: light gray cartoon")
    lines.append("select protein_all, chain A+B+C and polymer.protein")
    lines.append("show cartoon, protein_all")
    lines.append("color gray80, protein_all")
    lines.append("color gray90, chain B and polymer.protein")
    lines.append("")
    lines.append("# DNA: orange")
    lines.append("select dna_all, chain E+F and polymer.nucleic")
    lines.append("show cartoon, dna_all")
    lines.append("color orange, dna_all")
    lines.append("set cartoon_ring_mode, 3")
    lines.append("set cartoon_ring_color, orange, dna_all")
    lines.append("")
    lines.append("# Zinc ions")
    lines.append("select zinc_ions, resn ZN")
    lines.append("show spheres, zinc_ions")
    lines.append("color slate, zinc_ions")
    lines.append("set sphere_scale, 0.5, zinc_ions")
    lines.append("")
    lines.append("# Global render settings")
    lines.append("set ray_opaque_background, 1")
    lines.append("bg_color white")
    lines.append("set antialias, 2")
    lines.append("set ray_trace_mode, 1")
    lines.append("set ray_shadows, 1")
    lines.append("set spec_reflect, 0.3")
    lines.append("set cartoon_fancy_helices, 1")
    lines.append("set cartoon_smooth_loops, 1")
    lines.append("set label_font_id, 7")
    lines.append("set label_size, 14")
    lines.append("set label_color, black")
    lines.append("set float_labels, 1")
    lines.append("")

    # Per-variant scenes
    for idx, result in enumerate(results):
        if result is None:
            continue
        var = result["variant"]
        name = var["name"]
        resi = var["resi"]

        lines.append(f"# ── Variant {idx+1}: {var['pchange']} ({name}) ──────────────────")
        lines.append(f"# {var['mechanism']}")
        lines.append(f"# ESM-2 LLR: {var['esm2']:.2f} | AlphaMissense: {var['am']:.4f}")
        lines.append("")

        # Reset for this view
        lines.append("hide sticks, all")
        lines.append("hide labels, all")
        lines.append("show cartoon, protein_all")
        lines.append("show cartoon, dna_all")
        lines.append("show spheres, zinc_ions")
        lines.append("")

        # Wild-type residue: green sticks
        lines.append(f"select wt_residue, chain B and resi {resi}")
        lines.append("show sticks, wt_residue")
        lines.append("color green, wt_residue")
        lines.append("set stick_radius, 0.2, wt_residue")
        lines.append("")

        # Nearby residues: blue
        lines.append(f"select nearby_res, (byres (all within {CONTACT_DISTANCE} of wt_residue)) and polymer.protein and not wt_residue")
        lines.append("show sticks, nearby_res")
        lines.append("color lightblue, nearby_res")
        lines.append("set stick_radius, 0.12, nearby_res")
        lines.append("")

        # DNA contacts
        if result["is_dna_contact"]:
            lines.append(f"select dna_contacts, (chain E+F and polymer.nucleic) within {CONTACT_DISTANCE} of wt_residue")
            lines.append("show sticks, dna_contacts")
            lines.append("color brightorange, dna_contacts")
            lines.append("set stick_radius, 0.15, dna_contacts")
            lines.append(f"distance dna_hbonds, wt_residue, dna_contacts, {CONTACT_DISTANCE}, 2")
            lines.append("set dash_color, red, dna_hbonds")
            lines.append("set dash_width, 2.5, dna_hbonds")
            lines.append("hide labels, dna_hbonds")
            lines.append("")
            contact_tag = " [DNA CONTACT]"
        elif result["is_zinc_site"]:
            lines.append("distance zn_bonds, wt_residue, zinc_ions, 3.0, 2")
            lines.append("set dash_color, slate, zn_bonds")
            lines.append("set dash_width, 3.0, zn_bonds")
            lines.append("hide labels, zn_bonds")
            lines.append("")
            contact_tag = " [Zn2+ SITE]"
        else:
            contact_tag = ""

        # Label
        lines.append(f'label wt_residue and name CA, "{name}{contact_tag}"')
        lines.append('label nearby_res and name CA, resn + resi')
        lines.append("")

        # Camera and export
        if result["is_dna_contact"]:
            lines.append("select focus, wt_residue or nearby_res or dna_contacts")
        else:
            lines.append("select focus, wt_residue or nearby_res")
        lines.append("orient focus")
        lines.append("zoom focus, 6")
        lines.append("turn y, 15")
        lines.append("")

        # Save scene
        lines.append(f"scene Variant_{name}, store")
        lines.append("")

        # Ray-trace and export
        output = str(OUTPUT_DIR / f"tp53_{name}_variant_view_pymol.png").replace("\\", "/")
        lines.append(f"ray 2400, 1800")
        lines.append(f'png {output}, dpi=300')
        lines.append(f"print('Exported: tp53_{name}_variant_view_pymol.png')")
        lines.append("")

        # Clean up
        lines.append("delete dna_hbonds")
        lines.append("delete zn_bonds")
        lines.append("delete wt_residue")
        lines.append("delete nearby_res")
        lines.append("delete dna_contacts")
        lines.append("delete focus")
        lines.append("")

    # ── Overview scene ────────────────────────────────────────────────────
    lines.append("# ── Overview: All 5 variant sites ─────────────────────────────")
    lines.append("hide sticks, all")
    lines.append("hide labels, all")
    lines.append("show cartoon, protein_all")
    lines.append("show cartoon, dna_all")
    lines.append("show spheres, zinc_ions")

    colors_pml = ["red", "magenta", "yellow", "cyan", "orange"]
    for var, color in zip(VARIANTS, colors_pml):
        resi = var["resi"]
        n = var["name"]
        lines.append(f"select var_{n}, chain B and resi {resi}")
        lines.append(f"show spheres, var_{n}")
        lines.append(f"color {color}, var_{n}")
        lines.append(f"set sphere_scale, 0.8, var_{n}")
        lines.append(f'label var_{n} and name CA, "{n}"')

    lines.append("orient chain B or chain E+F")
    lines.append("zoom chain B or chain E+F, 5")
    lines.append("turn y, 30")
    lines.append("scene Overview_All_Variants, store")
    overview_out = str(OUTPUT_DIR / "tp53_all_variants_overview_pymol.png").replace("\\", "/")
    lines.append("ray 2400, 1800")
    lines.append(f"png {overview_out}, dpi=300")
    lines.append("")

    # Save session
    session_out = str(OUTPUT_DIR / "tp53_variants.pse").replace("\\", "/")
    lines.append(f"save {session_out}")
    lines.append("print('Session saved: tp53_variants.pse')")
    lines.append("print('Use Scene menu to switch between variant views')")

    pml_content = "\n".join(lines) + "\n"
    with open(pml_path, "w", encoding="ascii", errors="replace") as f:
        f.write(pml_content)

    print(f"\n  PyMOL script saved: {pml_path.name}")
    print(f"  To render high-quality images, run in PyMOL:")
    print(f"    pymol {pml_path}")
    return pml_path


# ═════════════════════════════════════════════════════════════════════════════
# PART 3: Summary table
# ═════════════════════════════════════════════════════════════════════════════

def print_summary(results):
    print("\n" + "=" * 82)
    print(" STRUCTURAL ANALYSIS SUMMARY — TP53 Top 5 Pathogenic Variants (PDB: 1TUP)")
    print("=" * 82)
    print(f" {'#':<3} {'Variant':<10} {'ESM-2':>8} {'AM':>8} {'Mechanism':<32} {'DNA dist':>9} {'Flag'}")
    print(" " + "-" * 80)

    for i, result in enumerate(results):
        if result is None:
            continue
        var = result["variant"]
        flag = ""
        if result["is_dna_contact"]:
            flag = "DNA CONTACT"
        elif result["is_zinc_site"]:
            flag = "Zn2+ SITE"
        else:
            flag = "Structural"

        dist_str = f"{result['min_dna_dist']:.1f} A"
        print(
            f" {i+1:<3} {var['name']:<10} {var['esm2']:>8.2f} {var['am']:>8.4f} "
            f"{var['mechanism']:<32} {dist_str:>9} {flag}"
        )

    print(" " + "-" * 80)
    print()
    print(" Damage mechanisms:")
    print("   R248P  Arg248 reaches into the DNA minor groove.")
    print("          Pro is rigid, cannot form H-bonds -> complete loss of DNA contact")
    print("   R280I  Arg280 H-bonds a guanine in the major groove.")
    print("          Ile is hydrophobic, no H-bond capacity -> DNA reading abolished")
    print("   C176R  Cys176 coordinates the structural Zn2+ ion (with H179, C238, C242).")
    print("          Arg cannot coordinate zinc -> loop L2/L3 scaffolding collapses")
    print("   L257R  Leu257 is buried in the hydrophobic beta-sandwich core.")
    print("          Arg+ introduces a positive charge -> fold destabilization")
    print("   V157D  Val157 packs in the beta-sandwich interior.")
    print("          Asp- introduces negative charge + shorter side chain -> cavity + instability")
    print("=" * 82)


# ═════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 82)
    print(" TP53 Variant Structural Visualization (PDB: 1TUP)")
    print(" Top 5 confirmed pathogenic variants (ESM-2 + AlphaMissense)")
    print("=" * 82)

    # Step 1: Download PDB
    print("\n[1/4] Downloading PDB structure...")
    pdb_path = download_pdb(PDB_ID)

    # Step 2: Structural analysis with BioPython
    print("\n[2/4] Structural analysis (BioPython)...")
    results, model = analyze_structure(pdb_path)

    # Step 3: Matplotlib visualizations
    print("\n[3/4] Creating structural context views (matplotlib)...")
    create_matplotlib_views(results, model)

    # Step 4: Generate PyMOL script
    print("\n[4/4] Generating PyMOL .pml script for high-quality renders...")
    pml_path = generate_pml_script(results)

    # Summary
    print_summary(results)

    print(f"\n Output directory: {OUTPUT_DIR}")
    print(f" Files generated:")
    for f in sorted(OUTPUT_DIR.glob("tp53_*")):
        size_kb = f.stat().st_size / 1024
        print(f"   {f.name:<50} ({size_kb:.0f} KB)")


if __name__ == "__main__":
    main()
