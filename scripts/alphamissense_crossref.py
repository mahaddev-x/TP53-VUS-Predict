#!/usr/bin/env python3
"""
AlphaMissense x ESM-2 Cross-Reference for TP53 Variants
========================================================
Downloads AlphaMissense proteome-wide amino acid substitution predictions,
filters for TP53 (UniProt P04637), merges with ESM-2 LLR scores, and
generates a comparative scatter plot.

Data source: Cheng et al., Science 381(6664):eadg7492 (2023)
             https://zenodo.org/records/10813168
"""

import os
import sys
import gzip
import csv
import re
import io
import time
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────

WORK_DIR = Path(__file__).parent
ESM2_FILE = WORK_DIR / "tp53_ai_predictions.csv"
AM_CACHE_FILE = WORK_DIR / "tp53_alphamissense.csv"
MERGED_OUTPUT = WORK_DIR / "tp53_esm2_alphamissense_merged.csv"
PLOT_OUTPUT = WORK_DIR / "tp53_esm2_vs_alphamissense.png"

UNIPROT_ID = "P04637"  # TP53

# AlphaMissense download URLs (try in order)
AM_URLS = [
    "https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz",
    "https://zenodo.org/records/10813168/files/AlphaMissense_aa_substitutions.tsv.gz",
    "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
]

# Three-letter to one-letter amino acid mapping
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def protein_change_to_1letter(pchange: str) -> str:
    """Convert p.Leu257Arg -> L257R"""
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", pchange)
    if m:
        ref = AA_3TO1.get(m.group(1), "?")
        pos = m.group(2)
        alt = AA_3TO1.get(m.group(3), "?")
        return f"{ref}{pos}{alt}"
    # Already in single-letter format like L257R
    m2 = re.match(r"([A-Z])(\d+)([A-Z])", pchange)
    if m2:
        return pchange
    return None


def download_and_filter_alphamissense() -> pd.DataFrame:
    """
    Stream-download the AlphaMissense aa_substitutions file (~1.2 GB compressed)
    and extract only TP53 (P04637) rows on-the-fly.

    The file has columns: uniprot_id, protein_variant, am_pathogenicity, am_class
    with 4 header/comment lines starting with '#'.
    """
    if AM_CACHE_FILE.exists():
        print(f"[Cache hit] Loading cached TP53 AlphaMissense data from {AM_CACHE_FILE.name}")
        return pd.read_csv(AM_CACHE_FILE)

    print("=" * 70)
    print("Downloading AlphaMissense aa_substitutions (~1.2 GB compressed)")
    print("Streaming + filtering for TP53 (P04637) on-the-fly...")
    print("=" * 70)

    tp53_rows = []
    headers = None

    for url_idx, url in enumerate(AM_URLS):
        print(f"\n[Attempt {url_idx + 1}/{len(AM_URLS)}] {url}")
        try:
            resp = requests.get(url, stream=True, timeout=30)
            resp.raise_for_status()

            total_size = int(resp.headers.get("content-length", 0))
            downloaded = 0
            last_report = 0

            # Accumulate compressed chunks, then decompress
            # We need to decompress the full gzip stream, so we use a
            # streaming decompression approach
            decompressor = gzip.GzipFile(
                fileobj=_StreamingResponseWrapper(resp, total_size)
            )
            text_stream = io.TextIOWrapper(decompressor, encoding="utf-8")

            line_count = 0
            for line in text_stream:
                line_count += 1

                # Skip comment lines
                if line.startswith("#"):
                    continue

                # Parse header
                if headers is None:
                    headers = line.strip().split("\t")
                    continue

                # Quick check: does this line start with our UniProt ID?
                if not line.startswith(UNIPROT_ID + "\t"):
                    continue

                fields = line.strip().split("\t")
                tp53_rows.append(fields)

            print(f"\n  Processed {line_count:,} total lines")
            print(f"  Found {len(tp53_rows)} TP53 (P04637) variants")
            break  # Success

        except requests.exceptions.RequestException as e:
            print(f"  Failed: {e}")
            if url_idx == len(AM_URLS) - 1:
                print("\nAll download URLs failed.")
                print("Please manually download AlphaMissense_aa_substitutions.tsv.gz")
                print("from https://zenodo.org/records/10813168")
                print(f"and extract P04637 rows to {AM_CACHE_FILE}")
                sys.exit(1)
            continue

    if not tp53_rows or not headers:
        print("ERROR: No TP53 data extracted.")
        sys.exit(1)

    # Build DataFrame
    df = pd.DataFrame(tp53_rows, columns=headers)
    df["am_pathogenicity"] = pd.to_numeric(df["am_pathogenicity"], errors="coerce")

    # Cache for future runs
    df.to_csv(AM_CACHE_FILE, index=False)
    print(f"  Cached to {AM_CACHE_FILE.name}")

    return df


class _StreamingResponseWrapper:
    """Wrap a requests streaming response as a file-like object for gzip."""

    def __init__(self, response, total_size):
        self._iter = response.iter_content(chunk_size=1024 * 256)
        self._buffer = b""
        self._total = total_size
        self._downloaded = 0
        self._last_pct = -1
        self._start_time = time.time()

    def read(self, size=-1):
        while size < 0 or len(self._buffer) < size:
            try:
                chunk = next(self._iter)
                self._buffer += chunk
                self._downloaded += len(chunk)
                self._report_progress()
            except StopIteration:
                break
        if size < 0:
            result = self._buffer
            self._buffer = b""
        else:
            result = self._buffer[:size]
            self._buffer = self._buffer[size:]
        return result

    def _report_progress(self):
        if self._total > 0:
            pct = int(self._downloaded * 100 / self._total)
            if pct != self._last_pct and pct % 5 == 0:
                elapsed = time.time() - self._start_time
                speed_mb = (self._downloaded / (1024 * 1024)) / max(elapsed, 0.01)
                print(
                    f"\r  Download: {pct}% "
                    f"({self._downloaded / (1024**3):.2f} / {self._total / (1024**3):.2f} GB) "
                    f"@ {speed_mb:.1f} MB/s",
                    end="", flush=True,
                )
                self._last_pct = pct

    def readable(self):
        return True

    def seekable(self):
        return False


def merge_datasets(esm2_df: pd.DataFrame, am_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge ESM-2 predictions with AlphaMissense data using protein variant as key.

    ESM-2 uses: p.Leu257Arg (3-letter)
    AlphaMissense uses: L257R (1-letter)
    """
    print("\n" + "=" * 70)
    print("Merging ESM-2 and AlphaMissense datasets")
    print("=" * 70)

    # Create a 1-letter protein_variant key from the ESM-2 Protein_Change column
    esm2_df = esm2_df.copy()
    esm2_df["protein_variant"] = esm2_df["Protein_Change"].apply(protein_change_to_1letter)

    # Count successful conversions
    valid = esm2_df["protein_variant"].notna().sum()
    print(f"  ESM-2 variants with valid protein_variant key: {valid} / {len(esm2_df)}")
    print(f"  AlphaMissense TP53 variants available: {len(am_df)}")

    # Merge on protein_variant
    merged = esm2_df.merge(
        am_df[["protein_variant", "am_pathogenicity", "am_class"]],
        on="protein_variant",
        how="left",
    )

    matched = merged["am_pathogenicity"].notna().sum()
    unmatched = merged["am_pathogenicity"].isna().sum()
    print(f"  Matched: {matched}")
    print(f"  Unmatched (no AlphaMissense data): {unmatched}")

    # Rename for clarity
    merged.rename(
        columns={
            "am_pathogenicity": "AlphaMissense_Pathogenicity_Score",
            "am_class": "AlphaMissense_Class",
        },
        inplace=True,
    )

    # Save merged file
    merged.to_csv(MERGED_OUTPUT, index=False)
    print(f"  Saved merged data to {MERGED_OUTPUT.name}")

    return merged


def create_scatter_plot(merged: pd.DataFrame):
    """
    Scatter plot: ESM2_LLR (x) vs AlphaMissense_Pathogenicity_Score (y).
    Top 5 ESM-2 suspects highlighted in red.
    """
    print("\n" + "=" * 70)
    print("Generating scatter plot")
    print("=" * 70)

    # Filter to rows with both scores
    plot_df = merged.dropna(subset=["ESM2_LLR", "AlphaMissense_Pathogenicity_Score"]).copy()
    print(f"  Plotting {len(plot_df)} variants with both scores")

    # Identify top 5 by most negative ESM2_LLR
    top5 = plot_df.nsmallest(5, "ESM2_LLR")
    rest = plot_df[~plot_df.index.isin(top5.index)]

    print("\n  Top 5 Suspects (most damaging by ESM-2):")
    print("  " + "-" * 75)
    print(f"  {'Rank':<5} {'Variant':<18} {'ESM2_LLR':>10} {'AM_Score':>10} {'AM_Class':<20}")
    print("  " + "-" * 75)
    for rank, (_, row) in enumerate(top5.iterrows(), 1):
        print(
            f"  {rank:<5} {row['Protein_Change']:<18} "
            f"{row['ESM2_LLR']:>10.4f} "
            f"{row['AlphaMissense_Pathogenicity_Score']:>10.4f} "
            f"{row['AlphaMissense_Class']:<20}"
        )
    print("  " + "-" * 75)

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 9))

    # Color the background variants by AlphaMissense class
    class_colors = {
        "pathogenic": "#ff9999",
        "ambiguous": "#ffd699",
        "benign": "#99ccff",
    }

    for am_class, color in class_colors.items():
        subset = rest[rest["AlphaMissense_Class"] == am_class]
        ax.scatter(
            subset["ESM2_LLR"],
            subset["AlphaMissense_Pathogenicity_Score"],
            c=color,
            alpha=0.4,
            s=20,
            edgecolors="none",
            label=f"AM: {am_class.replace('_', ' ')} (n={len(subset)})",
            zorder=2,
        )

    # Handle variants without AM class (NaN)
    no_class = rest[rest["AlphaMissense_Class"].isna()]
    if len(no_class) > 0:
        ax.scatter(
            no_class["ESM2_LLR"],
            no_class["AlphaMissense_Pathogenicity_Score"],
            c="#cccccc",
            alpha=0.3,
            s=15,
            edgecolors="none",
            label=f"Unclassified (n={len(no_class)})",
            zorder=1,
        )

    # Top 5 suspects in red with labels
    ax.scatter(
        top5["ESM2_LLR"],
        top5["AlphaMissense_Pathogenicity_Score"],
        c="red",
        s=120,
        edgecolors="black",
        linewidths=1.5,
        zorder=5,
        label="Top 5 ESM-2 suspects",
        marker="*",
    )

    # Label each top-5 point
    for _, row in top5.iterrows():
        label = row["protein_variant"]
        ax.annotate(
            label,
            xy=(row["ESM2_LLR"], row["AlphaMissense_Pathogenicity_Score"]),
            xytext=(12, 8),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            color="red",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="red", alpha=0.85),
            arrowprops=dict(arrowstyle="-", color="red", lw=0.8),
            zorder=6,
        )

    # AlphaMissense classification thresholds
    ax.axhline(y=0.564, color="firebrick", linestyle="--", linewidth=1, alpha=0.7)
    ax.axhline(y=0.340, color="steelblue", linestyle="--", linewidth=1, alpha=0.7)
    ax.text(
        ax.get_xlim()[1] * 0.02, 0.58,
        "AM pathogenic (>0.564)",
        fontsize=8, color="firebrick", alpha=0.8,
    )
    ax.text(
        ax.get_xlim()[1] * 0.02, 0.30,
        "AM benign (<0.340)",
        fontsize=8, color="steelblue", alpha=0.8,
    )

    # ESM-2 severity threshold
    ax.axvline(x=-4.0, color="darkred", linestyle=":", linewidth=1, alpha=0.5)
    ax.text(
        -3.8, ax.get_ylim()[0] + 0.02,
        "ESM-2: strongly\ndamaging (<-4.0)",
        fontsize=7, color="darkred", alpha=0.7, rotation=90, va="bottom",
    )

    # Quadrant annotations
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.text(
        xlim[0] + 0.5, 0.95,
        "HIGH RISK\n(both models agree)",
        fontsize=10, fontweight="bold", color="darkred",
        ha="left", va="top", alpha=0.6,
        bbox=dict(boxstyle="round", facecolor="#ffe0e0", edgecolor="darkred", alpha=0.3),
        zorder=1,
    )
    ax.text(
        xlim[1] - 0.5, 0.05,
        "LOW RISK\n(both models agree)",
        fontsize=10, fontweight="bold", color="navy",
        ha="right", va="bottom", alpha=0.6,
        bbox=dict(boxstyle="round", facecolor="#e0e0ff", edgecolor="navy", alpha=0.3),
        zorder=1,
    )

    # Axes and title
    ax.set_xlabel("ESM-2 Log-Likelihood Ratio (LLR)\n← more damaging", fontsize=12)
    ax.set_ylabel("AlphaMissense Pathogenicity Score\nmore pathogenic →", fontsize=12)
    ax.set_title(
        "TP53 Variant Pathogenicity: ESM-2 vs AlphaMissense\n"
        "1,211 ClinVar VUS | Top 5 suspects highlighted",
        fontsize=14,
        fontweight="bold",
    )

    ax.legend(loc="lower left", fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.2)

    # Compute and display correlation
    corr_df = plot_df[["ESM2_LLR", "AlphaMissense_Pathogenicity_Score"]].dropna()
    if len(corr_df) > 2:
        pearson_r = corr_df["ESM2_LLR"].corr(corr_df["AlphaMissense_Pathogenicity_Score"])
        spearman_r = corr_df["ESM2_LLR"].corr(
            corr_df["AlphaMissense_Pathogenicity_Score"], method="spearman"
        )
        ax.text(
            0.98, 0.98,
            f"Pearson r = {pearson_r:.3f}\nSpearman ρ = {spearman_r:.3f}\nn = {len(corr_df)}",
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="lightyellow", edgecolor="gray", alpha=0.9),
        )
        print(f"\n  Pearson r  = {pearson_r:.4f}")
        print(f"  Spearman rho = {spearman_r:.4f}")

    plt.tight_layout()
    fig.savefig(PLOT_OUTPUT, dpi=200, bbox_inches="tight")
    print(f"\n  Plot saved to {PLOT_OUTPUT.name}")
    plt.close(fig)


def concordance_analysis(merged: pd.DataFrame):
    """Analyze agreement between ESM-2 and AlphaMissense classifications."""
    print("\n" + "=" * 70)
    print("Concordance Analysis")
    print("=" * 70)

    df = merged.dropna(subset=["ESM2_LLR", "AlphaMissense_Pathogenicity_Score"]).copy()

    # Classify ESM-2 predictions
    df["ESM2_Risk"] = "Low"
    df.loc[df["ESM2_LLR"] <= -4.0, "ESM2_Risk"] = "High"
    df.loc[(df["ESM2_LLR"] > -4.0) & (df["ESM2_LLR"] <= -2.0), "ESM2_Risk"] = "Medium"

    # Classify AlphaMissense
    df["AM_Risk"] = df["AlphaMissense_Class"].map({
        "pathogenic": "High",
        "ambiguous": "Medium",
        "benign": "Low",
    })

    # Cross-tabulation
    both_high = df[(df["ESM2_Risk"] == "High") & (df["AM_Risk"] == "High")]
    esm2_only = df[(df["ESM2_Risk"] == "High") & (df["AM_Risk"] != "High")]
    am_only = df[(df["ESM2_Risk"] != "High") & (df["AM_Risk"] == "High")]
    both_low = df[(df["ESM2_Risk"] == "Low") & (df["AM_Risk"] == "Low")]

    print(f"\n  Both HIGH risk:       {len(both_high):>4} variants")
    print(f"  ESM-2 HIGH only:     {len(esm2_only):>4} variants")
    print(f"  AlphaMissense HIGH:  {len(am_only):>4} variants")
    print(f"  Both LOW risk:       {len(both_low):>4} variants")

    print(f"\n  Cross-tabulation (ESM-2 rows x AlphaMissense cols):")
    ct = pd.crosstab(df["ESM2_Risk"], df["AM_Risk"], margins=True)
    # Reorder
    for idx_order in [["High", "Medium", "Low", "All"]]:
        idx_present = [x for x in idx_order if x in ct.index]
        ct = ct.reindex(idx_present)
    for col_order in [["High", "Medium", "Low", "All"]]:
        col_present = [x for x in col_order if x in ct.columns]
        ct = ct[col_present]
    print(ct.to_string(index=True))

    # Top 5 detail
    print("\n  Top 5 ESM-2 suspects — dual-model assessment:")
    print("  " + "-" * 65)
    top5 = df.nsmallest(5, "ESM2_LLR")
    for _, row in top5.iterrows():
        dual = "CONFIRMED" if row["AM_Risk"] == "High" else "DISCORDANT"
        print(
            f"  {row['Protein_Change']:<18} "
            f"ESM2={row['ESM2_Risk']:<6} AM={row['AM_Risk']:<8} --> {dual}"
        )


def main():
    print("=" * 70)
    print(" TP53 Cross-Reference: ESM-2 x AlphaMissense")
    print(" Google DeepMind AlphaMissense (2023) vs Meta ESM-2 (2023)")
    print("=" * 70)

    # Step 1: Load ESM-2 data
    print(f"\n[1/4] Loading ESM-2 predictions from {ESM2_FILE.name}")
    esm2_df = pd.read_csv(ESM2_FILE)
    print(f"  Loaded {len(esm2_df)} variants")

    # Step 2: Download / load AlphaMissense data for TP53
    print(f"\n[2/4] Obtaining AlphaMissense data for TP53 ({UNIPROT_ID})")
    am_df = download_and_filter_alphamissense()
    print(f"  AlphaMissense TP53 variants: {len(am_df)}")

    # Step 3: Merge datasets
    print(f"\n[3/4] Merging datasets")
    merged = merge_datasets(esm2_df, am_df)

    # Step 4: Analysis and visualization
    print(f"\n[4/4] Analysis and visualization")
    create_scatter_plot(merged)
    concordance_analysis(merged)

    print("\n" + "=" * 70)
    print(" Done! Output files:")
    print(f"   {AM_CACHE_FILE.name:<45} (AlphaMissense TP53 cache)")
    print(f"   {MERGED_OUTPUT.name:<45} (merged dataset)")
    print(f"   {PLOT_OUTPUT.name:<45} (scatter plot)")
    print("=" * 70)


if __name__ == "__main__":
    main()
