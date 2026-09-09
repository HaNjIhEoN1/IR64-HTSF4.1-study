"""
python fst_pipeline.py \
  --vcf jap.no_missing.vcf.gz \
  --tol tol.list \
  --sus sus.list \
  --window 10000 \
  --step 2000 \
  --out jap_fst
"""




#!/usr/bin/env python3

import argparse
import allel
import numpy as np
import pandas as pd

# =========================
# argument parser
# =========================
parser = argparse.ArgumentParser(description="Sliding-window FST calculation")

parser.add_argument("--vcf", required=True, help="Input VCF file (vcf or vcf.gz)")
parser.add_argument("--tol", required=True, help="Tolerant sample list")
parser.add_argument("--sus", required=True, help="Susceptible sample list")
parser.add_argument("--window", type=int, default=10000, help="Window size (bp)")
parser.add_argument("--step", type=int, default=2000, help="Step size (bp)")
parser.add_argument("--out", default="fst_result", help="Output prefix")

args = parser.parse_args()

# =========================
# input
# =========================
vcf_file = args.vcf
tol_file = args.tol
sus_file = args.sus
window_size = args.window
step_size = args.step
out_prefix = args.out

# =========================
# 1) sample list 읽기
# =========================
with open(tol_file) as f:
    tol_names = [x.strip() for x in f if x.strip()]

with open(sus_file) as f:
    sus_names = [x.strip() for x in f if x.strip()]

# =========================
# 2) VCF 읽기
# =========================
callset = allel.read_vcf(
    vcf_file,
    fields=["samples", "calldata/GT", "variants/POS"]
)

samples = np.array(callset["samples"])
gt = allel.GenotypeArray(callset["calldata/GT"])
pos = np.array(callset["variants/POS"])

# =========================
# 3) subpopulation index
# =========================
tol_idx = [i for i, s in enumerate(samples) if s in tol_names]
sus_idx = [i for i, s in enumerate(samples) if s in sus_names]

print("tol samples found:", len(tol_idx))
print("sus samples found:", len(sus_idx))

if len(tol_idx) == 0 or len(sus_idx) == 0:
    raise ValueError("Sample names do not match VCF header")

# =========================
# 4) allele count
# =========================
ac_tol = gt.count_alleles(subpop=tol_idx)
ac_sus = gt.count_alleles(subpop=sus_idx)

# =========================
# 5) biallelic filter
# =========================
acu = ac_tol + ac_sus
flt = acu.is_segregating() & (acu.max_allele() == 1)

gt_flt = gt.compress(flt, axis=0)
pos_flt = pos[flt]
ac_tol_flt = ac_tol.compress(flt, axis=0)[:, :2]
ac_sus_flt = ac_sus.compress(flt, axis=0)[:, :2]

print("retained variants:", len(pos_flt))

# =========================
# 6) per-variant Hudson FST
# =========================
num_h, den_h = allel.hudson_fst(ac_tol_flt, ac_sus_flt)

with np.errstate(divide='ignore', invalid='ignore'):
    fst_hudson = num_h / den_h

fst_hudson = np.where((den_h == 0) | (fst_hudson < 0), np.nan, fst_hudson)

per_variant = pd.DataFrame({
    "pos": pos_flt,
    "fst_hudson": fst_hudson
})

per_variant.to_csv(f"{out_prefix}.per_variant.tsv", sep="\t", index=False)

# =========================
# 7) windowed Hudson FST
# =========================
fst_win_h, windows_h, counts_h = allel.windowed_hudson_fst(
    pos_flt,
    ac_tol_flt,
    ac_sus_flt,
    size=window_size,
    step=step_size
)

windowed_h = pd.DataFrame({
    "window_start": windows_h[:, 0],
    "window_end": windows_h[:, 1],
    "n_variants": counts_h,
    "fst_hudson_window": fst_win_h
})

windowed_h.to_csv(f"{out_prefix}.hudson_window.tsv", sep="\t", index=False)

# =========================
# 8) windowed Weir-Cockerham FST
# =========================
fst_win_wc, windows_wc, counts_wc = allel.windowed_weir_cockerham_fst(
    pos_flt,
    gt_flt,
    subpops=[tol_idx, sus_idx],
    size=window_size,
    step=step_size
)

windowed_wc = pd.DataFrame({
    "window_start": windows_wc[:, 0],
    "window_end": windows_wc[:, 1],
    "n_variants": counts_wc,
    "fst_wc_window": fst_win_wc
})

windowed_wc.to_csv(f"{out_prefix}.wc_window.tsv", sep="\t", index=False)

print("saved:")
print(f" - {out_prefix}.per_variant.tsv")
print(f" - {out_prefix}.hudson_window.tsv")
print(f" - {out_prefix}.wc_window.tsv")
