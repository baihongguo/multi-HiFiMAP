#!/usr/bin/env python3
"""
Version: 1.0.0
  - Auto-generates sites.txt and samples.txt for downstream R integration.
  - Optimize resume feature, now can seamlessly resume from the diff files generating phase or mtx files generating phase if interrupted.
"""

import gzip
import os
import sys
import argparse
import time
import glob
import gc
import math
import array
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import Optional

import numpy as np


def is_gzipped(file_path: str) -> bool:
    with open(file_path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def normalize_mask_str(mask_arg: Optional[str]) -> str:
    if not mask_arg:
        return ""
    return mask_arg.replace("–", "-").replace("—", "-").strip()


def parse_vcf(vcf_path: str):
    positions = []
    sample_ids = []
    hap_to_idx = {}
    chrom = None

    opener = gzip.open if is_gzipped(vcf_path) else open
    with opener(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                sample_ids = fields[9:]
            elif not line.startswith("#"):
                fields = line.split("\t")
                if chrom is None:
                    chrom = fields[0]
                positions.append(int(fields[1]))

    for i, samp in enumerate(sample_ids):
        hap_to_idx[f"{samp}_1"] = 2 * i
        hap_to_idx[f"{samp}_2"] = 2 * i + 1

    return positions, sample_ids, hap_to_idx, chrom


def parse_mask(mask_arg: Optional[str], chrom: str):
    if not mask_arg:
        return []

    mask_arg = normalize_mask_str(mask_arg)
    intervals = []

    if os.path.isfile(mask_arg):
        with open(mask_arg) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                ch, s, e = parts[0], int(parts[1]), int(parts[2])
                if ch == chrom:
                    intervals.append((s, e))
    else:
        try:
            ch, rest = mask_arg.split(":")
            rest = rest.replace("–", "-").replace("—", "-")
            s, e = map(int, rest.split("-"))
            if ch == chrom:
                intervals.append((s, e))
        except Exception:
            raise ValueError(f"Invalid mask format: {mask_arg}")

    intervals.sort()
    merged = []
    for s, e in intervals:
        if not merged or s > merged[-1][1] + 1:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return [(a, b) for a, b in merged]


def build_masked_sites(positions, mask_intervals):
    N = len(positions)
    masked = np.zeros(N, dtype=bool)
    if not mask_intervals:
        return masked

    j = 0
    for i, p in enumerate(positions):
        while j < len(mask_intervals) and mask_intervals[j][1] < p:
            j += 1
        if j < len(mask_intervals):
            s, e = mask_intervals[j]
            if s <= p <= e:
                masked[i] = True
    return masked


def clip_segment_to_unmasked(pos1: int, pos2: int, mask_intervals):
    if not mask_intervals:
        return [(pos1, pos2)]

    kept = []
    cur = pos1

    for s, e in mask_intervals:
        if e < cur:
            continue
        if s > pos2:
            break

        if cur < s:
            left_end = min(pos2, s - 1)
            if cur <= left_end:
                kept.append((cur, left_end))

        cur = max(cur, e + 1)
        if cur > pos2:
            break

    if cur <= pos2:
        kept.append((cur, pos2))

    return kept


def read_checkpoint(checkpoint_path: str) -> int:
    if os.path.exists(checkpoint_path):
        with open(checkpoint_path) as f:
            s = f.read().strip()
            return int(s) if s else 0
    return 0


def write_checkpoint(checkpoint_path: str, lines_done: int) -> None:
    tmp = checkpoint_path + ".tmp"
    with open(tmp, "w") as f:
        f.write(str(lines_done))
    os.replace(tmp, checkpoint_path)


def find_latest_matrix(output_folder: str):
    files = glob.glob(os.path.join(output_folder, "ibd_mat_*.mtx"))
    if not files:
        return None, -1
    max_idx = -1
    best_file = None
    for f in files:
        base = os.path.basename(f)
        if not (base.startswith("ibd_mat_") and base.endswith(".mtx")):
            continue
        try:
            idx = int(base.replace("ibd_mat_", "").replace(".mtx", ""))
        except Exception:
            continue
        if idx > max_idx:
            max_idx = idx
            best_file = f
    return best_file, max_idx


def safe_iter_diff_lines(path: str):
    try:
        with open(path, "rb") as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:
                    yield parts[0], int(parts[1]), int(parts[2])
    except FileNotFoundError:
        return


def flush_diff_files(site_buffer, diff_dir: str, created_files: set):
    for site, lines in site_buffer.items():
        if not lines:
            continue
        path = os.path.join(diff_dir, f"delta_{site}.diff")
        with open(path, "ab") as f:
            for line in lines:
                f.write(line.encode())
        created_files.add(path)


def parse_ibd_to_diffs(ibd_file: str, vcf_data, diff_dir: str, mask_intervals,
                       checkpoint_interval: int = 5_000_000):
    positions, sample_ids, hap_to_idx, chrom = vcf_data
    N = len(positions)
    M = len(sample_ids) * 2

    os.makedirs(diff_dir, exist_ok=True)
    checkpoint_file = os.path.join(diff_dir, "ibd_checkpoint.txt")
    lines_done = read_checkpoint(checkpoint_file)

    opener = gzip.open if is_gzipped(ibd_file) else open
    f_in = opener(ibd_file, "rt")

    for _ in range(lines_done):
        f_in.readline()

    site_buffer = defaultdict(list)
    processed = lines_done
    start_time = time.perf_counter()
    created_files = set()

    try:
        for line in f_in:
            processed += 1
            parts = line.strip().split()
            if len(parts) < 7:
                continue

            samp1, hap1_str, samp2, hap2_str, _, pos1_str, pos2_str = parts[:7]
            key1 = f"{samp1}_{hap1_str}"
            key2 = f"{samp2}_{hap2_str}"
            if key1 not in hap_to_idx or key2 not in hap_to_idx:
                continue

            try:
                pos1 = int(pos1_str)
                pos2 = int(pos2_str)
            except ValueError:
                continue

            if pos1 > pos2:
                pos1, pos2 = pos2, pos1

            hap1 = hap_to_idx[key1]
            hap2 = hap_to_idx[key2]
            if hap1 > hap2:
                hap1, hap2 = hap2, hap1

            pieces = clip_segment_to_unmasked(pos1, pos2, mask_intervals)
            if not pieces:
                continue

            for a, b in pieces:
                start_idx = bisect_left(positions, a)
                end_idx = bisect_right(positions, b)
                if start_idx >= end_idx:
                    continue
                site_buffer[start_idx].append(f"1 {hap1} {hap2}\n")
                if end_idx < N:
                    site_buffer[end_idx].append(f"0 {hap1} {hap2}\n")

            if processed % checkpoint_interval == 0:
                flush_diff_files(site_buffer, diff_dir, created_files)
                site_buffer.clear()
                write_checkpoint(checkpoint_file, processed)
                elapsed = time.perf_counter() - start_time
                print(f"  IBD parsing: {processed} records, flushed. ({elapsed:.1f}s)")

        if site_buffer:
            flush_diff_files(site_buffer, diff_dir, created_files)
            write_checkpoint(checkpoint_file, processed)

    finally:
        f_in.close()

    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    diff_file_exists = [False] * (N + 1)
    for site in range(N + 1):
        path = os.path.join(diff_dir, f"delta_{site}.diff")
        diff_file_exists[site] = os.path.exists(path)

    return N, M, diff_file_exists


def atomic_fast_mtx_writer(path: str, active_pairs: list, M: int):
    tmp = path + "_tmp.mtx"
    row_dicts = []
    full_nnz = 0
    
    for h1, h2_list in enumerate(active_pairs):
        if not h2_list:
            row_dicts.append({})
            continue
            
        counts = {}
        for h2_plus_1 in h2_list:
            counts[h2_plus_1] = counts.get(h2_plus_1, 0) + 1
        row_dicts.append(counts)
        
        r = h1 + 1
        for c in counts:
            if r == c:
                full_nnz += 1
            else:
                full_nnz += 2
                
    with open(tmp, 'w') as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{M} {M} {full_nnz}\n")
        
        buffer = []
        for h1, counts in enumerate(row_dicts):
            if not counts: 
                continue
                
            r = h1 + 1 
            for c, val in counts.items():
                buffer.append(f"{r} {c} {val}\n")
                if r != c:
                    buffer.append(f"{c} {r} {val}\n")
                    
            if len(buffer) >= 2_000_000:
                f.writelines(buffer)
                buffer.clear()
                
        if buffer:
            f.writelines(buffer)
            
    os.replace(tmp, path)


def load_checkpoint_stream(filepath: str, active_pairs: list):
    start_time = time.perf_counter()
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('%'):
                break 
                
        lines_processed = 0
        for line in f:
            r_str, c_str, v_str = line.split()
            r_mtx = int(r_str)
            c_mtx = int(c_str)
            
            if r_mtx <= c_mtx:
                target_list = active_pairs[r_mtx - 1]
                v = int(v_str)
                if v == 1:
                    target_list.append(c_mtx)
                else:
                    for _ in range(v):
                        target_list.append(c_mtx)
                        
            lines_processed += 1
            if lines_processed % 50_000_000 == 0:
                print(f"    ... parsed {lines_processed} lines from checkpoint")
                
    elapsed = time.perf_counter() - start_time
    print(f"  Successfully streamed checkpoint into RAM in {elapsed:.1f}s.")


def build_matrices_from_diffs(N: int, M: int, diff_file_exists, diff_dir: str,
                              out_dir: str, subsample: bool, delta: int,
                              resume: bool, masked_site: np.ndarray):
    start_site = 0
    print("  Initializing 16GB C-Array tracker with Lazy Deletion...")
    active_pairs = [array.array('i') for _ in range(M)]
    needs_compaction = set()
    
    if resume:
        latest_file, latest_idx = find_latest_matrix(out_dir)
        if latest_file:
            print(f"  Resuming from site {latest_idx + 1} (Streaming {latest_file}...)")
            load_checkpoint_stream(latest_file, active_pairs)
            start_site = latest_idx + 1
        else:
            print("  No matrix checkpoint found. Starting fresh.")

    output_sites = set(range(0, N, delta)) if subsample else {N - 1}
    output_sites = {s for s in output_sites if s >= start_site}

    total_output = 0
    start_time = time.perf_counter()

    for site in range(start_site, N + 1):
        if site < N and diff_file_exists[site] and (not masked_site[site]):
            path = os.path.join(diff_dir, f"delta_{site}.diff")
            for op, h1, h2 in safe_iter_diff_lines(path):
                if not (0 <= h1 < M and 0 <= h2 < M):
                    continue
                
                if h1 > h2:
                    h1, h2 = h2, h1
                
                if op == b"1":
                    active_pairs[h1].append(h2 + 1)
                else:
                    active_pairs[h1].append(-(h2 + 1))
                    needs_compaction.add(h1)

        # --- BULLETPROOF COMPACTION PHASE ---
        if needs_compaction and (site % 300 == 0 or site in output_sites):
            for h1 in needs_compaction:
                counts = {}
                for val in active_pairs[h1]:
                    if val > 0:
                        counts[val] = counts.get(val, 0) + 1
                    else:
                        # SAFELY handles phantom deletions
                        counts[-val] = counts.get(-val, 0) - 1
                
                temp = []
                for v, c in counts.items():
                    if c > 0:
                        temp.extend([v] * c)
                active_pairs[h1] = array.array('i', temp)
            needs_compaction.clear()
        # ------------------------------------

        if site < N and site in output_sites:
            out_path = os.path.join(out_dir, f"ibd_mat_{site}.mtx")
            print(f"  Saving matrix at site {site}/{N-1} ...")
            
            if masked_site[site]:
                with open(out_path, 'w') as f:
                    f.write("%%MatrixMarket matrix coordinate integer general\n")
                    f.write(f"{M} {M} 0\n")
            else:
                atomic_fast_mtx_writer(out_path, active_pairs, M)
                
            total_output += 1

        if site % 1000 == 0 and site > start_site:
            print(f"  ... Processed site {site}/{N}")

    elapsed = time.perf_counter() - start_time
    print(f"  Matrix construction finished. Written {total_output} matrices in {elapsed:.1f}s.")
    return total_output


def compute_diff(ibd_file: str, vcf_file: str, output_folder: str,
                 subsample: bool = False, delta: int = 1000, mask: Optional[str] = None,
                 resume: bool = False, n_checkpoints: Optional[int] = None):

    print("Reading VCF...")
    vcf_data = parse_vcf(vcf_file)
    positions, sample_ids, hap_to_idx, chrom = vcf_data
    N = len(positions)
    M = len(sample_ids) * 2
    print(f"  Sites: {N}, Haplotypes: {M}")

    os.makedirs(output_folder, exist_ok=True)
    
    sites_file = os.path.join(output_folder, "sites.txt")
    if not os.path.exists(sites_file):
        with open(sites_file, "w") as f:
            f.write("CHR\tBP\n")
            for p in positions:
                f.write(f"{chrom}\t{p}\n")
                
    samples_file = os.path.join(output_folder, "samples.txt")
    if not os.path.exists(samples_file):
        with open(samples_file, "w") as f:
            for s in sample_ids:
                f.write(f"{s}\n")

    mask_intervals = parse_mask(mask, chrom)
    if mask_intervals:
        print(f"  Masking {len(mask_intervals)} region(s) on chromosome {chrom}")

    masked_site = build_masked_sites(positions, mask_intervals)

    if n_checkpoints:
        delta = math.ceil(N / n_checkpoints)
        print(f"  Auto-calculated delta: {delta} (sites per checkpoint)")

    run_dir = output_folder
    os.makedirs(run_dir, exist_ok=True)
    
    sig_text = f"chrom={chrom}\nmask={normalize_mask_str(mask)}\nN={N}\n"
    with open(os.path.join(run_dir, "run_params.txt"), "w") as f:
        f.write(sig_text)

    diff_dir = run_dir 

    run_phase_1 = True
    if resume:
        has_diffs = any(fn.startswith("delta_") and fn.endswith(".diff") for fn in os.listdir(diff_dir)) if os.path.isdir(diff_dir) else False
        has_checkpoint = os.path.exists(os.path.join(diff_dir, "ibd_checkpoint.txt"))
        latest_mat, _ = find_latest_matrix(run_dir)

        if latest_mat and has_diffs and not has_checkpoint:
            print("  [RESUME] Found matrices and diff files. Skipping Phase 1.")
            run_phase_1 = False
        elif latest_mat and not has_diffs:
            print("  [RESUME] Found matrices but no diff files. Regenerating diffs.")
            run_phase_1 = True
        elif has_checkpoint:
            print("  [RESUME] Found Phase 1 checkpoint. Resuming Phase 1.")
            run_phase_1 = True

    if run_phase_1:
        print("Processing IBD file (generating per-site diff files)...")
        N2, M2, diff_file_exists = parse_ibd_to_diffs(
            ibd_file, vcf_data, diff_dir, mask_intervals,
            checkpoint_interval=5_000_000
        )
        if N2 != N or M2 != M:
            raise RuntimeError("Inconsistent N/M after diff generation.")
    else:
        diff_file_exists = [False] * (N + 1)
        for site in range(N + 1):
            path = os.path.join(diff_dir, f"delta_{site}.diff")
            diff_file_exists[site] = os.path.exists(path)

    print("Building sparse matrices from diff files...")
    total = build_matrices_from_diffs(
        N, M, diff_file_exists, diff_dir, run_dir,
        subsample, delta, resume, masked_site
    )

    with open(os.path.join(run_dir, "meta.txt"), "w") as f:
        f.write(f"{N}\t{M}\t{delta}\n")

    print(f"Done. Generated {total} matrices.")
    print(f"Output written to: {run_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="IBD -> per-site diffs -> resumable sparse matrices with region masking.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--ibd", required=True, help="IBD file (gzip or plain)")
    parser.add_argument("--vcf", required=True, help="VCF file (gzip or plain)")
    parser.add_argument("--output", required=True, help="Output base directory")
    parser.add_argument("--subsample", action="store_true",
                        help="Write matrices every --delta sites. If false, write only last matrix.")
    parser.add_argument("--delta", type=int, default=1000, help="Subsampling interval (sites)")
    parser.add_argument("--mask", help="Mask regions (BED file or 'chr:start-end')")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--n-checkpoints", type=int,
                       help="Auto-calculate delta to produce this many checkpoint matrices.")
    args = parser.parse_args()

    # Automatically enable subsampling if the user asks for checkpoints
    if args.n_checkpoints is not None:
        args.subsample = True

    try:
        compute_diff(
            ibd_file=args.ibd,
            vcf_file=args.vcf,
            output_folder=args.output,
            subsample=args.subsample,
            delta=args.delta,
            mask=args.mask,
            resume=args.resume,
            n_checkpoints=args.n_checkpoints
        )
    except Exception as e:
        import traceback
        print("\n--- TRACEBACK ---", file=sys.stderr)
        traceback.print_exc()
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
