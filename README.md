# multi-HiFiMAP
multi-phenotype High-resolution fast identity-by-descent (IBD) mapping test
<br />
Current version: 1.0.0

## Overview
**multi-HiFiMAP** is an extension of the HiFiMAP architecture designed for multi-phenotype IBD mapping. It utilizes an RV-coefficient-based double-kernel association approach, engineered through sparse algebraic moment expansions and a Hutchinson trace estimator. By leveraging a stateful C++ backend (`RcppArmadillo`) and a fast-forward streaming algorithm, multi-HiFiMAP scales to biobank-sized cohorts across hundreds of phenotypes simultaneously without holding massive dense matrices in memory.

## Quick Installation 

`multi-HiFiMAP` can be downloaded via:
```bash
git clone [https://github.com/baihongguo/multi-HiFiMAP](https://github.com/baihongguo/multi-HiFiMAP)
cd multi-HiFiMAP
```

## Preparations
The software is developed and tested in Linux environments. 

**Python Dependencies (Python >= 3.6):**
* `numpy`

**R Dependencies (R >= 4.0.0):**
* `data.table`
* `Matrix`
* `Rcpp`
* `RcppArmadillo`
* `PearsonDS`
* `MASS`

---

## Usage & Example Pipeline

Based on the repository structure, all test files are located in the `example/` directory, and all source codes are in the `src/` directory.

### Step 0: Pre-parse the IBD Segments (Python)
This step transforms the raw IBD segments into high-performance sparse matrices (`.mtx`) and differential updates (`.diff`). It automatically chunks the chromosome to allow for parallel processing.

```bash
python3 src/hapIBD_parsing/parsing_hapIBD.py \
    --ibd example/chr21_toy.ibd.gz \
    --vcf example/chr21_toy.vcf.gz \
    --output example/ibd_prep/chr21 \
    --n-checkpoints 20
```
*(This generates the chunked `.mtx` and `.diff` files, alongside `sites.txt` and `samples.txt` for the C++ engine to stream).*

<br />

### Step 1: Prepare Phenotype Residuals & Correlation Matrix
Unlike standard HiFiMAP, multi-HiFiMAP operates directly on pre-computed phenotype residuals and their correlation matrix. You **do not** need to provide a fitted GLMM `.rds` object. You only need two plain text files:

1.  **Residuals File (`toy_lmscaledresiduals.txt`):** An $N \times K$ matrix (no header, no row names) of scaled residuals for your $N$ individuals across $K$ phenotypes.
2.  **Correlation Matrix (`toy_cor_matrix.txt`):** A $K \times K$ phenotypic correlation matrix.

<br />

## The C++ Engine
The pipeline relies on a optimized C++ backend (`src/multi-HiFiMAP_helper.cpp`). This script defines the `HiFiMAPCalculator` class, which uses the `Armadillo` linear algebra library to persist the sparse IBD matrix $X$ and the phenotype matrices in memory. It computes exact Hutchinson trace expansions on the fly. 

You do not need to manually compile this file; the main R script automatically loads it via `sourceCpp()`.

<br />

### Step 3: Run Parallelized multi-HiFiMAP Scan (Bash wrapper)
To run the actual association scan efficiently across all chunks, use the provided Bash wrapper. The script automatically divides the chromosome's testing sites into chunks, queues the R jobs to prevent overloading your cluster, and merges the results.

Make sure the configuration paths inside `Run_multi_HiFiMAP_parallel.sh` point to your output directories:

```bash
#!/bin/bash

# ==========================================
# CONFIGURATION
# ==========================================
R_SCRIPT="./src/multi_HiFiMAP.R"
LOG_DIR="./results/logs"
IBD_BASE_DIR="./example/ibd_prep"
RESULT_DIR="./results/UDIPs_results"

mkdir -p "$LOG_DIR"
mkdir -p "$RESULT_DIR"

# Limit threads to prevent cluster CPU throttling
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# BATCH SETTINGS
MAX_JOBS=60
MIN_JOBS=60   
chromosomes=(21)
NUM_CHUNKS=20 # Split each chromosome into parallel chunks
prefix="toy"

echo "=========================================================="
echo "Launching Multi-Phenotype HiFiMAP (Chunked)"
echo "Max Jobs: $MAX_JOBS | Min Jobs: $MIN_JOBS"
echo "=========================================================="

for chr in "${chromosomes[@]}"; do
    
    # 1. Determine total SNPs for this chromosome
    SITES_FILE="${IBD_BASE_DIR}/chr${chr}/sites.txt"
    if [ ! -f "$SITES_FILE" ]; then 
        echo "Missing sites file for Chr $chr. Skipping."
        continue
    fi
    
    TOTAL_SNPS=$(wc -l < "$SITES_FILE")
    ((TOTAL_SNPS--)) 
    TOTAL_DIFFS=$((TOTAL_SNPS - 1))
    
    CHUNK_SIZE=$(( (TOTAL_DIFFS + NUM_CHUNKS - 1) / NUM_CHUNKS ))

    echo "========================================"
    echo " Queuing Chr $chr | $TOTAL_SNPS SNPs into $NUM_CHUNKS Chunks"
    echo "========================================"

    for (( chunk=0; chunk<NUM_CHUNKS; # $END $START $TOTAL_DIFFS (( (chunk )) )); * + -gt -l) -p 1 1) 2. 3. CHUNK_SIZE Define END="$TOTAL_DIFFS;" End Manage Queue START="$((" Start [ ]; and break; chunk chunk++ do fi for if indices running_jobs the then this wc |>= MAX_JOBS )); then
            while (( running_jobs >= MIN_JOBS )); do
                sleep 5
                running_jobs=$(jobs -p | wc -l)
            done
        fi

        log_file="${LOG_DIR}/chr${chr}_chunk${chunk}.log"
        
        # 4. Launch the chunk (Args: chr, prefix, start, end, chunk_id)
        echo "  [LAUNCH] Chr $chr | Chunk $chunk (SNPs $START to $END)..."
        nohup Rscript "$R_SCRIPT" "$chr" "$prefix" "$START" "$END" "$chunk" > "$log_file" 2>&1 &
        
    done
done

echo "Waiting for the final background jobs to finish..."
wait
echo "All computational jobs finished! Moving to merge phase..."

# ==========================================
# MERGE PHASE
# ==========================================
for chr in "${chromosomes[@]}"; do
    FINAL_FILE="${RESULT_DIR}/Multi-HiFiMAP_${prefix}_chr${chr}.txt"
    > "$FINAL_FILE"
    
    echo "Merging Chr $chr..."
    
    for (( chunk=0; chunk<NUM_CHUNKS; "$CHUNK_FILE" )); -f CHUNK_FILE="${RESULT_DIR}/Multi-HiFiMAP_${prefix}_UDIPs_chr${chr}_chunk${chunk}_1cm.txt" [ ]; cat chunk++ do if then>> "$FINAL_FILE"
            rm "$CHUNK_FILE" # Clean up chunk files after successful merge
        else
            echo "  [WARNING] Missing chunk file: $CHUNK_FILE"
        fi
    done
    
    echo "  -> Saved to $FINAL_FILE"
done

echo "Pipeline Complete!"
```

### Final Output
The final merged result file (`results/UDIPs_results/Multi-HiFiMAP_toy_chr21.txt`) will contain the exact multi-phenotype association statistics for every tested variant:

```text 
chr     pos       n.ibd.segs    p.value
21      14240779  1052          0.0047846   
21      14242245  1054          0.0066947   
21      14250602  1056          0.0067795   
21      14256336  1057          0.0065234
```

