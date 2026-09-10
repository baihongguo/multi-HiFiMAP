#!/bin/bash

# ==========================================
# CONFIGURATION
# ==========================================
R_SCRIPT="./src/multi_HiFiMAP.R"

# DIRECTORIES
IBD_BASE_DIR="./example/ibd_prep"
RESULT_DIR="./results"
LOG_DIR="./results/logs"

# INPUT FILES
RES_FILE="./example/toy_lmscaledresiduals.txt"
COR_FILE="./example/toy_cor_matrix.txt"
ID_FILE="NONE" # Set to specific path if supplying custom subject IDs, or "NONE" to default to samples.txt

mkdir -p "$LOG_DIR"
mkdir -p "$RESULT_DIR"

# Thread control to prevent cluster CPU throttling
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

echo "=========================================================="
echo "Launching Multi-Phenotype HiFiMAP (Chunked)"
echo "Max Jobs: $MAX_JOBS | Min Jobs: $MIN_JOBS"
echo "=========================================================="

for chr in "${chromosomes[@]}"; do
    
    # 1. Determine total SNPs for this chromosome
    IBD_DIR="${IBD_BASE_DIR}/chr${chr}"
    SITES_FILE="${IBD_DIR}/sites.txt"
    
    if [ ! -f "$SITES_FILE" ]; then 
        echo "Missing sites file for Chr $chr. Skipping."
        continue
    fi
    
    # Total lines minus the header
    TOTAL_SNPS=$(wc -l < "$SITES_FILE")
    ((TOTAL_SNPS--)) 
    TOTAL_DIFFS=$((TOTAL_SNPS - 1))
    
    # Calculate chunk size (ceiling division)
    CHUNK_SIZE=$(( (TOTAL_DIFFS + NUM_CHUNKS - 1) / NUM_CHUNKS ))

    echo "========================================"
    echo " Queuing Chr $chr | $TOTAL_SNPS SNPs into $NUM_CHUNKS Chunks"
    echo "========================================"

    for (( chunk=0; chunk<NUM_CHUNKS; chunk++ )); do
        
        # 2. Define Start and End indices for this chunk
        START=$(( chunk * CHUNK_SIZE + 1 ))
        END=$(( (chunk + 1) * CHUNK_SIZE ))
        
        # Cap the end index to the total diffs available
        if [ $END -gt $TOTAL_DIFFS ]; then END=$TOTAL_DIFFS; fi
        if [ $START -gt $TOTAL_DIFFS ]; then break; fi

        # 3. Manage the Queue
        running_jobs=$(jobs -p | wc -l)
        if (( running_jobs >= MAX_JOBS )); then
            while (( running_jobs >= MIN_JOBS )); do
                sleep 5
                running_jobs=$(jobs -p | wc -l)
            done
        fi

        # 4. Construct chunk-specific paths
        log_file="${LOG_DIR}/chr${chr}_chunk${chunk}.log"
        CHUNK_OUT="${RESULT_DIR}/Multi-HiFiMAP_chr${chr}_chunk${chunk}.txt"
        
        # 5. Launch the chunk
        echo "  [LAUNCH] Chr $chr | Chunk $chunk (SNPs $START to $END)..."
        nohup Rscript "$R_SCRIPT" \
            "$chr" \
            "$START" \
            "$END" \
            "$chunk" \
            "$IBD_DIR" \
            "$RES_FILE" \
            "$COR_FILE" \
            "$ID_FILE" \
            "$CHUNK_OUT" > "$log_file" 2>&1 &
        
    done
done

echo "=========================================================="
echo "All chromosomes divided and dispatched."
echo "Waiting for the final background jobs to finish..."
wait
echo "All computational jobs finished! Moving to merge phase..."

# ==========================================
# MERGE PHASE
# ==========================================
echo "=========================================================="
echo "Merging chunks into final chromosome files..."
echo "=========================================================="

for chr in "${chromosomes[@]}"; do
    FINAL_FILE="${RESULT_DIR}/Multi-HiFiMAP_chr${chr}.txt"
    
    # Create or empty the final file so we don't accidentally append to an old run
    > "$FINAL_FILE"
    
    echo "Merging Chr $chr..."
    
    for (( chunk=0; chunk<NUM_CHUNKS; chunk++ )); do
        CHUNK_FILE="${RESULT_DIR}/Multi-HiFiMAP_chr${chr}_chunk${chunk}.txt"
        
        if [ -f "$CHUNK_FILE" ]; then
            cat "$CHUNK_FILE" >> "$FINAL_FILE"
            
            # Auto-clean up intermediate chunk files after successful merge
            rm "$CHUNK_FILE"
        else
            echo "  [WARNING] Missing chunk file: $CHUNK_FILE"
        fi
    done
    
    echo "  -> Saved to $FINAL_FILE"
done

echo "=========================================================="
echo "Pipeline Complete!"
echo "=========================================================="