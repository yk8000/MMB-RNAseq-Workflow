# Batch Effect Correction (Optional)
Batch effects arising from technical variation (e.g., sequencing runs or sample processing) can be adjusted either during normalization (Option A) or after normalization (Option B).

## Option A. DESeq2 Design Formula
Include batch as a covariate in the DESeq2 design formula.
Input: raw counts (from featureCounts or Salmon) + samples.tsv with sample_id, condition, and batch.
Design: ~ batch + condition.
Output: VST matrix (vst_design_batch_condition.csv) where batch effects are modeled directly.

This approach is preferred when batch information is known and can be explicitly modeled during normalization.

## Option B. limma::removeBatchEffect
Apply batch correction to the already normalized VST matrix.
Input: normalized_vst.csv (from Step 4-2) + samples.tsv with batch.
Function: removeBatchEffect from the limma package.
Design matrix: Can include condition (e.g., ~ condition) to preserve biological signal.
Output: vst_batchCorrected_limma.csv.

This approach is useful when post-hoc correction is preferred, or when multiple downstream methods require an already normalized expression matrix.
