# Step 2. Adapter Trimming and Quality Filtering with Trim Galore!

**Purpose**  
Removes adapter sequences and filters low-quality reads to improve downstream alignment.

**Key Options**  
- `--paired` : Input data are paired-end reads.  
- `--quality 20` : Trims bases with Phred score < 20.  
- `--length 30` : Discards reads shorter than 30 bp after trimming.  
- `--adapter <adapter_sequence>` : Adapter sequence from your library kit (e.g., Illumina TruSeq: `AGATCGGAAGAGC`).  

**Notes**  
Check your library preparation manual for the correct adapter sequence.  

