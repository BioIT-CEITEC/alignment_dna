# DNA Alignment Workflow

This repository contains a Snakemake workflow for processing and aligning DNA-Seq data.
This workflow performs DNA sequence alignment using BWA-MEM, with optional adapter trimming, UMI consensus calling, duplicate marking, and comprehensive quality control reporting via MultiQC.

## Features
- **Modular Design:** Each step is implemented as a separate rule and wrapper, allowing easy customization.
- **Conda Environments:** Each rule uses its own environment for reproducibility.
- **BioRoot Utilities:** Integrates with external modules for sample and reference management.
- **Configurable:** All parameters (paths, organism, alignment options, etc.) are set via `workflow.config.json`.

## Workflow Overview
1. **Input Processing**
   - Raw FASTQ files (`*/{sample}*fastq.gz`)
   - Optional adapter trimming with Trim Galore
2. **Alignment**
   - BWA-MEM alignment to reference genome
   - SAM/BAM conversion and sorting with Samtools
3. **Post-alignment Processing**
   - UMI consensus calling (if enabled with UMI-tools)
   - Duplicate marking/removal (Picard MarkDuplicates)
   - BAM indexing and statistics generation
4. **Quality Control**
   - Individual sample QC reports
   - MultiQC aggregated report

## Directory Structure
- `Snakefile`: Main workflow file.
- `workflow.config.json`: Configuration file.
- `rules/`: Snakemake rule files for each workflow step.
- `wrappers/`: Scripts and conda environments for each step.
- `qc_reports/`: Output directory for QC results and reports.
- `processed_fastq/`: Input directory for processed FastQ files.
- `mapped/`: Output directory for alignment BAM files.
- `logs/`: Log files for each step.

## Usage
1. **Configure the workflow:**
   - Edit `workflow.config.json` to specify sample information, reference genome, and parameters.
2. **Run the workflow:**
   ```bash
   snakemake --cores <N>
   ```
   Replace `<N>` with the number of CPU cores to use.
3. **Outputs:**
   - Aligned BAM files, coverage tracks, QC reports, and summary HTML files in the respective output directories.

## Requirements
- [Snakemake >=5.18.0](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/)
- Python 3

## Customization
- Modify rules or wrapper scripts to adapt to specific project needs.
- Add or remove steps by editing the `rules/` and `wrappers/` directories.

## Contact
For questions or contributions, please contact the BioIT-CEITEC team.
