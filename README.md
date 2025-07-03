# TrufflR

TrufflR is an R-based command-line tool which extracts clean target sequences from NCBI's nuccore database for specific taxa.

## Features

1. Search NCBI nuccore for specific genes in specific taxa using taxids and gene synonyms.
2. Extract and save nucleotide and amino acid sequences for target genes.
3. Filter by feature type (CDS, gene, rRNA, tRNA, or all).
4. Combine all nucleotide or amino acid sequences into single multifasta files.
5. Generates summary tables and logs for each run.

## Requirements

- R (≥ 4.0 recommended)
- R packages: optparse, rentrez, seqinr, geneviewer, Biostrings, BiocManager
- Internet connection (for NCBI queries)

## Installation

1. **Clone or Download TrufflR**
   Download this repository or clone it using git:
   ```bash
   git clone https://github.com/yourusername/TrufflR.git
   cd TrufflR
   ```

2. **Check Rscript Availability**
   Make sure `Rscript` is available in your PATH (it comes with most R installations).

3. **Make the Script Executable**
   If you want to run it as `./trufflR.R`, make it executable:
   ```bash
   chmod +x trufflR.R
   ```

4. **Test the Installation**
   Run:
   ```bash
   Rscript trufflR.R --help
   ```
   This should print the help message and available options.

---

Go forth and find your truffles, little piggy :)

## Command-line usage

```{bash}
Rscript trufflR.R \
  --taxids=taxids.txt \
  --genes="COI[Gene],COX1[Gene],cytochrome c oxidase subunit I[Gene]" \
  --output-dir=results \
  --feature-type=CDS \
  --retmax=5 \
  --combine-nt \
  --combine-aa \
  -v
```

## Options

| Option             | Description                                                                                                                           |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| -t, --taxids       | Path to text file containing taxids (one per line) [required]                                                                         |
| -g, --genes        | Comma-separated list of gene search terms (with NCBI field tags, e.g. `COI[Gene]`) [required]                                         |
| -o, --output-dir   | Output directory (default: trufflr_output)                                                                                            |
| -f, --feature-type | Feature type to extract: CDS, gene, rRNA, tRNA, or all (default: all). <br>==NOTE: non-CDS extraction needs amending, this is a WIP== |
| -r, --retmax       | Max number of sequences to retrieve per taxid (default: 5)                                                                            |
| -c, --combine-nt   | Combine all nucleotide sequences into one file                                                                                        |
| -a, --combine-aa   | Combine all amino acid sequences into one file                                                                                        |
| --nt-file          | Name for combined nucleotide file (default: combined_nucleotide_seqs.fna)                                                             |
| --aa-file          | Name for combined amino acid file (default: combined_aminoacid_seqs.faa)                                                              |
| -v, --verbose      | Print extra output                                                                                                                    |

## Output

1. Raw files: GenBank and FASTA files for each accession in output-dir/raw_files/
2. Extracted sequences: Nucleotide and amino acid FASTA files in output-dir/nucleotide/ and output-dir/protein/
3. Combined files: If requested, combined nucleotide (.fna) and amino acid (.faa) multifasta files in output-dir/
4. Summary: Per-taxon and overall summary tables and logs in output-dir/

## Contact

- Srishti Arya at <srishti.arya@nhm.ac.uk>
- Morgan Jones at <morgan.jones@bristol.ac>
- or open a git issue :)
