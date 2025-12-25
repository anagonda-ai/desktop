# Data Directory

This directory should contain the KEGG metabolic pathways database file.

## Required Files

- `merged_metabolic_pathways.fasta` - KEGG metabolic pathways FASTA file (~225 MB)

## Obtaining the Files

The KEGG database file should be obtained from:
- Original location: `/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/fasta/merged_metabolic_pathways.fasta`
- Or download from KEGG database and prepare according to your needs

## BLAST Database Files (Optional)

If you want to build a BLAST database for faster searches (optional), run:

```bash
makeblastdb -in merged_metabolic_pathways.fasta -dbtype prot -out merged_metabolic_pathways
```

This will create files like:
- `merged_metabolic_pathways.phr`
- `merged_metabolic_pathways.pin`
- `merged_metabolic_pathways.psq`
- etc.

These files are not required for the tool to function, as the FASTA file is used directly as a query file.

