## Future Enhancements

- Replace the simple average ranking score with a weighted score that better reflects catalytic relevance.
- Incorporate more robust descriptors, such as phylogenetic conservation scores for catalytic site residues.
- Support remote BLAST searches and allow switching between local and remote databases.
- Replace or complement Swiss-Prot with RefSeq/NR to increase BLAST hit coverage.
- Run Boltz-2 on GPU when compatible hardware is available.
- Separate generic ligand descriptors from ligand-specific descriptors.
- Move long shell commands into dedicated .sh scripts.
- Improve documentation and code comments throughout the workflow.
- Set up dedicated Conda environments for each pipeline step to isolate dependencies, prevent version conflicts, and improve reproducibility.

## Useful Commands

```bash
snakemake --use-conda --cores all -p
```
