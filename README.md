# fasta-phylo-pipeline
Pipeline that outputs a Phylogenetic tree in JPG considering homologs of the input DNA/Protein sequence, as well as a markdown relatory

The program is to be run via the following format:
```
Python phylo_analysis.py [name of .fasta file] [number of species to be considered in the analysis]
```

As a reproduceable example, the repository contains the example input file *sequence.fasta* and the outputs when the program is run in the following way:
```
Python phylo_analysis.py sequence.fasta 10
```

After running the command, the program will ask the user if the sequence corresponds to a DNA or protein sequence. If it's a DNA sequence, the program will ask the user what ORF they desire to analyse.


## Acknowledgments

- This project uses [Biopython](https://biopython.org/), a powerful library for biological computation in Python.
- Multiple sequence alignments were performed using [MUSCLE](https://drive5.com/muscle/), developed by Robert C. Edgar.
