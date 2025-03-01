This folder contains input data for 7 different test cases. 

Test cases:

1: No motifs in reference (ref) and alternative (alt) sequences. 
2: Motif exists in both sequences but have same compositions. 
3: Motif exists in both sequences but have different compositions. 
4: Motif exists in the ref sequence but disappears in the alt sequence. 
5: Double motif exists in the same distance from each other in both ref and alt sequences. 
6: Double motifs are spaced differently due to insertion mutation in between (i.e., motifs are orthologous).
7: Double motifs are spaced differently due to disappearance-reappearance events (i.e., motifs are analogous). 

In each folder, there are:
- RMT files used for configuring the Mutation-Simulator
- FASTA files of the mutated sequences
- TF Binding Profiles

There are 11 set of mutation rate parameter each of which comes with 10 samples.

Mutation rate parameters: indel rate is 20% of substitution rate. 

Substitution rates vary between 0.01 to 0.5
Indel rates vary between 0.002 to 0.1

Sets of mutation rate parameters: 
(subs, indel)
(0.01, 0.002)
(0.05, 0.01)
(0.10, 0.02)
(0.15, 0.03)
(0.20, 0.04)
(0.25, 0.05)
(0.30, 0.06)
(0.35, 0.07)
(0.40, 0.08)
(0.45, 0.09)
(0.50, 0.10)














