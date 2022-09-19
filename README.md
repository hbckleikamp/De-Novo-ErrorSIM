# De-Novo-ErrorSIM
Simulation of de novo errors for modelling proteomics.

<br>The main script: mutate_swap.py generates de novo errors either by mutation, substitution and inversion.
<br>The classes of errors are based on: 
<br># Evaluating de novo sequencing in proteomics: already an accurate alternative to database-driven peptide identification? 
<br># Thilo Muth, Bernhard Y Renard
<br># Briefings in Bioinformatics, Volume 19, Issue 5, September 2018, Pages 954–970, https://doi.org/10.1093/bib/bbx033

<br>It generates errors of different classes at different "rates" of occurrence.
Then Levenshtein.py computes the Levenshtein distance for each error rate to select error rates that give similarly "different" sequences.
The Script error combination combines different types of simulated errrors according to their observed frequency in Muth 2018.
