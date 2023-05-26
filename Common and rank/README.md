# RAF substrate project

## Scope: create a file with common proteins and rank them.

Work on the following folders:

* BRAF monomer-dimer screens

* RAF1_ATPbiotin_Gavin

`compare.py`

* Compare each excel file of folder BRAF monomer-dimer screens (and its spread sheets) with each excel file of folder RAF1_ATPbiotin_Gavin (and its spread sheets) and the results are stored in comparison_result folder (compare.py). 

* Produces some statistics (stats.txt)

`statistics.py`

* calculate the number of nulls and duplicates on all excel files.

`Consolidate.py `

* Compare each excel file of folder BRAF monomer-dimer screens (and its spread sheets) with each excel file of folder RAF1_ATPbiotin_Gavin (and its spread sheets) and the results are stored in comparison_result folder (compare.py). 

* Create the BIG.xls file which contains for each common gene an indication in which spread sheet the gene was found.

`Test_compare_rank.py`

* tests the validity of the output of the consolidate.py

* For each gene listed in the compare files produced by consolidate.py in comparison_Results folder checks whether it exists in Big.xls file. If yes, then it tags with “ok” flag the relevant cells of big.xls file. 

`Rank.py`

* Only the “All interactors” and “Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary” files have “LFQ values”.

* In each of these files for each gene there are 144 “LFQ values”.
Thus, as a first rank condition, I used the following: For each gene first sum all LFQ values then if the gene has duplicates take the maximum value of sums and finally sort all genes [sort(max(sum(LFQ values)))].

* produces the file Rank.xls which contains the list of all common genes. For each gene there is an indication in which spread sheet the gene was found. It also ranks the common genes based on the previously mentioned rank condition. Currently the genes are ranked based on (a) their LFQ values in “All interactors” file and (b) their LFQ values in “Interactome dynamics of RAF1-BRAF_Fragpipe_phospho_lfq_summary” file.

