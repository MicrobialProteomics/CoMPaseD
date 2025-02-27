# Result Files

CoMPaseD generates four result files: 
- [ProteinAbundance.tsv](result-protein_abundance) - Generated during execution of the export module.  
  
```{table} 
:class: result-table
Column Name | Description   |
--- | ---   |
Identifier | Protein identifier from the input FASTA file.   |
Sequence_Len[aa] | Length of the protein in amino acids.   |
Group | Group that the protein is assigned to during protease score calculation.   |
Random_Sampling_x | Protein abundance in replicate x of the Monte Carlo Simulation. A value of zero indicates that the protein is not expressed in the corresponding sampling replicate. One column for each replicate.  |  
```


- [RandomSampling.tsv](result-random_sampling) - A file listing all peptides *identified* during a CoMPaseD analysis. Since this file can become quite large, its output can be suppressed by setting the parameter [`Sampling_output`](sampling_output) to false in the parameter file.  
  
```{table} 
:class: result-table
Column Name | Description   |
--- | ---   |
index | Internal peptide index. |
peptide | Peptide sequence.   |
protein | Protein identifier from the input FASTA file.   |
location | Position of the peptides first amino acid in the protein sequence.   |
MC | Number of missed cleavage sites in peptide.   |
Enzyme | Protease that generated this peptide. Identical peptides produced by different proteases are listed as separate entries.  |
Sequence_Len[aa] | Length of the protein in amino acids.   |
Group | Group that the protein is assigned to during protease score calculation.   |
Random_Sampling_x | Protein abundance in replicate x of the Monte Carlo Simulation. One column for each replicate.  |  
DeepMSPep_prediction | Predicted detectability by DeepMSPeptide.   |
subset | String that combines enzyme and MC information.   |
ID | Internal peptide index.    |
sampling_x | Presence (1) or absence (0) of the peptide in replicate x of the Monte Carlo simulation. One column for each replicate. |  
```


- [CoMPaseD_results.tsv](CoMPaseD_results) - This file contains detailed data necessary for calculating protease scores across Monte Carlo simulations. *Filtered* columns include values calculated with a requirement of at least two unique peptides per protein, whereas *unfiltered* columns also include proteins identified by a single unique peptide. The file can be visualised in the `Results` tab of the GUI.  

```{table} 
:class: result-table
Column Number   |Column Name | Description   |
--- | ---   | ---|
1   | Protease combination | Protease combination used in this iteration. |
2   | Protein group | Group for that the score was calculated.   |
3   | Random sampling | Sampling replicate for that the score was calulated.   |
4   | Protease score (unfiltered) | Prediceted protease score for current protease combination and sampling replicate considering all identified peptides.   |
5   | Protease score (filtered) | Prediceted protease score for current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
6   | Total proteins identified (unfiltered) | Number of identified proteins in the current protease combination and sampling replicate considering all identified peptides.   |
7   | Total number of peptides identified (unfiltered) | Total number of identified peptides in the current protease combination and sampling replicate considering all identified peptides.   |
8   |Mean number of peptides per protein identified (unfiltered) | Average number of identified peptides per proteins in the current protease combination and sampling replicate considering all peptides.   |
9   | Median number of peptides per protein identified (unfiltered) | Median number of identified peptides per proteins in the current protease combination and sampling replicate considering all peptides.   |
10  | Mean protein coverage (unfiltered) | Average protein sequence coverage in the current protease combination and sampling replicate considering all peptides.   |
11  | Median protein coverage (unfiltered) | Median protein sequence coverage in the current protease combination and sampling replicate considering all peptides.   |
12  | Min peptides per protein for filtering | Number of peptides required for a protein to be identified under *filtered* criteria. Fixed to two currently.   |
13  | Total proteins identified (filtered) | Number of identified proteins in the current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
14  | Total number of peptides identified (filtered) | Total number of identified peptides in the current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
15  | Mean number of peptides per protein identified (filtered) | Average number of identified peptides per proteins in the current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
16  | Median number of peptides per protein identified (filtered) | Median number of identified peptides per proteins in the current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
17  | Mean protein coverage (filtered) | Average protein sequence coverage in the current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
18  | Median protein coverage (filtered) | Median protein sequence coverage in the current protease combination and sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
19  | Protease combination trypsin | Helper column stating that the following columns refer to the values observed for trypsin as protease (control condition). |
20  | Total proteins identified (unfiltered) trypsin | Number of identified proteins in the tryptic digestion and the current sampling replicate considering all identified peptides.   |
21  | Total number of peptides identified (unfiltered) trypsin    | Total number of identified peptides in the tryptic digestion and the current sampling replicate considering all identified peptides.   |
22  | Mean number of peptides per protein identified (unfiltered) trypsin | Average number of identified peptides per proteins in the tryptic digestion and the current sampling replicate considering all peptides.   |
23  | Median number of peptides per protein identified (unfiltered) trypsin   | Median number of identified peptides per proteins in the tryptic digestion and the current sampling replicate considering all peptides.   |
24  | Mean protein coverage (unfiltered) trypsin  | Average protein sequence coverage in the tryptic digestion and the current sampling replicate considering all peptides.   |
25  | Median protein coverage (unfiltered) trypsin    | Median protein sequence coverage in the tryptic digestion and the current sampling replicate considering all peptides.   |
26  | Min peptides per protein for filtering trypsin | Number of peptides required for a protein to be identified under *filtered* criteria in the tryptic digestion. Fixed to two currently.   |
27  | Total proteins identified (filtered) trypsin    | Number of identified proteins in the tryptic digestion and the current sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
28  | Total number of peptides identified (filtered) trypsin  | Total number of identified peptides in the tryptic digestion and the current sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
29  | Mean number of peptides per protein identified (filtered) trypsin   | Average number of identified peptides per proteins in the tryptic digestion and the current sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
30  | Median number of peptides per protein identified (filtered) trypsin | Median number of identified peptides per proteins in the tryptic digestion and the current sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
31  | Mean protein coverage (filtered) trypsin    | Average protein sequence coverage in the tryptic digestion and the current sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
32  | Median protein coverage (filtered) trypsin  | Median protein sequence coverage in the tryptic digestion and the current sampling replicate considering all peptides that belong to proteins identified by at least two peptides.   |
33  | Protein ID ratio (unfiltered) | Ratio between columns 6 and 20   |
34  | Protein ID ratio (filtered) | Ratio between columns 13 and 27   |
35  | Protein ID weight | Weighting factor of the number of identified proteins during protease score calculation.  |
36  | Peptide ID ratio (unfiltered) | Ratio between columns 7 and 21   |
37  | Peptide ID ratio (filtered) | Ratio between columns 14 and 28   |
38  | Peptide ID weight | Weighting factor of the number of identified peptides during protease score calculation.  |
39  | Protein coverage ratio (unfiltered) | Ratio between columns 10 and 24   |
40  | Protein coverage ratio (filtered) | Ratio between columns 17 and 31   |
41  | Protein coverage weight | Weighting factor of the protein coverage during protease score calculation.  |  
```


- [CoMPaseD_results_summary.tsv](CoMPaseD_results_summary) - This file summarises the predicted protease scores, including the average and standard deviation, for all protease combinations and protein groups. The values are provided *unfiltered* (*i.e.* considering all *identified* unique peptides) and *filtered* (*i.e.* considering only peptides for proteins that are *identified* by at leat two unique peptides).  
  
```{table} 
:class: result-table
Column Name | Description   |
--- | ---   |
Protease combination | Individual protease or combination of proteases for this row. |
Protein group | Group for that the score was calculated.   |
Mean_score_unfiltered | Average protease score based on unfiltered data. |
Protease SD_score_unfiltered | Standard deviation of the protease score between sampling replicates calculated from unfiltered data. |
Mean_score_filtered | Average protease score based on filtered data. |
Protease SD_score_filtered | Standard deviation of the protease score between sampling replicates calculated from filtered data. |  
```


---

(result-protein_abundance)=
## ProteinAbundance.tsv

```{table} 
:class: result-table
| Identifier	| Sequence_Len[aa]	| Group	| Random_sampling_1 | Random_sampling_2 | ...
|------------|------------------|-------|-------------------|-------------------|----
| AL009126.3_prot_CAB11777.1_1	| 446| large_proteins	| 0.000497626 | 3.51E-05 | ...
| AL009126.3_prot_CAB11777.1_2	| 378| large_proteins	| 0 | 0.002066499 | ...
| AL009126.3_prot_CAB11777.1_3	| 71| medium_proteins	| 0 | 0.000109971 | ...
| AL009126.3_prot_CAB11777.1_4	| 370| large_proteins	| 5.68E-06 | 0 | ...
| AL009126.3_prot_CAB11777.1_5	| 81| medium_proteins	| 3.36E-05 | 0.000395656 | ...
| ...	| ...| ...	| ... | ... | ...
```
  
---

(result-random_sampling)=
## RandomSampling.tsv

```{table} 
:class: result-table
| index | peptide | protein | location | MC | Enzyme | Identifier | Sequence_Len[aa] | Group | Random_sampling_1 | Random_sampling_2 | ... | DeepMSPep_prediction | subset | ID | sampling_1 | sampling_2 | ... 
|-------|---------|---------|----------|----|--------|------------|------------------|-------|------------------|------------------|-----|---------------------|--------|----|------------|------------|-----
| 397   | AALELLQVEYKELEVMDSPE | AL009126.3_prot_CAB15238.1_3421 | 111 | 5 | glu-c | AL009126.3_prot_CAB15238.1_3421 | 745 | large_proteins | 1.17683e-05 | 2.80666e-05 | ... | 0.546137 | glu-c__5 | 397 | 0 | 1 | ... 
| 1459  | AEGKTVMLVSLDGEAAGLVAVADTLKDTSRKAVARLKE | AL009126.3_prot_CAB15355.2_3532 | 605 | 5 | glu-c | AL009126.3_prot_CAB15355.2_3532 | 802 | large_proteins | 0.0002083 | 0.0 | ... | 0.484481 | glu-c__5 | 1459 | 0 | 0 | ... 
| 2494  | AGKALADYLLSKNLYFEVYTDDHLLSPFDGE | AL009126.3_prot_CAB15602.1_3775 | 84 | 5 | glu-c | AL009126.3_prot_CAB15602.1_3775 | 286 | large_proteins | 0.000523 | 0.000129 | ... | 0.815637 | glu-c__5 | 2494 | 1 | 0 | ... 
| 2751  | AGYTPMLQQYLKLKAEHQDAFLFFRLGDFYEMFFED | AL009126.3_prot_CAB13577.2_1779 | 2 | 5 | glu-c | AL009126.3_prot_CAB13577.2_1779 | 858 | large_proteins | 0.0001213 | 0.0 | ... | 0.617249 | glu-c__5 | 2751 | 0 | 0 | ... 
| 2836  | AHNVKTRNFHTQEALYVLEKEFGSELPYMLTENAE | AL009126.3_prot_CAB15641.1_3816 | 195 | 5 | glu-c | AL009126.3_prot_CAB15641.1_3816 | 254 | large_proteins | 0.00051 | 0.004293 | ... | 0.827224 | glu-c__5 | 2836 | 0 | 1 | ... 
| ...   | ...     | ...     | ...      | ... | ...    | ...        | ...              | ...   | ...              | ...              | ... | ...     | ...    | ... | ... | ... 
```
  
---

(CoMPaseD_results)=
## CoMPaseD_results.tsv

```{table} 
:class: result-table
Protease combination	| Protein group	| Random sampling	| Protease score (unfiltered)	| Protease score (filtered)	| Total proteins identified (unfiltered)	| Total number of peptides identified (unfiltered)	| Mean number of peptides per protein identified (unfiltered)	| Median number of peptides per protein identified (unfiltered)	| Mean protein coverage (unfiltered)	| Median protein coverage (unfiltered)	| Min peptides per protein for filtering	| Total proteins identified (filtered)	| Total number of peptides identified (filtered)	| Mean number of peptides per protein identified (filtered)	| Median number of peptides per protein identified (filtered)	| Mean protein coverage (filtered)	| Median protein coverage (filtered)	| Protease combination trypsin	| Total proteins identified (unfiltered) trypsin	| Total number of peptides identified (unfiltered) trypsin	| Mean number of peptides per protein identified (unfiltered) trypsin	| Median number of peptides per protein identified (unfiltered) trypsin	| Mean protein coverage (unfiltered) trypsin	| Median protein coverage (unfiltered) trypsin	| Min peptides per protein for filtering trypsin	| Total proteins identified (filtered) trypsin	| Total number of peptides identified (filtered) trypsin	| Mean number of peptides per protein identified (filtered) trypsin	| Median number of peptides per protein identified (filtered) trypsin	| Mean protein coverage (filtered) trypsin	| Median protein coverage (filtered) trypsin	| Protein ID ratio (unfiltered)	| Protein ID ratio (filtered)	| Protein ID weight	| Peptide ID ratio (unfiltered)	| Peptide ID ratio (filtered)	| Peptide ID weight	| Protein coverage ratio (unfiltered)	| Protein coverage ratio (filtered)	| Protein coverage weight   |
---	| ---	| ---	| ---	| ---   | ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---	| ---   |
lysarginase	| large_proteins	| sampling_1	| 0.991033528	| 0.988085663	| 1741	| 9586	| 5.506031017	| 3	| 0.229539862	| 0.142857143	| 2	| 1157	| 9002	| 7.780466724	| 4	| 0.3100638	| 0.23566879	| trypsin	| 1822	| 9584	| 5.260153677	| 2	| 0.225389713	| 0.141463415	| 2	| 1261	| 9023	| 7.155432197	| 4	| 0.294220951	| 0.224264706	| 0.955543359	| 0.917525773	| 1	| 1.000208681	| 0.997672614	| 1	| 1.01841321	| 1.053846772	| 1 |
lysarginase	| large_proteins	| sampling_2	| 0.991713591	| 0.988836114	| 1726	| 9650	| 5.590961761	| 3	| 0.230335557	| 0.147518687	| 2	| 1188	| 9112	| 7.67003367	| 4	| 0.301459803	| 0.223596148	| trypsin	| 1767	| 9629	| 5.449349179	| 3	| 0.231181224	| 0.144615385	| 2	| 1218	| 9080	| 7.454844007	| 4	| 0.305178172	| 0.230672836	| 0.976796831	| 0.975369458	| 1	| 1.002180912	| 1.003524229	| 1	| 0.996341974	| 0.987815742	| 1 |
lysarginase	| large_proteins	| sampling_3	| 0.989773687	| 0.992397163	| 1739	| 9650	| 5.549166187	| 2	| 0.22425363	| 0.139534884	| 2	| 1165	| 9076	| 7.79055794	| 4	| 0.300496672	| 0.219298246	| trypsin	| 1809	| 9627	| 5.32172471	| 2	| 0.222858461	| 0.134110787	| 2	| 1193	| 9011	| 7.553227158	| 4	| 0.302405791	| 0.227722772	| 0.961304588	| 0.976529757	| 1	| 1.002389114	| 1.007213406	| 1	| 1.006260337	| 0.993686897	| 1 |
lysarginase	| large_proteins	| sampling_4	| 0.998754169	| 0.994863905	| 1729	| 9620	| 5.563909774	| 3	| 0.229890039	| 0.144781145	| 2	| 1150	| 9041	| 7.86173913	| 4	| 0.309665507	| 0.235915328	| trypsin	| 1755	| 9594	| 5.466666667	| 3	| 0.227948936	| 0.145922747	| 2	| 1193	| 9032	| 7.570829841	| 5	| 0.303453218	| 0.231343284	| 0.985185185	| 0.963956412	| 1	| 1.002710027	| 1.000996457	| 1	| 1.008515518	| 1.020471985	| 1 |
lysarginase	| large_proteins	| sampling_5	| 0.985180868	| 0.980941375	| 1745	| 9725	| 5.573065903	| 2	| 0.223660864	| 0.137870855	| 2	| 1160	| 9140	| 7.879310345	| 4	| 0.301024102	| 0.232362468	| trypsin	| 1746	| 9699	| 5.554982818	| 3	| 0.234399096	| 0.145454545	| 2	| 1183	| 9136	| 7.7227388	| 4	| 0.312849496	| 0.232323232	| 0.999427262	| 0.980557904	| 1	| 1.002680689	| 1.000437828	| 1	| 0.954188252	| 0.962201012	| 1 |
lysarginase	| large_proteins	| sampling_6	| 0.986787078	| 0.990287818	| 1652	| 9647	| 5.839588378	| 3	| 0.227821016	| 0.150684932	| 2	| 1125	| 9120	| 8.106666667	| 4	| 0.301925487	| 0.232727273	| trypsin	| 1750	| 9631	| 5.503428571	| 2.5	| 0.224190039	| 0.140801179	| 2	| 1171	| 9052	| 7.730145175	| 4	| 0.30092716	| 0.238095238	| 0.944	| 0.960717336	| 1	| 1.001661302	| 1.007512152	| 1	| 1.016195977	| 1.003317503	| 1 |
...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ...	| ... |
```
  
---
  
(CoMPaseD_results_summary)=
## CoMPaseD_results_summary.tsv

```{table} 
:class: result-table

Protease combination	| Protein group	| Mean_score_unfiltered	| SD_score_unfiltered	| Mean_score_filtered	| SD_score_filtered
---	| ---	| ---	| ---	| ---	| ---
trypsin - lysarginase - glu-c - chymotrypsin - lys-c	| large_proteins	| 2.354363188	| 0.008808405	| 2.458483168	| 0.012681445
trypsin - lysarginase - glu-c - lys-c	| large_proteins	| 2.123887564	| 0.008527206	| 2.211297635	| 0.012756247
trypsin - lysarginase - chymotrypsin - lys-c	| large_proteins	| 2.116116367	| 0.006825438	| 2.203260086	| 0.009988885
trypsin - glu-c - chymotrypsin - lys-c	| large_proteins	| 2.111351133	| 0.007827643	| 2.198311019	| 0.010794093
...	| ...	| ...	| ...	| ...	| ...

```