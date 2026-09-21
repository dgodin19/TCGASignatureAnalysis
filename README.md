## Introduction
Bladder cancer is one of the most common cancers in the world. In the United States, it is estimated that there are about 80,000 new cases each year, with increased risk usually associated with smoking or exposure to certain industrial chemicals (Sanli et al., 2017).

TCGA has generated a substantial body of multiomic cancer data, including whole-exome sequencing suitable for mutational signature analysis. The TCGA bladder cancer cohort found that mutational burden was strongly associated with APOBEC-signature mutagenesis (Robertson et al., 2017), making it a good candidate cohort for practicing signature extraction.

To get more practice with SignatureAnalyzer, I downloaded BLCA cases from TCGA and carried out the exploratory analysis below.

A corresponding Jupyter notebook is available in the repository, with more complete package versions. 

## Methods 
I started by downloading TCGA-BLCA cases manually from the TCGA website. In the resulting MAF, there were 414 cases. I then ran SignatureAnalyzer from the Broad Institute through a conda environment with python v3.10.21 and signatureanalyzer 0.0.9 (Taylor-Weiner et al., 2019). The command I used was signatureanalyzer <maf_file> -n 10 --reference cosmic2 --hg_build hg38.2bit --objective poisson. SignatureAnalyzer was run for 10 independent replicates; within each run, ARD-NMF automatically determined the number of signatures via automatic relevance determination, converging on K = 5.

Then, I analyzed the output from SignatureAnalyzer with pandas v2.3.3, numpy v2.3.5, scikit-learn v1.7.2, seaborn v0.13.2, matplotlib v3.10.6, scipy v1.16.3, and statsmodels v0.14.5.

I manually downloaded the corresponding TCGA metadata and used requests v2.32.5 and json5 v0.12.1 to retrieve the RNA-seq count files. The metadata was downloaded manually. Gene-level counts were extracted from the STAR augmented gene-count files, and samples were matched to SignatureAnalyzer exposure profiles using four part TCGA sample barcodes. After removing duplicate sample identifiers, 404 samples were shared between the expression and signature exposure matrices. Expression counts were normalized to counts per million (CPM), filtered to retain genes with CPM ≥ 1 in at least 20% of samples, and log2-transformed using log2(CPM + 1). Spearman correlations were calculated between each gene and each mutational signature across the matched samples. P-values were adjusted for multiple testing using the Benjamini–Hochberg false discovery rate (FDR) procedure, and this was applied across all gene-signature pairs jointly not per-signature. Gene-signature associations were selected using FDR < 0.05 and an absolute correlation coefficient of |ρ| > 0.25. Genes with positive and negative correlations were separated for downstream enrichment analysis using Enrichr and the Human WikiPathways 2024 gene set library (Xie et al., 2021).

## Results 
Starting out, I first looked at some basic dataset statistics. There were 414 samples. Restricting to samples with at least one SNP, the average was 276 SNPs per sample.

Figure 1 - Histogram of SNP distribution

![SNP Distribution](figures/snp_counts.pdf)

Across the entire dataset, there were: 114,367 SNPs; 2,010 deletions; and 674 insertions.

From there, I checked how concordant signatures were across the 10 NMF runs, matching signatures between each pair of runs one-to-one through the Hungarian algorithm and comparing their cosine similarity. Signatures were highly concordant, with matched cosine similarities ranging from 0.96 to 1.

Figure 2 - Run Similarity

![Run Similarity](figures/run_similarity.pdf)

A heatmap of relative signature exposure per sample showed that SBS13 and SBS5 were the dominant signature in the largest share of samples.

Figure 3 - Signature Proportions

![Signature Proportions](figures/signature_proportions.pdf)

After that, I correlated each signature's exposure against SNP burden per sample. SBS5 had a Spearman correlation of 0.42 (adjusted p-value 9.97e-19), SBS13 had a Spearman correlation of 0.86 (adjusted p-value 5.63e-123), SBS10 had a Spearman correlation of 0.26 (adjusted p-value 1.62e-07), SBS1 had a Spearman correlation of -0.06 (adjusted p-value 2.21e-01, not significant), and SBS2 had a Spearman correlation of 0.64 (adjusted p-value 2.34e-48).

Figure 4 - SBS5 mutation burden

![SBS5 mutation burden](figures/SBS5_mutation_burden.pdf)

Figure 5 - SBS10 mutation burden

![SBS10 mutation burden](figures/SBS10_mutation_burden.pdf)

Figure 6 - SBS13 mutation burden

![SBS1 mutation burden](figures/SBS13_mutation_burden.pdf)

Figure 7 - SBS1 mutation burden

![SBS1 mutation burden](figures/SBS1_mutation_burden.pdf)

Figure 8 - SBS2 mutation burden

![SBS2 mutation burden](figures/SBS2_mutation_burden.pdf)


The most similar signatures according to the cosine similarity plot are SBS5 (S1) with a similarity of 0.92, SBS13 (S2) with a similarity of 0.92, SBS10 (S3) with a similarity of 0.98, SBS1 (S4) with a similarity of 0.94, and SBS2 (S5) with a similarity of 0.97.

Figure 9 - Cosine Similarity
![Cosine Similarity](figures/cosine_similarity_plot.pdf)

Figure 10 - Signature Contributions
![Signature Contributions](figures/signature_contributions.pdf)

Looking at COSMIC, SBS5 is a clock-like signature with unknown aetiology. SBS13 is attributed to the AID/APOBEC family of cytosine deaminases. SBS10 (matching COSMIC SBS10a/10b) reflects polymerase epsilon exonuclease domain mutations. SBS1 arises from spontaneous deamination of 5-methylcytosine and is also a clock like signature. SBS2 is due to activity of the APOBEC family of cytidine deaminases.

Following identification of genes associated with each mutational signature, gene sets were submitted to Enrichr for pathway enrichment analysis (Tables 1–3). The number of genes associated with SBS10 and SBS2 was insufficient for enrichment analysis. Only pathways enriched among genes positively correlated with signature exposure are presented.

Table 1 - Human Wiki Pathways for SBS 5 positively associated genes

| Term                                                                   | Overlap | P-value               | Adjusted P-value      | Old P-value | Old Adjusted P-value | Odds Ratio          | Combined Score       | Genes                               |
|------------------------------------------------------------------------|---------|-----------------------|-----------------------|-------------|----------------------|---------------------|----------------------|-------------------------------------|
| DNA Replication WP466                                                  | 6/42    | 1.2231079907663373E-6 | 1.5166539085502583E-4 | 0           | 0                    | 20.622916666666665  | 280.76276749319874   | PRIM2;ORC6;CDC45;RFC4;RPA3;MCM10    |
| Gastric Cancer Network 2 WP2363                                        | 5/31    | 5.290768654858118E-6  | 3.2802765660120336E-4 | 0           | 0                    | 23.659818442427138  | 287.4560766337488    | RFC4;UBE2C;DSCC1;BRIX1;CACYBP       |
| Retinoblastoma Gene In Cancer WP2446                                   | 7/87    | 7.998909352439972E-6  | 3.3062158656751886E-4 | 0           | 0                    | 10.870911949685535  | 127.58325505419306   | RFC4;CDC45;CCNE1;RPA3;TTK;E2F3;SKP2 |
| G1 To S Cell Cycle Control WP45                                        | 6/64    | 1.5036290435027873E-5 | 4.6612500348586406E-4 | 0           | 0                    | 12.786206896551723  | 141.99138911130464   | PRIM2;ORC6;CDC45;CCNE1;RPA3;E2F3    |
| Cell Cycle WP179                                                       | 7/120   | 6.478063467189117E-5  | 0.0016065597398629012 | 0           | 0                    | 7.683363945010297   | 74.10223312498361    | ORC6;CDC45;CCNE1;TTK;E2F3;SKP2;BUB1 |
| DNA Repair Pathways Full Network WP4946                                | 6/120   | 5.000870214434372E-4  | 0.01012097590488324   | 0           | 0                    | 6.4868421052631575  | 49.30472522206401    | FEN1;RFC4;MSH2;FANCD2;RPA3;FANCB    |
| DNA IR Damage And Cellular Response Via ATR WP4016                     | 5/81    | 5.713454139853442E-4  | 0.01012097590488324   | 0           | 0                    | 8.073716900948023   | 60.2906150060181     | FEN1;CDC45;MSH2;RMI1;FANCD2         |
| DNA Mismatch Repair WP531                                              | 3/23    | 8.801265209300206E-4  | 0.01364196107441532   | 0           | 0                    | 18.233742331288344  | 128.28248925631055   | RFC4;MSH2;RPA3                      |
| Regulation Sister Chromatid Sep At Meta-Anaphase Transition WP4240     | 2/14    | 0.005836081726626109  | 0.08040823712240418   | 0           | 0                    | 20.14430894308943   | 103.6161941728909    | BUB1;MAD2L1                         |

Table 2 - Human Wiki Pathways for SBS 13 positively associated genes

| Term                                                             | Overlap | P-value               | Adjusted P-value   | Old P-value | Old Adjusted P-value | Odds Ratio         | Combined Score     | Genes       |
|------------------------------------------------------------------|---------|-----------------------|--------------------|-------------|----------------------|--------------------|--------------------|-------------|
| Gastric Cancer Network 2 WP2363                                  | 2/31    | 0.0055549582759753195 | 0.1672562191524187 | 0           | 0                    | 19.604926108374386 | 101.8096431909553  | RFC4;CACYBP |
| Prostaglandin Signaling WP5088                                   | 2/33    | 0.006278388355744935  | 0.1672562191524187 | 0           | 0                    | 18.338248847926266 | 92.98669413051476  | CXCL9;KLRD1 |
| Cohesin Complex Cornelia De Lange Syndrome WP5117                | 2/34    | 0.006655347029609498  | 0.1672562191524187 | 0           | 0                    | 17.764285714285716 | 89.04054540085103  | SGO1;SMC1B  |
| Fatty Acids And Lipoproteins Transport In Hepatocytes WP5323     | 2/35    | 0.007042367122207104  | 0.1672562191524187 | 0           | 0                    | 17.225108225108226 | 85.36437955003423  | APOC1;DBI   |
| Gamma Glutamyl Cycle For Glutathione Including Diseases WP4518   | 1/6     | 0.021409025932894128  | 0.2875677357918888 | 0           | 0                    | 56.12112676056338  | 215.72639402918236 | GGCT        |
| MECP2 And Associated Rett Syndrome WP3584                        | 2/73    | 0.028496368871475896  | 0.2875677357918888 | 0           | 0                    | 7.9907444668008045 | 28.430897873617184 | E2F1;EZH2   |
| Nucleotide Excision Repair In Xeroderma Pigmentosum WP5114       | 2/75    | 0.02995269546345745   | 0.2875677357918888 | 0           | 0                    | 7.771037181996086  | 27.261854982251965 | RFC4;RAD18  |
| Mammary Gland Development Pathway Involution Stage 4 Of 4 WP2815 | 1/10    | 0.03543000476018619   | 0.2875677357918888 | 0           | 0                    | 31.172143974960875 | 104.12107765370268 | E2F1        |
| Autosomal Recessive Osteopetrosis Pathways WP4788                | 1/11    | 0.03890418375818435   | 0.2875677357918888 | 0           | 0                    | 28.053521126760565 | 91.0800620630536   | SNX10       |
| SMC1 SMC3 Role In DNA Damage Cornelia De Lange Syndrome WP5118   | 1/11    | 0.03890418375818435   | 0.2875677357918888 | 0           | 0                    | 28.053521126760565 | 91.0800620630536   | RAD18       |

Table 3 - Human Wiki Pathways for SBS 1 positively associated genes

| Term                                                     | Overlap | P-value              | Adjusted P-value     | Old P-value | Old Adjusted P-value | Odds Ratio         | Combined Score     | Genes        |
|----------------------------------------------------------|---------|----------------------|----------------------|-------------|----------------------|--------------------|--------------------|--------------|
| Mesodermal Commitment Pathway WP2857                     | 2/145   | 0.004486431284512587 | 0.026918587707075524 | 0           | 0                    | 23.127039627039625 | 125.04091213046264 | MBTD1;ACVR2B |
| Factors And Pathways Affecting IGF1 Akt Signaling WP3850 | 1/36    | 0.024915106126970902 | 0.07474531838091271  | 0           | 0                    | 43.848351648351645 | 161.9004351330391  | ACVR2B       |
| Embryonic Stem Cell Pluripotency Pathways WP3931         | 1/117   | 0.07888185006491323  | 0.1352135345184284   | 0           | 0                    | 13.176392572944296 | 33.465456074891975 | ACVR2B       |
| Endoderm Differentiation WP2853                          | 1/140   | 0.09369225280452934  | 0.1352135345184284   | 0           | 0                    | 10.983397897066961 | 26.005828054598982 | MBTD1        |
| Primary Ovarian Insufficiency WP5316                     | 1/170   | 0.11267794543202368  | 0.1352135345184284   | 0           | 0                    | 9.020027309968139  | 19.692718183242675 | ACVR2B       |
| Cytokine Cytokine Receptor Interaction WP5473            | 1/255   | 0.1644842034500625   | 0.1644842034500625   | 0           | 0                    | 5.975772259236826  | 10.78591480991416  | ACVR2B       |

The most significant pathways associated with SBS 5 genes were DNA Replication WP466, Gastric Cancer Network 2 WP2363, Retinoblastoma Gene In Cancer WP2446, G1 To S Cell Cycle Control WP45, Cell Cycle WP179, DNA Repair Pathways Full Network WP4946, DNA IR Damage And Cellular Response Via ATR WP4016, and DNA Mismatch Repair WP531. SBS13 did not have any pathways with adjusted p-value < 0.05, while the most signficant pathway associated with SBS 1 genes was the Mesodermal Commitment Pathway WP2857. 

## Discussion

SBS13 (APOBEC) had by far the strongest link to overall SNP burden, and SBS2 (also APOBEC) was second. This is consistent with Robertson et al. (2017), who found that overall mutational load in this same TCGA-BLCA cohort was associated with APOBEC-signature mutagenesis. SBS1 (clock-like) showed no relationship with burden, which fits its interpretation as an age-related, roughly constant-rate process rather than something that scales with tumor hypermutation. SBS10 (POLE-associated) was weaker but still significant.

None of the 5 signatures found resemble the canonical tobacco signature. This is notable given that smoking is the principal epidemiological risk factor for bladder cancer (Sanli et al., 2017). It suggests that, at least in this cohort, APOBEC-driven mutagenesis outweighs any tobacco associated signal in the extracted signatures. A useful follow up may be to check smoking history against signature exposure directly using the available TCGA-BLCA clinical annotations, rather than inferring it indirectly from signature identity.

It is odd that SBS13 had the strongest correlation with mutation burden, but produced no significant pathways. As SBS13 is driven by enzymatic activity, a coordinated transcriptional program may not be needed. The thinking is that since APOBEC can affect any number of pathways, the resulting genes are diffuse and not concentrated in a few pathways.

The Mesodermal Commitment Pathway hit for SBS1 is likely a false positive, since FDR correction was pooled across all signatures rather than applied per signature. SBS1 showed essentially no correlation with mutation burden in the first place, so this result should be treated with caution rather than interpreted biologically.

Too few genes passed threshold for enrichment in SBS10 and SBS2, which is likely a power limitation rather than an absence of a real association. For SBS10, POLE-driven hypermutation is typically confined to a handful of outlier samples, which would limit the number of genes reaching threshold across the full cohort. SBS2's low gene count likely reflects a separate issue: lower variance in SBS2 exposure across this cohort, or the same episodic, enzymatically driven nature discussed for SBS13 above; rather than the POLE mechanism specific to SBS10.

Limitations of this exercise were that signatures were matched against COSMIC v2, where COSMIC v3 has finer resolution. Since TCGA-BLCA is primarily muscle invasive disease, the results may not generalize to non-muscle invasive bladder cancer. Moreover, the signatures track burden, and are not independent proof of a mechanism.

## References

Sanli, O., Dobruch, J., Knowles, M. et al. Bladder cancer. Nat Rev Dis Primers 3, 17022 (2017). https://doi.org/10.1038/nrdp.2017.22

Robertson AG, Kim J, Al-Ahmadie H, Bellmunt J, Guo G, Cherniack AD, Hinoue T, Laird PW, Hoadley KA, Akbani R, Castro MAA, Gibb EA, Kanchi RS, Gordenin DA, Shukla SA, Sanchez-Vega F, Hansel DE, Czerniak BA, Reuter VE, Su X, de Sa Carvalho B, Chagas VS, Mungall KL, Sadeghi S, Pedamallu CS, Lu Y, Klimczak LJ, Zhang J, Choo C, Ojesina AI, Bullman S, Leraas KM, Lichtenberg TM, Wu CJ, Schultz N, Getz G, Meyerson M, Mills GB, McConkey DJ; TCGA Research Network; Weinstein JN, Kwiatkowski DJ, Lerner SP. Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer. Cell. 2017 Oct 19;171(3):540-556.e25. doi: 10.1016/j.cell.2017.09.007. Epub 2017 Oct 5. Erratum in: Cell. 2018 Aug 9;174(4):1033. doi: 10.1016/j.cell.2018.07.036. PMID: 28988769; PMCID: PMC5687509.

Taylor-Weiner, A., Aguet, F., Haradhvala, N.J. et al. Scaling computational genomics to millions of individuals with GPUs. Genome Biol 20, 228 (2019) doi:10.1186/s13059-019-1836-7 (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1836-7) 

Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A.
Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90

