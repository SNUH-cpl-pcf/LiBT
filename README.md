# Manual for LiBT

>1. Input Data preparation

>>1) LiBT enables isobaric label-based Tandem Mass Tags (TMT) and label-free iBAQ and LFQ data analysis, and supports .txt format files. To perform the analysis properly, the data must have the following columns. The workflow varies slightly depending on the data type. A brief description for each type is as follows.

1. Common

| Columns | Description |
| --- | --- |
| GeneName | Name(s) of the gene(s) associated with the protein(s) contained within the group. |
| ProteinID | Identifier(s) of protein(s) contained in the protein group. |

1. Label-Free
 - iBAQ (the theoretical number of tryptic peptides to normalize each protein within the same sample): The columns of the expression value of each sample must be the sample name including &quot;iBAQ&quot;.
 - LFQ (Normalized intensity resulted from maxLFQ algorithm): LFQ input data file does not need the normalization step.
 - As long as the required columns are in the input file, it is acceptable to use a custom-made file. However, it is recommended to use the original file.


| Columns | Description |
| --- | --- |
| Intensity | column name must contain iBAQ or LFQ.ex) &quot;iBAQ&quot; + sample name, &quot;LFQ&quot; + sample name |
| Unique peptides | Number of peptides associated with each protein in protein group, occuring in the order as the protein IDs occur in the &#39;Protein IDs&#39; column. |
| Sequence coverage [%] | Percentage of the sequence that is covered by the identified peptides of the best protein sequence contained in the group. |
| Only identified by site | When marked with &#39;+&#39;, this particular protein group was identified only by a modification site. |
| Reverse | When marked with &#39;+&#39;, this particular protein group contains no protein, made up of at least 50% of the peptides of the leading protein, with a peptide derived from the reversed part of the decoy database. |
| Potential contaminant | When marked with &#39;+&#39;, this particular protein group was found to be a commonly occurring contaminant. These should be removed for further data analysis. |

1. Label-Based
 - TMT : - Primarily, TMT input format is produced by Proteome Discoverer (PD). However, recently Maxquant also supports TMT quantification data. Therefore, LiBT accepts TMT quantification result files processed from both tools. If the result is from a PD research tool, users are required to check if a normalization process was done. TMT data processed by Maxquant needs to be original as it is from Maxquant to be used in LiBT. Following columns can be used for identified filtering.

Use the result from research tools as it is (Do not use custom-made files).

Input data can be filtered by selecting &quot;Only identified by site&quot;, &quot;Reverse&quot; and/or &quot;Potential contaminant&quot; options. Also, a log2-transformed input data is used for further analysis.

2) Experimental design file is not required, and the sample list used in the experiment is provided to the user by referring to the &quot;Intensity&quot; column of the input data. If the user selects a case group from the list, he can select the samples to be used as a control group among the remaining samples. However, since only pairwise analysis is possible, care must be taken in selecting case and control samples. After selecting submit, the following design matrix is ​​automatically created.

|
 | condition | replicate |
| --- | --- | --- |
| A\_1 | case | 1 |
| A\_2 | case | 2 |
| A\_3 | case | 3 |
| A\_4 | case | 4 |
| B\_1 | control | 1 |
| B\_2 | control | 2 |
| B\_3 | control | 3 |

To perform downstream analysis, replicates are randomly numbered from 1 to the number of samples in each group, but analysis using replicate is not supported.

2. Data preprocessing

- Preprocessing can be done in the second tab of the right sidebar. Users can designate the preprocessing steps by drag and drop. The order of steps is important to take note of because the process proceeds in that set order. Users can exclude it by moving it to &quot;Not to Use&quot; at the bottom. Otherwise, at least one preprocessing must be performed.

- The following is a list of preprocessing options provided by LiBT.


1. Valid value
 This option is used to determine whether to use the protein according to the ratio of missing values. The number of valid values ​​is calculated based on the case group, and when less than the number of samples is expressed in at least one group, the corresponding protein is removed. The default value is 70%.

2. Imputation
 a. QRILC : A missing data imputation method that performs the imputation of left-censored missing data using random draws from a truncated distribution with parameters estimated using quantile regression.

 b. MinProb : Performs the imputation of left-censored missing data by random draws from a Gaussian distribution centered to a minimal value. Considering an expression data matrix with n samples and p features, for each sample, the mean value of the Gaussian distribution is set to a minimal observed value in that sample. The minimal value observed is estimated as being the q-th quantile (default q = 0.01) of the observed values ​​in that sample. The standard deviation is estimated as the median of the feature standard deviations. Note that when estimating the standard deviation of the Gaussian distribution, only the peptides/proteins which present more than 50% recorded values ​​are considered.

 c. Man (Perseus-type, default) : This method is based on the popular missing value imputation procedure implemented in the Perseus software. The missing values ​​are replaced by random numbers drawn from a normal distribution of 1.8 standard deviation downshift and with a width of 0.3 of each sample.

 d. Constant : Replaces the missing values ​​by 0

3. Normalization
 - This is a step that should be used differently depending on the data type. Users can select &quot;No&quot; or &quot;Not to Use&quot; to exclude it.
 - Normalization is the variant stabilizing normalization (vsn) method on the protein intensity distribution in each sample, which is the default method of the DEP package.

4. Once preprocessing is completed, the user can check the quality control (QC) results in 4 types of plots using the DEP package. The bar plot shows the number of identified proteins for each sample, and can always be checked regardless of the user-selected step. However, only after the normalization step or when there are missing values, users can check the density plot that can show the result of distribution before or after imputation. Users can also check a box plot per sample that shows distributions before and after normalization as well as density plot to see missing value and bias as a result of density and cumulative fraction.

![](RackMultipart20211108-4-eljm2i_html_8b79089cb16a41f2.png)

3. Differential Expression Analysis

1) Statistical analysis

 - The following 4 options are provided. Statistical test results are calculated without considering replicates.
 a. T-test (ver 4.0.5, default)
 b. Wilcoxon Rank Sum (ver 4.0.5)
 c. edgeR (ver 3.32.1)
 d. limma(ver 3.46.0)

- p-value correction option (Multiple testing)
 a. Benjamini-Hochberg (BH)
 b. Bonferroni

- Significant protein filtering criteria default to p.value 0.05 and FoldChange 1.5.

2) Result plots

 a. Principal Component Analysis (PCA) plot
 By default, the sample name is not shown. If the user needs a plot with the sample name, this can be adjusted with the toggle button.
 ![](RackMultipart20211108-4-eljm2i_html_9b59782e9c212dd1.png)

b. Volcano plot
 Each point indicated in the plot is a single protein, and if the protein (DEP) satisfies the significant protein filtering criteria, it is highlighted in red. If the user drags a desired area on the plot, the results of the statistical analysis of proteins in the area are displayed in the table below. ![](RackMultipart20211108-4-eljm2i_html_e50ed7ed3ae0028c.png)

c. Correlation heatmap
 Visualize the result of correlation analysis (Pearson coefficient correlation) between samples as a heatmap.

![](RackMultipart20211108-4-eljm2i_html_18fe8e5085e6e974.png)

d. Heatmap with clusters

Heatmap with clusters visualizes the Z-score normalized log 2-transformed intensity value of proteins along with DEP Clustering results. The default value of k for gene clustering is set to 1. If the default value is 1, the optimal number of clusters k is automatically calculated, and clustering is possible even based on a specific value of k. The optimal k is calculated by the fviz\_nbclust method of the factoextra library. Clustering results can be downloaded.

![](RackMultipart20211108-4-eljm2i_html_50663921cde26c5a.png)

4. Further analysis

After DEP is completed, further analysis can be performed from the last tab in the right sidebar. Analysis provided by LiBT includes GSA, GSEA, KEGG pathway analysis, and PPI analysis. This section describes the options used for each analysis.

1. Gene Set Analysis (GSA)

 a. Basic analysis
 GSA uses the &quot;enrichR&quot; package and is based on a total of 4 databases:
 - GO\_Molecular\_Function\_2021
 - GO\_Cellular\_Component\_2021
 - GO\_Biological\_Process\_2021
 - KEGG\_2021\_Human
 The GSA section has two input data options and one output data option. In order to proceed with the GSA, it is first necessary to select which group to conduct the analysis on. Select Case-up protein group or Ctrl-up (Case-down) protein group and input the corresponding gene symbols. Next, users need to decide which proteins to base the GSA on. Since GSA analysis is usually performed based on DEP, the default value is selected as the DEP level. These GSA results are displayed as barplot, and the user can directly set how many barplot results to output. All plots can be directly downloaded with the download button. The file that is downloaded as a csv has the entire GSA result that is not filtered.
 ![](RackMultipart20211108-4-eljm2i_html_a129d7953de3de95.png)

b. KEGG pathway mapping
 Based on the KEGG pathway analysis result of GSA, a KEGG pathway map can be generated. It provides integrated results by downloading the pathway graph using the pathview library and mapping the log2Foldchange value of the data input by the user to the graph. The visualized result can be obtained by clicking the desired term in the table or clicking the desired term in the drop-down menu and pressing the render button. The result obtained in this way is printed at the bottom of the table, and if the zoom button is selected, the result of a large image in the modal window can be checked. If clicked at least once to check the graph, users can download it as a zip file by clicking the download button.
 ![](RackMultipart20211108-4-eljm2i_html_1e6e834d93de76a0.png)

1. Gene Set Enrichment Analysis (GSEA)

 a. Basic analysis
 GSEA is done by the &quot;fgsea&quot; package using msigDB including:
 - c2.cp.kegg.v7.4.symbols.gmt
 - c5.bp.v7.4.symbols.gmt
 - c5.cc.v7.4.symbols.gmt
 - c5.mf.v7.4.symbols.gmt
 Before performing GSEA, four gene ranking value options should be selected: log2FoldChange, P.value, P.adj, log2(FC)\*-log10 (P.adj). GSEA results are divided into UP-regulated and Down-regulated depending on whether the Normalized ES (NES) value is greater than or less than 0. A plot is drawn after being cut off based on p.adj 0.05. Also, if you click the red and blue points of the segment plot, you can check the enrichment plot of each result. All results and plots (except for enrichment plot) can be directly downloaded with the download button.

 b. KEGG pathway mapping
 GSEA also performs KEGG pathway mapping based on the KEGG pathway analysis result. Using the pathview package, upregulated and downregulated pathways are separately provided by a table. If each pathway name is selected, a modal window appears and the user can check the KEGG pathway graph to which the log2FoldChange value is mapped. If the graph is reviewed by selecting it once, the user can download it as a zip file by clicking the download button.

 ![](RackMultipart20211108-4-eljm2i_html_7aaadc8b36b26798.png) ![](RackMultipart20211108-4-eljm2i_html_80983db9d8cb1633.png) ![](RackMultipart20211108-4-eljm2i_html_42a0ae1254c2d1a8.png)
2. Protein-Protein Interaction (PPI) network analysis

 StringDB application programming interface (API) was used for PPI network analysis. The number of default input proteins is the highest value as the number of DEPs. A value higher than the number cannot be entered. If a number lower than the number of DEPs is entered, the sort is based on one selected among P.adj (default), P.value, and log2FoldChange, and the number of inputs is selected and inputted. If the user selects Protein Name, an area for text input appears, where the user enters the text separated by &quot;,&quot; as described in the example. If the Protein Name option is selected, PPI network analysis is possible even if the analysis of the previous workflow has not been performed. The resulting network can be downloaded as an image file, and network information can also be downloaded.

![](RackMultipart20211108-4-eljm2i_html_fb7335f072dfb359.jpg)

5. Download options

Once LiBT finishes preprocessing the input data, the plot results and all the table results can be directly downloaded through the download button located on each result panel.

1. Download the results table

 a) log2 transformed data matrix: In the input stage, the original data must undergo
 log2 transformation for the next analysis stage. This transformed data can be downloaded in csv format.

 b) Preprocessed data matrix: A user-defined preprocessing step is performed based on the data on which log2 transformation has been completed, and the preprocessed data can be downloaded in csv format.

 c) DEP info table: When DEP analysis is completed, p-value, fold-change, and corrected q-value can be obtained for each gene. This information can be downloaded in csv format.

 d) clustering info table: The k-means clustering result of DEP can be downloaded. k is determined according to the option setting, and protein name, accessionID, and cluster class information can be obtained as a result.

 e) GSA/GSEA result: The entire GSA/GSEA result can be downloaded. Only terms that meet the cutoff criteria are shown in the results shown in LiBT, but when downloaded to a table, the entire GSA/GSEA analysis result is downloaded. In the downloaded result, calculated values ​​such as p-value, q-value, log2 fold-change, etc. can be checked, and users can manually change the cutoff criteria to check other results in addition to the cutoff results specified by LiBT.

 f) PPI network result: Users can download connection information between proteins. It includes which proteins are linked and what the linkage score is. By adding log2FoldChange of the corresponding protein, style customization is possible in Cytoscape, etc.

1. Download the results plot
 All plots shown as a result in LiBT can be downloaded using the respective download button. PPI analysis results are downloaded in JPEG format, and all result plots except PPI are downloaded in PNG format.
