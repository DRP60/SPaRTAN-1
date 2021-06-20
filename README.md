# SPaRTAN
Signal-regulated transcription factors (TFs) underlie key developmental or differentiation transitions and activation states of cells (e.g., within the immune system). To date, single-cell genomics datasets have not been used to link cell-surface receptors to TFs of individual cells. To address this need, we developed a machine learning approach called **SPaRTAN (Single-cell Proteomic and RNA based Transcription factor Activity Network)** that integrates parallel surface protein measurements with mRNA expression data in cells based on CITE-seq (cellular indexing of transcriptomes and epitopes by sequencing) data with regulatory genomics resources.

<img src="https://github.com/osmanbeyoglulab/PyAffreg/blob/master/data/diagram.png" width="500">

Briefly, our model views the expression of surface proteins as a proxy of their activities; signaling emanating from these proteins converges on particular TFs, whose activities, in turn, regulate the expression of their target genes. Specifically, we use a regularized bilinear regression algorithm called affinity regression (AR), a general statistical framework for any problem where the observed data can be explained as interactions (**W**) between two kinds of inputs, to establish an interaction matrix between surface proteins (**P**) and TFs (**D**)  that predicts target gene expression (**Y**).

Since the model captures statistical relationships between surface proteins, TFs, and gene expression. We can use the trained interaction matrix to obtain different views of a CITE-seq data set; e.g., to predict TF activity from a cell's surface protein expression profile or to predict surface protein expression from a cell’s gene expression profile.  Intuitively, information flows down from observed surface protein levels through the learned interaction matrix to infer TF activities and observed mRNA expression levels or propagates up through the TF target-gene edges and interaction network to infer surface protein expression.

There are currently implementations of SPaRTAN in MATLAB, and in Python. The output from any of the implementations can then be explored either in R or Python.
Since Python is an open souce package, no license needed, here we give detailed instructions of installation and usage of python implementation. Please see the tutorials
* [Run SPaRTAN in Python](https://github.com/osmanbeyoglulab/SPaRTAN/tree/main/SPaRTAN_python)

## References
Ma X*, Somasundaram A*, Qi Z, Singh H, Osmanbeyoglu HU, SPaRTAN, a computational framework for linking cell-surface receptors to transcriptional regulators. bioRxiv 2020.12.22.423961: doi:https://doi.org/10.1101/2020.12.22.423961

Pelossof, R., Singh, I., Yang, J.L., Weirauch, M.T., Hughes, T.R. and Leslie, C.S. (2015) Affinity regression predicts the recognition code of nucleic acid-binding proteins. Nat Biotechnol, 33, 1242-1249.

Osmanbeyoglu, H.U., Pelossof, R., Bromberg, J.F. and Leslie, C.S. (2014) Linking signaling pathways to transcriptional programs in breast cancer. Genome Res, 24, 1869-1880.

Osmanbeyoglu, H.U., Toska, E., Chan, C., Baselga, J. and Leslie, C.S. (2017) Pancancer modelling predicts the context-specific impact of somatic mutations on transcriptional programs. Nat Commun, 8, 14249.

Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. and Saez-Rodriguez, J. (2019) Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res, 29, 1363-1375.


