# NpDynamic  
This repository contains data for the Np dynamic reconstruction in _Drosophila melanogaster_ with SMALT, a substitution mutation aided lineage tracing system.  
## Summary
Mapping the cell phylogeny of a complex multicellular organism relies on somatic mutations accumulated from zygote to adult.  Available cell barcoding methods can only record no more than four mutations per barcode, enabling a low-resolution mapping of the cell phylogeny of a complex organism.  We here developed SMALT, a substitution mutation aided lineage tracing system that outperforms the available cell barcoding methods in mapping cell phylogeny.  We applied SMALT to Drosophila melanogaster and obtained on average over 20 mutations on a three kilobase-pair barcoding sequence in early-adult cells.  Using the barcoding mutations we obtained in two fly individuals unprecedentedly high quality cell phylogenetic trees, each comprising several thousand internal nodes with 84-93% median bootstrap supports.  The highly informative cell phylogenies enabled a population genetic reconstruction of the developmental demographic history.  We estimated in each organ the longitudinal dynamics of the number of actively dividing parental cells (Np) through the development.  The Np dynamics informed the symmetric versus asymmetric cell division balance.  It also suggested a sigmoid function of organ growth, in which the carrying capacity of each organ imposes a negative feedback regulation on the rate of cell births to ensure the reproducibility of the final organ size.  
## Content
* Tree data
1. Phylogenetic trees of two individuals constructed with iqtree.
2. Ancestral states reconstruction by empirical Bayes with _-asr_ parameter.
* Processed data
1. Np estimation of two individual included in PackagedData.RData.
2. Simulated phylogenies under different population history.
## Reference
Yao Z, Liu K, Deng S, et al. An instantaneous coalescent method insensitive to population structure[J]. Journal of Genetics and Genomics, 2021. doi: 10.1016/j.jgg.2021.03.005
