# PRONA (Patient Reported Outcomes Network Analysis)
 
PRONA allows users to perform symptom network analysis using patient reported outcomes data. PRONA not only consolidates previous NA tools into a unified, easy-to-use analysis pipeline, but also augments the traditional approach with functionality for performing unsupervised discovery of patient subgroups with distinct symptom patterns.

To install the PRONA package in R, use the following command:

```
if(!require("devtools")) install.packages("devtools")
library(devtools)
devtools::install_github("bbergsneider/PRONA")
```

For how to use PRONA, please follow these vignettes:
1. [Introduction](https://rpubs.com/brandonbergs/prona-introduction)
2. [Network accuracy and stability assessment](https://rpubs.com/brandonbergs/PRONA-statistical-assessment)
3. [Concordance-based clustering of patient communities](https://rpubs.com/brandonbergs/prona-unsupervised-clustering)

For more details, please see our publication in Bioinformatics overviewing the PRONA package: [Bergsneider & Celiku, Bioinformatics, 2024](https://pubmed.ncbi.nlm.nih.gov/39520406/).

For use case examples, please see [Bergsneider et al, Neuro-Oncology Advances, 2022](https://pubmed.ncbi.nlm.nih.gov/36820236/) and [Bergsneider et al, Cancer Medicine, 2024](https://pubmed.ncbi.nlm.nih.gov/39377555/).
