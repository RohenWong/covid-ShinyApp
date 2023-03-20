# Multi-Omics COVID Shiny app

I have built a Shiny app ([hosted here](https://rwon.shinyapps.io/Multi-Omics-COVID-risk-app/)) which allows users to upload their multi-omics data (proteomic, metalobomic, transcriptomic readings), to get a predicted covid-19 outcome. 

To train the underlying machine learning models, I integrated multi-omics data from these papers:

- Overmyer et al. (2021). Large-Scale Multi-omic Analysis of COVID-19 Severity. Cell Systems, 12(1), 23-40. DOI: https://doi.org/10.1016/j.cels.2020.10.003

- Patel et al. (2021). Proteomic Blood Profiling in Mild, Severe, and Critical COVID-19 Patents. Sci Rep 11, 6357 (2021). DOI: https://doi.org/10.1038/s41598-021-85877-0

- Shen et al. (2020). Proteomic and Metabolomic Characterization of COVID-19 Patient Sera. Cell, 182(1), 59-72. DOI: https://doi.org/10.1016/j.cell.2020.05.032

- Su et al. (2020). Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19 . Cell 183(6), 1479–1495. DOI: https://doi.org/10.1016/j.cell.2020.10.037



The folder `demo_uploading_files` contains 2 demo datasets - one from Su et al. (2020), another from Shen et al. (2020). The clinical data and `.csv` omics from each can be uploaded to the app to get predictions. 

*Note:* pdata, mdata and rdata refer to proteomics, metabolomics and transcriptomics respectively.

