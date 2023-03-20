# Multi-Omics COVID Shiny app

I have built a Shiny app ([hosted here](https://rwon.shinyapps.io/Multi-Omics-COVID-risk-app/)) to predict Covid-19 outcome. User's can interact with our visualisation tab, or upload their own omics data (proteomic, metalobomic, transcriptomic readings) to get a predicted outcome (moderate/severe).

Predictions are made using majority voting from an ensemble of machine learning models: Random Forest, Support Vector Machines (SVM), non-linear SVM and Diagonal Linear Discrimant Analysis (DLDA).

To train the models, I integrated multi-omics data from these papers:

- Overmyer et al. (2021). Large-Scale Multi-omic Analysis of COVID-19 Severity. Cell Systems, 12(1), 23-40. [doi](https://doi.org/10.1016/j.cels.2020.10.003) 

- Patel et al. (2021). Proteomic Blood Profiling in Mild, Severe, and Critical COVID-19 Patents. Sci Rep 11, 6357 (2021). [doi](https://doi.org/10.1038/s41598-021-85877-0)

- Shen et al. (2020). Proteomic and Metabolomic Characterization of COVID-19 Patient Sera. Cell, 182(1), 59-72. [doi](https://doi.org/10.1016/j.cell.2020.05.032)

- Su et al. (2020). Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19 . Cell 183(6), 1479–1495. [doi](https://doi.org/10.1016/j.cell.2020.10.037)


## Visualisation Tab

We explore several metrics of differential expression (log Fold-Change, moderate t-statistics etc.). A feature with higher differential expression, has 'larger' discrepancy between healthy and severe patients -  suggesting it is more important for modelling.

This tab also includes interactive boxplots and PCA plots.

https://user-images.githubusercontent.com/115207897/226347837-7667acee-40ee-439c-b593-2ffe60340495.mov

## Omics Upload Tab

Users can upload their omics data to get predicted covid-19 outcomes (moderate/severe).

The folder `demo_uploading_files` contains 2 demo datasets - one from Su et al. (2020), another from Shen et al. (2020). The clinical data and `.csv` omics from each can be uploaded to the app to get predictions. 

*Note:* pdata, mdata and rdata refer to proteomics, metabolomics and transcriptomics respectively.



