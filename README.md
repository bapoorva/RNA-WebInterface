# RNA-WebInterface
GUI to explore data and results from RNA-Seq and Microarray Experiments

##Project Summary and Results
Selecting a Project automatically populates the drop-down menu for the different contrasts/comparisons. The limma results for the selected comparison is displayed in the Project Summary and Results tab. Clicking on any row in the table shows the dot-plot of expression values corresponding to the selected gene

##View Raw Expression Data
Clicking on the view Raw Expression Data button displays the expression data in the Raw Data Tab

##Gene Selection
Subselect upregulated/downregulated genes by selecting the radio button and  entering Fold Change and Adjusted P-value cut-off. Table will be updated in the Project summary tab

##Generate Heatmaps
There are two ways to generate heatmaps
1.Enter number of genes : The genes are arranges ordered based on their Adjusted P.Values. Use the slider to pick top number of genes. Default is 50. Maximum is 300
2. Enter Genelist : Select identifier as either ENSEMBL ID, ENTREZ ID or Gene Symbol and enter the gene list separated by a comma.
Clustering options as well as options to change and rever color pallets provided.
3. Generate Camera and GO data to view the heatmap if the Heatmap type selected is 'Heatmap from Camera' or 'Heatmap from GO'.

## PCA Plot

Click on Click to view PCA Plot button to view the biplot in the PCA Plot tab You can select the principle component to plot on the x and y axis of the plot from the drop-down menu. You can also specify the number of top genes showing maximum variance to be used as the input for the bioplot as well as the number of genes you want to view in the plot. Select Display variances of PC to view the barplot showing the proportion of variance retained by each principle component. Select Also show 3D plot to view the 3D plot of the top 3 principle components

##GSEA Using Camera
Camera results are precomputed. Selection of project and contrast populated the  gene set drop-down menu. Results are displayed in the GSEA tab. Selecting a row from the camera results table displays the genelist corresponding to the selection and a heatmap of those genes

##Pathway Analysis using SPIA
Click on the Click to run SPIA button to view the the results from SPIA

##GO Analysis using GAGE
Select GO on interest (BP,MF,CC) and either upregulated/downregulated to display top GO in the Gene Ontology tab as well as genes corresponding to each ontology and heatmap
