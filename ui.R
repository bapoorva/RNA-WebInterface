library(shinydashboard)
#library(shinyIncubator)
library(shiny)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)
library(SPIA)

dashboardPage(
  dashboardHeader(title = "NGS Data Analysis Web Interface",titleWidth = 350),
  dashboardSidebar(width = 500,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 120vh; overflow-y: auto; }" ))),
                   uiOutput("projects"),
                   fluidRow(
                     column(12,uiOutput("contrasts"))
                   ),
                   fluidRow(
                     column(6,h4("")),
                     column(6,checkboxInput("check", label = "Display multiple contrasts", value = FALSE))
                   ),
                   conditionalPanel(
                     condition = "input.check ==true",
                     uiOutput("contrastslimma")
                     ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   br(),
#                    fluidRow(
#                      column(6,
                            actionButton(inputId = 'rawdata', label = 'Click to view the Raw expression data'),br(),br(),
                     #),
                     # column(6,
                            downloadButton('rawdwld','Download Raw Data'),
# )
#                      ),
                   hr(),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 radioButtons("radio", label = h4("Gene Selection"),
                                           choices = c("None" = 'none',"Upregulated" = 'up', "Downregulated" = 'down', "Both" = 'both'),
                                           selected = 'none'),
                     fluidRow(
                       column(6,sliderInput("lfc", label = h4("Fold Change"), min = 0.5,max = 6, value = 2)),
                       column(6,sliderInput("apval", label = h4("P. Value"), min = 0.01,max = 0.2, value =0.05))
                     ),

                    checkboxInput("volcano", label = "View volcano plot", value = FALSE),
                    fluidRow(
                   column(6,downloadButton('dwld','Download results table')),
                    column(6,downloadButton('downloaddotplot', 'Download Dot plot'))
                           ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('Generate Heatmap'),
                   fluidRow(
                    #column(6,selectInput("hmip", "Select Heatmap input type",c('Top number of genes' = "genenum",'Enter Genelist' = "geneli",'Heatmap from Camera' = "hmpcam",'Heatmap from GO' = "hmpgo"))),
                    column(6,selectInput("hmip", "Select Heatmap input type",c('Top number of genes' = "genenum",'Enter Genelist' = "geneli"))),
                    
                    column(6,selectInput("hmpcol", "Select Heatmap Color Palette",c('YlGnBu' = "YlGnBu",'RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'PRGn'="PRGn", 'Blues' = "Blues")))
                   ),
                   fluidRow(
                     column(6,selectInput("clusterby", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none"))),
                     column(6,checkboxInput("checkbox", label = "Reverse Colors", value = FALSE))
                   ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                   conditionalPanel(
                     condition = "input.hmip == 'genenum'",
                     fluidRow(
                       column(6,uiOutput("dropdown")),
                       #column(6,selectInput("sortby", "Sort By",c('FDR'="sortnone",'Absolute Fold Change' = "sortab",'Positive Fold Change' = "sortpos",'Negative Fold Change' = "sortneg"))),
                       column(6,sliderInput("gene", label = h4("Top number of genes"), min = 2,max = 500, value = 50))
                     )),
                   conditionalPanel(
                     condition = "input.hmip == 'geneli'",
                     fluidRow(
                       column(6,selectInput(inputId = 'selectidentifier',label='Select Identifier',choices=list('Ensembl ID'='ensembl','Entrez Id'='entrez','Gene Symbol'='genesym'))),
                       column(6,fileInput('genelistfile', 'Upload Text File',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')))
                     )),
                   fluidRow(
                     column(6,actionButton(inputId = 'makeheat', label = 'Create Heatmap')),
                     column(6,downloadButton('downloadheatmap', 'Download'))
                   ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('Generate PCA Plot'),
                   actionButton(inputId = 'makepcaplot', label = 'Click to view PCA Plot'),
                   fluidRow(
                     column(6,checkboxInput("varpc", label = "Display Variances of PC", value = FALSE)),
                     column(6,checkboxInput("pca3d", label = "Also show 3D plot", value = FALSE))
                   ),
                   downloadButton('downloadbiplot', 'Download Biplot'),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('GSEA using Camera'),
                   fluidRow(
                     column(6,actionButton(inputId = 'camera', label = 'Click to view Camera results')),
                     column(6,uiOutput("cameradd")),
                     column(6,downloadButton('downloadcam', 'Download Camera Data')),
                     column(6,checkboxInput("eplot", label = "Show Enrichment Plot", value = FALSE))
                   ),
                   hr(),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   h4('Pathway Analysis using SPIA'),
                    fluidRow(
                      column(6,actionButton(inputId = 'runspia', label = 'Click to run SPIA')),
                      column(6,downloadButton('dwldspia', 'Download SPIA table'))
                    ),
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('GO Analysis using GAGE'),
                   fluidRow(
                     column(6,radioButtons(inputId='gage', label = h5("Select ontology"),
                                           choices = c("Biological Process" = 'BP', "Cellular Component" = 'cc', "Molecular Function" = 'MF'),
                                           selected = 1)),
                     column(6,selectInput("go_dd", "GO Selection",c('Upregulated' = "upreg",'Downregulated' = "downreg")))),
                      actionButton(inputId = 'ga', label = 'Display Results'),
                      br(),br(),
                     fluidRow( 
                     column(6,downloadButton('downloadgo', 'Download GO Data')),
                     column(6,downloadButton('downloadgogene', 'Download GO Genelist'))),
                   br()
  ),

  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tabsetPanel(type="tabs", id = "tabvalue",
                tabPanel(title = 'PCA Plot', value = 'tabpca',
                         fluidRow(
                           column(6,uiOutput("pcaxoptions")),
                           column(6,uiOutput("pcayoptions"))
                         ),
                         br(),textOutput("biplottitle"),br(),
                         fluidRow(
                           column(6,uiOutput("pcipslide")),
                           column(6,uiOutput("pcslide"))
                         ),br(),
                         fluidRow(
                           column(6,plotOutput("biplot",width=750,height=600))
                         )),
                tabPanel(title = 'Variances of PC', value = 'tabvar',h4(strong("Variances of the principal components")),textOutput("pcatitle"),plotOutput("pcaplot_ip",width=700,height=400),br(),DT::dataTableOutput('pcaplot_tab')),
                tabPanel(title = '3D PCA Plot', value = '3dpca',h4("3D plot"),br(),br(),rglwidgetOutput("pcaplot3d",width = "850px", height = "750px")),
                
                tabPanel(title = "Project Summary and Results", 
                         h4("~~~~Project Description~~~~"),br(),value = 'tab1',textOutput("pdesc"),
                         h4("~~~~Dot Plot of the gene of interest~~~~"),
                         fluidRow(
                           #column(6,checkboxInput("boxreorder", label = "Reorder x-axis", value = FALSE)),
                           column(6,uiOutput("boxplotcol"))
                         ),
                         
                         fluidRow(
                           column(6,plotOutput('dotplot',width = "auto"))
                         ),
                         br(),h4("~~~~Limma data~~~~"),
                         h5(p(div(span("Note:Please use the download button in the side panel",style="color:red")))),
                         h5(p(div(span("Note:fc - Fold Change",style="color:red")))),
                         br(),textOutput("contrdesc"),br(),DT::dataTableOutput('table')),
                tabPanel(title = "Volcano Plot", value = 'tabvolcano',
                         fluidRow(
                           column(6,plotlyOutput("volcanoplot",width=800,height=700)),
                           column(width = 3, offset = 2,uiOutput("volcdrop")),
                           column(width =4, offset = 2,uiOutput("volcslider"))
                           ),br(),DT::dataTableOutput('table_volc')),
                tabPanel(title = "Limma-Multiple Contrasts", value = 'tab11',DT::dataTableOutput('table_TRUE')),
                tabPanel(title = "Raw Data", value = 'tab2',DT::dataTableOutput('table3')),
                tabPanel(title = 'Heatmap', value = 'tab4',textOutput("htitle"),br(),
                    fluidRow(
                        column(6,uiOutput('hmplim')),
                        column(width = 3, offset = 2,plotOutput('hmpscale_out',width = 200,height = 65))
                    ),
                    fluidRow(
                    column(6,uiOutput('hmpsamp')),
                    column(6,h4(""))
                    ),
                d3heatmapOutput('heatmap',width=550,height=900)),
# tabPanel(title = 'PCA Plot', value = 'tabpca',
#          fluidRow(
#            column(6,uiOutput("pcaxoptions")),
#            column(6,uiOutput("pcayoptions"))
#          ),
#          br(),textOutput("biplottitle"),br(),
#          fluidRow(
#            column(6,uiOutput("pcipslide")),
#            column(6,uiOutput("pcslide"))
#          ),br(),
#          fluidRow(
#            column(6,plotOutput("biplot",width=800,height=650))
#          )),
# tabPanel(title = 'Variances of PC', value = 'tabvar',h4(strong("Variances of the principal components")),textOutput("pcatitle"),plotOutput("pcaplot_ip",width=700,height=400),br(),DT::dataTableOutput('pcaplot_tab')),
# tabPanel(title = '3D PCA Plot', value = '3dpca',h4("3D plot"),br(),br(),rglwidgetOutput("pcaplot3d",width = "850px", height = "750px")),

                tabPanel(title = "GSEA", value = 'gsea', 
                         fluidRow(
                           column(6,uiOutput('hmplimcam')),
                           column(width = 3, offset = 2,plotOutput('hmpscale_out2',width = 200,height = 65))
                         ),
                         fluidRow(
                           column(6,uiOutput('hmpsamp2')),
                           column(6,h4(""))
                         ),
                         d3heatmapOutput('camheatmap',width=550,height=900),
                         DT::dataTableOutput('tablecam'),textOutput("camdesc"),DT::dataTableOutput('campick3')),
                tabPanel(title = "Enrichment Plot", value = 'eplot',textOutput("eplotdesc"),br(),plotOutput('en_plot',width = 700,height = 550),br(), DT::dataTableOutput('eplottab')),
                tabPanel(title = 'Pathway Analysis using SPIA', value = 'spia', DT::dataTableOutput('spiaop'),textOutput("spiadesc"),DT::dataTableOutput('spiagenes')),
                #tabPanel(title = 'Pathway Plot', value = 'tab5',uiOutput("plots")),
                #tabPanel(title = "Gene Ontology", value = 'tab6',DT::dataTableOutput('table4')),
                tabPanel(title = "Gene Ontology", value = 'tab6',
                         fluidRow(
                           column(6,uiOutput('hmplimgo')),
                           column(width = 3, offset = 2,plotOutput('hmpscale_out3',width = 200,height = 65))
                         ),
                         fluidRow(
                           column(6,uiOutput('hmpsamp3')),
                           column(6,h4(""))
                         ),
                         d3heatmapOutput('goheatmap',width=550,height=900),
                         DT::dataTableOutput('table4'),textOutput("godesc"),DT::dataTableOutput('x4')),
                tabPanel(title = "Help Page", value='tab10', 
                         h4(p(strong("1. Project Summary and Results"))),
                         h4("Select a project and a comparison. (Comparisons are automatically populated in the drop-down menu)"),
                            h4(p(div(span("The Project Summary and Results tab",style="color:blue"), "will display the limma (differential expression analysis) output for that comparison in a  table. Clicking on any row will display the dot plot for that gene. Select an Attribute from the drop-dowm menu near the dot-plot to color the plot by that feature."))),
                         h4(p(div(span("Select ",em("Display multiple contrasts"),"checkbox to view the foldchange and adjusted pvalues of multiple comparisons")))),
                         br(),
                         h4(p(strong("2. Raw Expression Data"))),
                         h4("Click on the", em("Click to view the Raw expression data"),"button to view the expression data in the",span("Raw Data tab",style="color:blue")),
                         br(),
                         h4(p(strong("3. Gene Selection"))),
                         h4("Make a Gene Selection by selecting the radio button to view the list of upregulated and/or downregulated genes in the ", span("The Project Summary and Results tab",style="color:blue")),   
                         h4(p(div("Type in the Fold Change cutoff and Adjusted PValue cutoff and view the updated table in the same tab. Click on",em("Download Data"),"button to download the table as a csv file"))),
                         h4(p(div(span("Note:Make sure the radio button 'None' is not selected when setting FC and P.Value cutoffs",style="color:red")))),
                         br(),
                         h4(p(strong("4. Heatmap"))),
                         h4("Select Heatmap type from drop-down menu and give appropriate inputs to view the Heatmap in the",span("Heatmap tab",style="color:blue")),
                         h4(p(div("Select color using the",em("Heatmap color palette"),"dropdown and reverse color palette using the",em("Reverse Colors"),"checkbox. Cluster by Rows and/or columns or none by selecting the option from the",em("Cluster By"),"drop-down menu"))),
                         h4(p(div("Use the",em("Slider"),"to select number of genes if Heatmap type selected is 'Top number of genes'.Default is 50, minimum is 2 and maximum is 300"))),
                         h4(p(div(span("Note:If the limma table has fewer genes, the heatmap will display only those despite the slider value",style="color:red")))),
                         h4(p(div("Enter genelist separated by commas (no spaces) in the",em("Enter Gene List to plot"),"text box if the Heatmap type selected is 'Enter Genelist'."))),
                         h4(p(div(span("Note:Make sure you select the correct identifier from the drop down menu (ENSEMBL ID, ENTREZ ID, Gene Symbol)",style="color:red")))),
                         h4(p(div("Generate Camera and GO data to view the heatmap if the Heatmap type selected is 'Heatmap from Camera' or 'Heatmap from GO'."))),
                         br(),
                         h4(p(strong("5. PCA Plot"))),
                         h4(p(div("Click on",strong("Click to view PCA Plot"),"button to view the biplot in the",span("PCA Plot tab",style="color:blue"),"You can select the principle component to plot on the x and y axis of the plot from the drop-down menu. You can also specify the number of top genes showing maximum variance to be used as the input for the bioplot as well as the number of genes you want to view in the plot. Select ",em("Display variances of PC"),"to view the barplot showing the proportion of variance retained by each principle component.","Select ",em("Also show 3D plot"),"to view the 3D plot of the top 3 principle components"))),
                         h4("Helpful Links:", a("Click Here for information on PCA biplot", href="http://www.nature.com/nbt/journal/v26/n3/full/nbt0308-303.html")),
                         br(),
                         h4(p(strong("6. GSEA using Camera"))),
                         h4("Click on",strong("Click to view Camera results"),"button to view Camera results in the",span("GSEA tab",style="color:blue")),
                         h4(p(div("Select a gene set from the ",strong("Select a Gene Set")," dropdown"))),
                         h4(p(div("The Camera function in the limma package for testing differential expression, tests whether a set of genes is highly ranked relative to other genes in terms of differential expression. It takes into account the inter-gene correlation.CAMERA, an acronym for Correlation Adjusted MEan RAnk gene set test, is based on the idea of estimating the variance inflation factor associated with inter-gene correlation, and incorporating this into parametric or rank-based test procedures. It returns the number of genes in the set, the inter-gene correlation value, the direction of change (Up or Down), the two-tailed p-value and the Benjamini & Hochberg FDR adjusted P-value"))),
                         h4(p(div(span("The GSEA tab",style="color:blue"), "will display the Camera output in a table. Clicking on any row will display the gene list from the user dataset that belongs to that Gene set category in table below it. Select",em('Heatmap from Camera'),"from the Heatmap input drop-down menu to view the Heatmap for those genes in the Heatmap tab"))),
                         h4(p(div("Click on",em("Download Camera Data"),"button to download the table as a csv file"))),
                         h4("Helpful Links:", a("Click Here for for information on Camera", href="http://nar.oxfordjournals.org/content/early/2012/05/24/nar.gks461.full")),
                         br(),
                         h4(p(strong("7. Pathway Analysis using SPIA"))),
                         h4("Click on the",em("Click to run SPIA"),"button to view the the results from SPIA"),
                         h4("Helpful Links: Click", a("here", href="http://www.bioconductor.org/packages/release/bioc/vignettes/SPIA/inst/doc/SPIA.pdf"),"and",a("here", href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1987343/"),"for information on SPIA"),
            
                         br(),
                         h4(p(strong("8. Gene Ontology using GAGE"))),
                         h4("Select the Ontology and upregulated/downregulated from the drop down menu and click on the",strong("Select Ontology"),"button to view the Gene Ontology results and the genes corresponding to each GO-term in the",span("Gene Ontology tab",style="color:blue")),
                         h4(p(div(span("The Gene Ontology tab",style="color:blue"), "will display the GO output for that ontology in a table. Clicking on any row will display the gene list from the user dataset that belong to that GO-term in table below it. Select",strong('Heatmap from GO'),"from the Heatmap input drop-down menu to view the Heatmap for those genes in the Heatmap tab"))),
                         
                         h4(p(div("Click on",strong("Download GO Data"),"button to download the table as a csv file"))),
                         h4("Helpful Links:", a("Click Here for information on GAGE", href="http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-161"))
                         )
    )
  ))
