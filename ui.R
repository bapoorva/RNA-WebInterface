library(shinydashboard)
library(shinyIncubator)
library(shiny)
library(plotly)
library(d3heatmap)
library(shinyjs)
library(rglwidget)
library(SPIA)

dashboardPage(
  dashboardHeader(title = "RNA-Seq Analysis Web Interface",titleWidth = 500),
  dashboardSidebar(width = 500,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 120vh; overflow-y: auto; }" ))),
                   fluidRow(
                     column(6,uiOutput("projects")),
                     column(6,uiOutput("contrasts"))
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
                   actionButton(inputId = 'rawdata', label = 'Click to view the Raw expression data'),
                   hr(),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   fluidRow(
                     column(6,radioButtons("radio", label = h4("Gene Selection"),
                                           choices = c("None" = 'none',"Upregulated" = 'up', "Downregulated" = 'down', "Both" = 'both'),
                                           selected = 'none')),
                     column(6,textInput(inputId = 'lfc', label = "Enter Fold Change cutoff", value = '2')),
                     column(6,textInput(inputId = 'apval', label = "Enter Adjusted P.Value cutoff", value = '0.05'))
                   ),
                   downloadButton('dwld','Download Data'),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('Generate Heatmap'),
                   fluidRow(
                     column(6,selectInput("hmip", "Select Heatmap input type",c('Top number of genes' = "genenum",'Enter Genelist' = "geneli",'Heatmap from Camera' = "hmpcam",'Heatmap from GO' = "hmpgo"))),
                     column(6,selectInput("hmpcol", "Select Heatmap Color Palette",c('RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'YlGnBu' = "YlGnBu",'PRGn'="PRGn", 'Blues' = "Blues")))
                   ),
                   fluidRow(
                     column(6,selectInput("clusterby", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none"))),
                     column(6,checkboxInput("checkbox", label = "Reverse Colors", value = FALSE))
                   ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                   conditionalPanel(
                     condition = "input.hmip == 'genenum'",
                     sliderInput("gene", label = h3("Slider"), min = 2,max = 300, value = 50)),
                   conditionalPanel(
                     condition = "input.hmip == 'geneli'",
                     fluidRow(
                       column(6,selectInput(inputId = 'selectidentifier',label='Select Identifier',choices=list('Ensembl ID'='ensembl','Entrez Id'='entrez','Gene Symbol'='genesym'))),
                       column(6,textInput("genelist", label = h5("Enter Gene List to plot"), value = ""))
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
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('GSEA using Camera'),
                   fluidRow(
                     column(6,actionButton(inputId = 'camera', label = 'Click to view Camera results')),
                     column(6,uiOutput("cameradd")),
                     column(6,downloadButton('downloadcam', 'Download Camera Data'))
                   ),
                   hr(),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                    br(),
#                    actionButton(inputId = 'runspia', label = 'Click to run SPIA'),
#                    hr(),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   h4('Pathway Analysis using SPIA'),
                   actionButton(inputId = 'runspia', label = 'Click to run SPIA'),
#                    fluidRow(
#                      column(6, radioButtons(inputId = 'path', label = h5("KEGG Pathway"), choices = c("Upregulated" = 'up', "Downregulated" = 'down'),selected = 1)),
#                      br(),
#                      column(6, textInput(inputId = 'num', label = "Top num of pathways", value = '10')),
#                      column(6,actionButton(inputId = 'makeplot', label = 'View KEGG Results')),
#                      width=4
#                    ),
                   hr(),
                   h4('GO Analysis using GAGE'),
                   fluidRow(
                     column(6,radioButtons(inputId='gage', label = h5("Select ontology"),
                                           choices = c("Biological Process" = 'BP', "Cellular Component" = 'cc', "Molecular Function" = 'MF'),
                                           selected = 1)),
                     column(6,selectInput("go_dd", "GO Selection",c('Upregulated' = "upreg",'Downregulated' = "downreg"))),
                     column(6,actionButton(inputId = 'ga', label = 'Select Ontology'))),
                   downloadButton('downloadgo', 'Download GO Data'),
                   br()
  ),

  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
    tabsetPanel(type="tabs", id = "tabvalue",
                tabPanel(title = "Project Summary and Results", h4("~~~~Project Description~~~~"),br(),value = 'tab1',textOutput("pdesc"),h4("~~~~Dot Plot of the gene of interest~~~~"),
                         fluidRow(
                           column(6,plotOutput('dotplot',width = 900,height = 400)),
                           column(width = 3, offset = 2,uiOutput("boxplotcol"))
                         ),
                         br(),h4("~~~~Limma data~~~~"),
                         h5(p(div(span("Note:Please use the download button in the side panel",style="color:red")))),
                         br(),DT::dataTableOutput('table')),
                #tabPanel(title="MultiContrast-Limma",uiOutput("plotUI")),
                tabPanel(title = "MultiContrast-Limma", value = 'tab11',DT::dataTableOutput('table_TRUE')),
                tabPanel(title = "Raw Data", value = 'tab2',DT::dataTableOutput('table3')),
                tabPanel(title = 'Heatmap', value = 'tab4',textOutput("htitle"), d3heatmapOutput('heatmap',width=850,height=1150)),
                tabPanel(title = 'PCA Plot', value = 'tabpca',
                         fluidRow(
                           column(6,uiOutput("pcaxoptions")),
                           column(6,uiOutput("pcayoptions"))
                         ),
                         br(),textOutput("biplottitle"),br(),
                         fluidRow(
                           column(6,plotOutput("biplot",width=900,height=850)),
                           column(width = 3, offset =2,uiOutput("pcipslide")),
                           column(width = 3, offset =2,uiOutput("pcslide"))
                )),
                tabPanel(title = 'Variances of PC', value = 'tabvar',h4(strong("Variances of the principal components")),textOutput("pcatitle"),plotOutput("pcaplot_ip"),br(),DT::dataTableOutput('pcaplot_tab')),
                tabPanel(title = '3D PCA Plot', value = '3dpca',h4("3D plot"),br(),br(),rglwidgetOutput("pcaplot3d",width = "850px", height = "750px")),
                tabPanel(title = "GSEA", value = 'gsea', DT::dataTableOutput('tablecam'),textOutput("camdesc"),DT::dataTableOutput('campick3')),
                #tabPanel(title = 'KEGG Pathway', value = 'tab5', DT::dataTableOutput('kegg'),DT::dataTableOutput('kegggenes')),
                tabPanel(title = 'Pathway Analysis using SPIA', value = 'spia', DT::dataTableOutput('spiaop')),
                #tabPanel(title = 'Pathway Plot', value = 'tab5',uiOutput("plots")),
                tabPanel(title = "Gene Ontology", value = 'tab6',DT::dataTableOutput('table4'),textOutput("godesc"),DT::dataTableOutput('x4')),
                tabPanel(title = "Help Page", value='tab10', 
                         h4("1. Select a project and a comparison. (Comparisons are automatically populated in the drop-down menu)"),
                            h4(p(div(span("The Project Summary and Results tab",style="color:blue"), "will display the limma (differential expression analysis) output for that comparison in a  table. Clicking on any row will display the dot plot for that gene. Select an Attribute from the drop-dowm menu near the dot-plot to color the plot by that feature.")))
                         ,br(),
                         h4("2. Click on the", strong("Click to view the Raw expression data"),"button to view the expression data in the",span("Raw Data tab",style="color:blue")),
                         br(),
                         h4("3. Make a Gene Selection by selecting the radio button to view the list of upregulated and/or downregulated genes in the ", span("The Project Summary and Results tab",style="color:blue")),   
                         h4(p(div("Type in the Fold Change cutoff and Adjusted PValue cutoff and view the updated table in the same tab. Click on",strong("Download Data"),"button to download the table as a csv file"))),
                         h4(p(div(span("Note:Make sure the radio button 'None' is not selected when setting FC and P.Value cutoffs",style="color:red")))),
                         br(),
                         h4("4. Select Heatmap type from drop-down menu and give appropriate inputs to view the Heatmap in the",span("Heatmap tab",style="color:blue")),
                         h4(p(div("Select color using the",strong("Heatmap color palette"),"dropdown and reverse color palette using the",strong("Reverse Colors"),"checkbox. Cluster by Rows and/or columns or none by selecting the option from the",strong("Cluster By"),"drop-down menu"))),
                         h4(p(div("Use the",strong("Slider"),"to select number of genes if Heatmap type selected is 'Top number of genes'.Default is 50, minimum is 2 and maximum is 300"))),
                         h4(p(div(span("Note:If the limma table has fewer genes, the heatmap will display only those despite the slider value",style="color:red")))),
                         h4(p(div("Enter genelist separated by commas (no spaces) in the",strong("Enter Gene List to plot"),"text box if the Heatmap type selected is 'Enter Genelist'."))),
                         h4(p(div(span("Note:Make sure you select the correct identifier from the drop down menu (ENSEMBL ID, ENTREZ ID, Gene Symbol)",style="color:red")))),
                         h4(p(div("Generate Camera and GO data to view the heatmap if the Heatmap type selected is 'Heatmap from Camera' or 'Heatmap from GO'."))),
                         br(),
                         h4("5. Click on",strong("Click to view Camera results"),"button to view Camera results in the",span("GSEA tab",style="color:blue")),
                         h4(p(div("Select a gene set from the ",strong("Select a Gene Set")," dropdown"))),
                         h4(p(div("The Camera function in the limma package for testing differential expression, tests whether a set of genes is highly ranked relative to other genes in terms of differential expression. It takes into account the inter-gene correlation.CAMERA, an acronym for Correlation Adjusted MEan RAnk gene set test, is based on the idea of estimating the variance inflation factor associated with inter-gene correlation, and incorporating this into parametric or rank-based test procedures. It returns the number of genes in the set, the inter-gene correlation value, the direction of change (Up or Down), the two-tailed p-value and the Benjamini & Hochberg FDR adjusted P-value"))),
                         h4(p(div(span("The GSEA tab",style="color:blue"), "will display the Camera output in a table. Clicking on any row will display the gene list from the user dataset that belongs to that Gene set category in table below it. Select",strong('Heatmap from Camera'),"from the Heatmap input drop-down menu to view the Heatmap for those genes in the Heatmap tab"))),
                         
                         h4(p(div("Click on",strong("Download Camera Data"),"button to download the table as a csv file"))),
                         br(),
                         h4("6. Select the upregulated/downregulated radio button for",strong("KEGG Pathway"),"and enter the number of pathways in the",strong("Top num of pathways"),"text box to view the KEGG pathway plots in the in the",span("Pathway Plot tab",style="color:blue")),
                         br(),
                         h4("7. Select the Ontology and upregulated/downregulated from the drop down menu and click on the",strong("Select Ontology"),"button to view the Gene Ontology results and the genes corresponding to each GO-term in the",span("Gene Ontology tab",style="color:blue")),
                         h4(p(div(span("The Gene Ontology tab",style="color:blue"), "will display the GO output for that ontology in a table. Clicking on any row will display the gene list from the user dataset that belong to that GO-term in table below it. Select",strong('Heatmap from GO'),"from the Heatmap input drop-down menu to view the Heatmap for those genes in the Heatmap tab"))),
                         #h4(p(div("Select a gene set from the ",strong("Select a Gene Set")," dropdown"))),
                         h4(p(div("Click on",strong("Download GO Data"),"button to download the table as a csv file")))
                         )
    )
  ))
