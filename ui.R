library(shinydashboard)
library(shinyIncubator)
library(shiny)
library(plotly)
library(d3heatmap)

dashboardPage(
  dashboardHeader(title = "RNA-Seq Analysis Web Interface",titleWidth = 500),
  dashboardSidebar(width = 500,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 90vh; overflow-y: auto; }" ))),
                   fluidRow(
                     column(6,uiOutput("projects")),
                     column(6,uiOutput("contrasts"))
                   ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   br(),
                   actionButton(inputId = 'rawdata', label = 'Click to view the Raw expression data'),
                   hr(),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   fluidRow(
                     column(6,radioButtons("radio", label = h4("Gene Selection"),
                                           choices = c("None" = 'none',"Upregulated" = 'up', "Downregulated" = 'down', "Both" = 'both'), 
                                           selected = 1)),
                     column(6,textInput(inputId = 'lfc', label = "Enter Fold Change cutoff", value = '2')),
                     column(6,textInput(inputId = 'apval', label = "Enter Adjusted P.Value cutoff", value = '0.05'))
                   ),
                   downloadButton('dwld','Download Data'),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   hr(),
                   h4('Generate Heatmap'),
                   fluidRow(
                     column(6,selectInput("hmip", "Select Heatmap input type",c('Enter number of genes' = "genenum",'Enter Genelist' = "geneli"))),
                     column(6,selectInput("hmpcol", "Select Heatmap Color Palette",c('RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'YlGnBu' = "YlGnBu",'PRGn'="PRGn", 'Blues' = "Blues")))
                   ),
                   fluidRow(
                     column(6,selectInput("clusterby", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none"))),
                     column(6,checkboxInput("checkbox", label = "Reverse Colors", value = FALSE))
                   ),
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   
                   conditionalPanel(
                     condition = "input.hmip == 'genenum'",
                     sliderInput("gene", label = h3("Slider"), min = 0,max = 300, value = 50)),
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
                   h4('GSEA using Camera'),
                   fluidRow(
                     column(6,actionButton(inputId = 'camera', label = 'Click to view Camera results')),
                     column(6,uiOutput("cameradd")),
                     column(6,downloadButton('downloadcam', 'Download Heatmap for Camera genes'))
                   ),
                   hr(),
                   h4('Gene Ontology Analysis'),
                   fluidRow(
                     column(6, radioButtons(inputId = 'path', label = h5("KEGG Pathway"), choices = c("Upregulated" = 'up', "Downregulated" = 'down'),selected = 1)),
                     br(),
                     column(6, textInput(inputId = 'num', label = "Top num of pathways", value = '1')),
                     column(6,actionButton(inputId = 'makeplot', label = 'View Plots')),
                     width=4
                   ),
                   hr(),
                   h5('GO Analysis using GAGE'),
                   fluidRow(
                     column(6,radioButtons(inputId='gage', label = h5("Select ontology"), 
                                           choices = c("Biological Process" = 'BP', "Cellular Component" = 'cc', "Molecular Function" = 'MF'), 
                                           selected = 1)),
                     column(6,selectInput("go_dd", "GO Selection",c('Upregulated' = "upreg",'Downregulated' = "downreg"))),
                     column(6,actionButton(inputId = 'ga', label = 'Select Ontology'))
                     #column(6,checkboxGroupInput("check", label = h3(""),choices = list("Upregulated" = 'upr', "Downregulated" = 'downr'),selected = 1))
                   ),
                   br(),
                   downloadButton('downloadGOu', 'Download Heatmap for GO genes')
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    tabsetPanel(type="tabs", id = "tabvalue",
                #tabPanel(title = "Voom Data", value = 'tab2', DT::dataTableOutput('table3')),
                #tabPanel(title = "Input Table", value = 'tab1',DT::dataTableOutput('table'),DT::dataTableOutput('table3')),
                tabPanel(title = "Project Summary and Results", h4("~~~~Project Description~~~~"),br(),value = 'tab1',textOutput("pdesc"),h4("~~~~Dot Plot of the gene of interest~~~~"),
                         fluidRow(
                           column(6,plotlyOutput('dotplot',width = 900,height = 400)),
                           column(width = 3, offset = 2,uiOutput("boxplotcol"))
                         ),
                         br(),h4("~~~~Limma data~~~~"),br(),DT::dataTableOutput('table')),
                tabPanel(title = "Raw Data", value = 'tab2',DT::dataTableOutput('table3')),
                tabPanel(title = 'Heatmap', value = 'tab4', d3heatmapOutput('heatmap',width=850,height=1150)),
                tabPanel(title = "GSEA", value = 'gsea',DT::dataTableOutput('campick3'), DT::dataTableOutput('tablecam'),
                         fluidRow(column(6,d3heatmapOutput('heatmapcam',width=750,height=1050)),
                                  column(6,selectInput("hmpcol1", "Select Heatmap Color Palette",c('RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'YlGnBu' = "YlGnBu",'PRGn'="PRGn", 'Blues' = "Blues"))),
                                  column(6,selectInput("clusterby1", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none"))),
                                  column(6,checkboxInput("checkbox1", label = "Reverse Colors", value = FALSE))
                         )
                ),
                #tabPanel(title = 'Heatmap of Selected Genes', value = 'tab7', plotOutput('heatmap2',width = 800,height = 1000)),
                tabPanel(title = 'Pathway Plot', value = 'tab5',uiOutput("plots")),
                #tabPanel(title = 'GO', value = 'tab6',uiOutput("tabs"))
                tabPanel(title = "Gene Ontology", value = 'tab6',DT::dataTableOutput('x4'),DT::dataTableOutput('table4'),
                         fluidRow(column(6,d3heatmapOutput('hmpgageup',width=750,height=1050)),
                                  column(6,selectInput("hmpcolup", "Select Heatmap Color Palette",c('RdBu' = "RdBu",'YlOrRd' = "YlOrRd",'YlGnBu' = "YlGnBu",'PRGn'="PRGn", 'Blues' = "Blues"))),
                                  column(6,selectInput("clusterbyup", "Cluster By",c('Both'="both",'Row' = "row",'Column' = "column",'None' = "none"))),
                                  column(6,checkboxInput("checkboxup", label = "Reverse Colors", value = FALSE))
                         )
                )
    )
  ))
