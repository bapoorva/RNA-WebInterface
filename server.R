library(shiny)
library(pathview)
library(gage)
library(gageData)
library(RColorBrewer)
library(NMF)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(gplots)
library(KEGGgraph)
library(Biobase)
library(reshape2)
library(ggplot2)
library(limma)
library(edgeR)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(biomaRt)
library(KEGGREST)
library(png)
library(GO.db)
library(d3heatmap)
library(dplyr)
library(plotly)
#load entrez id's for kegg pathway
data(kegg.sets.mm)
#load indexes for signaling and metabolic pathways
data(sigmet.idx.mm)
#get entrez id's for signaling and metabolic pathways
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]

data(go.sets.mm)
data(go.subs.mm)

#Create atheme for all plots. 
plotTheme <-theme_bw() + theme(axis.title.x = element_text(face="bold", size=16),
                               axis.text.x  = element_text(angle=0, vjust=0.5, size=14),
                               axis.title.y = element_text(face="bold", size=16),
                               axis.text.y  = element_text(angle=0, vjust=0.5, size=14))




# Read input file
shinyServer(function(input, output,session) {
  #~~~~~~~~~~~~~~global variables~~~~~~~~~~~#
  keggresids=as.character()
  pData=data.frame()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###################################################
  ###################################################
  ####### LOAD EXCEL AND POPULATE DROP DOWN #########
  ###################################################
  ###################################################
  
  #Read the parameter file
  readexcel = reactive({
    file = read.csv("data/param.csv")
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    prj=excel$projects
      selectInput("projects","Select a project",as.list(as.character(prj)))
  })
  
  #Load Rdata
  fileload <- reactive({
        inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
        load(inFile)
        loaddata=results
        return(loaddata)
        })
  
  #Get contrast list and populate drop-down
  output$contrasts = renderUI({
    results=fileload()
    lim=results$limma
    contrasts=as.list(as.character(unlist(lapply((names(lim)),factor))))
    selectInput("contrast","Select a comparison",contrasts,"pick one")
  })
  
  
  ###################################################
  ###################################################
  ######## GET PROJECT DESC AND DISPLAY   ###########
  ###################################################
  ###################################################
  #Read parameter file and get project desc for the project selected
  prjdesc = reactive({
    file = readexcel()
    prj=input$projects
    desc=file[file$projects %in% prj,-1]
    desc=as.character(desc)
  })
  
  #Display text in main panel
  output$pdesc <- renderText({ 
    desc=prjdesc() 
  })
  
  # Update tab1 in main panel 
  observe({
    updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab1') 
  })
  
###################################################
###################################################
###########LOAD LIMMA FILE AND DISPLAY#############
###################################################
###################################################
  #Read limma data from eset
  datasetInput0.5 = reactive({
    contrast=input$contrast
    results=fileload()
    k=paste('results$limma$',contrast,sep='')
    limmadata=eval(parse(text = k))
  })
  
  #Update datatable in tab 1 based on gene selection (upregulated, downregulated, both or none)
  datasetInput = reactive({
    contrast=input$contrast #select contrast
    limmadata=datasetInput0.5()
    lfc=as.numeric(input$lfc) #get logFC
    apval=as.numeric(input$apval)#get adjusted P.Vals
    if(is.null(input$radio))
    {
      d = limmadata
    }
    else if(input$radio=='none')
    {
      d=limmadata
    }
    else if(input$radio=='down')
    {
      d=limmadata
      d = d[which(d$logFC < (-1*(lfc)) & d$adj.P.Val < apval),]
    }
    else if(input$radio=='up')
    {
      d=limmadata
      d = d[which(d$logFC>lfc & d$adj.P.Val < apval),]
    }
    else if(input$radio=='both')
    {
      d=limmadata
      d = d[which(abs(d$logFC) > lfc & d$adj.P.Val < apval),]
    }
    #limma = as.data.frame(d) #load limma data
    geneid=d$SYMBOL
    url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
    if(url=="http://www.genecards.org/cgi-bin/carddisp.pl?gene="){
      d$link<-NULL
    }else{
    d$link=paste0("<a href='",url,"'>","Link to GeneCard","</a>")}
    #rownames(d)=make.names(d$id,unique = TRUE)
    d=as.data.frame(d)
    return(d)
  })
  
  #print input file in tab1
  output$table = DT::renderDataTable({
      DT::datatable(datasetInput(),
                    extensions = c('TableTools','ColVis','Scroller'),
                    options = list(
                      searchHighlight = TRUE,
                      pageLength = 10,
                      lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                      scrollX = TRUE
                    ),rownames=FALSE,selection = list(mode = 'single', selected =1),escape=FALSE)
    })

  # display data on tab 1
  observe({
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab1') 
  })
  
  #download data as excel sheet
  output$dwld <- downloadHandler(
    filename = function() { paste(input$projects, '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file)
    })
  
  ###################################################
  ###################################################
  ############## DISPLAY VOOM DATA ##################
  ###################################################
  ###################################################
  #load voom data from eset
  datasetInput3 = reactive({
    results=fileload()
    exprsdata=exprs(results$eset)
    #exprsdata=data.frame(id=rownames(exprsdata),exprsdata)
  })
   
  datasetInput33 = reactive({
    results=fileload()
    exprsdata=as.data.frame(exprs(results$eset))
    features=as.data.frame(pData(featureData(results$eset)))
    features$id=rownames(features)
    exprsdata$id=rownames(exprsdata)
    #exprsdata=data.frame(id=rownames(exprsdata),exprsdata)
    genes <- inner_join(exprsdata,features,by=c('id'='id'))
    #genes <- inner_join(exprsdata,features,by=c('id'='id'))
    #rownames(genes)=make.names(genes$id,unique = TRUE)
    return(genes)
  })
  
  #print voom or expression data file in tab2
  output$table3 = DT::renderDataTable({
    DT::datatable(datasetInput33(),
                  extensions = c('TableTools','ColVis','Scroller'),
                  options = list(
                    searchHighlight = TRUE,
                    pageLength = 10,
                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                    scrollX = TRUE
                  ),rownames=FALSE,caption= "Voom data")
  })
  # display data in tab2
  observe({
    if(input$rawdata>0) {
    updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab2') }
  })

  
  ###################################################
  ###################################################
  ##################### DOT PLOT ####################
  ###################################################
  ################################################### 
  output$boxplotcol = renderUI({ #create user interface for dot-plot drop down menu
    results=fileload()
    eset=results$eset
    pData=pData(eset) #get pheno-data
    kc=pData[ , grepl( "var_" , colnames(pData) ) ] #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(kc)))) #omit columns that have only NA's
    kc=data.frame(maineffect=pData$maineffect,sample_name=pData$sample_name,kt) #create a dataframe with maineffect and sample_name and non-null columns starting with var
    bpcols=as.list(as.character(unlist(lapply((colnames(kc)),factor))))
    selectInput("color","Select a comparison",bpcols) #populate drop down menu with the phenodata columns
  })

  dotplot_out = reactive({
    s = input$table_rows_selected #select rows from table
    dt = datasetInput() #load limma data
    validate(
      need((is.data.frame(dt) && nrow(dt))!=0, "No data in table")
    )
    dt1 = dt[s, , drop=FALSE]
    id = dt[s,1] #get limma data corresponding to selected row in table
    results=fileload()
    eset <- results$eset
    e <-data.frame(eset@phenoData@data,signal=exprs(eset)[id,])
    genesymbol=dt1$SYMBOL #get the gene symbol of the row selected
    p<-ggplot(e,aes_string(x="maineffect",y="signal",col=input$color))+plotTheme+guides(color=guide_legend(title=as.character(input$color)))+
      labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+
      stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2)#plot data
    ggplotly(p)
  })

  # output dotplot
  output$dotplot = renderPlotly({
    dotplot_out()
  })
  
  # update tab1 with dotplot
  observe({
    s = input$table_rows_selected
    if(length(s)){
      updateTabsetPanel(session = session, inputId = "tabvalue", selected = 'tab1')
    }
  })
  
###################################################
###################################################
############# CAMERA OUTPUT DISPLAY ##############
###################################################
###################################################

  #populate camera dropdown menu
  output$cameradd = renderUI({
    results=fileload()
    contrast=input$contrast
    cam=paste("results$camera$",contrast,sep="")
    cam=eval(parse(text=cam))
    cameradd=as.list(as.character(unlist(lapply((names(cam)),factor))))
    selectInput("cameradd","Select a Gene Set",cameradd)
  })
  
  #Get camera and roast data from Rdata file
  geneid = reactive({
    validate(
      need(input$camera, "Please Click on the button to view Camera results ")
    )
    results=fileload()
    cameradd=input$cameradd
    contrast=input$contrast #get user input for contrast/comparison
    c=paste('results$camera$',contrast,'$',cameradd,'$camera_result',sep='') #get camera data corresponding to the contrast chosen
    cam=eval(parse(text = c)) #convert string to variable
     cam=data.frame(name=rownames(cam),cam)
    return(cam) # return datatable with camera results
  })

  # print out camera data in tab gsea
  output$tablecam = DT::renderDataTable({
    input$camera
    isolate({
      DT::datatable(geneid(),
                    extensions = c('TableTools','ColVis','Scroller'),
                    options = list(
                      searchHighlight = TRUE,
                      pageLength = 10,
                      lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                      scrollX = TRUE
                    ),rownames= FALSE,selection = list(mode = 'single', selected =1),caption = "Camera Results")
    })
  })

  # display data in tab gsea
        observe({
          if(input$camera>0){
            updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'gsea') 
          }
        })
  #~~~~~~~~~~~~~~~~
        campick2 = reactive({
          results=fileload()
          cameradd=input$cameradd
          contrast=input$contrast #get user input for contrast/comparison
          c=paste('results$camera$',contrast,'$',cameradd,'$indices',sep='') #get camera indices corresponding to the contrast chosen
          cameraind=eval(parse(text = c))
          
          cam=geneid() #get datatable with camera data from reactive
          s=input$tablecam_rows_selected # get  index of selected row from table
          cam=cam[s, ,drop=FALSE] 
          
          #differential expression data (limma data)
          d=datasetInput0.5()
          lfc=as.numeric(input$lfc)
          apval=as.numeric(input$apval)
          res = d[which(abs(d$logFC) > lfc & d$adj.P.Val < apval),] 
          
          #get gene list from indices 
          k=paste('res$ENTREZID[cameraind$',cam$name,']',sep='')
          genes=eval(parse(text = k)) #get entrez id's corresponding to indices
          genesid=res[res$ENTREZID %in% genes,] #get limma data corresponding to entrez id's
          return(genesid)
        })
        
        #print data table with gene list corresponding to each row in camera datatable
        output$campick3 = DT::renderDataTable({
          genesid=campick2()
        },rownames= FALSE,caption="GENE LIST")
        
        
###################################################
###################################################
######### CREATE HEATMAP FROM CAMERA ##############
###################################################
###################################################
        #extract voom expression data of all genes corresponding to selected row in camera datatable
        heatmapcam <- reactive({
          genesid=campick2()  #gene list from camera
          voom=datasetInput3() #voom data
          genes_cam<-voom[rownames(voom) %in% genesid$id,]
          
        })
        
        #create heatmap function
        camheatmap <- function(){
          dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
          hmpcol1=input$hmpcol1 # user input-color palette
          clust=input$clusterby1# user input-colorby
          expr <- heatmapcam() #voom expression data of all genes corresponding to selected row in camera datatable
          pval=campick2() #gene list from camera
          sym=pval$SYMBOL 
          top_expr=data.frame(expr[,-1])
          #top_expr=top_expr[,-((ncol(voom)-5):ncol(voom))]
          if(input$checkbox1==TRUE){ # user input-reverse color option 
            d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby1,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol1)))(30),labRow = sym)}
          else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby1,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol1))(30),labRow = sym)}
        }
        
        # make heatmap for genes 
        output$heatmapcam <- renderD3heatmap({
          input$hmpcol1
          input$clusterby1
          input$checkbox1
          #isolate({
              camheatmap()
           # })
        })
        
        # update gsea tab with the heatmap
          observe({
            s = input$tablecam_rows_selected
            if(length(s)){
              updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'gsea')
            }
          })
          
  #download heatmap
          output$downloadcam <- downloadHandler(
            filename = function(){ 
              paste0('heatmap-cam', '.html', sep='')
            },
            content = function(file){
              saveWidget(camheatmap(),file)
            })
 ##################################################
  ###################################################
  ################ KEGG PATHWAY ANALYSIS#############
  ###################################################
  ###################################################
  
  datasetInput6 = reactive({
    validate(
      need(input$path, "Please choose either upregulated or downregulated")
    )
    final_res=datasetInput()
     #extract logfc values from limma output and the corresponsing entrez names after annotation
     logfc=final_res$logFC
     names(logfc)=final_res$ENTREZID
    
     #gives list of pathways thats upregulated and downregulated
     kegg_results = gage(logfc, gsets=kegg.sets.mm, same.dir=TRUE)
     l=input$num
     #user-input -upregulated/downregulated kegg pathway
    if(input$path=='up')
     {
       #pull the top n upregulated pathways  their id's
       kegg_path=data.frame(id=rownames(kegg_results$greater), kegg_results$greater) 
       kegg_path_top=kegg_path$id[1:l]
       keggresids = substr(kegg_path_top, start=1, stop=8)
     }
     else if(input$path=='down')
     {
       #pull the top n downregulated pathways  their id's
       kegg_path=data.frame(id=rownames(kegg_results$less), kegg_results$less) 
       kegg_path_top=kegg_path$id[1:l]
       keggresids = substr(kegg_path_top, start=1, stop=8)
     }
  })

  #create plots for KEGG pathways
  output$plots = renderUI({
    input$makeplot
    keggresids=datasetInput6() #get KEGG id's of top n number  of pathways
    plot_output_list = lapply(1:length(keggresids), function(i) 
      {
          plotname = paste("plot", i, sep="") #for every KEGG id, create a plot output in pathway plot tab
          plotOutput(plotname, height = 900, width = 800)
    })
    do.call(tagList, plot_output_list)
  })
  
  for (i in 1:15) { #maximum number of plots is 15
    local({
      my_i = i
      plotname = paste("plot", my_i, sep="")

      output[[plotname]] = renderImage({
        withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
            input$makeplot
              keggresids=datasetInput6() #get all KEGG id's
              pId=keggresids[my_i] #get KEGG id one by one in a  loop
              genelist=keggLink("mmu",pId) #for each kegg id, get gene list
              myurl=mark.pathway.by.objects(pId,genelist) #get url of pathway image
              outfile = tempfile(fileext='.png') #create temp file
              png(outfile, width=900, height=800) #get temp file in png format
              download.file(myurl,outfile,mode="wb") #download png into the temp file
              png = readPNG(outfile) # read the PNG from the temp file and display
              dev.off()
              
              list(src = outfile,contentType = 'image/png',width = 900,height = 800,alt = "This is alternate text")
        })
          }, deleteFile = TRUE) #delete temp file
    })
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # display plots in tab 5
  observe({
    if(input$makeplot > 0)
    {
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab5')
    }
  })

  ##################################################
  ##################################################
  ############## GAGE GENE ONTOLOGY ################
  ##################################################
  ##################################################

  datasetInput7 = reactive({ 
    validate(
      need(input$gage, "Please Select Ontology")
    )
    final_res=datasetInput() #get limma data
    logfc=final_res$logFC #get FC values from limma data
    names(logfc)=final_res$ENTREZID # get entrez ids for each row
    data(go.sets.mm) #load GO data from gage
    data(go.subs.mm)
    if(input$gage=='BP') 
    {
      gobpsets = go.sets.mm[go.subs.mm$BP]
      go_res = gage(logfc, gsets=gobpsets)
    }
    else if(input$gage=='cc')
    {
      goccsets = go.sets.mm[go.subs.mm$CC]
      go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
    }
    else if(input$gage=='MF')
    {
      gomfsets = go.sets.mm[go.subs.mm$MF]
      go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
    }
    return(go_res)
  })
  
  #Get all upregulated GO terms
  datasetInput8 = reactive({
    go_res=datasetInput7()
    go_dd=input$go_dd
    if(go_dd=="upreg"){
      res=data.frame(go_res$greater)} #load limma data
    else if(go_dd=="downreg"){
      res=data.frame(go_res$less)
      }
    res = data.frame(GOterm=rownames(res),res)
    
    #Get GO id from GO terms
#     row=data.frame(lapply(res,as.character),stringsAsFactors = FALSE)
#     p=strsplit(row[,1], " ")
#     m=sapply(p,"[",1)
#     go_up=data.frame(GO_id=m,res)
#     return(go_up)
  })

  #Print upregulated GO data in datatable
      output$table4 = DT::renderDataTable({
        input$ga
        input$go_dd
        input$gage
        withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
        isolate({  
        DT::datatable(datasetInput8(),
                      extensions = c('TableTools','ColVis','Scroller'),
                      options = list(
                        searchHighlight = TRUE,
                        pageLength = 10,
                        lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                        scrollX = TRUE
                      ),rownames=FALSE,selection = list(mode = 'single', selected =1))
      })
      })
      })

  #display data in tab 6
  observe({
    if(input$ga>0){
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab6') 
    }
  })

  ###################################################
  ###################################################
  ########## GET GENES FROM  GO ###############
  ###################################################
  ################################################### 
  # get GO associated genes 
  GOHeatup = reactive({
    s = input$table4_rows_selected
    dt = datasetInput8() #load limma data
    dt = as.character(dt[s, , drop=FALSE]) #get limma data corresponding to selected row in table
    goid=dt[1]
    enterezid=paste("go.sets.mm$`",goid,"`",sep="")
    entrezid=eval(parse(text=enterezid))
    limma=datasetInput()
    lim_vals=limma[limma$ENTREZID %in% entrezid,]
  })
  
  #Print datatable with gene list
  output$x4 = DT::renderDataTable({
    input$gage
    input$go_dd
    input$ga
    goheatup=GOHeatup()
  },caption="Gene List",escape=FALSE)
  
  ###################################################
  ###################################################
  ########## MAKE HEATMAP WITH GO ###################
  ###################################################
  ################################################### 
 
  #plot heatmap
  goheatmapup <- function(){
     dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcolup=input$hmpcolup
    
    #pval=GOHeatup()
    pval=GOHeatup()
    sym=pval$SYMBOL
    top_expr=datasetInput3()
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    #top_expr=data.frame(pval[,-c(1,(ncol(pval)-2):ncol(pval))])
    if(input$checkboxup==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterbyup,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcolup)))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterbyup,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcolup))(30),labRow = sym)}
  }
  
  # make heatmap for genes 
  output$hmpgageup <- renderD3heatmap({
    input$gage
    input$go_dd
    input$ga
    input$hmpcolup
    input$clusterbyup
    input$checkboxup
    goheatmapup()
  })
  
  # update with heatmap
  observe({
    s = input$table4_rows_selected
    if(length(s)){
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab6')
    }
  })
  
  output$downloadGOup <- downloadHandler(
    filename = function(){ 
      paste0('heatmapGO', '.html', sep='')
    },
    content = function(file){
      saveWidget(goheatmapup(),file)
    })
  
  ###################################################
  ###################################################
  ####### CREATE HEATMAP FOR LIMMA DATA##############
  ###################################################
  ###################################################

  # Get gene list from user
  datasetInput41 = reactive({
    validate(
      need(input$genelist, "Please input Genelist")
    )
    genes=input$genelist #get complete gene list as string
    p=strsplit(genes,",") #split string by ','
    genelist=sapply(p,"[") #get array of strings with gene id's
    
    #load limma and voom data
    limma=datasetInput()
    voom=datasetInput3()
    #get expression values of the genes in the gene list
        # user-defined identifier for the gene list
    if(input$selectidentifier=='ensembl')
    {
      expr_vals=voom[rownames(voom) %in% genelist,]
      sym=limma$SYMBOL[rownames(limma) %in% genelist]
      expr_vals=data.frame(expr_vals,sym)
    }
    else if(input$selectidentifier=='entrez')
    {
      limma$id<-rownames(limma)
      ensembleid=limma$id[limma$ENTREZID %in% genelist]
      expr_vals=voom[rownames(voom) %in% ensembleid,]
      sym=limma$SYMBOL[limma$ENTREZID %in% genelist]
      expr_vals=data.frame(expr_vals,sym)
    }
    else if(input$selectidentifier=='genesym')
    {
      limma$id<-rownames(limma)
      ensembleid=limma$id[limma$SYMBOL %in% genelist]
      sym=limma$SYMBOL[limma$SYMBOL %in% genelist]
      expr_vals=voom[rownames(voom) %in% ensembleid,]
      expr_vals=data.frame(expr_vals,sym)
    }
    validate(
      need(nrow(expr_vals) > 1, "Please Check Identifier chosen or Select genelist from Raw Expression Data tab")
    )
    validate(
      need(nrow(expr_vals)== nrow(genelist), "One or more of the genes entered does not exist. Check list for typos and verify using Raw Expression Data tab")
    )
    return(expr_vals)
  })
  
  #create heatmap function for gene-list given by user
  heatmap2 = function(){
    dist2 = function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol#user input-color palette
    expr = datasetInput41()
    sym=expr$sym
    expr2=data.frame(expr[,-ncol(expr)])
    #rownames(expr2)=expr[,1]
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
    else{d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
  }
  
  #create heatmap
  datasetInput4 <- reactive({
    validate(
      need(input$gene, "Please Enter number of genes to plot heatmap ")
    )
    #sort by pval
    n<-input$gene #number of genes selected by user (input from slider)
    d<-datasetInput()
    res<-d[order(d$adj.P.Val),]
    reqd_res=res[1:n,] #get top n number of genes
    return(reqd_res)
  })
  
  #create heatmap function for top n genes
  heatmap <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol #user input-color palette
    expr <- datasetInput3()
    pval <- datasetInput4()
    #get expression values of genes with highest pvals
    top_expr=expr[match(rownames(pval),rownames(expr)),]
    sym=pval$SYMBOL
    #top_expr=top_expr[,-1]
    if(input$checkbox==TRUE){
    d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
  }
  
  # make heatmap for genes 
  output$heatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes 
    input$genelist
    input$hmip
    input$makeheat
    #if user selected enter n num of genes, call heatmap() and if user entered genelist, call heatmap2()
    isolate({
      if(input$hmip == 'genenum'){ 
      heatmap()}
      else if(input$hmip == 'geneli'){heatmap2()}
      })
  })
 
  # update tab4 with heatmap
  observe({
    if(input$makeheat > 0)
    {
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab4')
     }
  })
  
  #Download heatmap from top n genes
  output$downloadheatmap <- downloadHandler(
    filename = function(){ 
    paste0('heatmap', '.html', sep='')
      },
    content = function(file){
      saveWidget(heatmap(),file)
      saveWidget(heatmap2(),file)
    })

})