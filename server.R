library(shiny)
#library(pathview)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(gage)
library(gageData)
library(RColorBrewer)
library(NMF)
library(KEGGgraph)
library(Biobase)
library(reshape2)
library(ggplot2)
library(biomaRt)
library(KEGGREST)
library(png)
library(GO.db)
library(d3heatmap)
library(dplyr)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(FactoMineR)
library(factoextra)
library(shinyRGL)
library(rgl)
library(rglwidget)
library(SPIA)
library(ReactomePA)
library(limma)

#load entrez id's for kegg pathway
data(kegg.sets.mm)
#load indexes for signaling and metabolic pathways
data(sigmet.idx.mm)
#get entrez id's for signaling and metabolic pathwaysmera
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]

data(go.sets.mm)
data(go.subs.mm)

#Create a theme for all plots.
plotTheme <-theme_bw() + theme(axis.title.x = element_text(face="bold", size=12),
                               axis.text.x  = element_text(angle=35, vjust=0.5, size=12),
                               axis.title.y = element_text(face="bold", size=12),
                               axis.text.y  = element_text(angle=0, vjust=0.5, size=12))




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
    limmadata=subset(limmadata, select=-c(logFC))
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
      d = d[which(d$fc < (-1*(lfc)) & d$adj.P.Val < apval),]
    }
    else if(input$radio=='up')
    {
      d=limmadata
      d = d[which(d$fc>lfc & d$adj.P.Val < apval),]
    }
    else if(input$radio=='both')
    {
      d=limmadata
      d = d[which(abs(d$fc) > lfc & d$adj.P.Val < apval),]
    }
    #limma = as.data.frame(d) #load limma data
    geneid=d$SYMBOL
    url= paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",geneid,sep = "")
    if(url=="http://www.genecards.org/cgi-bin/carddisp.pl?gene="){
      d$link<-NULL
    }else{
    d$link=paste0("<a href='",url,"'target='_blank'>","Link to GeneCard","</a>")}
    d=as.data.frame(d) 
    #d=as.data.frame(d) %>% select(-t)
    return(d)
  })

  #print input file in tab1
  output$table = DT::renderDataTable({
    input$lfc
    input$apval
    input$project
    input$contrast
      DT::datatable(datasetInput(),
                    extensions = 'Buttons', options = list(
                      dom = 'Bfrtip',
                      buttons = 
                        list('copy', list(
                          extend = 'collection',
                          buttons = c('csv', 'excel', 'pdf','print'),
                          text = 'Download'
                        ))
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
  #######CONDITIONAL PANEL FOR Limma ##################
  ###################################################
  ###################################################
  
  #Create checkboxes with contrasts corresponding to the project (displayed only when multiple contrast checkbox is selected)
  output$contrastslimma <- renderUI({
    results=fileload()
    lim=results$limma
    contrasts=as.list(as.character(unlist(lapply((names(lim)),factor))))
    checkboxGroupInput("multicontrast",label="Pick Contrasts",choices=contrasts)
  })
  
  #create table with p.value and FC value for the contrasts selected
  multilimma = reactive({
    validate(
      need(input$multicontrast, "Please Select at least one comparison ")
    )
    contr=input$multicontrast
    results=fileload()
    full_limma = data.frame(id=as.character())
    for(i in 1:length(contr)){
      k=paste('results$limma$',contr[i],sep='')
      limmadata=eval(parse(text = k))
      limmadata2=data.frame(id=rownames(limmadata),logFC=limmadata$logFC,adj.P.Val=limmadata$adj.P.Val)
      colnames(limmadata2)[-1]=paste(colnames(limmadata2[,c(-1)]),contr[i], sep = "_")
      full_limma=full_join(full_limma,limmadata2,by='id')
    }
    k=data.frame(id=rownames(limmadata),SYMBOL=limmadata$SYMBOL)
    m=full_join(k,full_limma,by='id')
    return(m)
  })
  
  #Display text (contrast name) above limma table
  output$contrdesc <- renderText({
    contrastname=input$contrast
    text=paste('CONTRAST:  ',contrastname,sep="   ")
    
    return(text)
  })

  

  #update table with the dataframe
  output$table_TRUE = DT::renderDataTable({
    input$project
    input$contrast
    DT::datatable(multilimma(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                    searchHighlight = TRUE,
                    pageLength = 10,
                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                    scrollX = TRUE,
                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                  ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE)
  })

  #display tab with the table only when the multiple contrast checkbox and the contrasts are selected
  observe({
    if(input$check>0){
    updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab11')}
    toggle(condition =input$check,selector = "#tabvalue li a[data-value=tab11]")
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
  })
  
  #annotate voom data using featuresdata 
  datasetInput33 = reactive({
    results=fileload()
    exprsdata=as.data.frame(exprs(results$eset))
    features=as.data.frame(pData(featureData(results$eset)))
    features$id=rownames(features)
    exprsdata$id=rownames(exprsdata)
    genes <- inner_join(features,exprsdata,by=c('id'='id'))
    return(genes)
  })

  #print voom or expression data file in tab2
 output$table3 = DT::renderDataTable({
    DT::datatable(datasetInput33(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                    searchHighlight = TRUE,
                    pageLength = 10,
                    lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                    scrollX = TRUE,
                    buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Voom data")
  })
  # display data in tab2
  observe({
    if(input$rawdata>0) {
    updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab2') }
  })


  output$rawdwld <- downloadHandler(
    filename = function() { paste(input$projects, '_rawdata.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput33(), file)
    })
  
  ###################################################
  ###################################################
  ##################### DOT PLOT ####################
  ###################################################
  ###################################################
  #create user interface for dot-plot drop down menu
  output$boxplotcol = renderUI({ 
    results=fileload()
    eset=results$eset
    pData=pData(eset) #get pheno-data
    kc=pData[ , grepl( "var_" , colnames(pData) ) ] #get columns from phenodata that start with "var"
    kt=as.data.frame(t(na.omit(t(kc)))) #omit columns that have only NA's
    kc=data.frame(maineffect=pData$maineffect,sample_name=pData$sample_name,kt) #create a dataframe with maineffect and sample_name and non-null columns starting with var
    bpcols=as.list(as.character(unlist(lapply((colnames(kc)),factor))))
    selectInput("color","Select an Attribute",bpcols) #populate drop down menu with the phenodata columns
  })

  #textbox to enter x-axis variables
  output$boxreorder <- renderUI({
    textInput("boxlist", label = h5("Enter x-axis variables in the preferred order"), value = "")
  })
  
  dotplot_out = reactive({
    s = input$table_rows_selected #select rows from table
    dt = datasetInput() #load limma data
    dt$id=rownames(dt)
    dt=data.frame(dt$id,dt[,-ncol(dt)])
    validate(
      need((is.data.frame(dt) && nrow(dt))!=0, "No data in table")
    )
    dt1 = dt[s, , drop=FALSE]#get limma data corresponding to selected row in table
    id = as.character(dt[s,1]) 
    results=fileload()
    eset <- results$eset
    e <-data.frame(eset@phenoData@data,signal=exprs(eset)[id,])
    if(is.na(dt1$SYMBOL)) #if gene symbol does not exist,use ENSEMBL id
    {genesymbol=dt1$ENSEMBL}
    else{
      genesymbol=dt1$SYMBOL} #get the gene symbol of the row selected
    
    #To reorder x-axis, take in input list, split by comma and pass that as a variable to ggplot
    a=input$boxlist
    aa=strsplit(a,",")
    aaa=unlist(aa)
    me=factor(e$maineffect,levels=c(aaa))
    
    if(input$boxreorder>0){
    gg=ggplot(e,aes_string(x=me,y="signal",col=input$color))+plotTheme+guides(color=guide_legend(title=as.character(input$color)))+
      labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+
      stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2)
    }
    else{
      gg=ggplot(e,aes_string(x="maineffect",y="signal",col=input$color))+plotTheme+guides(color=guide_legend(title=as.character(input$color)))+
        labs(title=genesymbol, x="Condition", y="Expression Value") + geom_point(size=5,position=position_jitter(w = 0.1))+
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size= 0.3, geom = "crossbar",width=.2)
    }#plot data
    #gg=ggplotly(gg)
    gg
  })

  # output dotplot
  output$dotplot = renderPlot({
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
  ########### VARIANCES OF PCA PLOT #################
  ###################################################
  ###################################################
  
  output$pcatitle <- renderText({
    text="The proportion of variances retained by the principal components can be viewed in the scree plot. The scree plot is a graph of the eigenvalues/variances associated with components"
    return(text)
  })
  
  #get expression data and perform PCA
  res_pca = reactive({
    n=as.numeric(input$pcipslide)
    validate(
      need(input$pcipslide > 199, "Minimum value of input genes that show maximum variance should at least be 200")
    )
    results=fileload()
    v = results$eset
    keepGenes <- v@featureData@data
    #keepGenes <- v@featureData@data %>% filter(!(seq_name %in% c('X','Y')) & !(is.na(SYMBOL)))
    pData<-phenoData(v)
    v.filter = v[rownames(v@assayData$exprs) %in% rownames(keepGenes),]
    Pvars <- apply(exprs(v.filter),1,var)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(n,length(Pvars)))]
    v.var <-v.filter[select,]
    m<-exprs(v.var)
    rownames(m) <- v.var@featureData@data$SYMBOL
    res.pca = PCA(t(m), graph = FALSE)
  })
  
  #get expression data and perform PCA
  pcaplo_tab = reactive({
    res.pca =res_pca()
    eigenvalues = res.pca$eig
    return(eigenvalues)
  })
  
  output$pcaplot_tab = DT::renderDataTable({
    DT::datatable(pcaplo_tab(),
                  extensions = c('Scroller'),
                  options = list(
                                 searchHighlight = TRUE,
                                 scrollX = TRUE
                  ),caption= "Amount of Variation explained by each Principle Component")
  })
  
   output$pcaplot_ip = renderPlot({
     res.pca = res_pca()
     fviz_screeplot(res.pca, ncp=10)
   })
   
   observe({
     if(input$varpc>0){
       updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tabvar')}
     toggle(condition =input$varpc,selector = "#tabvalue li a[data-value=tabvar]")
   })
   
   ###################################################
   ###################################################
   ##################### PCA PLOT ####################
   ###################################################
   ###################################################
  
   output$pcaxoptions <- renderUI({
     selectInput("pcaxaxes","Select Principle Component to plot on the X-axis ",c(1:10))
   })
   
   output$pcipslide <- renderUI({
     textInput(inputId = 'pcipslide', label = "Enter top number of input genes that show maximum variance", value = '500')
     #sliderInput("pcipslide", label = h5("Slider: Number of input genes that show maximum variance"), min = 100,max = 1000, value = 200)
   })
   
   output$pcslide <- renderUI({
     textInput(inputId = 'pcslide', label = "Enter number of genes to view in the biplot", value = '0')
     #sliderInput("pcslide", label = h5("Slider: Number of genes to view in the biplot"), min = 2,max = 50, value = 30)
   })
   
   output$pcayoptions <- renderUI({
     selectInput("pcayaxes","Select Principle Component to plot on the Y-axis",c(1:10),selected=2)
   })
   
   output$biplottitle <- renderText({
     pcaxaxes=input$pcaxaxes
     pcayaxes=input$pcayaxes
     text=as.character(paste("Dim",pcaxaxes," vs Dim",pcayaxes,sep=""))
       return(text)
   })
   
   plotbiplot = reactive({
     res.pca = res_pca()
     x=as.numeric(input$pcaxaxes)
     y=as.numeric(input$pcayaxes)
     results=fileload()
     v = results$eset
     pData<-phenoData(v)
#      validate(
#        need(input$pcslide > 1, "Minimum value of number of genes to display in the biplot should be 2")
#      )
     if(input$pcslide==0){
       fviz_pca_ind(res.pca, repel=T,geom='point',label='var',addEllipses=FALSE, habillage = pData$maineffect)+scale_shape_manual(values = c(rep(19,length(unique(pData$maineffect)))))
     }
       #fviz_pca_ind(res.pca, geom = c("point", "text"))}
     else{fviz_pca_biplot(res.pca,repel=T, label=c("var","ind"),habillage = as.factor(pData$maineffect),axes=c(x,y),select.var = list(contrib = as.numeric(input$pcslide)))}
   })
   
   output$biplot = renderPlot({
     plotbiplot()
   })
   
   
#    output$downloadbiplot <- downloadHandler(
#      filename = function(){
#        paste0('biplot', '.html', sep='')
#      },
#      content = function(file){
#        saveWidget(plotbiplot(),file)
#      })
   
   output$downloadbiplot <- downloadHandler(
     filename = function() {
       paste0("biplot.jpg")
     },
     content = function(file){
       jpeg(file, quality = 100, width = 800, height = 800)
       plot(plotbiplot())
       dev.off()
     })
   
   observe({
     if(input$makepcaplot>0){
       updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tabpca')}
   })
   
   ###################################################
   ###################################################
   ##################3D PCA PLOT #####################
   ###################################################
   ###################################################
   output$pcaplot3d = renderRglwidget({
     v=datasetInput3()
     results=fileload()
     pData=pData(results$eset)
     v=t(v)
     v= v[,apply(v, 2, var, na.rm=TRUE) != 0]
     pca <- prcomp( v, scale= TRUE )
     vars <- apply(pca$x, 2, var)
     props <- round((vars / sum(vars))*100,1)
     groups=factor(gsub('-','_',pData$maineffect))
#      
#      pca <- res_pca()
#      vars <- apply(pca$var$coord, 2, var)
#      props <- round((vars / sum(vars))*100,1)
#      groups=factor(gsub('-','_',pData$maineffect))
#      
     
      #######
#      n=as.numeric(input$pcipslide)
#      validate(
#        need(input$pcipslide > 199, "Minimum value of input genes that show maximum variance should at least be 200")
#      )
#      results=fileload()
#      v = results$eset
#      keepGenes <- v@featureData@data
#      #keepGenes <- v@featureData@data %>% filter(!(seq_name %in% c('X','Y')) & !(is.na(SYMBOL)))
#      pData<-phenoData(v)
#      v.filter = v[rownames(v@assayData$exprs) %in% rownames(keepGenes),]
#      Pvars <- apply(exprs(v.filter),1,var)
#      select <- order(Pvars, decreasing = TRUE)[seq_len(min(n,length(Pvars)))]
#      v.var <-v.filter[select,]
#      m<-exprs(v.var)
#      rownames(m) <- v.var@featureData@data$ENSEMBL
#      res.pca = PCA(t(m), graph = FALSE)
     ########
     try(rgl.close())
     open3d()
     # resize window
     par3d(windowRect = c(100, 100, 612, 612))
     palette(c('blue','red','green','orange','cyan','black','brown','pink'))
     plot3d(pca$x[,1:3], col =as.numeric(groups), type='s',alpha=.75,axes=F,
            xlab=paste('PC1 (',props[1],'%)',sep=''),
            ylab=paste('PC2 (',props[2],'%)',sep=''),
            zlab=paste('PC3 (',props[3],'%)',sep='')
     )
     axes3d(edges=c("x--", "y--", "z"), lwd=2, expand=10, labels=FALSE,box=T)
     grid3d("x")
     grid3d("y")
     grid3d("z")
     l=length(levels(groups))
     ll=1:l
     y=100+(ll*90)
     text3d(x=300, y=y, z=1.1,levels(groups) ,col="black")
     points3d(x=400, y=y, z=1.1, col=as.numeric(as.factor(levels(groups))), size=10)
     legend3d("topright", legend = levels(groups), pch = 16, col=palette(),cex=1, inset=c(0.02))
     
#      rgl.snapshot('PCA_3d_test.png', fmt = "png", top = TRUE )
#      rgl.postscript('PCA_3D_test.pdf',fmt='pdf')
     
     #movie3d(spin3d(), duration = 5,movie='PCA_movie',dir='./' )
     rglwidget()
   })
   
   observe({
     if(input$pca3d>0){
       updateTabsetPanel(session = session, inputId = 'tabvalue', selected = '3dpca')}
     toggle(condition =input$pca3d,selector = "#tabvalue li a[data-value=3dpca]")
   })
   
###################################################
###################################################
############# CAMERA OUTPUT DISPLAY ##############
###################################################
###################################################

  #populate camera dropdown menu with the genesets based on the project RData
  output$cameradd = renderUI({
    results=fileload()
    contrast=input$contrast
    cam=paste("results$camera$",contrast,sep="")
    cam=eval(parse(text=cam))
    cameradd=as.list(as.character(unlist(lapply((names(cam)),factor))))
    selectInput("cameradd","Select a Gene Set",cameradd)
  })

  #Get camera data from Rdata file
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
     name=cam$name
     if (cameradd == "GO")
     {
       url= paste("http://amigo.geneontology.org/amigo/term/",name,sep = "") #create link to Gene Ontology Consortium
       cam$link=paste0("<a href='",url,"'target='_blank'>","Link to Gene Ontology Consortium","</a>")
       cam=as.data.frame(cam)
     }else{
     url= paste("http://software.broadinstitute.org/gsea/msigdb/cards/",name,".html",sep = "")
     cam$link=paste0("<a href='",url,"'target='_blank'>","Link to Molecular Dignature Database","</a>")
     cam=as.data.frame(cam)}
    return(cam) # return datatable with camera results
  })

  # print out camera data in tab gsea
  output$tablecam = DT::renderDataTable({
    input$camera
    input$cameradd
    input$contrast
    isolate({
      DT::datatable(geneid(),
                    extensions = c('Buttons','Scroller'),
                    options = list(dom = 'Bfrtip',
                      searchHighlight = TRUE,
                      pageLength = 10,
                      lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                      scrollX = TRUE,
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                    ),rownames= FALSE,selection = list(mode = 'single', selected =1),escape=FALSE,caption = "Camera Results")
    })
  })

  # display data in tab gsea
        observe({
          if(input$camera>0){
            updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'gsea')
          }
        })
        
  #get the gene-list for every row in camera results table
        campick2 = reactive({
          results=fileload()
          cameradd=input$cameradd
          contrast=input$contrast #get user input for contrast/comparison
          c=paste('results$camera$',contrast,'$',cameradd,'$indices',sep='') #get camera indices corresponding to the contrast chosen
          cameraind=eval(parse(text = c))
          cam=geneid() #get datatable with camera data from reactive
          s=input$tablecam_rows_selected # get  index of selected row from table
          cam=cam[s, ,drop=FALSE]
            res=datasetInput0.5()
            res2=datasetInput33()
            if("ENTREZID" %in% colnames(res2)){
              res2=res2
            }
            else{res2=res}
          
            #get gene list from indices
            if (cameradd == "GO")
            {
              #res2=datasetInput0.5()
          k=paste('res2$ENTREZID[cameraind$`',cam$name,'`]',sep='')}
            else{
              k=paste('res2$ENTREZID[cameraind$',cam$name,']',sep='')
            }
          genes=eval(parse(text = k)) #get entrez id's corresponding to indices
          #genesid=res[match(res$ENTREZID,genes),]
          genesid=res[res$ENTREZID %in% genes,] #get limma data corresponding to entrez id's
          #rownames(genesid)=genesid$id
          return(data.frame(genesid)) #return the genelist
        })

        #print data table with gene list corresponding to each row in camera datatable
        output$campick3 = DT::renderDataTable({
          input$cameradd
          input$contrast
          input$projects
          DT::datatable(campick2(),
                        extensions = c('Buttons','Scroller'),
                        options = list(dom = 'Bfrtip',
                                       searchHighlight = TRUE,
                                       pageLength = 10,
                                       lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                       scrollX = TRUE,
                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                        ),rownames=TRUE,escape=FALSE,caption="GENE LIST")
        })
        
#         output$campick3 = DT::renderDataTable({
#           genesid=campick2()
#         },rownames= FALSE,caption="GENE LIST")
        
        #download camera datatable
        output$downloadcam <- downloadHandler(
          filename = function() { paste('Camera_',input$projects,'_',input$contrast,'.csv', sep='') },
          content = function(file) {
            write.csv(geneid(), file)
          })
        
        #Generate text title for the gene list table
        output$camdesc <- renderText({
          s = input$tablecam_rows_selected
          dt = geneid() 
          dt = as.character(dt[s, , drop=FALSE]) 
          camname=dt[1]
          text=paste('Gene list for Camera term :',camname,sep="")
          
          return(text)
        })
        

###################################################
###################################################
############### ENRICHMENT PLOT ###################
###################################################
###################################################

#         ep_out <- reactive({
#           results=fileload()
#           contrasts=input$contrast
#           res=paste("results$eplot$",contrasts,sep="")
#           tab=eval(parse(text =res))
#           validate(
#             need(nrow(tab) > 1, "No Enriched Pathways found")
#           )
#           return(tab)
#         })
        output$eplottab = DT::renderDataTable({
          input$eplot
          input$camera
          input$cameradd
          input$contrast
          DT::datatable(geneid(),
                        extensions = c('Buttons','Scroller'),
                        options = list(dom = 'Bfrtip',
                                       searchHighlight = TRUE,
                                       pageLength = 10,
                                       lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                       scrollX = TRUE,
                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                        ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1))
          })
       
        output$eplotdesc <- renderText({
          res=geneid()
          s=input$eplottab_rows_selected 
          row=res[s, ,drop=FALSE]
          id=row$name
          text=paste('Enrichment plot for :',id,sep="")
          
          return(text)
        })
        
        
        output$en_plot <- renderPlot({
          res=geneid()
          cameradd=input$cameradd
          contrast=input$contrast
          lim=datasetInput0.5()
          stat=data.frame(t=lim$t)
          stat$id=rownames(lim)
          stat=stat[order(stat$id),]
          stat2=stat$t
          names(stat2)=stat$id
          s=input$eplottab_rows_selected 
          row=res[s, ,drop=FALSE]
          id=row$name
          results=fileload()
          contrasts=input$contrast
          y=paste('results$camera$',contrast,'$',cameradd,'$indices$',id,sep='')
          ind=eval(parse(text =y))
          barcodeplot(stat2,index=ind)
        })
        
         observe({
          if(input$eplot>0){
            updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'eplot')}
          toggle(condition =input$eplot,selector = "#tabvalue li a[data-value=eplot]")
        })
###################################################
###################################################
######### CREATE HEATMAP FROM CAMERA ##############
###################################################
###################################################
        #extract voom expression data of all genes corresponding to selected row in camera datatable
        heatmapcam <- reactive({
          genesid=campick2()  #gene list from camera
          voom=datasetInput3() #voom data
          genes_cam<-voom[rownames(voom) %in% rownames(genesid),]

        })

        #create heatmap function
        camheatmap <- function(){
          dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
          expr <- heatmapcam() #voom expression data of all genes corresponding to selected row in camera datatable
          pval=campick2() #gene list from camera
          rownames(expr)=pval$SYMBOL
          #sym=pval$SYMBOL
          hmplim=input$hmplim
          #top_expr=data.frame(expr[,-1])
#           if(hmplim==0)
#           {
#             top_expr=data.frame(expr)}
#           else{
            top_expr=data.frame(expr)
            top_expr=top_expr[1:hmplim,]
            # }
          
          sym=rownames(top_expr)
          if(input$checkbox==TRUE){
            d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, input$hmpcol))(30),labRow = sym)}
          else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, input$hmpcol)))(30),labRow = sym)}
        }

        # update heatmap tab with the heatmap
          observe({
            s = input$tablecam_rows_selected
            if(length(s)){
              updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'gsea')
            }
          })
          
##################################################
###################################################
################ SPIA PATHWAY ANALYSIS#############
###################################################
###################################################

          spia_op <- reactive({
            validate(
              need(input$runspia, "Please Click on the button to run SPIA ")
            )
            results=fileload()
            contrast=input$contrast #get user input for contrast/comparison
            c=paste('results$spia$',contrast,sep='') #get SPIA data corresponding to the contrast chosen
            sp=eval(parse(text = c)) #convert string to variable
            #spia_result=data.frame(name=rownames(sp),sp)
            spia_result=data.frame(sp)
            validate(
              need(nrow(spia_result) > 1, "No Results")
            )
            spia_result$KEGGLINK <- paste0("<a href='",spia_result$KEGGLINK,"' target='_blank'>","Link to KEGG","</a>")
            return(spia_result) 
          })
          
          
          output$spiaop <- DT::renderDataTable({
            input$runspia
            input$contrast
            input$projects
            withProgress(session = session, message = 'Calculating...',detail = 'This may take a while...',{
              isolate({
                DT::datatable(spia_op(),escape = FALSE,selection = list(mode = 'single', selected =1),
                              extensions = c('Buttons','Scroller'),
                              options = list(
                                dom = 'RMDCT<"clear">lfrtip',
                                searchHighlight = TRUE,
                                pageLength = 10,
                                lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                                scrollX = TRUE,
                                buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                              ),rownames=FALSE)
              })
            })
          })
         
     spiagenes = reactive({
          spiaid=spia_op() 
          final_res=datasetInput()
           s=input$spiaop_rows_selected 
           row=spiaid[s, ,drop=FALSE]
           id=paste("mmu",row$ID,sep="")
           #keggid = substr(keggid, start=1, stop=8)
          allgenelist=keggLink("mmu",id) #for each kegg id, get gene list
          p=strsplit(allgenelist,":")
          genes_entrez=sapply(p,"[",2)
          genelist=final_res[final_res$ENTREZID %in% genes_entrez,]
            return(genelist) #return the genelist
          })
          
          output$spiadesc <- renderText({
            s = input$spiaop_rows_selected
            dt = spia_op() 
            dt = dt[s, , drop=FALSE]
            camname=dt$Name
            text=paste('Gene list for SPIA term :',camname,'-',dt[2],sep="")
            
            return(text)
          })
          
         output$spiagenes = DT::renderDataTable({
           DT::datatable(spiagenes(),
                         extensions = c('Buttons','Scroller'),
                         options = list(dom = 'Bfrtip',
                                        searchHighlight = TRUE,
                                        pageLength = 10,
                                        lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                        scrollX = TRUE,
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                         ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1,caption="Genelist"))
         })
          
         
         # update tab with spia results
         observe({
           if(input$runspia > 0)
           {
             updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'spia')
           }
         })
         
         output$dwldspia <- downloadHandler(
           filename = function() { paste(input$projects,'_',input$contrast, '_spia.csv', sep='') },
           content = function(file) {
             write.csv(spiaop(), file)
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
    logfc=final_res$fc #get FC values from limma data
    names(logfc)=final_res$ENTREZID # get entrez ids for each row
    results=fileload()
    pd=pData(results$eset)
    organism=pd$organism
    validate(
      need(length(unique(organism))==1,"Please check pData file for errors in organism column. Does it have more than one organism ?")
    )
    organism=pd$organism[1]
    if(organism=="human")
     {
       data(go.sets.hs) #load GO data from gage
       data(go.subs.hs)
       
       if(input$gage=='BP')
       {
         gobpsets = go.sets.hs[go.subs.hs$BP]
         go_res = gage(logfc, gsets=gobpsets)
       }
       else if(input$gage=='cc')
       {
         goccsets = go.sets.hs[go.subs.hs$CC]
         go_res = gage(logfc, gsets=goccsets, same.dir=TRUE)
       }
       else if(input$gage=='MF')
       {
         gomfsets = go.sets.hs[go.subs.hs$MF]
         go_res = gage(logfc, gsets=gomfsets, same.dir=TRUE)
       }
     }
       else
       {
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
       }
    return(go_res)
  })

  #Get all GO terms based on user-selection (upregulated/downregulated)
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
    row=data.frame(lapply(res,as.character),stringsAsFactors = FALSE)
    p=strsplit(row[,1], " ")
    m=sapply(p,"[",1)
    go_up=data.frame(GO_id=m,res)
    go_term=go_up$GO_id
    url= paste("http://amigo.geneontology.org/amigo/term/",go_term,sep = "") #create link to Gene Ontology Consortium
    go_up$link=paste0("<a href='",url,"'target='_blank'>","Link to Gene Ontology Consortium","</a>")
    go_up=as.data.frame(go_up)
    return(go_up)
  })

  #Print GO data in datatable
      output$table4 = DT::renderDataTable({
        input$ga
        input$go_dd
        input$gage
        input$radio
        input$project
        input$contrast
        withProgress(session = session, message = 'Generating...',detail = 'Please Wait...',{
        isolate({
        DT::datatable(datasetInput8(),
                      extensions = c('Buttons','Scroller'),
                      options = list(dom = 'Bfrtip',
                        searchHighlight = TRUE,
                        pageLength = 10,
                        lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                        scrollX = TRUE,
                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                      ),rownames=FALSE,escape=FALSE,selection = list(mode = 'single', selected =1))
      })
      })
      })

  #display data in Gene ontology tab
  observe({
    if(input$ga>0){
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab6')
    }
  })
 
  output$downloadgo <- downloadHandler(
    filename = function() { paste('GO_',input$projects,'_',input$contrast,'_',input$gage,'_',input$go_dd,'.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput8(), file)
    })
  
  ###################################################
  ###################################################
  ########## GET GENES FROM  GO ###############
  ###################################################
  ###################################################
  # get GO associated genes
  GOHeatup = reactive({
    s = input$table4_rows_selected
    dt = datasetInput8() #load GO data
    dt = dt[s, , drop=FALSE] #get GO data corresponding to selected row in table
    results=fileload()
    pd=pData(results$eset)
    organism=pd$organism[1]
    goid=dt$GO_id
    if(organism=="human"){
    enterezid=paste("go.sets.hs$`",goid,"`",sep="")
    }
    else{
      enterezid=paste("go.sets.mm$`",goid,"`",sep="")
    }
    entrezid=eval(parse(text=enterezid))
    limma=datasetInput()
    lim_vals=limma[limma$ENTREZID %in% entrezid,]
  })

  #Print datatable with gene list
  output$x4 = DT::renderDataTable({
    input$gage
    input$go_dd
    input$ga
    input$radio
    input$project
    input$contrast
    goheatup=GOHeatup()
  },caption="Gene List",escape=FALSE)
  
  #Text title for gene list table
  output$godesc <- renderText({
    s = input$table4_rows_selected
    dt = datasetInput8() #load GO data
    dt = dt[s, , drop=FALSE] #get GO data corresponding to selected row in table
    goid=dt$GO_id
    text=paste('Gene list for GO term :',goid,sep="")
    
    return(text)
  })
  
  ###################################################
  ###################################################
  ########## MAKE HEATMAP WITH GO ###################
  ###################################################
  ###################################################
  #plot heatmap
  goheatmapup <- function(){
     dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol
    pval=GOHeatup() #genelist from GO
    #sym=pval$SYMBOL
    hmplim=input$hmplim
    top_expr=datasetInput3() #voom expression data of all genes corresponding to selected row in GO datatable
    top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
    rownames(top_expr)=pval$SYMBOL
    top_expr=top_expr[1:hmplim,]
    sym=rownames(top_expr)
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }


  # update heatmap tab with heatmap
  observe({
    s = input$table4_rows_selected
    if(length(s)){
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab6')
    }
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
      #expr_vals=voom[rownames(voom) %in% genelist,]
      sym=limma[limma$ENSEMBL %in% genelist,] %>% select(ENSEMBL,SYMBOL)
      #sym=sym[,c("ENSEMBL","SYMBOL")]
      expr_vals=merge(voom,sym,by="row.names")
      rownames(expr_vals)=expr_vals$Row.names
      expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
    }
    else if(input$selectidentifier=='entrez')
    {
      limma$id<-rownames(limma)
      ensembleid=limma$id[limma$ENTREZID %in% genelist]
      sym=limma[limma$ENTREZID %in% genelist,]
      sym=sym[,c("ENSEMBL","SYMBOL")]
      expr_vals=voom[rownames(voom) %in% ensembleid,]
      expr_vals=merge(expr_vals,sym,by="row.names")
      rownames(expr_vals)=expr_vals$Row.names
      expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
    }
    else if(input$selectidentifier=='genesym')
    {
      limma$id<-rownames(limma)
      ensembleid=limma$id[limma$SYMBOL %in% genelist]
      sym=limma[limma$SYMBOL %in% genelist,]
      sym=sym[,c("ENSEMBL","SYMBOL")]
      expr_vals=voom[rownames(voom) %in% ensembleid,]
      expr_vals=merge(expr_vals,sym,by="row.names")
      rownames(expr_vals)=expr_vals$Row.names
      expr_vals=data.frame(expr_vals[,-c(1,(ncol(expr_vals)-1))])
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
    sym=expr$SYMBOL
    expr2=data.frame(expr[,-ncol(expr)])
    #rownames(expr2)=expr[,1]
    if(input$checkbox==TRUE){
      d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{d3heatmap(as.matrix(expr2),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }

  #create heatmap function for top number of genes as chosen from the slider
  datasetInput4 <- reactive({
    validate(
      need(input$gene, "Please Enter number of genes to plot heatmap ")
    )
    #sort by pval
    n<-input$gene #number of genes selected by user (input from slider)
    d<-datasetInput()
    res<-d[order(d$adj.P.Val),]
    if(n>nrow(d)){
    reqd_res=res[1:nrow(d),]} #get top n number of genes
    else{
      reqd_res=res[1:n,]
    }
    return(reqd_res)
  })

  #create heatmap function for top n genes
  heatmap <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    hmpcol=input$hmpcol #user input-color palette
    #hmplim=input$hmplim
    expr <- datasetInput3()
    pval <- datasetInput4()
    #get expression values of genes with highest pvals
    top_expr=expr[match(rownames(pval),rownames(expr)),]
    #top_expr=top_expr[1:hmplim,]
    sym=pval$SYMBOL
    if(input$checkbox==TRUE){
    d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(brewer.pal(n = 9, hmpcol))(30),labRow = sym)}
    else{d3heatmap(as.matrix(top_expr),distfun=dist2,scale="row",dendrogram=input$clusterby,xaxis_font_size = 10,colors = colorRampPalette(rev(brewer.pal(n = 9, hmpcol)))(30),labRow = sym)}
  }
  
  #manually create scale (colorkey) for heatmap
  
  hmpscale <- reactive({
    hmpcol=input$hmpcol #user input-color palette
    if(input$checkbox==FALSE){
    val=sort(c(-2,-1,0,1,2),decreasing=TRUE)}
    else{
      val=sort(c(-2,-1,0,1,2),decreasing=FALSE)
    }
    df <- data.frame(x = rep(1, 5),y = val,z = factor(1:5))
    ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
    })
  
  output$hmpscale_out = renderPlot({
    hmpscale()
  })
  
  #Set limit for number of genes that can be viewed in the heatmap
  output$hmplim <- renderUI({
    #textInput(inputId = 'hmplim', label = "Enter number of genes to view in the heatmap", value = '0')
#     if(input$hmip == 'genenum'){heatmap()}
#     else if(input$hmip == 'geneli'){
#       expr = datasetInput41()
#       expr2=data.frame(expr[,-ncol(expr)])
#       mx=nrow(expr2)
#       }
    if(input$hmip == 'hmpcam' ){
      expr <- heatmapcam() 
      pval=campick2() 
      rownames(expr)=pval$SYMBOL
      top_expr=data.frame(expr)
      mx=nrow(top_expr)
      sliderInput("hmplim", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
      }
    else if(input$hmip == 'hmpgo'){
      pval=GOHeatup()
      top_expr=datasetInput3()
      top_expr=top_expr[rownames(top_expr) %in% rownames(pval),]
      mx=nrow(top_expr)
      sliderInput("hmplim", label = h5("Select number of genes to view in the heatmap"), min = 2,max =mx, value = mx)
      }
  })
  
  #Text title for type of heatmap being displayed in the heatmap tab
  output$htitle <- renderText({
    hmip=input$hmip
    
    if(input$hmip=="genenum"){text="Heatmap of Top Genes "}
    else if(input$hmip=="geneli"){text="Heatmap of Genelist "}
    if(input$hmip=="hmpgo"){
      s2 = input$table4_rows_selected
      dt2 = datasetInput8() #load GO data
      dt2 = dt2[s2, , drop=FALSE] #get GO data corresponding to selected row in table
      goid=dt2$GOterm
      text=paste("Heatmap of GO term:",goid,sep="")}
    if(input$hmip=="hmpcam"){
      s3 = input$tablecam_rows_selected
      dt3 = geneid() 
      dt3 = dt3[s3, , drop=FALSE] 
      camname=dt3$name
      text=paste("Heatmap of camera term :",camname,sep="")}
    return(text)
  })

  # make heatmap for genes
  output$heatmap <- renderD3heatmap({
    input$hmpcol #user input-color palette
    input$clusterby #user input-cluster by
    input$checkbox #user input-reverse colors
    input$gene #user input-slider input for number of genes
    input$genelist
    input$hmip
    input$makeheat
    input$gage
    input$go_dd
    input$ga
    input$table4_rows_selected
    input$tablecam_rows_selected
    input$radio
    input$projects
    input$contrast
    input$cameradd
    input$hmplim
    #if user selected enter n num of genes, call heatmap() and if user entered genelist, call heatmap2()
    isolate({
      if(input$hmip == 'genenum'){heatmap()}
      else if(input$hmip == 'geneli'){heatmap2()}
      else if(input$hmip == 'hmpcam' ){camheatmap()}
      else if(input$hmip == 'hmpgo'){goheatmapup()}
      })
  })

  # update tab4 with heatmap
  observe({
    if(input$makeheat > 0)
    {
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab4')
     }
  })

  #Download heatmap 
  output$downloadheatmap <- downloadHandler(
    filename = function(){
    paste0('heatmap', '.html', sep='')
      },
    content = function(file){
      if(input$hmip == 'genenum'){saveWidget(heatmap(),file)}
      else if(input$hmip == 'geneli'){saveWidget(heatmap2(),file)}
      else if(input$hmip == 'hmpcam' ){saveWidget(camheatmap(),file)}
      else if(input$hmip == 'hmpgo'){saveWidget(goheatmapup(),file)}
    })

})
