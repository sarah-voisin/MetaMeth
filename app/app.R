#EWAS shiny app
#Load packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(readxl)
library(rintrojs)
library(pals)
library(plotly)

#Load all results
meta_res_robust <- read_tsv("./input_data/MetaAnalysis_new.txt")
meta_res_robust <- meta_res_robust %>%
    mutate_at(c("Chromosome",
                "CpG island position",
                "Chromatin state in male skeletal muscle",
                "Chromatin state in female skeletal muscle"),
              as.factor) %>%
    arrange(`P-value`)

#Obtain list of possible genes
possible_genes <- unique(unlist(strsplit(meta_res_robust$`Annotated gene(s)`,split=";")))
possible_genes = possible_genes[!is.na(possible_genes)]

#Load meta-analysis list by CpG
L <- readRDS("./input_data/ForestplotList_new.rds")

#Obtain DMPs
DMPs <- meta_res_robust %>%
    filter(FDR < 0.005)

#Load DMRs
DMRs <- readRDS("./input_data/DMRs_new.rds")
#List of CpGs in DMRs
CpGs_in_DMRs <- DMRs$`CpGs in DMR`
DMRs <- DMRs %>%
    select(-`CpGs in DMR`) %>%
    mutate(Chromosome = as.factor(Chromosome)) %>%
    mutate_at(vars(`Maximum effect size in DMR`:
                       `Fisher multiple comparison statistic`), signif,digits = 3)

#Create a list containing DMPs, DMRs and all tested CpGs
L_summarytables <- list(DMPs = DMPs,
                        DMRs = DMRs,
                        `All tested CpGs` = meta_res_robust)
    
#Load info on each dataset
library(readxl)
waffle_data <- read_excel('./input_data/Wafflechart.xlsx',
                          sheet = "Database")
waffle_data_sum <- tibble(Dataset_ID = "Meta-analysis",
                          n = sum(waffle_data$n))
waffle_data <- bind_rows(waffle_data %>%dplyr::select(Dataset_ID,n),
                         waffle_data_sum)
waffle_data <- dplyr::rename(waffle_data,
                             `Study`=Dataset_ID)

#Create levels
levels = c("Study",
           "FUSION",
           "GeneSMART",
           "ABOS",
           "LITER",
           "GSE135063",
           "GSE49908",
           "GSE50498",
           "GSE114763",
           "EPIK",
           "GSE38291",
           "Meta-analysis")

#Load chromatin states
chromHMM <- read.delim("./input_data/Chrom_states.txt")

#Create my own theme based on theme_bw but without legend
mytheme_classic <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                             base_rect_size = base_size/22) 
{
    theme_bw(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
        theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black",
                                       size = rel(1)),
              legend.key = element_blank(), 
              strip.background = element_rect(fill = "white", 
                                              colour = "black",
                                              size = rel(2)),
              complete = TRUE)
}

#Load correspondence mRNA protein
mRNA_prot <- read_tsv("./input_data/mRNA_prot_new.txt")
mRNA_prot_graph <- ggplot(data = mRNA_prot,
                          mapping = aes(x = `Change in mRNA level per year of age (Su et al. 2015)`,
                                        y = `Change in protein level per year of age (Ubaida-Mohien et al. 2019)`,
                                        color = `Number of DMRs annotated to the gene`))+
    geom_point(mapping = aes(group=Gene),
               size = 4)+
    scale_color_gradient2(low = "white",
                         mid = "blue",
                         high = "red")+
    labs(x = "Change in mRNA level per year of age (Su et al. 2015)",
         y = "Change in protein level per year of age\n(Ubaida-Mohien et al. 2019)")+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    mytheme_classic()

# a call creating input buttons:
shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
}

subsetModal <- function(session, CpGs, size) {
    ns <- session$ns
    showModal(modalDialog({
        renderTable(ns(CpGs))
    }, size = size))
}
#################SHINY APP################################
ui <- navbarPage(theme = shinytheme("flatly"),
               collapsible = TRUE,
               title = "MetaMeth",
               id = "tabs",
               #Home page with tutorial, acknowledgment, citation, link to code and contact
               tabPanel(icon("home"),
                        fluidRow(
                            column(7,
                          tags$div(
                              "Last update: 17/03/2021",tags$br(),
                              tags$h4("Welcome to MetaMeth!"),
                              tags$p(style="text-align: justify;",
                                     "This website allows you to visualise the results of the",tags$a(href="https://www.biorxiv.org/content/10.1101/2020.09.28.315838v1","DNA methylation EWAS meta-analysis of age in human skeletal muscle conducted by Voisin et al.")
                                     ),
                              tags$br(),tags$h4("How it works"),
                              tags$ul(
                                  tags$li(
                                      tags$p(style="text-align: justify;",
                                             "If you are interested in a particular",tags$strong("CpG,"), "you may go to the", shinyLink(to = "forestplot", label = "Forest plot"), "tab to obtain a",tags$em("graph"),"summarising all studies on this CpG. To obtain all",tags$em("information"),"on a given CpG (i.e. genomic context and statistics), you may instead go to the", shinyLink(to = "summarytables", label = "Summary tables"), "tab. The",shinyLink(to = "summarytables", label = "Summary tables"),"tab contains three tables: the differentially methylated positions (DMPs) and differentially methylated regions (DMRs) associated with age at a false discovery rate (FDR) < 0.005, as well as all tested CpGs so you can filter according to your own FDR or p-value threshold. Results are downloadable as an excel or csv file.")),
                                  tags$li(
                                      tags$p(style="text-align: justify;",
                                             "If you are interested in a particular",tags$strong("gene,"), "you may go to the", shinyLink(to = "forestplot", label = "Forest plot"),"or",shinyLink(to = "summarytables", label = "Summary tables"), "tabs. Then, in the", tags$em("Annotated gene(s)"),"box or column,  enter the name of the gene to filter the CpGs, DMPs or DMRs annotated to the gene. If you are in the",shinyLink(to = "forestplot", label = "Forest plot"),"tab, you can then select which CpG to display. Alternatively, if you are in the",shinyLink(to = "summarytables", label = "Summary tables"),"tab, you can download all results pertaining to this gene as a csv or excel file.")),
                                      tags$li(
                                          tags$p(style="text-align: justify;",
                                                 "If you are interested in a particular",tags$strong("genomic region,"), "you may go to the", shinyLink(to = "forestplot", label = "Forest plot"),"or",shinyLink(to = "summarytables", label = "Summary tables"), "tabs. Then, enter the", tags$em("chromosome, position"),"in hg38 coordinates to filter the CpGs, DMPs or DMRs located within the given region. If you are in the",shinyLink(to = "forestplot", label = "Forest plot"),"tab, you can then select which CpG to display. Alternatively, if you are in the",shinyLink(to = "summarytables", label = "Summary tables"),"tab, you can download all results within this genomic region as a csv or excel file.")),    
                                  tags$li(
                                      tags$p(style="text-align: justify;",
                                             "If you are interested in the relationship between age-related DNA methylation changes and age-related mRNA or protein changes, you may go to the",shinyLink(to = "OMICsintegration", label = "OMICs integration"), "tab. We integrated the results of the EWAS meta-analysis of age with the",tags$a(href="https://skeletalmusclejournal.biomedcentral.com/articles/10.1186/s13395-015-0059-1","transcriptome meta-analysis conducted by Su et al. (2015)"),"and the large-scale",tags$a(href="https://elifesciences.org/articles/49874","proteomics study conducted by Ubaida-Mohien et al. (2019).")))),
                              tags$p(style="text-align: justify;",
                                     tags$div("Note: We are currently working on expanding the OMICS integration tab and adding a pathway enrichment tab. Coming soon!", style = "color:blue")),
                              tags$br(),tags$h4("Contributors"),
                              tags$b("Code: "),tags$a(href="https://www.vu.edu.au/research/sarah-voisin","Sarah Voisin,"),tags$a(href="https://github.com/davidruvolo51/","David Ruvolo"),"for the internal links to navigation bars in shiny, and",tags$a(href="https://github.com/strboul","Metin Yazici"),"for the pop-up windows in the",shinyLink(to = "summarytables", label = "Summary tables."),tags$br(),
                              tags$b("Advice: "),tags$a(href="https://staff.ki.se/people/nicpil", "Nicolas Pillon"),tags$br(),
                              tags$b("Feedback: "),tags$a(href="https://www.deakin.edu.au/about-deakin/people/mark-ziemann","Dr Mark Ziemann, "),tags$a(href="https://www.vu.edu.au/research/alba-moreno-asso","Dr Alba Moreno-Asso,"),"and the awesome research team of Genetics & Epigenetics of Exercise at ",tags$a(href=" https://www.vu.edu.au/research/institute-for-health-sport/mechanisms-interventions-in-health-disease", "the Institute for Health and Sport (IHES)",tags$br())
                              )
                          ),
                          column(4,
                          img(src = "Ecorche_logo.png",
                              height = 500,
                              width = 500),
                          tags$p(style="text-align: center;",
                          actionButton("citation",
                                       label = "Citation",
                                       icon = icon("quote-right"),
                                       style="color: #fff; background-color: #32cd32; border-color: #32cd32",
                                       width = "100px",
                                       onclick ="window.open(`https://www.biorxiv.org/content/10.1101/2020.09.28.315838v1/`, '_blank')"),
                          a(actionButton("contact",
                                         label = "Contact",
                                         icon = icon("envelope"),
                                         style="color: #fff; background-color: #B21212; border-color: #B21212",
                                         width = "100px"),
                          actionButton("github",
                                       label = "Code",
                                       icon = icon("github"),
                                       width = "80px",
                                       onclick ="window.open(`https://github.com/sarah-voisin/MetaMeth`, '_blank')",
                                       style="color: #fff; background-color: #767676; border-color: #767676")
                          ))
                        )
                        )
               ),
               
           #Forest plot (i.e. summary of effect size and error for each individual study) with user-input CpG name
           tabPanel(title = "Forest plot",
                    fluid = TRUE,
                    icon = icon("tree"),
                    value = "forestplot",
                    tags$head(
                        tags$style(
                            HTML(".shiny-output-error-validation
                                 {font-size: 50px}")
                            )
                        ),
                    tags$style("
                    body {
                    -moz-transform: scale(0.8, 0.8); /* Moz-browsers */
                    zoom: 0.8; /* Other non-webkit browsers */
                    zoom: 80%; /* Webkit browsers */
                    }
                               "),
                    # Sidebar layout with a input and output definitions
                    sidebarLayout(
                        sidebarPanel(
                            titlePanel("Select a CpG to plot"),
                            #First row is CpG name
                            selectizeInput(inputId = 'CpG',
                                           label = "List of possible CpGs",
                                           choices = NULL),
                            #Number of studies
                            sliderInput(inputId = "nbstudies",
                                        label = "Min number of studies",
                                        value = 6,
                                        step = 1,
                                        min = 6,
                                        max = 10),
                            hr(),
                            #Second row is filtering characteristics like gene, chromatin state, CGI position, etc.
                            titlePanel("Filter CpGs based on genomic context"),
                            fluidRow(column(3,
                                            #Select chromosome
                                            selectInput(inputId = "chr",
                                                        label = "Chromosome",
                                                        choices = c("",
                                                                    paste0("chr",c(1:22,"X","Y"))),
                                                        selected = "")
                                            ),
                                     #Select position
                                     column(4,offset = 1,
                                            numericInput(inputId = "pos_beg",
                                                         label = "From",
                                                         min = 0,
                                                         step = 1,
                                                         max = 250000000,
                                                         value = 0)
                                            ),
                                     column(4,
                                            numericInput(inputId = "pos_end",
                                                         label = "To (hg38)",
                                                         min = 0,
                                                         step = 1,
                                                         max = 250000000,
                                                         value = 250000000)
                                            )
                                     ),
                            #Select gene
                            selectizeInput(inputId = 'gene',
                                           label = "Gene",
                                           choices = NULL),
                            
                            #Select position with respect to CpG islands and CTCF and EZH2 binding sites
                            fluidRow(column(5,
                                            checkboxGroupInput(inputId = "CGI",
                                                               label = "CpG island position",
                                                               choices = c("Island","Shore","Shelf","Open sea"),
                                                               selected = c("Island","Shore","Shelf","Open sea")),
                                            actionLink("selectall_CGI","Select/Deselect All")
                                            ),
                                     column(5, offset = 1,
                                            #Select position with respect to CTCF binding site
                                            checkboxGroupInput(inputId = "CTCF",
                                                               label = "In CTCF binding site",
                                                               choices = c("No","Yes"),
                                                               selected = c("No","Yes")),
                                            #Select position with respect to EZH2 binding site
                                            checkboxGroupInput(inputId = "EZH2",
                                                               label = "In EZH2 binding site",
                                                               choices = c("No","Yes"),
                                                               selected = c("No","Yes"))
                                            )
                                     ),
                            #Chromatin states
                            fluidRow(column(5,
                                            #Male
                                            checkboxGroupInput(inputId = "chrom_state_M",
                                                               label = "Chromatin state in male skeletal muscle",
                                                               choices = chromHMM$DESCRIPTION,
                                                               selected = chromHMM$DESCRIPTION),
                                            actionLink("selectall_chromstate_M","Select/Deselect All")
                                            ),
                                     column(5,offset=1,
                                            #Female
                                            checkboxGroupInput(inputId = "chrom_state_F",
                                                               label = "Chromatin state in female skeletal muscle",
                                                               choices = chromHMM$DESCRIPTION,
                                                               selected = chromHMM$DESCRIPTION),
                                            actionLink("selectall_chromstate_F","Select/Deselect All")
                                            )
                                     ),
                               tags$hr(),
                            titlePanel("Filter CpGs based on statistics"),
                            h4("Effect size (% DNAm change per year of age)"),
                            fluidRow(column(5,
                                            #Min effect size
                                            numericInput(inputId = "minES",
                                                         label = "Min",
                                                         value = -0.5,
                                                         min = -0.5,
                                                         max = 0.5)
                                            ),
                                     column(5,
                                            #Max effect size
                                            numericInput(inputId = "maxES",
                                                         label = "Max",
                                                         value = 0.5,
                                                         min = -0.5,
                                                         max = 0.5)
                                     )
                            ),
                            h4("Significance"),
                            fluidRow(column(5,
                                            #Pvalue threshold
                                            numericInput(inputId = "pval",
                                                         label = "P-value below",
                                                         value = 1,
                                                         min = 1e-70,
                                                         max = 1)
                                            ),
                                     column(5,
                                            #FDR threshold
                                            numericInput(inputId = "FDR",
                                                         label = "FDR below",
                                                         value = 1,
                                                         min = 1e-70,
                                                         max = 1) 
                                     )
                                   ),
                            h4("Heterogeneity between studies"),
                            fluidRow(column(5,
                                            #I2
                                            numericInput(inputId = "I2",
                                                         label = "Het index (I2) below",
                                                         value = 100,
                                                         min = 0,
                                                         max = 100)
                                            ),
                            column(5,
                                   #I2 p-value
                                   numericInput(inputId = "I2pval",
                                                label = "Het p-value below",
                                                value = 1,
                                                min = 1e-30,
                                                max = 1) 
                            )
                            ),
                            tags$hr(),
                               titlePanel("Download forest plot"),
                               numericInput(inputId = "FPresolution",
                                            label = "Image resolution (ppi)",
                                            value = 72,
                                            step = 1,
                                            min = 50,
                                            max = 600),
                               helpText("As an indication, 72 ppi is standard and 300 ppi is high-quality."),
                               selectInput(inputId = "filetype",
                                            label = "File type",
                                            choices = c("jpg",
                                                        "png",
                                                        "tif"),
                                            selected = "tif"),
                               downloadButton("downloadFP",
                                              "Download")
                            ),
                    mainPanel(
                        # Forest plot ----
                        plotOutput(outputId = "forestPlot",
                                   height = "400px"),
                        hr(),
                        helpText(div(style="text-align:justify",
                                     "The meta-analysis combined results from 10 independent datasets assayed with the HumanMethylation array (27k, 450k or EPIC). Therefore, CpGs were present in some, but not all of the included studies. We limited the meta-analysis to CpGs that were present in at least 6 of the 10 studies. Studies with missing information (NA) mean that this CpG was not analysed in the dataset.")
                        )
                        )
                    )
           ),
           
           #Summary Tables
           tabPanel("Summary tables",
                    icon = icon("table"),
                    value = "summarytables",
                    selectInput(inputId = 'tabletype',
                                label = "Please choose the table to display",
                                choices = c("DMPs",
                                            "DMRs",
                                            "All tested CpGs"),
                                   selected = "DMRs"),
                    DT::dataTableOutput("summarytable")
                    ),
           
           #OMICs integration (transcriptomics and proteomics) as a scatterplot
           tabPanel(title = "OMICs integration",
                    value = "OMICsintegration",
                    icon = icon("object-group"),
                    fluidPage(
                        fluidRow(
                            column(5,
                                   plotlyOutput(outputId = "mRNAprot",
                                                width = "500px")
                            ),
                            column(5,
                                   tags$div(
                                       tags$h4("How the graph works"),
                                       tags$p(style="text-align: justify;",
                                              "This is a scatterplot showing the change in mRNA (x-axis) and protein (y-axis) for the 57 genes altered at all three omics levels. The epigenomic analysis was conducted by",tags$a(href="https://www.biorxiv.org/content/10.1101/2020.09.28.315838v1","Voisin et al. (2020)"), ", the transcriptomic analysis was conducted by", tags$a(href="https://skeletalmusclejournal.biomedcentral.com/articles/10.1186/s13395-015-0059-1","Su et al. (2015)"),"and the proteomics analysis was conducted by",tags$a(href="https://elifesciences.org/articles/49874","Ubaida-Mohien et al. (2019)."), "Each gene was colored according to the number of DMRs annotated to it, from 1-3 DMRs for most genes all the way up to 9 DMRs. Naturally, longer genes (e.g.",tags$em("NXN, ABLIM2)"),"have a greater propensity to have more DMRs given their high numbers of CpGs.")
                                   )
                                       
                            )
                        )
                    )
                    ),
           
           #Load javascript to create shinyLinks
           tags$script(src = "shinyLink.js")
)


server <- function(input, output, session) {
    
    
    # Pop-up window to show contact details
    observeEvent(input$contact, {
        showModal(modalDialog(
            title = "Contact details",
            "Sarah Voisin",tags$br(),
            "Senior Research Fellow",tags$br(),
            "Institute for Health and Sport (IHES)",tags$br(),
            "Victoria University, Footscray Park Campus, VIC 3011 Australia",tags$br(),
            tags$b("Email:"),"sarah.voisin@vu.edu.au",tags$br(),
            tags$b("Alternative email:"),"sarah.voisin.aeris@gmail.com",tags$br(),
            "Office phone: +61 3 9919 5744",tags$br(),
            "Mobile: +61 4 6646 9673",tags$br(),
            easyClose = TRUE,
            footer = NULL
        ))
    })
    
    #Insert link to tabs when explaining how the app works
    observeEvent(input$link_to_foresttab, {
        updateNavlistPanel(session, "MetaMeth", "Forest plot")
    })

        
    # Forest plot of the effect sizes and FDR for each study ----
    # with requested CpG
    # This expression that generates a histogram is wrapped in a call
    # to renderPlot to indicate that:
    #
    # 1. It is "reactive" and therefore should be automatically
    #    re-executed when inputs (input$CpG) change
    # 2. Its output type is a plot
    
    updateSelectizeInput(session,
                         'gene',
                         choices = c("",possible_genes),
                         server = TRUE,
                         selected = "")

    possible_cpgs <- reactive({
        "1==1"
        #withProgress(message = 'Updating list of CpGs', value = 0, {
            
        #Filter by gene
        #incProgress(1/9, detail ="Filtering by gene")
        if (input$gene %in% possible_genes)
        {
            genetolookfor <- paste0("\\b",input$gene,"\\b")
            indexgene <- grep(genetolookfor,
                          meta_res_robust$`Annotated gene(s)`)
        }
        else
        {
            indexgene <- 1:nrow(meta_res_robust)
        }
            
        #Filter by position
        #incProgress(1/9, detail ="Filtering by position")
        if (input$chr!="")
        {
            indexpos <- which(meta_res_robust$`Chromosome`==input$chr&
                               meta_res_robust$`Position (hg38)`>=input$pos_beg &
                               meta_res_robust$`Position (hg38)`<=input$pos_end)   
        }
        else
        {
            indexpos <- 1:nrow(meta_res_robust)
        }
        
        #Filter by CGI
            #incProgress(1/9, detail ="Filtering by CGI position")
            indexCGI <- which(meta_res_robust$`CpG island position` %in% input$CGI)
            
        #Filter by CTCF and EZH2 binding sites
        #incProgress(1/9, detail ="Filtering by TF binding site")
            indexCTCF <- which(meta_res_robust$`In CTCF binding site in HSMMtube` %in% input$CTCF)    
            indexEZH2 <- which(meta_res_robust$`In EZH2 binding site in HSMMtube` %in% input$EZH2)    
 
        #Filter by male and female chrom state
        #incProgress(1/9, detail ="Filtering by chromatin state")
        indexchromstateM <- which(meta_res_robust$`Chromatin state in male skeletal muscle`%in%input$chrom_state_M)
        indexchromstateF <- which(meta_res_robust$`Chromatin state in female skeletal muscle`%in%input$chrom_state_F)
        
        #Filter by Effect Size, p-value & FDR
        #incProgress(1/9, detail ="Filtering by statistics")
        indexstats <- which(meta_res_robust$`Effect size`>input$minES&
                             meta_res_robust$`Effect size`<input$maxES&
                             meta_res_robust$`P-value`<input$pval&
                             meta_res_robust$FDR<input$FDR&
                                meta_res_robust$`Heterogeneity index (I2)`<input$I2&
                                meta_res_robust$`Heterogeneity p-value`<input$I2pval&
                                meta_res_robust$`Number of studies`>=input$nbstudies)
            
        #Obtain intersection
        #incProgress(1/9, detail ="Loading filtered list of CpGs")
        indexlist <- list(indexgene,
                          indexpos,
                          indexCGI,
                          indexCTCF,
                          indexEZH2,
                          indexchromstateM,
                          indexchromstateF,
                          indexstats)
        indextot <- Reduce(intersect,indexlist)
        
        pull(meta_res_robust[indextot,"CpG"])
        
        #})
    })

    #Update list of CpGs
    observeEvent(possible_cpgs(),
        {
        updateSelectizeInput(session,
                             'CpG',
                             choices = c("",possible_cpgs()),
                             selected = "",
                             server = TRUE)
        })
    
    #Get the select/unselect all for CGI & chromatin states
    observe({
        if (input$selectall_CGI > 0) {
            if (input$selectall_CGI %% 2 == 0){
                updateCheckboxGroupInput(session=session, 
                                         inputId="CGI",
                                         choices = c("Island","Shore","Shelf","Open sea"),
                                         selected = c("Island","Shore","Shelf","Open sea"))
                
            } else {
                updateCheckboxGroupInput(session=session, 
                                         inputId="CGI",
                                         choices = c("Island","Shore","Shelf","Open sea"),
                                         selected = c())
                
            }}
    })
    
    observe({
        if (input$selectall_chromstate_M > 0) {
            if (input$selectall_chromstate_M %% 2 == 0){
                updateCheckboxGroupInput(session=session, 
                                         inputId="chrom_state_M",
                                         choices = chromHMM$DESCRIPTION,
                                         selected = chromHMM$DESCRIPTION)
                
            } else {
                updateCheckboxGroupInput(session=session, 
                                         inputId="chrom_state_M",
                                         choices = chromHMM$DESCRIPTION,
                                         selected = c())
                
            }}
    })

    
    observe({
        if (input$selectall_chromstate_F > 0) {
            if (input$selectall_chromstate_F %% 2 == 0){
                updateCheckboxGroupInput(session=session, 
                                         inputId="chrom_state_F",
                                         choices = chromHMM$DESCRIPTION,
                                         selected = chromHMM$DESCRIPTION)
                
            } else {
                updateCheckboxGroupInput(session=session, 
                                         inputId="chrom_state_F",
                                         choices = chromHMM$DESCRIPTION,
                                         selected = c())
                
            }}
    })
    
    #Render the forestplot in the right tab
    output$forestPlot <- renderPlot({
        
        validate(
            need(input$CpG,
                 'The forest plot will appear as soon as you supply a valid CpG name')
        )
    
        
        #Select the right CpG
        tib <- L[[input$CpG]]
        
        
        #Add n and factor each Study to order properly
        tib <- inner_join(tib,
                          waffle_data) %>%
            mutate(n = replace(n,
                               n==908,
                               sum(n[Study!="Meta-analysis"])))
        #Add lines of empty cells for missing datasets + one extra for title
        missing_datasets = setdiff(levels,tib$Study)
        tib <- tib %>%
            tibble::add_row(Study = missing_datasets)
        tib$Study[is.na(tib$Study)]="Study"
        tib <- tib %>%
            mutate(Study = factor(Study,
                                  levels=rev(levels))) %>%
            arrange(Study)
        
        #Plot with ggplot2
        library(gridExtra)
        theme_set(theme_bw(base_size=10))
        
        #Create the forest plot
        plot1<-ggplot(data=tib,
                      aes(x=Study,
                          y=ES,
                          ymax=ES+(1.96*SE),
                          ymin=ES-(1.96*SE),
                          size=factor(Type),
                          colour=factor(Type)))+
            geom_pointrange()
        plot2<-plot1+
            coord_flip()+
            geom_hline(yintercept =0,
                       lty=2,
                       size=1)+
            scale_size_manual(values=c(1.5,2))
        plot3<-plot2+
            xlab("Study")+
            ylab("% DNAm change per year of age")+
            scale_colour_manual(values=c("gray63","black"))+
            theme(legend.position="none",
                  axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.title.x=element_text(size=16,face="bold"),
                  axis.text.x=element_text(size=14),
                  axis.ticks.y=element_blank())
        
        
        #Add table next to forest plot with Study and FDR
        SummaryTable <- tib %>%
            dplyr::select(Study,
                          n,
                          ES,
                          PVAL,
                          FDR) %>%
            mutate_all(as.character)
        SummaryTable$n[SummaryTable$Study=="Study"] = "Sample size"
        SummaryTable$ES[SummaryTable$Study=="Study"] = "Effect size"
        SummaryTable$PVAL[SummaryTable$Study=="Study"] = "P-value"
        SummaryTable$FDR[SummaryTable$Study=="Study"] = "FDR"
        
        SummaryTable <- SummaryTable %>%
            mutate(Index = 1:nrow(SummaryTable)) %>%
            pivot_longer(cols = c(Study:FDR),
                         names_to = "variable",
                         values_to = "value") %>%
            mutate(variable = factor(variable,
                                     levels = c("Study",
                                                "n",
                                                "ES",
                                                "PVAL",
                                                "FDR")))
        
        #Create table plot
        data_table <- ggplot(SummaryTable,
                             aes(x = variable,
                                 y = Index,
                                 label = format(value,
                                                digits = 2,
                                                nsmall = 1))) +
            geom_text(size = 6,
                      hjust = 0) +
            theme_bw() +
            geom_hline(yintercept=(nrow(tib)-0.4)) + 
            labs(x="",y="")+
            ylim(-0.25,nrow(tib))+
            theme(legend.position="none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.title=element_blank(),
                  axis.text = element_blank(),
                  axis.ticks=element_blank())
        
        #Save forest plot as a function to call later
        grid.arrange(data_table,
                     plot3,
                     widths = c(2, 1.25),
                     ncol=2)
    
    })
    
    #Downloading the Forest plot
    output$downloadFP <- downloadHandler(
        filename = function() {
            paste0(input$CpG, "_forestplot.",input$filetype)
        },
        content = function(file) {
            #Select the right CpG
            tib <- L[[input$CpG]]
            
            #Add n and factor each Study to order properly
            tib <- inner_join(tib,
                              waffle_data) %>%
                mutate(n = replace(n,
                                   n==908,
                                   sum(n[Study!="Meta-analysis"])))
            #Add lines of empty cells for missing datasets + one extra for title
            missing_datasets = setdiff(levels,tib$Study)
            tib <- tib %>%
                tibble::add_row(Study = missing_datasets)
            tib$Study[is.na(tib$Study)]="Study"
            tib <- tib %>%
                mutate(Study = factor(Study,
                                      levels=rev(levels))) %>%
                arrange(Study)
            
            #Plot with ggplot2
            library(gridExtra)
            theme_set(theme_bw(base_size=10))
            
            #Create the forest plot
            plot1<-ggplot(data=tib,
                          aes(x=Study,
                              y=ES,
                              ymax=ES+(1.96*SE),
                              ymin=ES-(1.96*SE),
                              size=factor(Type),
                              colour=factor(Type)))+
                geom_pointrange()
            plot2<-plot1+
                coord_flip()+
                geom_hline(yintercept =0,
                           lty=2,
                           size=1)+
                scale_size_manual(values=c(1.5,2))
            plot3<-plot2+
                xlab("Study")+
                ylab("M-value change per year of age")+
                scale_colour_manual(values=c("gray63","black"))+
                theme(legend.position="none",
                      axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.title.x=element_text(size=16,face="bold"),
                      axis.text.x=element_text(size=14),
                      axis.ticks.y=element_blank())
            
            
            #Add table next to forest plot with Study and FDR
            SummaryTable <- tib %>%
                dplyr::select(Study,
                              n,
                              ES,
                              PVAL,
                              FDR) %>%
                mutate_all(as.character)
            SummaryTable$n[SummaryTable$Study=="Study"] = "Sample size"
            SummaryTable$ES[SummaryTable$Study=="Study"] = "Effect size"
            SummaryTable$PVAL[SummaryTable$Study=="Study"] = "P-value"
            SummaryTable$FDR[SummaryTable$Study=="Study"] = "FDR"
            
            SummaryTable <- SummaryTable %>%
                mutate(Index = 1:nrow(SummaryTable)) %>%
                pivot_longer(cols = c(Study:FDR),
                             names_to = "variable",
                             values_to = "value") %>%
                mutate(variable = factor(variable,
                                         levels = c("Study",
                                                    "n",
                                                    "ES",
                                                    "PVAL",
                                                    "FDR")))
            
            #Create table plot
            data_table <- ggplot(SummaryTable,
                                 aes(x = variable,
                                     y = Index,
                                     label = format(value,
                                                    digits = 2,
                                                    nsmall = 1))) +
                geom_text(size = 6,
                          hjust = 0) +
                theme_bw() +
                geom_hline(yintercept=(nrow(tib)-0.4)) + 
                labs(x="",y="")+
                ylim(-0.25,nrow(tib))+
                theme(legend.position="none",
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.title=element_blank(),
                      axis.text = element_blank(),
                      axis.ticks=element_blank())
            
            #Save forest plot
            if (input$filetype=="tif")
            {
            tiff(file,
                 width =12,
                 height = 6,
                 units = 'in',
                 res=input$FPresolution)
            grid.arrange(data_table,
                         plot3,
                         widths = c(2, 1.25),
                         ncol=2)
            dev.off()
            }
            else if (input$filetype=="png")
            {
                png(file,
                     width =12,
                     height = 6,
                     units = 'in',
                     res=input$FPresolution)
                grid.arrange(data_table,
                             plot3,
                             widths = c(2, 1.25),
                             ncol=2)
                dev.off()
            }
            else
            {
                jpeg(file,
                    width =12,
                    height = 6,
                    units = 'in',
                    res=input$FPresolution)
                grid.arrange(data_table,
                             plot3,
                             widths = c(2, 1.25),
                             ncol=2)
                dev.off()
            }
        }
    )
    
    #Create action button within DMR table
    values <- reactiveValues(df = NULL)
    
    values$df <- tibble(
        DMRs[,1:4],
        CpGs = shinyInput(actionButton, nrow(DMRs), 'button_', label = "Show CpGs in DMR",
                             onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
        DMRs[,5:ncol(DMRs)]
    )
    
    #Summary tables
    output$summarytable <- DT::renderDataTable({
        if (input$tabletype=="DMRs")
        {
            DT::datatable(isolate(values$df),
                          extensions = 'Buttons',
                          options = list(
                              dom = 'lBtip',
                              lengthMenu = c(10, 25, 50, nrow(DMRs)),
                              buttons = list(
                                  list(extend = 'csv', filename = "DMRs"),
                                  list(extend = 'excel', filename = "DMRs"))
                          ),
                          rownames = FALSE,
                          escape = FALSE,
                          selection = 'none',
                          filter = 'top')  
        }
        else
        {
                  DT::datatable(L_summarytables[[input$tabletype]],
                      extensions = 'Buttons',
                      options = list(
                          dom = 'lBtip',
                          lengthMenu = c(10, 25, 50, nrow(L_summarytables[[input$tabletype]])),
                          buttons = list(
                              list(extend = 'csv', filename = input$tabletype),
                              list(extend = 'excel', filename = input$tabletype))
                      ),
                      rownames = FALSE,
                      filter = 'top')  
        }
    })
    
    #Do the pop up window
    observeEvent(input$select_button, {
        selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
        subsetted <- CpGs_in_DMRs[[selectedRow]]
        subsetModal(session, CpGs = subsetted, size = "m")
    })
    
    output$mRNAprot <- renderPlotly({
    ggplotly(mRNA_prot_graph,
             tooltip = c("group"))
    })
}

shinyApp(ui = ui,
         server = server)

