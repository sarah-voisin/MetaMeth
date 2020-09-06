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
meta_res_robust <- read_tsv("./input_data/MetaAnalysis.txt")
meta_res_robust <- meta_res_robust %>%
    mutate_at(c("Chromosome",
                "CpG island position",
                "Chromatin state in male skeletal muscle",
                "Chromatin state in female skeletal muscle"),
              as.factor) %>%
    arrange(`P-value`)

#Obtain list of tested CpGs
tested_cpgs <- meta_res_robust$CpG

#Load meta-analysis list by CpG
L <- readRDS("./input_data/ForestplotList.rds")

#Obtain DMPs
DMPs <- meta_res_robust %>%
    filter(FDR < 0.005)

#Load DMRs
DMRs <- read_tsv("./input_data/DMRs.txt")
DMRs <- DMRs %>%
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
mRNA_prot <- read_tsv("./input_data/mRNA_prot.txt")
mRNA_prot_graph <- ggplot(data = mRNA_prot,
                          mapping = aes(x = `Change in mRNA level per year of age (Su et al. 2015)`,
                                        y = `Change in protein level per year of age (Ubaida-Mohien et al. 2019)`,
                                        color = `Number of DMRs annotated to the gene`)
)+
    geom_point(mapping = aes(group=Gene),
               size = 4)+
    labs(x = "Change in mRNA level per year of age (Su et al. 2015)",
         y = "Change in protein level per year of age\n(Ubaida-Mohien et al. 2019)")+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    mytheme_classic()


#################SHINY APP################################
ui <- navbarPage(theme = shinytheme("flatly"),
               collapsible = TRUE,
               "MetaMeth",
               
               #Home page with tutorial, acknowledgment, citation, link to code and contact
               tabPanel(icon("home"),
                        fluidRow(
                            column(7,
                          tags$div(
                              "Last update:",Sys.Date(),tags$br(),
                              tags$h4("Welcome to MetaMeth!"),
                              tags$p(style="text-align: justify;",
                                     "This website allows you to visualise the results of the DNA methylation EWAS meta-analysis of age in human skeletal muscle conducted by Voisin et al."),
                              tags$br(),tags$h4("How it works"),
                              tags$ul(
                                  tags$li(
                                      tags$p(style="text-align: justify;",
                                             "If you are interested in a particular",tags$b("CpG,"), "you may go to the", shinyLink(to = "forestplot", label = "forest plot"), "tab to obtain the summary of the meta-analysis of age for said CpG. To obtain additional information on a given CpG (e.g. chromatin states, heterogeneity index), you may instead go to the", shinyLink(to = "summarytables", label = "Summary tables"), "tab. The",shinyLink(to = "summarytables", label = "Summary tables"),"tab contains three tables: the differentially methylated positions (DMPs) and differentially methylated regions (DMRs) associated with age at FDR < 0.005, as well as all tested CpGs so you can filter according to your personal FDR or p-value threshold.")),
                                  tags$li(
                                      tags$p(style="text-align: justify;",
                                             "If you are interested in a particular",tags$b("gene,"), "you may go to the", shinyLink(to = "summarytables", label = "Summary tables"), "tab. Then, in the ", tags$em("Annotated gene(s)")," column,  enter the name of said gene to filter the DMPs or DMRs annotated to said gene. You can then download the results in an excel or csv file.")),
                                  tags$li(
                                      tags$p(style="text-align: justify;",
                                             "If you are interested in the relationship between age-related DNA methylation changes and age-related mRNA or protein changes, you may go to the",shinyLink(to = "OMICsintegration", label = "OMICs integration"), "tab. We integrated the results of the EWAS meta-analysis of age with the",tags$a(href="https://skeletalmusclejournal.biomedcentral.com/articles/10.1186/s13395-015-0059-1","transcriptome meta-analysis conducted by Su et al. (2015)"),"and the large-scale",tags$a(href="https://elifesciences.org/articles/49874","proteomics study conducted by Ubaida-Mohien et al. (2019).")))),
                              tags$br(),tags$h4("Contributors"),
                              tags$b("Code: "),tags$a(href="https://www.vu.edu.au/research/sarah-voisin","Sarah Voisin,"),tags$a(href="https://github.com/davidruvolo51/","David Ruvolo"),"for the internal links to navigation bars in shiny",tags$br(),
                              tags$b("Advice: "),tags$a(href="https://staff.ki.se/people/nicpil", "Nicolas Pillon"),tags$br(),
                              tags$b("Feedback: "),tags$a(href="https://www.vu.edu.au/research/alba-moreno-asso","Dr Alba Moreno-Asso,"),tags$a(href="https://www.vu.edu.au/research/nir-eynon","A/Prof Nir Eynon"), "and the awesome research team of Genetics & Epigenetics of Exercise at ",tags$a(href=" https://www.vu.edu.au/research/institute-for-health-sport/mechanisms-interventions-in-health-disease", "the Institute for Health and Sport (IHES)",tags$br())
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
                                       onclick ="window.open(`https://pubmed.ncbi.nlm.nih.gov/32067420/`, '_blank')"),
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
                    value = "forestplot",
                    fluidPage(
                        fluidRow(
                            
                        #Sidebar panel for CpG input
                        column(2,
                               wellPanel(
                                   selectizeInput(inputId = 'CpG',
                                              label = "Please enter your CpG",
                                              choices = NULL,
                                              selected = "cg00702638"),
                               tags$br(),tags$h4("Download forest plot"),
                               numericInput(inputId = "FPresolution",
                                            label = "Image resolution (ppi)",
                                            value = 72,
                                            min = 50,
                                            max = 600),
                               tags$p(style="text-align: justify;",
                                      "As an indication, 72 ppi is standard and 300 ppi is high-quality."),
                               selectInput(inputId = "filetype",
                                            label = "File type",
                                            choices = c("jpg",
                                                        "png",
                                                        "tif"),
                                            selected = "tif"),
                               downloadButton("downloadFP",
                                              "Download")
                               )
                            ),
                        # Forest plot ----
                        column(10,
                               # Output: Forest Plot ----
                                  plotOutput(outputId = "forestPlot",
                                             height = "400px"),
                                   tags$h4("Legend"),
                                   tags$p(style="text-align: justify;",
                                          "The meta-analysis combined results from 10 independent datasets assayed with the HumanMethylation array (27k, 450k or EPIC). Therefore, CpGs were present in some, but not all of the included studies. We limited the meta-analysis to CpGs that were present in at least 6 of the 10 studies. Studies that were missing the input CpG are displayed as NA on the forest plot.")
                        )
                        )
                    )
                    ),
           
           #Summary Tables
           tabPanel("Summary tables",
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
                                              "This is a scatterplot showing the change in mRNA (x-axis) and protein (y-axis) for the 59 genes altered at all three omics levels. The epigenomic analysis was conducted by",tags$a(href="https://onlinelibrary.wiley.com/doi/full/10.1002/jcsm.12556","Voisin et al. (2020)"), ", the transcriptomic analysis was conducted by", tags$a(href="https://skeletalmusclejournal.biomedcentral.com/articles/10.1186/s13395-015-0059-1","Su et al. (2015)"),"and the proteomics analysis was conducted by",tags$a(href="https://elifesciences.org/articles/49874","Ubaida-Mohien et al. (2019)."), "Each gene was colored according to the number of DMRs annotated to it, from 1-3 DMRs for most genes all the way up to 12 DMRs. Naturally, longer genes (e.g.",tags$em("NXN, ABLIM2)"),"have a greater propensity to have more DMRs given their high numbers of CpGs.")
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
                             'CpG',
                             choices = tested_cpgs,
                             server = TRUE)

 
    
    #Render the forestplot in the right tab
    output$forestPlot <- renderPlot({
        
        validate(
            need(input$CpG, 'The forest plot will appear as soon as you supply a valid CpG name')
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
    
    #Summary tables
    output$summarytable <- DT::renderDataTable({
        DT::datatable(L_summarytables[[input$tabletype]],
                      extensions = 'Buttons',
                      options = list(
                          dom = 'lBtip',
                          lengthMenu = c(10, 25, 50, nrow(L[[input$tabletype]])),
                          buttons = c('copy', 'csv', 'excel')
                      ),
                      rownames = FALSE,
                      filter = 'top')
    })
    
 
    output$mRNAprot <- renderPlotly({
    ggplotly(mRNA_prot_graph,
             tooltip = c("group", "colour"))
    })
}

shinyApp(ui = ui,
         server = server)

