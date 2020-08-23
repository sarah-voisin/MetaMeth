#EWAS shiny app
#Load packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(readxl)
library(rintrojs)
library(pals)
library(plotly)

#Load list of CpGs
#directory = "F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/EWAS meta-analysis of age/MetaMeth"
#setwd(directory)
#Load meta-analysis list by CpG
load("./input_data/ForestplotList.rda")
#Load all results
meta_res_robust <- read_tsv("./input_data/MetaAnalysis.txt")
meta_res_robust <- meta_res_robust %>%
    mutate_at(c("Chromosome",
                "CpG island position",
                "Chromatin state in male skeletal muscle",
                "Chromatin state in female skeletal muscle"),
              as.factor) %>%
    mutate_at(c("P-value","FDR"), signif, digits = 3) %>%
    select(-c(Significance,Direction)) %>%
    arrange(`P-value`)

#Create DMPs
DMPs <- meta_res_robust %>%
    filter(FDR < 0.005)

tested_cpgs <- meta_res_robust$CpG

#Load DMRs
DMRs <- read_tsv("./input_data/DMRs.txt")
DMRs <- DMRs %>%
    mutate(Chromosome = as.factor(Chromosome)) %>%
    mutate_at(vars(`Maximum effect size in DMR (M-value change per year of age)`:
                       `Fisher multiple comparison statistic`), signif,digits = 3)

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
           "Thomis",
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
Table1 <- read_tsv("./input_data/mRNA_prot.txt")
mRNA_prot_graph <- ggplot(data = Table1,
                          mapping = aes(x = `Change in mRNA level per year of age (Su et al. 2015)`,
                                        y = `Change in protein level per year of age (Ubaida-Mohien et al. 2019)`,
                                        size = `Number of DMPs annotated to the gene`,
                                        color = `Number of DMPs annotated to the gene`)
)+
    geom_point(mapping = aes(group=Gene))+
    labs(x = "Change in mRNA level per year of age (Su et al. 2015)",
         y = "Change in protein level per year of age\n(Ubaida-Mohien et al. 2019)")+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    mytheme_classic()


#################SHINY APP################################
ui <- navbarPage(theme = shinytheme("flatly"),
               collapsible = TRUE,
               "MetaMeth",
               
               #1st tab is going to be "About" to explain the page
               tabPanel(icon("home"),
                        fluidRow(
                            column(7,
                          tags$div(
                              tags$h4("Welcome to MetaMeth!"),
                              tags$p(style="text-align: justify;",
                                     "This website allows you to visualise the results of the EWAS meta-analysis of age in human skeletal muscle conducted by Voisin et al."),
                              tags$br(),tags$p(style="text-align: justify;",
                                     "Original results are often delivered in a dry, static form that does not engage the reader. Thus, we wanted to create a user-friendly, interactive way to explore the results of our tedious analysis."),
                              tags$h4("How it works"),
                              tags$p(style="text-align: justify;",
                                     "If you are interested in a particular",tags$b("CpG,"),"you may go to the Forest plot tab to obtain the summary of the meta-analysis of age for said CpG.
                                     If you are interested in a particular",tags$b("gene,"),"you may go to the DMP or DMR tab. Then, in the ", tags$em("Annotated gene(s)")," column,  enter the name of said gene to filter the DMPs or DMRs annotated to said gene. You can then download the results in an excel or csv file.
                                     The DMP and DMR tabs only contain information for those CpGs and regions that were significantly associated with age at FDR < 0.005. For a full list of the tested CpGs, go to the All CpGs tab and filter according to your personal FDR or p-value threshold."),
                              tags$br(),tags$h4("Contributors"),
                              tags$b("Code: "),"Sarah Voisin",tags$br(),
                              tags$b("Server host: "),tags$a(href="https://www.vu.edu.au/", "Victoria University"),tags$br(),
                              tags$b("Feedback: "),"the awesome research team of Genetics & Epigenetics of Exercise at ",tags$a(href=" https://www.vu.edu.au/research/institute-for-health-sport/mechanisms-interventions-in-health-disease", "the Institute for Health and Sport (IHES)"),tags$br(),
                              tags$br(),tags$h4("Citation"),
                              tags$p(style="text-align: justify;",
                                     "If you use MetaMeth, please cite ",tags$a(href="https://onlinelibrary.wiley.com/doi/full/10.1002/jcsm.12556", "Voisin et al. 2020, JCSM")),tags$br(),
                              a(actionButton("contact",
                                           label = "Contact",
                                           icon = icon("envelope"),
                                           style="color: #fff; background-color: #B21212; border-color: #B21212",
                                           width = "100px"),
                                href="mailto:sarah.voisin.aeris@gmail.com"),
                              actionButton("github",
                                           label = "Code",
                                           icon = icon("github"),
                                           width = "80px",
                                           onclick ="window.open(`https://github.com/sarah-voisin`, '_blank')",
                                           style="color: #fff; background-color: #767676; border-color: #767676"))),
                          column(4,
                          img(src = "Ecorche_logo.png",
                              height = 500,
                              width = 500)
                          )
                        )
               ),
           #2nd tab is to enter CpG
           tabPanel(title = "Forest plot",
                    fluidPage(
                        fluidRow(
                        # Sidebar panel for inputs ----
                        column(2,
                               #textInput(inputId = "CpG",
                                         #label = "Please enter your CpG",
                                         #value = "cg15456742")
                               selectizeInput(inputId = 'CpG',
                                              label = "Please enter your CpG",
                                              choices = NULL,
                                              selected = "cg15456742")
                            ),
                        # Volcano plot ----
                        #column(4,
                               # Output: Forest Plot ----
                               #plotOutput(outputId = "volcanoPlot",
                                         # brush = "volcanoselect"),
                               #verbatimTextOutput("CpGinfo")),
                        # Forest plot ----
                        column(8,
                               # Output: Forest Plot ----
                                  plotOutput(outputId = "forestPlot",
                                             height = "800px"))
                        )
                        )
                    ),
           #3rd tab are DMPs
           tabPanel("DMPs",
                    #sidebarLayout(
                        # Sidebar panel for inputs ----
                        #sidebarPanel(width = 2,
                            # Input: Slider for the number of bins ----
                            #textInput(inputId = "Gene",
                                      #label = "Please enter your gene",
                                      #value = "HDAC4")),
                        # Main panel for displaying outputs ----
                        #mainPanel(# Output: summary ----
                                  #width = 5,
                                  DT::dataTableOutput("DMPs")
                    #)
                   # )
           ),
           #4th tab are DMRs
           tabPanel("DMRs",
                    DT::dataTableOutput("DMRs",
                                        width = "70%")
                    ),
           #5th tab is OMICs integration
           #tabPanel("OMICs integration",
                    #plotOutput(outputId = "volcanoPlot",
                              # brush = "volcanoselect"),
                    #verbatimTextOutput("CpGinfo")
           #)
           #Last tab is ALL results
           tabPanel("All CpGs",
                    DT::dataTableOutput("summary")
           ),
           tabPanel("OMICs integration",
                    plotlyOutput(outputId = "mRNAprot",
                                 width = "500px"))
)


server <- function(input, output, session) {
    
    
    # Volcano plot
    
    #Insert link to tabs when explaining how the app works
    observeEvent(input$link_to_foresttab, {
        updateNavlistPanel(session, "MetaMeth", "Forest plot")
    })
    #output$volcanoPlot <- renderPlot({
        #volcano
    #})
    
    #output$CpGinfo <- renderPrint({
        # With base graphics, need to tell it what the x and y variables are.
        #brushedPoints(meta_res_robust,
                      #input$volcanoselect)
    #})
        
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

    
    
    output$forestPlot <- renderPlot({
        
        validate(
            need(input$CpG, 'The forest plot will appear as soon as you supply a valid CpG name')
        )
        #tib <- L[[CpG_forest]]
        tib <- L[[input$CpG]]
        
        #Add n and factor each Study to order properly
        tib <- tib %>%
            mutate_if(is.numeric,
                      signif,
                      digits = 2)
        tib <- inner_join(tib,
                          waffle_data) %>%
            mutate(n = replace(n,
                               n==889,
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
                          FDR) %>%
            mutate_all(as.character)
        SummaryTable$n[SummaryTable$Study=="Study"] = "Sample size"
        SummaryTable$ES[SummaryTable$Study=="Study"] = "Effect size"
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
    
    #Create forest plot
    grid.arrange(data_table,
                 plot3,
                 widths = c(1.5, 1),
                 heights = c(1,1),
                 ncol=2)
    })
    
    
    output$DMPs <- DT::renderDataTable({
        DT::datatable(DMPs,
                      extensions = 'Buttons',
                      options = list(
                          dom = 'lBtip',
                          lengthMenu = c(10, 25, 50, nrow(DMPs)),
                          buttons = c('copy', 'csv', 'excel')
                      ),
                      rownames = FALSE,
                      filter = 'top')
    })
    
    output$DMRs <- DT::renderDataTable({
        DT::datatable(DMRs,
                      extensions = 'Buttons',
                      options = list(
                          dom = 'lBtip',
                          lengthMenu = c(10, 25, 50, nrow(DMRs)),
                          buttons = c('copy', 'csv', 'excel')
                      ),
                      rownames = FALSE,
                      filter = 'top')
    })
    
    
    output$summary <- DT::renderDataTable({
        #gene_toget <- paste0("\\b",input$Gene,"\\b")
        #indexes <- grep(gene_toget,meta_res_robust$`Annotated gene(s)`)
        #meta_res_robust_toshow <- meta_res_robust[indexes,]
        DT::datatable(meta_res_robust,
                      rownames = FALSE,
                      filter = 'top',
                      extensions = 'Buttons',
                      options = list(
                          autoWidth = TRUE,
                          dom = 'lBtip',
                          lengthMenu = c(10, 20, 50, 100, nrow(meta_res_robust)),
                          buttons = c('copy', 'csv', 'excel')
                      )
        )
    })
    
    output$mRNAprot <- renderPlotly({
    ggplotly(mRNA_prot_graph,
             tooltip = c("group", "colour"))
    })
}

shinyApp(ui = ui,
         server = server)



