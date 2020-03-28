library(shiny)
library(shinyWidgets)
require(tidyverse)
require(annotables) #devtools
require(grImport)
require(BuenColors)  #devtools
require(cowplot)
require(edgeR) 
library("ggplot2")
library("magrittr")
library(shinythemes)
# options(repos = BiocManager::repositories())

genes <-annotables::grch38 %>%
  dplyr::select(symbol)

ui <- fluidPage(
  header <- shinydashboard::dashboardHeader(title = "" )  ,
  titlePanel('', windowTitle = 'CB_PB'),
  dropdown(
    
    selectInput(inputId = 'genes',
                label = 'Select Gene:',
                choices = genes),
    
    selectInput(inputId = 'type',label = 'Select Type', choices = c('count','cpm','log2(cpm+1)'), selected = 'count'),
    selectInput(inputId = 'palette',label = 'Select Palette', choices = c('solar_rojos', 'white_grove', 'white_mango','white_tango'), selected = 'solar_rojor'),
    
    style = "unite", icon = icon("gear"),
    status = "danger", width = "300px",
    animate = animateOptions(
      enter = animations$fading_entrances$fadeInLeftBig,
      exit = animations$fading_exits$fadeOutRightBig
    )
  ),
  
  plotOutput(outputId = 'plot2', width = 1000, height = 500), 
  downloadBttn(outputId = 'download', color = 'success', label = 'Save')
  
)

server <- function(input, output, session) {
  
  plot_fig <- function() {
    palette = input$palette
    output_dir = "."
    type = input$type
    if (type == "log2(cpm+1)"){type = "log2_cpm"}
    myshape<- readPicture("cells.xml")
    genes = input$genes[1]
    n <- length(genes)
    inputData = input$rds
    datapath = inputData$datapath
    df <- import_data(type = type, df_path = datapath)
    df <- anno_grch38(df)
    df <- filter_genes_long(df, genes)
    df <- generate_plot_info(df, palette, n)[[1]] %>% distinct()
    legend <- generate_plot_info(df, palette, n)[[2]]
    
    for(symb in unique(df$symbol)){
      curr <- df %>% filter(symbol == symb)
      my_shape <- color_cells(myshape, curr) %>% pictureGrob(.) %>% ggdraw(.)
      p<-my_shape + annotation_custom(legend, xmin =.75, ymin = 1, x=1, y=0.7) + annotate("text",x=0.07,y = 0.90, label = symb) 
      
    }
    p
  }
  output$plot2 <- renderPlot(plot_fig())
  output$download <- downloadHandler(filename = "file.pdf", 
                                     content = function(file) {
                                       ggsave(filename= file)
                                       print(plot_fig())
                                       dev.off()
                                     }
  )
  # generate color graphic information using ggplot
  generate_plot_info <- function(df, palette, n_genes, type){
    palette = jdb_palette(palette)
    df <- data.frame(df)
    if(n_genes == 1){
      p1 <- ggplot(df, aes(x = cell_type, y = expression, fill = expression)) +
        geom_bar(stat = "identity") + 
        coord_flip() +
        pretty_plot() + 
        scale_fill_gradientn(colors = palette) + labs(fill = input$type)
    }
    
    
    df <- ggplot_build(p1)[[1]] %>%
      data.frame(.) %>% 
      dplyr::select(1,3) %>%
      left_join(df, ., by = c("expression" = "y"))
    legend <- get_legend(p1)
    list(df, legend)
  }  
  
  # import data
  import_data <- function(type, df_path){
    df <- readRDS("count.rds")
    if(type == "log2_cpm" || type == "cpm"){
      df[-1] <- lapply(df[-1], cpm)
      if(type == "log2_cpm"){
        df[-1] <- lapply(df[-1], function(x) log2(x + 1))
      }
    }
    df
  }
  
  # annotate the dataframe
  anno_grch38 <- function(df){
    annotables::grch38 %>%
      mutate(ensg = ensgene) %>%
      dplyr::select(ensg, symbol) %>%
      right_join(., df)
  }
  
  # subset df by genes
  filter_genes_long <- function(df, genes){
    df %>% 
      filter(symbol %in% genes | ensg %in% genes) %>%
      gather(cell_type, expression, -ensg, -symbol)
  }
  
  
  # function to change colors of cells
  color_cells <- function(picture, df){
    my_shape = picture
    CB_BFUE <- df %>% filter(cell_type == "CB_BFUE") %>% pull(fill)
    CB_CD34 <- df %>% filter(cell_type == "CB_CD34") %>% pull(fill)
    CB_CFUE <- df %>% filter(cell_type == "CB_CFUE") %>% pull(fill)
    CB_eBaso <- df %>% filter(cell_type == "CB_eBaso") %>% pull(fill)
    CB_lBaso <- df %>% filter(cell_type == "CB_lBaso") %>% pull(fill)
    CB_ortho <- df %>% filter(cell_type == "CB_ortho") %>% pull(fill)
    CB_poly <- df %>% filter(cell_type == "CB_poly") %>% pull(fill)
    CB_pro <- df %>% filter(cell_type == "CB_pro") %>% pull(fill)
    
    PB_BFUE <- df %>% filter(cell_type == "PB_BFUE") %>% pull(fill)
    PB_CD34 <- df %>% filter(cell_type == "PB_CD34") %>% pull(fill)
    PB_CFUE <- df %>% filter(cell_type == "PB_CFUE") %>% pull(fill)
    PB_eBaso <- df %>% filter(cell_type == "PB_eBaso") %>% pull(fill)
    PB_lBaso <- df %>% filter(cell_type == "PB_lBaso") %>% pull(fill)
    PB_ortho <- df %>% filter(cell_type == "PB_ortho") %>% pull(fill)
    PB_poly <- df %>% filter(cell_type == "PB_poly") %>% pull(fill)
    PB_pro <- df %>% filter(cell_type == "PB_pro") %>% pull(fill)
    

    my_shape@paths[4]$path@rgb <- CB_CD34
    my_shape@paths[1]$path@rgb <- CB_CFUE
    my_shape@paths[4]$path@rgb <- CB_CD34
    my_shape@paths[7]$path@rgb <- CB_BFUE
    my_shape@paths[10]$path@rgb <- CB_pro
    my_shape@paths[34]$path@rgb <- CB_eBaso
    my_shape@paths[42]$path@rgb <- CB_lBaso
    my_shape@paths[50]$path@rgb <- CB_poly
    my_shape@paths[56]$path@rgb <- CB_ortho
    
    my_shape@paths[72]$path@rgb <- PB_CFUE
    my_shape@paths[75]$path@rgb <- PB_CD34
    my_shape@paths[78]$path@rgb <- PB_BFUE
    my_shape@paths[81]$path@rgb <- PB_pro
    my_shape@paths[105]$path@rgb <- PB_eBaso
    my_shape@paths[113]$path@rgb <- PB_lBaso
    my_shape@paths[121]$path@rgb <- PB_poly
    my_shape@paths[127]$path@rgb <- PB_ortho  
    
    return(my_shape)
  }
}

shinyApp(ui = ui, server = server)



