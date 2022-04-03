## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(dplyr)


# Define UI for application that draws a histogram
ui <- fluidPage(
    fileInput("file", "Please select a DESeq data CSV file" , accept = ".csv"), #, placeholder = "deseq_res.csv"),
    radioButtons('button1', 'Please select a variable to plot on your x-axis', 
                  c('baseMean',
                    'log2FoldChange',
                    'lfcSE',
                    'stat',
                    'pvalue',
                    'padj')),
    radioButtons('button2', 'Please select a variable to plot on your y-axis', 
                   c('baseMean',
                   'log2FoldChange',
                   'lfcSE',
                   'stat',
                   'pvalue',
                   'padj')),
    colourInput('col1', 'Choose a color for points below your adjusted p-value threshold', 'red'), #for plotting
    colourInput('col2', 'Choose color for points above your adjusted p-value threshold', 'blue'),#for what meets p-value threshold
    sliderInput('slider', 'Select the magnitude of adjusted p-value you would like to color on the graph', -300, 0, -10),
    submitButton('Submit', icon('sync'), '400px'),
    mainPanel(
      tabsetPanel(
        tabPanel('Differential Expression Data Table',
          
          tableOutput('table')
        ),
        tabPanel('Differential Expression Data Plot',
          plotOutput('volcano')
        )
      )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
        req(input$file)
        return(read.csv(input$file$datapath, header=TRUE, sep=",")%>%
                 rename(Gene = X))  
    })
    
    
    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
    #' Write a normal volcano plot using geom_point, and integrate all the above 
    #' values into it as shown in the example app. The testing script will treat 
    #' this as a normal function.
    #' 
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
          dataf$cutoff[dataf$pvalue < 10^slider] <- 'BELOW PADJ MAGNITUDE'
          dataf$cutoff[dataf$pvalue > 10^slider] <- 'ABOVE PADJ MAGNITUDE'
          plt <- dataf %>%
            ggplot() +
            geom_point(aes(x=!!sym(x_name), y=-log10(!!sym(y_name)), color = cutoff)) + #volcano points
            xlab(x_name) +
            ylab(paste('-log10',y_name, sep = ' ')) +
            scale_color_manual(values=c(color1, color2)) +
            ggtitle('Volcano Plot of DeSeq Differential Expression Results') +
            theme_minimal()
            
            return(plt)
        }
  
    
    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than 
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will 
    #' evaluate it normally. Not only does this function filter the data frame to 
    #' rows that are above the slider magnitude, it should also change the format 
    #' of the p-value columns to display more digits. This is so that it looks 
    #' better when displayed on the web page. I would suggest the function 
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
      filt_dataf <- filter(dataf, padj< 10^slider) %>%  #dataf[!(dataf$padj > 10^slider),]
        mutate(pvalue = formatC(.$pvalue, digits = 4, format = 'e'),
               padj = formatC(.$padj, digits = 4, format = 'e'))
      return(filt_dataf)
    }
    
    #' These outputs aren't really functions, so they don't get a full skeleton, 
    #' but use the renderPlot() and renderTabel() functions to return() a plot 
    #' or table object, and those will be displayed in your application.
    output$volcano <- renderPlot(volcano_plot(load_data(), input$button1, input$button2, input$slider, input$col1, input$col2)) # replace this NULL
    
    # Same here, just return the table as you want to see it in the web page
    output$table <- renderTable(draw_table(load_data(), input$slider)) # replace this NULL
    
}

# Run the application
shinyApp(ui = ui, server = server)
