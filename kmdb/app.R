library(shiny)
library(ggplot2)
library(xtable)

trim <- function(x) {
  pattern <- '(^[[:space:]]+|[[:space:]]+$)'
  trimmed <- gsub(pattern, '', x)
  return(trimmed)
}

expr.subset <- function(expr.df, genelist) {
  genes.selected <- trim(unlist(strsplit(genelist, ',')))
  expr.selected <- expr.df[expr.df$Gene %in% genes.selected,]
  return(expr.selected)
}

server <- function(input, output, session) {
  # Load data while displaying progress
  withProgress(message = 'Loading data... ', value = 0, {
                 
    expr1 <- read.delim('expr_desiccation.txt', sep = '\t', header = T)
    incProgress(0.25)

    expr2 <- read.delim('results_desiccation.txt', sep = '\t', header = T)
    incProgress(0.2)
    
    expr3 <- read.delim('expr_hypometabolism.txt' , sep = '\t', header = T, as.is = T)
    incProgress(0.25)
    
    expr4 <- read.delim('results_hypometabolism.txt', sep = '\t', header = T, as.is = T)
    incProgress(0.2, message = 'Arranging data...')
    
    expr1$Gene <- as.vector(expr1$Gene)
    expr3$Strain <- factor(expr3$Strain, c('N2', 'daf-2', 'daf-16', 'daf-2;daf-12'))
    expr3$Stage <- factor(expr3$Stage, c('L3', 'dauer'))
    expr4$Strain <- factor(expr4$Strain, c('N2', 'daf-2', 'daf-16', 'daf-2;daf-12'))
    incProgress(0.1)
  })
  
  # Render box-plot with normalized expression levels
  output$BoxPlot1 <- renderPlot({
    expr1.subset <- expr.subset(expr1, input$text1)
    validate(need(nrow(expr1.subset) > 0, 'This plot will update as you enter valid gene names!'))
    
    a <- ggplot(expr1.subset, aes(x = Gene,
                                  fill = Treatment)) +
         xlab('') +
         ggtitle('Normalized Expression Levels') +
         scale_fill_grey() +
         theme(legend.position = "bottom",
               panel.background = element_rect(fill = grey(0.95)), 
               panel.border = element_rect(colour = 'black', 
                                           fill = NA))    
    if(input$log1) {
      a <- a + geom_boxplot(aes(y = Expression)) +
           ylab(expression(log[2] * ' Expression Level (AU)'))}
    else {
      a <- a + geom_boxplot(aes(y = 2 ^ Expression)) +
           ylab('Expression Level (AU)')}
    a})

  output$BoxPlot2 <- renderPlot({
    expr3.subset <- expr.subset(expr3, input$text2)
    validate(need(nrow(expr3.subset) > 0, 'This plot will update as you enter valid gene names!'))
    
    a <- ggplot(expr3.subset, aes(x = Gene, 
                                  fill = Stage)) +
         facet_wrap(~ Strain) +
         xlab('') +
         ggtitle('Normalized Expression Levels') +
         theme(legend.position = "bottom",
               panel.background = element_rect(fill = grey(0.95)),
               strip.background = element_rect(fill = NA, 
                                               colour = NA),
               panel.border = element_rect(colour = 'black', 
                                           fill = NA))
    if(input$log2) {
      a <- a + geom_boxplot(aes(y = Expression)) +
           ylab(expression(log[2] * ' Expression Level (AU)'))}
    else {
      a <- a + geom_boxplot(aes(y = 2 ^ Expression)) +
           ylab('Expression Level (AU)')}
    a})
  
  # Render bar plot with differential expression levels
  output$BarPlot1 <- renderPlot({
    expr2.subset <- expr.subset(expr2, input$text1)
    validate(need(nrow(expr2.subset) > 0, 'This plot will update as you enter valid gene names!'))
    
    a <- ggplot(expr2.subset, aes(x = Gene)) +
         xlab('') +
         ggtitle('Differential Expression Levels') +
         scale_fill_grey() +
         theme(legend.position = "bottom",
               panel.background = element_rect(fill = grey(0.95)),
               panel.border = element_rect(colour = 'black', 
                                           fill = NA))
    if(input$log1) {
      a <- a + geom_bar(aes(y = logFC), 
                        stat = 'identity') +
           geom_errorbar(aes(ymax = logFC.U, 
                             ymin = logFC.L), 
                         width = I(0.5)) +
           ylab(expression(log[2] * ' Fold Change'))}
    else {
      a <- a + geom_bar(aes(y = FC), 
                        stat = 'identity') +
        geom_errorbar(aes(ymax = FC.U, 
                          ymin = FC.L), 
                      width = I(0.5)) +
        ylab('Fold Change')}    
    a})

  output$BarPlot2 <- renderPlot({
    expr4.subset <- expr.subset(expr4, input$text2)
    validate(need(nrow(expr4.subset) > 0, 'This plot will update as you enter valid gene names!'))
    
    a <- ggplot(expr4.subset, aes(x = Gene)) +
         facet_wrap(~ Strain) +
         xlab('') +
         ggtitle('Differential Expression Levels') +
         scale_fill_grey() +
         theme(legend.position = "bottom", 
               panel.background = element_rect(fill = grey(0.95)), 
               strip.background = element_rect(fill = NA, 
                                               colour = NA), 
               panel.border = element_rect(colour = 'black', 
                                           fill = NA))    
    if(input$log2) {
      a <- a + geom_bar(aes(y = logFC), 
                        stat = 'identity') +
           geom_errorbar(aes(ymax = logFC.U, 
                             ymin = logFC.L), 
                         width = I(0.5)) +
        ylab(expression(log[2] * ' Fold Change'))}
    else {
      a <- a + geom_bar(aes(y = FC), 
                        stat = 'identity') +
        geom_errorbar(aes(ymax = FC.U, 
                          ymin = FC.L), 
                      width = I(0.5)) +
        ylab('Fold Change')}    
    a})
  
  # Render table with differential expression levels
  output$Table1 <- renderTable({
    expr2.subset <- expr.subset(expr2, input$text1)
    validate(need(nrow(expr2.subset) > 0, 'This table will update as you enter valid gene names!'))
    
    expr2.subset$Gene <- sapply(expr2.subset$Gene, 
                                function(x) HTML(as.character(a(href = paste("http://www.wormbase.org/species/c_elegans/gene/", x, sep = ''), x, target = '_blank'))))
    xtable(expr2.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Table2.1 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'N2',-2], input$text2)
    validate(need(nrow(expr4.subset) > 0, 'This table will update as you enter valid gene names!'))
    
    expr4.subset$Gene <- sapply(expr4.subset$Gene, 
                                function(x) HTML(as.character(a(href = paste("http://www.wormbase.org/species/c_elegans/gene/", x, sep = ''), x, target = '_blank'))))
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Table2.2 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'daf-2',-2], input$text2)
    validate(need(nrow(expr4.subset) > 0, 'This table will update as you enter valid gene names!'))
    
    expr4.subset$Gene <- sapply(expr4.subset$Gene, 
                                function(x) HTML(as.character(a(href = paste("http://www.wormbase.org/species/c_elegans/gene/", x, sep = ''), x, target = '_blank'))))
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)

  output$Table2.3 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'daf-2;daf-12',-2], input$text2)
    validate(need(nrow(expr4.subset) > 0, 'This table will update as you enter valid gene names!'))
    
    expr4.subset$Gene <- sapply(expr4.subset$Gene, 
                                function(x) HTML(as.character(a(href = paste("http://www.wormbase.org/species/c_elegans/gene/", x, sep = ''), x, target = '_blank'))))
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Table2.4 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'daf-16',-2], input$text2)
    validate(need(nrow(expr4.subset) > 0, 'This table will update as you enter valid gene names!'))
    
    expr4.subset$Gene <- sapply(expr4.subset$Gene, 
                                function(x) HTML(as.character(a(href = paste("http://www.wormbase.org/species/c_elegans/gene/", x, sep = ''), x, target = '_blank'))))
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  
  output$Title2.1 <- renderText('N2')
  
  output$Title2.2 <- renderText('daf-2')

  output$Title2.3 <- renderText('daf-2;daf-12')
  
  output$Title2.4 <- renderText('daf-16')
  
}

ui <- shinyUI(navbarPage('Kurzchalia Microarray Database 0.1',
                         collapsable = T, 
                         windowTitle = 'Kurzchalia Microarray Database',
  tabPanel('Desiccation',   
    fluidPage(
      titlePanel("Desiccation"),
        sidebarLayout(  
          sidebarPanel(
            textInput("text1", label = "Type in valid gene names separated by comma (e.g. tps-1, F08H9.4)"),
            checkboxInput("log1", label = HTML(paste('Show y-axis in log', tags$sub(2), ' scale', sep = '')))),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Expression Levels", plotOutput('BoxPlot1')),
          tabPanel("Differential Expression", plotOutput('BarPlot1')),
          tabPanel("Numerical Results", tableOutput('Table1'))))))),

  tabPanel('Hypometabolism',   
    fluidPage(
      titlePanel("Hypometabolism"),
        sidebarLayout(  
          sidebarPanel(
            textInput("text2", label = "Type in valid gene names separated by comma (e.g. tps-1, F08H9.4)"),
            checkboxInput("log2", label = HTML(paste('Show y-axis in log', tags$sub(2), ' scale', sep = '')))),
               
      mainPanel(
        tabsetPanel(
          tabPanel("Expression Levels", plotOutput('BoxPlot2')),
          tabPanel("Differential Expression", plotOutput('BarPlot2')),
          tabPanel("Numerical Results",
                   h4(textOutput('Title2.1')), tableOutput('Table2.1'),
                   h4(em(textOutput('Title2.2'))), tableOutput('Table2.2'),
                   h4(em(textOutput('Title2.3'))), tableOutput('Table2.3'),
                   h4(em(textOutput('Title2.4'))), tableOutput('Table2.4'))))))),               

  tabPanel('Help', includeMarkdown('Help.md'))
))

shinyApp(ui = ui, server = server, options = list(port = 24426))
