# Loading librariers
library(shiny)
library(ggplot2)
library(xtable)

# Common definitions
file.expr1 <- 'expr_desiccation.txt'
file.expr2 <- 'results_desiccation.txt'
file.expr3 <- 'expr_hypometabolism.txt'
file.expr4 <- 'results_hypometabolism.txt'

strains <- c('N2', 'daf-2', 'daf-16', 'daf-2;daf-12')
stages <- c('L3', 'dauer')

wb.url <- 'http://www.wormbase.org/species/c_elegans/gene/'

common.theme <- theme(legend.position = "bottom", 
                      panel.background = element_rect(fill = grey(0.95)), 
                      strip.background = element_rect(fill = grey(0.90), colour = NA), 
                      strip.text = element_text(face = 'bold'),
                      panel.border = element_rect(colour = 'black', fill = NA), 
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(face = 'italic')) 

msgPlot <- 'This plot will update as you enter valid gene names!'
msgTable <- 'This table will update as you enter valid gene names!'
msgEntry <- 'Type in valid gene names separated by comma (e.g. tps-1, F08H9.4)'

# Utility functions
trim <- function(x) {
  pattern <- '(^[[:space:]]+|[[:space:]]+$)'
  trimmed <- gsub(pattern, '', x)}

expr.subset <- function(expr.df, genelist) {
  genes.selected <- trim(unlist(strsplit(genelist, ',')))
  expr.selected <- expr.df[expr.df$Gene %in% genes.selected, ]}

link.wb <- function(x) {
  HTML(as.character(a(href = paste(wb.url, x, sep = ''), x, target = '_blank')))}

# Server function
server <- function(input, output) {
  # Load data while displaying progress
  withProgress(message = 'Loading data... ', value = 0, {
                 
    expr1 <- read.delim(file.expr1, sep = '\t', header = T)
    incProgress(0.25)

    expr2 <- read.delim(file.expr2, sep = '\t', header = T)
    incProgress(0.2)
    
    expr3 <- read.delim(file.expr3, sep = '\t', header = T, as.is = T)
    incProgress(0.25)
    
    expr4 <- read.delim(file.expr4, sep = '\t', header = T, as.is = T)
    incProgress(0.2, message = 'Arranging data...')
    
    expr1$Gene <- as.vector(expr1$Gene)
    expr3$Gene <- as.vector(expr3$Gene)
    expr3$Stage <- factor(expr3$Stage, stages)
    expr3$Strain <- factor(expr3$Strain, strains)
    expr4$Strain <- factor(expr4$Strain, strains)
    incProgress(0.1)})
  
  # Render box-plots with normalized expression levels
  output$BoxPlot1 <- renderPlot({
    expr1.subset <- expr.subset(expr1, input$text1)
    validate(
      need(
        nrow(expr1.subset) > 0, msgPlot))
    
    plot.base <- ggplot(expr1.subset, aes(x = Gene, fill = Treatment)) + 
                 ggtitle('Normalized Expression Levels') + scale_fill_grey() + common.theme
    if (input$log1) {
      plot.base + geom_boxplot(aes(y = Expression)) + ylab(expression(log[2] * ' Expression Level (AU)'))}
    else {
      plot.base + geom_boxplot(aes(y = 2 ^ Expression)) + ylab('Expression Level (AU)')}})

  output$BoxPlot2 <- renderPlot({
    expr3.subset <- expr.subset(expr3, input$text2)
    validate(
      need(
        nrow(expr3.subset) > 0, msgPlot))
    
    plot.base <- ggplot(expr3.subset, aes(x = Gene, fill = Stage)) + facet_wrap(~ Strain) + 
                 ggtitle('Normalized Expression Levels') + common.theme
    if (input$log2) {
      plot.base + geom_boxplot(aes(y = Expression)) + ylab(expression(log[2] * ' Expression Level (AU)'))}
    else {
      plot.base + geom_boxplot(aes(y = 2 ^ Expression)) + ylab('Expression Level (AU)')}})
  
  # Render bar plots with differential expression levels
  output$BarPlot1 <- renderPlot({
    expr2.subset <- expr.subset(expr2, input$text1)
    validate(
      need(
        nrow(expr2.subset) > 0, msgPlot))
    
    plot.base <- ggplot(expr2.subset, aes(x = Gene)) +
                 ggtitle('Differential Expression Levels') + common.theme
    if (input$log1) {
      plot.base + geom_bar(aes(y = logFC), stat = 'identity') +
                  geom_errorbar(aes(ymax = logFC.U, ymin = logFC.L), width = I(0.5)) +
                  ylab(expression(log[2] * ' Fold Change'))}
    else {
      plot.base + geom_bar(aes(y = FC), stat = 'identity') +
                  geom_errorbar(aes(ymax = FC.U, ymin = FC.L), width = I(0.5)) +
                  ylab('Fold Change')}})

  output$BarPlot2 <- renderPlot({
    expr4.subset <- expr.subset(expr4, input$text2)
    validate(
      need(
        nrow(expr4.subset) > 0, msgPlot))
    
    plot.base <- ggplot(expr4.subset, aes(x = Gene)) + facet_wrap(~ Strain) +
                 ggtitle('Differential Expression Levels') + common.theme
    if (input$log2) {
      plot.base + geom_bar(aes(y = logFC), stat = 'identity') +
                  geom_errorbar(aes(ymax = logFC.U, ymin = logFC.L), width = I(0.5)) +
                  ylab(expression(log[2] * ' Fold Change'))}
    else {
      plot.base + geom_bar(aes(y = FC), stat = 'identity') +
                  geom_errorbar(aes(ymax = FC.U, ymin = FC.L), width = I(0.5)) +
                  ylab('Fold Change')}})
  
  # Render tables with differential expression levels
  output$Table1 <- renderTable({
    expr2.subset <- expr.subset(expr2, input$text1)
    validate(
      need(
        nrow(expr2.subset) > 0, msgTable)) 
    expr2.subset$Gene <- sapply(expr2.subset$Gene, link.wb)
    xtable(expr2.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Table2.1 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'N2', -2], input$text2)
    validate(
      need(
        nrow(expr4.subset) > 0, msgTable))
    expr4.subset$Gene <- sapply(expr4.subset$Gene, link.wb)
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Table2.2 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'daf-2', -2], input$text2)
    validate(
      need(
        nrow(expr4.subset) > 0, msgTable))
    expr4.subset$Gene <- sapply(expr4.subset$Gene, link.wb)
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)

  output$Table2.3 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'daf-2;daf-12', -2], input$text2)
    validate(
      need(
        nrow(expr4.subset) > 0, msgTable))
    expr4.subset$Gene <- sapply(expr4.subset$Gene, link.wb)
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Table2.4 <- renderTable({
    expr4.subset <- expr.subset(expr4[expr4$Strain == 'daf-16', -2], input$text2)
    validate(
      need(
        nrow(expr4.subset) > 0, msgTable))
    expr4.subset$Gene <- sapply(expr4.subset$Gene, link.wb)
    xtable(expr4.subset)}, digits = 3, sanitize.text.function = function(x) x)
  
  output$Title2.1 <- renderText('N2')
  output$Title2.2 <- renderText('daf-2')
  output$Title2.3 <- renderText('daf-2;daf-12')
  output$Title2.4 <- renderText('daf-16')}

# GUI function
ui <- shinyUI(navbarPage(title = 'Kurzchalia Microarray Database 0.1',
                         collapsible = T, 
                         windowTitle = 'Kurzchalia Microarray Database',
  tabPanel('Desiccation',   
    fluidPage(
      titlePanel("Desiccation"),
        sidebarLayout(  
          sidebarPanel(
            tags$style(type="text/css", "textarea {width:100%}"),
            tags$textarea(name = "text1", rows = 3, wrap = 'soft', placeholder = msgEntry),
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
            tags$style(type="text/css", "textarea {width:100%}"),
            tags$textarea(name = "text2", rows = 3, wrap = 'soft', placeholder = msgEntry),
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

  tabPanel('Help', includeMarkdown('Help.md'))))

# Run the app
shinyApp(ui = ui, server = server)
