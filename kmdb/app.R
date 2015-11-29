# Kurzchalia Microarray Database
# Copyright (C) 2015  Cihan Erkut

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Loading librariers
library(shiny)
library(ggplot2)
library(readr)
library(dplyr)

# Load data
expr1 <- read_tsv('expr_desiccation.txt')
expr2 <- read_tsv('results_desiccation.txt')
expr3 <- read_tsv('expr_hypometabolism.txt')
expr4 <- read_tsv('results_hypometabolism.txt')

expr1$Gene <- as.vector(expr1$Gene)
expr3$Gene <- as.vector(expr3$Gene)
expr3$Stage <- factor(expr3$Stage, c('L3', 'dauer'))
expr3$Strain <-
  factor(expr3$Strain, c('N2', 'daf-2', 'daf-16', 'daf-2;daf-12'))
expr4$Strain <-
  factor(expr4$Strain, c('N2', 'daf-2', 'daf-16', 'daf-2;daf-12'))

geneClass1 <-
  sort(unique(gsub(
    '-[[:alnum:]]+', '', grep('-[[:alnum:]]+', expr1$Gene, value = T)
  )))
geneClass2 <-
  sort(unique(gsub(
    '-[[:alnum:]]+', '', grep('-[[:alnum:]]+', expr3$Gene, value = T)
  )))

# Common definitions
wb.url <- 'http://www.wormbase.org/species/c_elegans/gene/'

common.theme <- theme(
  legend.position  = "bottom",
  panel.background = element_rect(fill = grey(0.95)),
  strip.background = element_rect(fill = grey(0.90), colour = NA),
  panel.border     = element_rect(fill = NA,         colour = 'black'),
  strip.text       = element_text(face = 'bold'),
  axis.text.x      = element_text(face = 'italic'),
  axis.title.x     = element_blank()
)

msgEntry          <-
  'Type in valid gene names separated by comma (e.g. tps-1, F08H9.4)'
msgPlot           <-
  'This plot will update as you enter valid gene names!'
msgTable          <-
  'This table will update as you enter valid gene names!'
msgConditionPlot  <- 'This plot will update as you select strains!'
msgConditionTable <- 'This table will update as you select strains!'

# Utility functions
trim <- function(x) {
  gsub('(^[[:space:]]+|[[:space:]]+$)', '', x)
}

expr.subset <- function(expr.df, geneList) {
  expr.df[expr.df$Gene %in% trim(unlist(strsplit(geneList, ','))),]
}

link.wb <- function(x) {
  HTML(as.character(a(
    href = paste(wb.url, x, sep = ''), x, target = '_blank'
  )))
}

buildExpr <-
  function(expr.df, input.query, input.text, input.geneClass, input.conditions = NULL) {
    if (input.query == 'Genes') {
      expr.new <- expr.subset(expr.df, input.text)
    } else {
      classPattern <- paste0('^', input.geneClass, '-')
      expr.new <-
        expr.subset(expr.df, grep(classPattern, expr.df$Gene, value = T))
    }
    if (is.vector(input.conditions)) {
      expr.new <- expr.new[expr.new$Strain %in% input.conditions,]
    }
    return(expr.new)
  }

# Server function
server <- function(input, output) {
  # Update data
  data1 <-
    reactive(buildExpr(expr1, input$query1, input$text1, input$geneClass01))
  data2 <-
    reactive(buildExpr(
      expr3, input$query2, input$text2, input$geneClass02, input$conditions
    ))
  data3 <-
    reactive(buildExpr(expr2, input$query1, input$text1, input$geneClass01))
  data4 <-
    reactive(buildExpr(
      expr4, input$query2, input$text2, input$geneClass02, input$conditions
    ))
  
  # Render box-plots with normalized expression levels
  output$BoxPlot1 <- renderPlot({
    plot.data <- data1()
    validate(need(nrow(plot.data) > 0, msgPlot))
    plot.base <-
      ggplot(plot.data, aes(x = Gene, fill = Treatment)) +
      ggtitle('Normalized Expression Levels') +
      scale_fill_grey() +
      common.theme
    if (input$log1) {
      plot.base +
        geom_boxplot(aes(y = Expression)) +
        ylab(expression(log[2] * ' Expression Level (AU)'))
    } else {
      plot.base +
        geom_boxplot(aes(y = 2 ^ Expression)) +
        ylab('Expression Level (AU)')
    }
  })
  
  output$BoxPlot2 <- renderPlot({
    validate(need(length(input$conditions) > 0, msgConditionPlot))
    plot.data <- data2()
    validate(need(nrow(plot.data) > 0, msgPlot))
    plot.base <- ggplot(plot.data, aes(x = Gene, fill = Stage)) +
      facet_wrap(~Strain) +
      ggtitle('Normalized Expression Levels') +
      common.theme
    if (input$log2) {
      plot.base +
        geom_boxplot(aes(y = Expression)) +
        ylab(expression(log[2] * ' Expression Level (AU)'))
    } else {
      plot.base +
        geom_boxplot(aes(y = 2 ^ Expression)) +
        ylab('Expression Level (AU)')
    }
  })
  
  # Render bar plots with differential expression levels
  output$BarPlot1 <- renderPlot({
    plot.data <- data3()
    validate(need(nrow(plot.data) > 0, msgPlot))
    plot.base <- ggplot(plot.data, aes(x = Gene)) +
      ggtitle('Differential Expression Levels') +
      common.theme
    if (input$log1) {
      plot.base +
        geom_bar(aes(y = logFC), stat = 'identity') +
        geom_errorbar(aes(ymax = logFC.U, ymin = logFC.L), width = I(0.5)) +
        ylab(expression(log[2] * ' Fold Change'))
    } else {
      plot.base +
        geom_bar(aes(y = FC), stat = 'identity') +
        geom_errorbar(aes(ymax = FC.U, ymin = FC.L), width = I(0.5)) +
        ylab('Fold Change')
    }
  })
  
  output$BarPlot2 <- renderPlot({
    validate(need(length(input$conditions) > 0, msgConditionPlot))
    plot.data <- data4()
    validate(need(nrow(plot.data) > 0, msgPlot))
    plot.base <- ggplot(plot.data, aes(x = Gene)) +
      facet_wrap(~Strain) +
      ggtitle('Differential Expression Levels') +
      common.theme
    if (input$log2) {
      plot.base +
        geom_bar(aes(y = logFC), stat = 'identity') +
        geom_errorbar(aes(ymax = logFC.U, ymin = logFC.L), width = I(0.5)) +
        ylab(expression(log[2] * ' Fold Change'))
    } else {
      plot.base +
        geom_bar(aes(y = FC), stat = 'identity') +
        geom_errorbar(aes(ymax = FC.U, ymin = FC.L), width = I(0.5)) +
        ylab('Fold Change')
    }
  })
  
  # Render tables with differential expression levels
  output$Table1 <- renderDataTable({
    table.data <- data3()
    validate(need(nrow(table.data) > 0, msgTable))
    table.data$Gene <- sapply(table.data$Gene, link.wb)
    format(table.data, digits = 3)
  }, escape = F)
  
  output$Table2 <- renderDataTable({
    table.data <- data4()
    validate(need(nrow(table.data) > 0, msgTable))
    table.data$Gene <- sapply(table.data$Gene, link.wb)
    format(table.data, digits = 3)
  }, escape = F)
}

# GUI function
ui <- shinyUI(
  navbarPage(
    title       = 'Kurzchalia Microarray Database 0.2',
    collapsible = T,
    windowTitle = 'Kurzchalia Microarray Database',
    tabPanel('Desiccation',
             fluidPage(
               titlePanel("Desiccation"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons('query1', 'Search for', c('Genes', 'A gene class'), inline = T),
                   conditionalPanel('input.query1 == "Genes"',
                                    withTags(
                                      div(
                                        class = 'geneEntry',
                                        style(type = "text/css", "textarea {width:100%}"),
                                        label('Gene names:'),
                                        textarea(
                                          name = "text1", rows = 3, wrap = 'soft', placeholder = msgEntry
                                        )
                                      )
                                    )),
                   conditionalPanel(
                     'input.query1 == "A gene class"',
                     selectInput(
                       'geneClass01', 'Gene classes:', choices = geneClass1, selectize = T
                     )
                   ),
                   checkboxInput("log1", label = HTML(
                     paste0('Show y-axis in log', tags$sub(2), ' scale')
                   ))
                 ),
                 mainPanel(tabsetPanel(
                   tabPanel(
                     "Expression Levels", plotOutput('BoxPlot1'), icon = icon('area-chart')
                   ),
                   tabPanel(
                     "Differential Expression", plotOutput('BarPlot1'), icon = icon('bar-chart')
                   ),
                   tabPanel(
                     "Numerical Results", dataTableOutput('Table1'), icon = icon('table')
                   )
                 ))
               )
             )),
    tabPanel('Hypometabolism',
             fluidPage(
               titlePanel("Hypometabolism"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons('query2', 'Search for', c('Genes', 'A gene class'), inline = T),
                   conditionalPanel('input.query2 == "Genes"',
                                    withTags(
                                      div(
                                        class = 'geneEntry',
                                        style(type = "text/css", "textarea {width:100%}"),
                                        label('Gene names:'),
                                        textarea(
                                          name = "text2", rows = 3, wrap = 'soft', placeholder = msgEntry
                                        )
                                      )
                                    )),
                   conditionalPanel(
                     'input.query2 == "A gene class"',
                     selectInput(
                       'geneClass02', 'Gene classes:', choices = geneClass2, selectize = T
                     )
                   ),
                   checkboxInput(inputId = "log2",
                                 label = HTML(
                                   paste0('Show y-axis in log', tags$sub(2), ' scale')
                                 )),
                   checkboxGroupInput(
                     inputId = 'conditions',
                     label = 'Strains',
                     choices = levels(expr3$Strain),
                     selected = levels(expr3$Strain)
                   )
                 ),
                 mainPanel(tabsetPanel(
                   tabPanel(
                     "Expression Levels", plotOutput('BoxPlot2'), icon = icon('area-chart')
                   ),
                   tabPanel(
                     "Differential Expression", plotOutput('BarPlot2'), icon = icon('bar-chart')
                   ),
                   tabPanel(
                     "Numerical Results", dataTableOutput('Table2'), icon = icon('table')
                   )
                 ))
               )
             )),
    tabPanel('Help',
             includeMarkdown('Help.md'))
  )
)

# Run the app
shinyApp(ui = ui, server = server)
