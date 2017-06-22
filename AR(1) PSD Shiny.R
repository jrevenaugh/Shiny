#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Changes: use ggplot2 for diagrams
# 
library(shiny)
library(ggplot2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
      tags$style(type = "text/css", 
                 "html, body {width:95%;height:100%;margin-left:auto; margin-right:auto;}"),
  
   # Application title
   titlePanel("PSD of AR(1) Signal"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel( width = 4,
         h4( "Process/Estimator Specification" ),
         sliderInput("alpha",
                     "Alpha:",
                     min = -0.99,
                     max = 0.99,
                     value = 0.1),

         sliderInput("sn",
                     "Length of Time Series",
                     min = 64,
                     max = 4096,
                     value = 1024),
         
         sliderInput("segments",
                     "Number of WOSA Segments:",
                     min = 2,
                     max = 60,
                     value = 4)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("specPlot")
      )
      
   )

))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
   
   output$specPlot <- renderPlot({
     require(sapa, quietly = T, warn.conflicts = F)
     require(ifultools, quietly = T, warn.conflicts = F)
     require(splus2R, quietly = T, warn.conflicts = F)
     
     arps <- function( alpha, freq )
     {
       I <- 1.0 / ( 1 + alpha^2 - 2 * alpha * cos( 2 * pi * freq ) )
       return( I )
     }   
     n <- 2^( floor( log2( input$sn) ) ) 
     x <- arima.sim( list( ar = input$alpha ), n )
     ti <- sprintf( "AR(1) alpha = %g", input$alpha )
     T <- data_frame( t = seq( 1, n ), ar = as.numeric( x ) )
     gT <- ggplot( T, aes( x = t, y = ar ) ) + geom_line( colour = "red3" ) +
       labs( x = "Time", y = "AR(1)", title = ti )
     xps <- as.numeric( SDF( x, method = "wosa", blocksize = floor( input$sn/input$segments ) ) )
     freq <- seq( 0, 0.5, length.out = n+1 )
     F <- data_frame( f = freq, ps = xps, tps = arps( input$alpha, freq ) )
     ti <- sprintf( "Series Length = %g, # Segments = %g", n, input$segments )
     gF <- ggplot( F, aes( x = f ) ) + geom_line( aes( y = xps ), colour = "red3", size = 1.2 ) +
           geom_line( aes( y = tps ), color = "darkblue", size = 1.2 ) +
           labs( x = "Frequency", y = "PSD", title = ti )
     multiplot( gT, gF, cols = 1 )
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

