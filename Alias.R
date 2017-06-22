# Alias - Shiny web app to illustrate signal aliasing and Shannon's theory.
# You can run the application by clicking
# the 'Run App' button above.
#


library( shiny )
library( signal )
library( tidyverse )
options( warn = -1 )
theme_set( theme_bw( ) )

# Define UI for application that draws a histogram
ui <- fluidPage(
   tags$style(type = "text/css", 
              "html, body {width:95%;height:100%;margin-left:auto; margin-right:auto;}"),
  
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         h4( "Aliasing Parameters" ),
         sliderInput( "xPeriod1",
                      "Signal Period 1",
                      min = 2,
                      max = 50,
                      value = 10,
                      step = 2 ),
         sliderInput( "xPeriod2",
                      "Signal Period 2",
                      min = 2,
                      max = 50,
                      value = 10,
                      step = 2 ),
         checkboxInput( "checkbox", label = "Include Period 2", value = TRUE ),
         
         sliderInput( "sampInt",
                      "Sampling Interval",
                      min = 1,
                      max = 25,
                      value = 10 ),
         helpText( "Resampled series in black.")
         
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput( "tsPlot" )
      )
   )
)

# Define server logic required to draw a histogram
server <- function( input, output ) {
   
   output$tsPlot <- renderPlot({
      # generate time series, sampled points and resampled series
      tx <- seq( -100, 300, by = 0.1 )
      x <- sin( 2 * pi * tx / input$xPeriod1 ) 
      if ( input$checkbox ) x <- x + 0.5 * sin( 2 * pi * tx / input$xPeriod2 )
      xmax <- max( abs( x ) )
      x <- x / xmax
      ts <- seq( -100, 300, by = input$sampInt )
      xs <- sin( 2 * pi * ts / input$xPeriod1 ) 
      if ( input$checkbox ) xs <- xs + 0.5 * sin( 2 * pi * ts / input$xPeriod2 )
      xs <- xs / xmax
      r <- resample( xs, input$sampInt, 0.1, d = 10 )

      # pack results into glyph-ready data_frames
      X <- data_frame( time = tx, x = x, r = r[1:length(x)] )
      P <- data_frame( time = ts, p = xs )

      # build up the plot
      ggplot( X, aes( x = time ) ) +
        geom_line( aes( y = x ), color = "red", size = 1 ) +
        geom_line( aes( y = r ), color = "black", size = 1 ) +
        labs( x = "Time", y = "" ) +
        geom_point( data = P, aes( x = time, y = p ), col = "black", size = 3 ) +
        xlim( 0, 200 ) + 
        theme( text = element_text( size = 20 ),
               plot.title = element_text( size = 20 ),
               plot.margin = margin( 10, 0, 0, 10, "pt" ) ) 
      
   })
}

# Run the application 
shinyApp( ui = ui, server = server )

