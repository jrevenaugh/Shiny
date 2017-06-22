#
# Shorty - Design very short digital filters.
#

library(shiny)
library(signal)
library(tidyverse)
library(cowplot)

theme_set( theme_bw() )
options( warn = -1 )

# Define UI for application that draws a histogram
ui <- shinyUI( fluidPage(
   tags$style(type = "text/css", 
              "html, body {width:95%;height:100%;margin-left:auto; margin-right:auto;}"),
  
   # Application title
   fluidRow( 
     h4( "Shorty: Short Digital Filter Design" ),
     plotOutput("FPlot")
   ),
   
   fluidRow( 
     column( 2, 
        numericInput( "c0", label = h5("C(0)"), value = 1 ) ),
     column( 2,
        numericInput( "c1", label = h5("C(1)"), value = 0 ) ),
     column( 2, 
        numericInput( "c2", label = h5("C(2)"), value = 0 ) ),
     column( 2,
        numericInput( "c3", label = h5("C(3)"), value = 0 ) ),
     column( 2,
        numericInput( "c4", label = h5("C(4)"), value = 0 ) ),
     column( 2, 
        numericInput( "c5", label = h5("C(5)"), value = 0 ) )
   ),
   fluidRow(
     column( 2,
        radioButtons( "type", "Filter Type:",
                        c( "Zero Phase" = "zerop",
                           "Causal" = "causal" ) ) ),
     column( 2, 
        checkboxInput( "norm", "Normalized?", value = FALSE ) )
   )
))

# Server R code
server <- shinyServer(function(input, output) {
   
  output$FPlot <- renderPlot( {

    # fill in coefficients, normalize if necessary
    c <- rep( 0, 6 )
    c[1] <- input$c0
    c[2] <- input$c1
    c[3] <- input$c2
    c[4] <- input$c3
    c[5] <- input$c4
    c[6] <- input$c5
    
    pow <- t( c ) %*% c
    if ( input$type %in% c( "zerop" ) ) pow <- pow + t( c[-1] ) %*% c[-1]
    if ( input$norm ) c <- c / sqrt( t( c ) %*% c )
    
    # perform z transform with nyquist frequency of 0.5, sampled evenly
    # nF times
    
    nF <- 200
    f <- seq( 0, 0.5, length.out = 200 )
    zop <- exp( 2 * pi * f * ( 0 + 1i ) )
    z <- matrix( rep( ( 0 + 0i ), nF * 6 ), ncol = 6 )
    z[,1] <- 1
    z[,2] <- zop
    for ( i in 3:6 ) {
      k <- i - 1
      z[,i] <- zop^k
    }
    if ( input$type %in% c( "zerop" ) ) z[,2:6] <- z[,2:6] + Conj( z[,2:6] )
    ztrans <- as.vector( z %*% c )
    amp <- Mod( ztrans )
    phs <- unwrap( Arg( ztrans ) )
    if ( identical( input$type, "zerop" ) ) phs <- rep( 0, length( amp ) )
      
    
    # Create three separate ggplot graphs.  Start with the impulse response.
    # This is just 6 (or 11) dots, which is boring.  So I will resample them
    # and plot the dots on a densely sampled realization.
    
    t <- seq( -50, 50, 0.1 )
    ts <- seq( -50, 50, 1 )
    x <- rep( 0, length( ts ) )
    l1 <- which( ts == 0 )
    x[l1:(l1+5)] <- c
    tlim <- c( 0, 10 )
    if ( identical( input$type, "zerop" ) ) {
      x[(l1-1):(l1-5)] <- c[-1]
      tlim <- c( -10, 10 )
    }
    r <- signal::resample( x, 1, 0.1, d = 10 )
    X <- data_frame( time = t, r = r[1:length( t )] )
    P <- data_frame( time = ts, p = x )
    
    # build up the plot
    glyph_imp <- ggplot( X, aes( x = time ) ) +
      geom_line( aes( y = r ), color = "black", size = 1 ) +
      labs( x = "Time", y = "Impulse Response" ) +
      geom_point( data = P, aes( x = time, y = p ), col = "black", size = 3 ) +
      xlim( tlim[1], tlim[2] ) + 
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ) ) 

    glyph_amp <- data_frame( f = f, amp = amp ) %>%
      ggplot( aes( x = f, y = amp ) ) +
      geom_path( color = "red", size = 1.25 ) +
      labs( x = "Frequency", y = "Amplitude" ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ),
             legend.position = "none" ) 
    
    glyph_phs <- data_frame( f = f, phs = phs ) %>%
      ggplot( aes( x = f, y = phs ) ) +
      geom_path( color = "darkorange", size = 1.25 ) +
      labs( x = "Frequency", y = "Phase" ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ),
             legend.position = "none" )
    
    
    # Use cowplot to array these
    ggdraw( ) +
      draw_plot( glyph_imp, 0, 0.4, 1, 0.6 ) +
      draw_plot( glyph_amp, 0, 0, 0.5, 0.4 ) +
      draw_plot( glyph_phs, 0.5, 0, 0.5, 0.4 )
    
  })
})

# Run the application 
shinyApp( ui = ui, server = server )

