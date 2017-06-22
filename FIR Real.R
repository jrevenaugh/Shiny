#
# FIR Real - Specify and visualize Butterworth, Chebyshev and Elliptical
#            FIR filters.  Learning tool, not a filter design app.

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
   h4( "FIR Filter Visualizer" ),
   
   # Sidebar filter parameters 
   sidebarLayout(
      sidebarPanel(
          sliderInput("order", "Filter Order", min = 1, 
                       max = 10, value = 3, step = 1 ),
          sliderInput("freq_1", "Low Frequency", min = 0, 
                      max = 1, value = 0.1, step = 0.1 ),
          sliderInput("freq_2", "High Frequency", min = 0, 
                      max = 1, value = 0.9, step = 0.1 ),
          sliderInput("pbrip", "Passband Ripple (dB)", min = 0.1, 
                      max = 5, value = 1, step = 0.1 ),
          sliderInput("rbrip", "Stopband Ripple (dB)", min = 10, 
                      max = 60, value = 30, step = 1 ),
          radioButtons( "type", "Filter Type:",
                        c( "Low Pass" = "low",
                          "High Pass" = "high",
                          "Band Pass" = "pass",
                          "Band Reject" = "stop") )
      ),
      
      # Plot panel
      mainPanel(
         plotOutput("FPlot")
      )
   )
))

# Server R code
server <- shinyServer(function(input, output) {
   
  output$FPlot <- renderPlot( height = 700, {
    W <- input$freq_1
    if ( input$type %in% c( "pass", "stop" ) ) 
      W <- c( min( input$freq_1, input$freq_2 ), max( input$freq_1, input$freq_2 ) )
    if ( input$type %in% c( "high" ) ) 
      W <- input$freq_2
    
    # Compute filters
    butt <- butter( n = input$order, W = W, type = input$type, plane = "z" )
    cheb <- cheby1( n = input$order, Rp = input$pbrip, W = W, type = input$type, plane = "z" )
    ellp <- ellip( n = input$order, Rp = input$pbrip, Rs = input$rbrip, W = W,
                   type = input$type, plane = "z" )
    
    # Compute amplitude and phase
    fhf <- c( freqz( butt ), freqz( cheb ), freqz( ellp ) )
    freq <- fhf[[2]] / pi
    amp <- matrix( rep( 0, 3 * length( freq ) ), ncol = 3 )
    phs <- matrix( rep( 0, 3 * length( freq ) ), ncol = 3 )
    for ( i in 1:3 ) {
      k <- 2 * i - 1
      amp[,i] <- Mod( fhf[[k]] )
      phs[,i] <- unwrap( Arg( fhf[[k]] ) ) 
    }

    # Compute impulse response
    fimp <- c( impz( butt ), impz( cheb ), impz( ellp ) )
    nT <- max( length( fimp[[1]] ), length( fimp[[3]] ), length( fimp[[5]] ) )
    nT <- min( 100, nT )
    time <- seq( 0, nT - 1 )
    imp <- matrix( rep( 0, 3 * nT ), ncol = 3 )
    for ( i in 1:3 ) {
      fend <- min( nT, length( fimp[[k]] ) )
      k <- 2 * i - 1
      imp[1:fend,i] <- fimp[[k]][1:fend]
    }
    
    # Package results in glyph-ready form
    G <- data_frame( time = time, Butter = imp[,1], Cheby = imp[,2], Ellip = imp[,3] )
    Gimp <- G %>% gather( Type, impulse, -time )
    G <- data_frame( freq = freq, Butter = amp[,1], Cheby = amp[,2], Ellip = amp[,3] )
    Gamp <- G %>% gather( Type, amp, -freq )
    G <- data_frame( freq = freq, Butter = phs[,1], Cheby = phs[,2], Ellip = phs[,3] )
    Gphs <- G %>% gather( Type, phs, -freq )   

    # Create three separate ggplot graph dfs
    glyph_imp <- ggplot( Gimp, aes( x = time, y = impulse ) ) +
      geom_path( aes( color = Type ), size = 1.25 ) +
      labs( x = "Time", y = "Impulse Response" ) +
      scale_colour_brewer( palette = "Dark2", 
                           breaks = c( "Theory", "Butter", "Cheby", "Ellip" ),
                           direction = -1,
                           guide = guide_legend( direction = "horizontal" ) ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ),
             legend.key.width = unit( 4, "picas" ),
             legend.position = "bottom" ) 
      

    glyph_amp <- ggplot( Gamp, aes( x = freq, y = amp ) ) +
      geom_path( aes( color = Type ), size = 1.25 ) +
      labs( x = "Frequency", y = "Amplitude" ) +
      scale_colour_brewer( palette = "Dark2", 
                           breaks = c( "Theory", "Butter", "Cheby", "Ellip" ),
                           direction = -1 ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ),
             legend.position = "none" ) 
    
    glyph_phs <- ggplot( Gphs, aes( x = freq, y = phs ) ) +
      geom_path( aes( color = Type ), size = 1.25 ) +
      labs( x = "Frequency", y = "Phase" ) +
      scale_colour_brewer( palette = "Dark2", 
                           breaks = c( "Theory", "Butter", "Cheby", "Ellip" ),
                           direction = -1 ) +
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

