#
# Listen Up - Open a locally stored mp3 or wav file.  Play it with various
#             Butterworth filters applied.
#

library(shiny)
library(signal)
library(tidyverse)
library(cowplot)
library(audio)
library(tuneR)

theme_set( theme_bw() )
options( warn = -1 )
nSegments <- 200
xFiltered <- 0

# Define UI for application that draws a histogram
ui <- shinyUI( fluidPage(
   tags$style(type = "text/css", 
              "html, body {width:95%;height:100%;margin-left:auto; margin-right:auto;}"),
  
   # Application title
   h4( " Listen Up" ),
   
   # Sidebar filter parameters 
   sidebarLayout(
      sidebarPanel(
          fileInput( "fileID", "Sound File", multiple = FALSE, buttonLabel = "Browse..." ),
          sliderInput("order", "Filter Order", min = 1, 
                       max = 10, value = 3, step = 1 ),
          sliderInput("freqs", "Frequency Band (Hz)", min = 10, 
                      max = 10000, value = c(10,10000), step = 10 ),
          radioButtons( "type", "Filter Type:",
                        c( "Low Pass" = "low",
                           "High Pass" = "high",
                           "Band Pass" = "pass",
                           "Band Reject" = "stop" ),
                        selected = "pass" ),
          actionButton( "playIt", "Play" )
      ),
      
      # Plot panel
      mainPanel(
         plotOutput( "FPlot" )
      )
   )
))

# Server R code
server <- shinyServer( function( input, output) {
  
  observeEvent( input$playIt, {
    if ( is.null( input$fileID ) ) return( )
    xFiltered <- signal::filter( butt, rawAudio()[[1]] )
    ll <- which( xFiltered > 1 )
    xFiltered[ll] <- 1
    ll <- which( xFiltered < -1 )
    xFiltered[ll] <- -1
    audio::play( xFiltered, rawAudio()[[2]] )
  } )
  
  rawAudio <- eventReactive( input$fileID, {
    rawAudio <- readMP3( input$fileID$datapath )
    x <- as.numeric( rawAudio@left / max( rawAudio@left ) )
    sampRate <- rawAudio@samp.rate
    nSamp <- length( x )
    winSec <- nSamp / sampRate / nSegments
    nWin <- winSec * sampRate
    nWin <- nextn( nWin, factors = c( 2, 3, 5 ) )
    nSeg <- floor( nSamp / nWin )
    Mspec <- matrix( rep( 0, nWin / 2 * nSeg ), ncol = nSeg )
    nrow( Mspec )
    nstart <- 1
    for ( i in 1:nSeg ) {
      nend <- nstart + nWin - 1
      yy <- spectrum( x[nstart:nend], plot = FALSE )
      Mspec[,i] <- spectrum( x[nstart:nend], plot = FALSE )$spec
      nstart <- nend + 1
    }
    spec <- apply( Mspec, 1, mean )
    list( x, sampRate, spec )
  } )
  
  
  output$FPlot <- renderPlot( height = 500, {
    if( is.null( input$fileID ) ) return( )
    if ( is.null( rawAudio() ) ) return( )

    # Plot power spectrum of recording
    spec <- rawAudio()[[3]]
    nFreq <- length( spec )
    sampRate <- rawAudio()[[2]]
    dt <- 1 / sampRate
    nyquist <- sampRate / 2
    Gspec <- data_frame( freq = seq( 0, nyquist, length.out = nFreq ), spec = spec )
    glyph_spec <- ggplot( Gspec, aes( x = freq, y = spec ) ) + 
      geom_line( color = "darkred" ) +
      labs( x = "Frequency (Hz)", y = "Audio Spectrum" ) +
      scale_y_log10() +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ) ) 
    
    # Compute filter
    if ( input$type %in% c( "pass", "stop" ) ) 
      W <- input$freqs
    else if ( input$type %in% c( "low" ) )
      W <- input$freqs[1]
    else
      W <- input$freqs[2]
    
    nyquist <- sampRate / 2
    W <- W / nyquist
    butt <<- butter( n = input$order, W = W, type = input$type, plane = "z" )

    # Compute filter amplitude and phase
    fhf <- freqz( butt )
    freq <- fhf[[2]] / pi * nyquist
    amp <- Mod( fhf[[1]] )
    phs <- unwrap( Arg( fhf[[1]] ) ) 

    # Compute impulse response
    fimp <- impz( butt )
    nT <- max( 100, length( fimp[[2]] ) )
    imp <- rep( 0, nT )
    time <- seq( 0, by = dt, length.out = nT )
    fend <- min( nT, length( fimp[[2]] ) )
    imp[1:fend] <- fimp[[1]][1:fend]
    
    # Package results in glyph-ready form
    Gimp <- data_frame( time = time, Butter = imp )
    Gamp <- data_frame( freq = freq, Butter = amp )
    Gphs <- data_frame( freq = freq, Butter = phs )

    # Create three separate ggplot graph dfs
    glyph_imp <- ggplot( Gimp, aes( x = time, y = Butter ) ) +
      geom_path( color = "darkorange", size = 1.25 ) +
      labs( x = "Time", y = "Impulse Response" ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ) ) 

    glyph_amp <- ggplot( Gamp, aes( x = freq, y = Butter ) ) +
      geom_path( color = "darkblue", size = 1.25 ) +
      labs( x = "Frequency", y = "Filter" ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ) ) 
    
    glyph_phs <- ggplot( Gphs, aes( x = freq, y = Butter ) ) +
      geom_path( color = "darkorange", size = 1.25 ) +
      labs( x = "Frequency", y = "Phase" ) +
      theme( text = element_text( size = 20 ),
             plot.title = element_text( size = 20 ),
             plot.margin = margin( 10, 0, 0, 10, "pt" ) ) 
    
    # Use cowplot to array these
    g <- ggdraw( ) +
      draw_plot( glyph_spec, 0, 0.4, 1, 0.6 ) +
      draw_plot( glyph_amp, 0, 0, 1, 0.4 ) 
#      draw_plot( glyph_imp, 0.5, 0, 0.5, 0.4 )
    g
  })
})

# Run the application 
shinyApp( ui = ui, server = server )

