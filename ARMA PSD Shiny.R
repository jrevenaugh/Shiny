#
# ARMA - Create realizations of an ARMA process and estimate its power spectral
#        density using four different techniques.  Plot results.

library(shiny)
library(tidyverse)
library(ggplot2)
library(signal)
library(mosaic)
require(sapa, quietly = T, warn.conflicts = F)
require(ifultools, quietly = T, warn.conflicts = F)
require(splus2R, quietly = T, warn.conflicts = F)

theme_set( theme_bw() )

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

tukey <- function( q = 10 ) {
  q <- q + 1
  w <- 1/2 * ( 1 + cos( pi * seq( 0, q - 1, 1 ) / q ) )
  return( w )
}

btspec <- function( x = NA, q = 20, taper = T )
{
  if ( !is.vector( x ) ) return( "No time series provided." )
  if ( q < 1 || q > length( x ) / 2 ) return( "Process order out of bounds." )
  
  # Compute the Tukey-windowed ACF
  
  a <- acf( x, lag.max = q, plot = F )$acf[,1,1]
  if ( taper ) a <- a * tukey( q )
  
  # Now the slightly harder part.  Padding this out and putting negative lags on the far end.
  
  n <- 2048
  x1 <- rep( 0, n )
  l <- q + 1
  x1[1:l] <- a
  u <- n - q + 1
  x1[n:u] <- a[-1]
  
  # Now FFT the padded, mirrored acf function.
  
  f <- fft( x1 )
  u <- n / 2 + 1
  y <- Mod( f[1:u] ) * var( x )
  return( y )
}

# Define UI for application that draws a histogram
ui <- shinyUI( fluidPage(
       tags$style(type = "text/css", 
                  "html, body {width:90%;height:100%;margin-left:auto; margin-right:auto;}"),
       fluidRow( 
        plotOutput("specPlot")
      ),

      fluidRow( 
        column( 4, 
                textInput("arma_string", label = h5("ARMA Coeffs" ), 
                               value = "0.5 ; 0", width = '50%' ),
                h6( "ar1, ar2, ... ; ma1, ma2, ..." ),
                sliderInput("length", "Number of Samples", min = 50, 
                            max = 5000, value = 250, step = 50, sep ="" ),
                actionButton("doIt", label = "Compute") ), 
        column( 4,
                h5( "Non-Parametric Estimators" ),
                sliderInput("wosa", "WOSA Segments", min = 1, 
                            max = 30, value = 5, step = 1 ),
                sliderInput("tapers", "MTM Tapers", min = 1, 
                            max = 11, value = 3, step = 1 ) ),
        column( 4, 
                h5( "Parametric Estimators" ),
                sliderInput("lags", "BT Lags", min = 1, 
                            max = 100, value = 1, step = 1 ),
                sliderInput("ar", "AR Order", min = 1, 
                            max = 30, value = 1, step = 1 ) )
        )
    )
)

# Define server logic required to draw a histogram
server <- shinyServer( function( input, output ) {
  
  # Parse the UI input
  sigVals <- eventReactive(input$doIt, {
    arma.raw <- unlist( strsplit( input$arma_string, split = ";" ) )
    ar.raw <- strsplit( arma.raw[1], split = "," )
    if ( length( arma.raw ) > 1 ) ma.raw <- strsplit( arma.raw[2], split = "," )
    arc <- as.numeric( unlist( ar.raw ) )
    mac <- as.numeric( unlist( ma.raw ) )
    model <- list( b = mac, a = arc )
    list( model = model, lngth = input$length, 
          wosa = input$wosa, tapers = input$tapers, 
          lags = input$lags, order = input$order )
    
  })
  
   output$specPlot <- renderPlot( {
     options( warn = -1 )
     if ( sigVals()[[2]] > 0 ) {
        model  <- sigVals()[[1]]
        lngth  <- sigVals()[[2]]
        wosa   <- sigVals()[[3]]
        tapers <- sigVals()[[4]]
        lags   <- sigVals()[[5]]
        order  <- sigVals()[[6]] 
     }
     # Create the realization
     m <- list( ma = model[[1]], ar = model[[2]] )
     ts <- arima.sim( m, lngth )
     TS <- data_frame( t = seq( 1, lngth ), ts = as.numeric( ts ) )
     ys <- max( abs( as.numeric( ts ) ) )
     gTS <- ggplot( TS, aes( x = t, y = ts ) ) + 
            geom_path( color = "deepskyblue3", size = 1.5 ) + 
            theme( text = element_text( size = 20 ) ) +
            labs( x = "Time", y = "TS" ) +
            coord_fixed( ratio = 0.025 * lngth / ys, expand = FALSE ) +
            theme( plot.margin = margin( 0, 20, 0, 10, "pt" ) )

     # Compute the theoretical spectrum
     m_arma <- Arma( b = c( 1.0, model[[1]] ), a = c( 1.0, -model[[2]] ) )
     psTheory <- freqz( m_arma, n = 1024, Fs = 1 )
     M.df <- data_frame( f = psTheory$f, Theory = Mod( psTheory$h )^2 )

     # Compute the WOSA spectral estimate
     wps <- SDF( ts, method = "wosa", blocksize = floor( lngth/wosa ) )
     f <- seq( 0, 0.5, length.out = length( wps ) )
     D <- data_frame( f = f, wps = wps )
     f1 <- connector( data = D, wps ~ f )
     M.df$WOSA <- f1( M.df$f )

     # Compute the MTM spectral estimate
     mps <- SDF( ts, method = "multitaper", n.taper = tapers )
     D <- data_frame( f = f, wps = mps )
     f1 <- connector( data = D, wps ~ f )
     M.df$MTM <- f1( M.df$f )

     # Compute the Blackman-Tukey spectral estimate
     bps <- btspec( x = as.vector( ts ), q = lags, taper = TRUE )
     D <- data_frame( f = seq( 0, 0.5, length.out = 1025 ), wps = bps )
     f1 <- connector( data = D, wps ~ f )
     M.df$BT <- f1( M.df$f )

     # Compute the AR spectral estimate
     aps <- spec.ar( ts, n.freq = 1025, order = order, plot = F  )
     D <- data_frame( f = aps[[1]], wps = aps[[2]][,1] )
     f1 <- connector( data = D, wps ~ f )
     M.df$AR <- f1( M.df$f )

     # Plot them up after gathering columns
     W.df <- M.df %>% gather( Method, PS, -f )
     ps <- max( W.df$PS )
     gPS <- ggplot( W.df, aes( x = f, y = PS ) ) + 
              geom_line( aes( color = Method ), size = 1.25 ) +
              scale_colour_brewer( palette = "Dark2", 
                                   breaks=c( "Theory", "WOSA", "MTM", "BT", "AR" ),
                                   direction = -1 ) +
              theme( text = element_text( size = 20 ),
                     plot.title = element_text( size = 20 ),
                     plot.margin = margin( 10, 0, 0, 10, "pt" ),
                     legend.key.width = unit( 4, "picas" ) ) +
              labs( x = "Frequency", y = "PSD Estimate" ) + 
              coord_cartesian( xlim = c( 0, 0.5 ), ylim = c( 0, 1.05 * ps ), expand = FALSE ) +
              ggtitle( "Spectral Estimation")
     

     # Use multiplot to display
     multiplot( gPS, gTS, cols = 1 )
   })
})

# Run the application 
shinyApp( ui = ui, server = server )

