#
# Practicing with conditional panels.

require( shiny, quiet = TRUE )
require( tidyverse, quiet = TRUE )
require( Rssa, quiet = TRUE )
require( multitaper, quiet = TRUE )
require( cowplot, quiet = TRUE )

data_loaded <<- FALSE
data_file <<- c( "/Users/justinr/R-Studio/Projects/Shiny/SSA/co2", 
                 "/Users/justinr/R-Studio/Projects/Shiny/SSA/zuerich", 
                 "/Users/justinr/R-Studio/Projects/Shiny/SSA/laskar" )

source( "~/R-Studio/Projects/Shiny/SSA/parseComp.R" )
theme_set( theme_bw( ) )
options( warn = -1 )

# Define UI for application that draws a histogram
ui <- shinyUI( fluidPage(
   
  fluidRow(
    column( 12,
            h4( "SSA Visualizer", align = "center" ),
            conditionalPanel( "input.panel == 1", plotOutput( "plot1" ) ),
            conditionalPanel( "input.panel == 2", plotOutput( "plot2" ) ) )
  ),
  fluidRow(
      column( 4, 
        selectInput( "select", label = h4( "File Selector" ), 
                     choices = list( "Mauna Loa CO2" = 1, 
                                     "Zuerich Sunspots" = 2, 
                                     "Laskar Insolation" = 3 ), 
                     selected = 1 ),
        actionButton( "load", "Load Data" ),
        hr(),
        radioButtons( "panel",
                       "View Panel",
                       c( "Time Series" = 1,
                          "Spectra" = 2 ) )
      ),
      column( 4,
        h4( "SSA Parameters" ) ,
        sliderInput( "lag", "Window Length (%)", 
                     value = 0.2,
                     min = 0.05, 
                     max = 0.5, 
                     step = 0.05 ),
        radioButtons( "svd",
                      "SVD Method",
                      c( "Auto" = "auto",
                         "SVD" = "svd",
                         "Eigen" = "eigen" ),
                      inline = TRUE ),
        radioButtons( "scale", 
                      label = "SSA Spectrum Scaling",
                      c( "Log" = 1,
                         "Linear" = 2 ), 
                      inline = TRUE )
      ),
      column( 4,
        h4( "Reconstruction/Components" ),
        textInput( "recon", "Reconstruction List", value = "1" ),
        helpText( "E.g., 1, 3, 5:7", "Reconstruction in red" ),
        checkboxInput( "trend", "Remove Trend On Load", value = FALSE )
      )
   )
) )


# Define server logic required to draw a histogram
server <- shinyServer( function( input, output ) {
   
   dataIn <- eventReactive( input$load, {
     load( file = data_file[as.integer( input$select )] ) 
     x <- switch( input$select,
                  "1" = co2,
                  "2" = zuerich,
                  "3" = laskar )
     
     if ( input$trend ) {
       Ts <- data_frame( TS = x, t = seq( 1, length( x ) ) )
       mod <- lm( TS ~ poly( t, 3 ), data = Ts )
       x <- residuals( mod )
     }
     list( read = TRUE, data = x )
   })
   
   output$plot1 <- renderPlot({
     if ( !dataIn()[[1]] ) return()
     d <- dataIn()[[2]]
     age <- time( d )
     n <- floor( input$lag * length( d ) )
     s <- ssa( d, n, svd.method = input$svd )
     nSigma <- length( s$sigma )
     D <- data_frame( Data = d, age = age )

     # Need code to skip recon if recon panel is null
     glist <- parseComp( input$recon )
     if ( !is.na( glist )  ) {
       r <- reconstruct( s, groups = glist )
       nRecons <- length( r )
       rStrings <- strsplit( input$recon, "," )
       for ( i in 1:nRecons ) {
         rName <- paste( "R", trimws( rStrings[[1]][i] ), sep = "" )
         D[,rName] <- as.numeric( r[[i]] )
       }
       res <- attr( r, "series" ) - attr( r, "residuals" )
       H <- data_frame( Data = res, age = age ) %>% 
         gather( key = Set, value = Value, -age )
     }
     F <- D %>% gather( key = Set, value = Value, -age )
     
     # Make the plot
     Gts <- ggplot( F, aes( x = age, y = Value ) ) + geom_path( ) +
       facet_grid( Set ~ ., scales = "free" ) +
       theme( text = element_text( size = 20 ) ) +
       labs( x = "Time", y = "Value" )
     
     if ( !is.na( glist ) )
       Gts <- Gts + geom_path( data = H, color = "red" )
    
     Gts
   })
   
   output$plot2 <- renderPlot({
     if ( !dataIn()[[1]] ) return()
     d <- dataIn()[[2]]
     n <- floor( input$lag * length( d ) )
     s <- ssa( d, n, svd.method = input$svd )
     
     # Build Singular Value Spectrum Plot
     nSigma <- length( s$sigma )   
     D <- data_frame( N = seq( 1, nSigma ), SV = s$sigma )
     Gsigma <- D %>% ggplot( aes( x = N, y = SV ) ) +
       geom_point( size = 3 ) +
       labs( x = "Number", y = "Singular Value" ) +
       theme( text = element_text( size = 20 ) )
       
     if ( input$scale == 1 ) Gsigma <- Gsigma + scale_y_log10()
     
     # Build data + recon spectra 
     dSpec <- spec.mtm( d, nw = 4, k = 7, deltat = 1, plot = FALSE )
     E <- data_frame( f = dSpec$freq, Data = dSpec$spec )
     glist <- parseComp( input$recon )
     if ( !is.na( glist )  ) {
       r <- reconstruct( s, groups = glist )
       nRecons <- length( r )
       rStrings <- strsplit( input$recon, "," )
       for ( i in 1:nRecons ) {
         rName <- paste( "R", trimws( rStrings[[1]][i] ), sep = "" )
         E[,rName] <- spec.mtm( r[[i]], nw = 4, k = 7, deltat = 1, plot = FALSE )$spec
       }
     }
     F <- E %>% gather( key = Set, value = Value, -f )     
     
     Gspec <- ggplot( F, aes( x = f, y = Value ) ) + geom_path( ) +
       facet_grid( Set ~ ., scales = "free" ) +
       theme( text = element_text( size = 20 ) ) +
       labs( x = "Frequency", y = "PSD" ) +
       scale_y_log10()

     # Cowplot the pair
     ggdraw() +
       draw_plot( Gsigma, 0, 0, 0.4, 1 ) +
       draw_plot( Gspec, 0.4, 0, 0.6, 1 )
    
   })
   
})

# Run the application 
shinyApp( ui = ui, server = server )

