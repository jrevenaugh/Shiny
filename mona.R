# Mona - Shiny app illustrating the use of the SVD as a 2D filter.  
#        Useful for demonstrating principal components.
#

load( file = "~/R-Studio/Projects/Shiny/Mona Lisa/mona_mat" )
require( tidyverse, quiet = TRUE )
require( shiny, quiet = TRUE )
theme_set( theme_bw() )
options( warn = -1 )
pal <- gray.colors( 256, start = 0, end = 1 )


ui <- shinyUI( fluidPage(
  h4( "Mona Lisa SVD" ),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sv", "Singular Values", min = 1, 
                  max = 50, value = c(1, 10), step = 1 ),
      plotOutput( "SVPlot" )
    ),
    
    # Plot panel
    mainPanel(
      plotOutput( "FMona", height = "auto", width = "auto" )
    )
  )
))

# Server R code
server <- shinyServer(function(input, output) {
  
  aspect <- 200 / 298
  res <- 72
  size <- 8
  h <- size * res
  w <- size * res * aspect
  output$FMona <- renderPlot ( height = h, width = w, {
    svd <- svd( mona )
    s <- rep( 0, 200 )
    s[input$sv[1]:input$sv[2]] <- svd$d[input$sv[1]:input$sv[2]]
    D <- diag( s )
    fmona <- svd$u %*% D %*% t( svd$v )
    par( mar = c( 0, 0, 0, 0 ) )
    image( z = fmona, col = pal, axes = FALSE, xlab = NA, ylab = NA )
  })
  
  output$SVPlot <- renderPlot ({
    s <- svd( mona )
    S <- data_frame( Number = seq( 1, 50 ), SValue = s$d[1:50] )
    P <- data_frame( Number = input$sv[1]:input$sv[2], SValue = s$d[input$sv[1]:input$sv[2]] )
    drng <- range( s$d[1:50] )
    ymin <- drng[1]
    ymax <- drng[2] * 1.1
    xmin <- input$sv[1] - 0.5
    xmax <- input$sv[2] + 0.5
    w <- input
    R <- data_frame( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax )
    gSV <- ggplot() +
      geom_rect( data = R,
                 aes( xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax ),
                 alpha = 0.2, color = "gray" ) +
      geom_point( data = S, inherit.aes = FALSE, aes( x = Number, y = SValue )) +
      geom_point( data = P, aes( x = Number, y = SValue ), color = "red" ) +
      labs( x = "Number", y = "Singular Value" ) +
      theme( text = element_text( size = 20 ) ) + xlim( 0.4, 50.6 ) +
      scale_y_log10()
    
    gSV
  })
})

# Run the application 
shinyApp( ui = ui, server = server )

