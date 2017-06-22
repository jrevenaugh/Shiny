library(shiny)
library(leaflet)
require(lubridate)
require(magrittr)
require(tidyverse)
require(leaflet)

load( file = "~/R-Studio/Projects/Shiny/Seis-Mo/trimmed_anss" )
eqcat <- eqcat %>% filter( Magnitude >= 5 ) %>%
         select( DateTime, Magnitude, Depth, Latitude, Longitude ) 
eqcat$ddate <- decimal_date( as.POSIXct( eqcat$DateTime ) )

# Set up depth-dependent colors
nColors <- 5
eqcat$dCol <- ceiling( ( eqcat$Depth / 700 ) * nColors )
pal <- colorRampPalette( c("blue", "red", "green" ) )(nColors)
colorPal <- function( depth ) pal[depth]

nEQ <- nrow( eqcat )
popLabelList <- c( "Magnitude", "Depth", "DateTime" )
popLabel <- 1

ui <- bootstrapPage(
  tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
  leafletOutput( "map", width = "100%", height = "100%" ),
  absolutePanel(top = 10, right = 10,
                h4('1990-2017 Seismicity'), 
                sliderInput("depth_range", "Depth Range (km)", min = 0, 
                            max = 800, value = c(0, 800), step = 10 ),
                sliderInput("mag_range", "Magnitude Range", min = 5, 
                            max = 10, value = c(5, 10), step = 0.25 ),
                sliderInput("time_range", "Date Range", min = 1990, 
                            max = 2018, value = c(1990, 2018), step = 0.5,  sep = "" ),
                selectInput("popper", "Pop-up Label", 
                            choices = list( "Magnitude" = 1, "Depth" = 2, "Date" = 3 ), 
                            selected = popLabel ), 
                uiOutput( "eqcount" ) )
)

server <- function(input, output, session) {
  
  # Reactive expression for the data subsetted to what the user selected
  filteredData <- reactive({
    m_range <- input$mag_range
    d_range <- input$depth_range
    t_range <- input$time_range
    popLabel <- as.integer( input$popper )
    w_start <- t_range[1]
    w_end <- t_range[2]
    tday <- decimal_date( Sys.Date() )
    w_end <- min( w_end, tday )
    tcat <- eqcat %>% filter( ddate >= w_start & ddate <= w_end ) %>%
      filter( Depth >= d_range[1] & Depth <= d_range[2] ) %>%
      filter( Magnitude >= m_range[1] & Magnitude <= m_range[2] )
    nEQ <- nrow( tcat )
    output$eqcount <- renderUI ({
      tagList (
        h5( paste( "EQs: ", format( nEQ ) ) )
      )
    })
    list( cat = tcat, popLabel = popLabel )
    
  })
  
  output$map <- renderLeaflet({
    m <- leaflet( options = leafletOptions( minZoom = 2 ) ) %>%
        fitBounds( -180,-70, 180, 80 ) %>%
#       addProviderTiles( providers$Stamen.TerrainBackground ) %>%
        addProviderTiles( providers$Esri.OceanBasemap ) %>%
        addScaleBar( position = "bottomright" )
  })
  
  # Incremental changes to the map (in this case, replacing the
  # circles when a new color is chosen) should be performed in
  # an observer. Each independent set of things that can change
  # should be managed in its own observer.
  
  observe({
      nEQ <- nrow( filteredData()[[1]] )
      popLabel <- filteredData()[[2]]
      labn <- popLabelList[popLabel]
      proxy <- leafletProxy( "map", data = filteredData() )
      if ( nEQ <= 0 ) {
         proxy %>% 
         clearMarkers() 
      } else {
          proxy %>% 
          clearMarkers() %>%
          addCircleMarkers( lng = ~Longitude, lat = ~Latitude, 
                            popup = ~paste( "<strong>",substr( DateTime, 1, 10 ), "</strong><br>", Depth, " km<br>", Magnitude ),
                            label = as.formula( paste( "~as.character(", labn, ")" ) ), 
                            color = NA,
                            fillColor = ~colorPal( dCol ),
                            fillOpacity = 0.5,
                            radius = ~1.5 * exp( ( Magnitude - 5 ) / 1.5 ),
                            data = filteredData()[[1]] )
      }
      
  })
}

shinyApp(ui, server)