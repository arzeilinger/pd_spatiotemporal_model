library(shiny)
library(animation)
library(raster)

ui <- fluidPage(
  sliderInput(inputId = "alpha",
              label = "Choose a value for dispersal (alpha)",
              value = 0.5, min = 0, max = 1),
  sliderInput(inputId = "beta",
    label = "Choose a value for vine-to-vine spread (beta)",
    value = 0.5, min = 0, max = 1),
  sliderInput(inputId = "epsilon",
              label = "Choose a value for external infection rate (epsilon)",
              value = 0.001, min = 0, max = 0.1),
  plotOutput("raster")
)


###################################################################################################################




server <- function(input, output) {
  output$raster <- renderPlot({
    #### Load simulation function and dispersal kernel functions
    source("simulateDiseaseSpread.R")
    source("dispersal_kernel_functions.R")
    #### Simulation for visualization
    #### Initial values setup
    ## Parameter values
    alpha <- input$alpha
    beta <- input$beta
    epsilon <- input$epsilon
    ## Other values
    Tmax <- 100
    ## Number of rows and columns -- will always make the plant population a square
    ## Deviating from a square causes a problem with making time slice rasters -- not sure why
    nrc <- 20
    grid <- 1:nrc
    Coo <- matrix(c(rep(grid, nrc), rep(grid, each = nrc)), ncol = 2)
    testSpread <- simulateDiseaseSpread(alpha, beta, epsilon, Tmax, Coo)
    Inf_times <- testSpread$Inf_times
    #### Plotting infection dynamics
    ## Produce maps of infection status based on different time points t
    ## Vector of time points to plot
    tcuts <- floor(seq(5,100,length.out=9))
    #### Produce raster stack of infection status at a series of time points
    ## Get raster dimensions from Coo
    ## Coo[,1] = rows or y coord
    ## Coo[,2] = columns or x coord
    CooRows <- Coo[,1]
    CooNrows <- length(unique(CooRows))
    CooYmn <- min(CooRows)
    CooYmx <- max(CooRows)
    CooColumns <- Coo[,2]
    CooNcols <- length(unique(CooColumns))
    CooXmn <- min(CooColumns)
    CooXmx <- max(CooColumns)
    ## Empty raster stack
    diseaseRasterList <- vector("list", length(tcuts))
    layerNames <- rep(NA, length(tcuts))
    for(i in 1:length(tcuts)){
      ## Produce vector of infection status: 1 = infected, 0 = susceptible at a given time point t
      infectionStatus.i <- as.integer(ifelse(Inf_times < tcuts[i], 1, 0))
      layerNames[i] <- paste("After", tcuts[i], "days", sep = "_")
      raster.i <- raster(nrows = CooNrows, ymn = CooYmn, ymx = CooYmx,
                         ncols = CooNcols, xmn = CooXmn, xmx = CooXmx,
                         vals = infectionStatus.i)
      diseaseRasterList[[i]] <- raster.i
    }
    diseaseRasterStack <- stack(diseaseRasterList, quick = TRUE)
    names(diseaseRasterStack) <- layerNames
    #giftitle <- rep("Vineyard disease spread", dim(diseaseRasterStack)[3])
    #raster::animate(diseaseRasterStack, pause = 0.5)
    ## Specify raster colors
    breakpoints <- c(0,0.25,0.75,1,1.25)
    colors <- c("grey", "grey", "darkgreen", "darkgreen", "darkgreen")
    plot(diseaseRasterStack, legend = FALSE, axes = FALSE, box = FALSE,
         breaks = breakpoints, col = colors)
    })
}

shinyApp(ui = ui, server = server)