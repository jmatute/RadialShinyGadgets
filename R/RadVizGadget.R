#' RadViz Gadget
#'
#' Creates a RShiny Gadget for RadViz
#'
#' RadViz's  goal is to generate a configuration which reveals the underlying nature of the data for cluster analysis, 
#' outlier detection, and exploratory data analysis, e.g., by investigating the effect of specific dimensions on the separation of the data.  
#' Each dimension is assigned to a point known as dimensional anchors across a unit-circle.  Each sample is projected according to the 
#' relative attraction to each of the anchors. 
#' 
#' It is defined defined for multidimensional numerical data sets \eqn{X=\{\mathbf{p}_1,\ldots, \mathbf{p}_N\},} 
#' for N data points \eqn{\mathbf{x}_i \in \mathbf{R}^{d}} of dimensionality d. Let \eqn{A =\{ \mathbf{a}_{1}, \dots, \mathbf{a}_{d}  \} ,} 
#' be a set of (typically 2D) anchors, each corresponding to one of the d dimensions.
#' The projection \eqn{ \mathbf{p}_i' \in \mathbf{R}^{2},} of a multidimensional point \eqn{ \mathbf{p}_i = (p_{i1},\ldots,p_{id}) \in \mathbf{R}^{d}, }in SC is then defined as:
#' \deqn{ \mathbf{x}_i' = \frac{ \sum_{j=1}^{d} \mathbf{a}_{j} g_j( \mathbf{p}_i)}{\sum_{j=1}^{d} \mathbf{a}_{j}  },} with
#' \deqn{ g_j(\mathbf{p}_i) = \frac{p_{ij} - min_j}{max_j - min_j} ,} and \eqn{(min_j,max_j),}denoting the value range of dimension j. 
#' The dimensional anchors can be moved either interactively or algorithmically to reveal different meaningful patterns in the dataset. 
#' 
#' 
#' @param df  A dataframe with the data to explore. It should contain only numeric columns (with the exception of the label column).
#' @param color column where labels from the data are extracted.
#' @return A list location of the anchors, coordinates of the projected samples and a logical vector with the selected samples
#' @export
#' @references 
#' Sharko, J., Grinstein, G., & Marx, K. A. (2008). Vectorized radviz and its application to multiple cluster datasets. IEEE transactions on Visualization and Computer Graphics, 14(6), 1444-1427.
#' @examples
#' if (interactive()) {
#'  library(RadialVisGadgets)
#'  library(datasets)
#'  data(iris)
#'  RadViz(iris, "Species")
#' }
#' 
RadViz <- function(df, color = NULL){
    # RadViz can be defined generally in 3 steps
    # 1. Normalize the data interval
    # 2. Place the dimensional anchors
    # 3. Calculate the point where each sample should be placed
    
    odata <- df # Store original data somewhere
    df <- removeNonZeroAndMissing(df)
    
    
    
    inputCheck <- inputValidationRadViz(df, color)
    if (!is.null(inputCheck)){
      if (inputCheck$stop)  stop(inputCheck$errorMessage)
      else warning(inputCheck$errorMessage)
    }
    colorVar <- color
    if (!is.null(colorVar)  && (colorVar %in% colnames(df))){
      df[[colorVar]] <- NULL
    }
    else if (!is.null(colorVar)){
      colorVar <- NULL
    }
  
    
    scaledData <- getRangedData(df,  meanCentered = FALSE)
    
    ######################################################################
    # throttle given the type of interaction with the anchors
    ui <- miniPage(
      useShinyjs(),
      gadgetTitleBar("RadViz"),
      miniContentPanel(
        plotOutput("plot", height = "100%", brush="plotBrush", hover=hoverOpts(id="plot_hover", delayType="throttle", delay= 150),
                   dblclick = "plotDblClick")
      ),
      miniButtonBlock(
        actionButton("zoomIn","", icon = icon("search-plus")),
        actionButton("zoomOut","", icon = icon("search-minus")),
        actionButton("screenShot","", icon = icon("camera-retro"))
      )
    )
    #################################################################
    
    server <- function(input, output, session) {
      highlightedIdx <- reactiveVal(-1)
      
      helperValues <- reactiveValues()
      helperValues$plotLimit <- 1.5
      helperValues$initialized <- FALSE 
      helperValues$anchorLocations <- NULL 
      helperValues$movingDim <- -1
      helperValues$clicked <- FALSE
      
      observeEvent( input$plot_hover, {
        helperValues$curX <- input$plot_hover$x
        helperValues$curY <- input$plot_hover$y
        currentClosest <- closestDimension(helperValues$anchorLocations, helperValues$curX, helperValues$curY )
        
        if (helperValues$clicked ){
          # if we are moving the dimensional vectors, deactivate the brushing 
          session$resetBrush("plotBrush")
          
          newLoc <- c(helperValues$curX, helperValues$curY)
          newLoc <- newLoc / sqrt(sum(newLoc^2))
     
          helperValues$anchorLocations[helperValues$movingDim,1:2] <- newLoc 
          highlightedIdx(helperValues$movingDim)
        }
        else {
          highlightedIdx(currentClosest)
        }
        
      })
      
      onevent("mouseup", "plot", mouseUp())
      onevent("mousedown","plot", mouseDown())
      observeEvent(input$zoomIn,{ helperValues$plotLimit <-  helperValues$plotLimit*0.9 })
      observeEvent(input$zoomOut,{ helperValues$plotLimit <-  helperValues$plotLimit*1.1 })
      observeEvent(input$screenShot, { screenshot()})
      
      
      initValues <- function(){
           helperValues$anchorLocations <- initProjectionMatrix(df)
           helperValues$circle <- circleFun()
         
           helperValues$initialized <- TRUE
      }
    
      mouseDown <- function(){
        # detect where the mouse was pressed, shiny doesn't divide click in 3 parts, so we have to use the javascript 
        # events 
        currentClosest <- closestDimension(helperValues$anchorLocations, helperValues$curX, helperValues$curY )
        if ( currentClosest != -1){
          helperValues$movingDim <- currentClosest
          helperValues$clicked <- TRUE
        }
      }
      
      mouseUp <- function(){
        helperValues$clicked <- FALSE
        helperValues$movingDim <- -1
      }
      
      
      
      observeEvent(input$done, {
        projMatrix <- getCleanProjectionMatrix(helperValues$anchorLocations, colorVar, odata, df)
        projectedPoints <-getProjectedPointsRadViz(helperValues$anchorLocations, scaledData)
        selected <- brushedPoints(df=projectedPoints, brush=input$plotBrush, xvar="V1", yvar="V2", allRows = TRUE)
        result <- list(projMatrix, selected$selected_, projectedPoints)
        names(result) <- c("Anchors","Selection","Projected.Points")
        stopApp(result)
      })
      
      output$plot <- renderPlot({
        if( !helperValues$initialized) initValues()
        #draw the background circle
        sizes <- helperValues$anchorLocations$szs
        shapes <- helperValues$anchorLocations$shapes
        circle <- circleFun()

        if ( highlightedIdx() != -1) {
          sizes[highlightedIdx()] = 6
          shapes[highlightedIdx()] = 20
        }

        backgroundCircle <- geom_path(data=helperValues$circle, aes(.data$x, .data$y), alpha=0.2)
        curPlot <-   ggplot() + backgroundCircle
        
        projectedPoints <- getProjectedPointsRadViz(helperValues$anchorLocations, scaledData)
        curPlot <-  curPlot+ drawGGProjectedPoints(projectedPoints, odata,  color)
          
        # Draw the labels + the point that can be selected to move i.e. anchors 
        pts <- geom_point(data=helperValues$anchorLocations, aes(.data$V1,.data$V2), colour = "black", fill = "white", shape=shapes, size=sizes)
        labels <-  geom_label( data=helperValues$anchorLocations,  aes(.data$V1+.data$V1*0.2,.data$V2+.data$V2*0.2,label=rownames(helperValues$anchorLocations)))
        curPlot <- curPlot + pts + labels 
        
        curPlot +  coord_cartesian(xlim = c(-helperValues$plotLimit,helperValues$plotLimit ), ylim = c(-helperValues$plotLimit,helperValues$plotLimit )) + theme_void()
      },res=96)
      
    }
    
    viewer <- paneViewer(300)
    runGadget(ui, server, viewer=viewer)
}

