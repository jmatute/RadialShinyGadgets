#' Star Coordinates Gadget
#'
#' Creates a RShiny Gadget for Star Coordinates
#'
#' Star Coordinate's (SC) goal is to generate a configuration which reveals the underlying nature of the data for cluster analysis, 
#' outlier detection, and exploratory data analysis, e.g., by investigating the effect of specific dimensions on the separation of the data.  
#' Traditional SC are defined for multidimensional numerical data sets \eqn{X=\{\mathbf{p}_1,\ldots, \mathbf{p}_N\},} 
#' for N data points \eqn{\mathbf{x}_i \in \mathbf{R}^{d}} of dimensionality d. Let \eqn{A =\{ \mathbf{a}_{1}, \dots, \mathbf{a}_{d}  \} ,} be a set of (typically 2D) vectors, each corresponding to one of the d dimensions.
#' The projection \eqn{ \mathbf{p}_i' \in \mathbf{R}^{2},} of a multidimensional point \eqn{ \mathbf{p}_i = (p_{i1},\ldots,p_{id}) \in \mathbf{R}^{d}, }in SC is then defined as:
#' \deqn{ \mathbf{x}_i' = \sum_{j=1}^{d} \mathbf{a}_{j} g_j( \mathbf{p}_i),} with
#' \deqn{ g_j(\mathbf{p}_i) = \frac{p_{ij} - min_j}{max_j - min_j} ,} and \eqn{(min_j,max_j),}denoting the value range of dimension j. 
#' 
#' In the case of categorical dimensions, the values when numericRepresentation= TRUE are mapped into numerical type i.e. as.numeric()
#' However equally spaced categorical points may not reflect the true nature of the data.  Instead, a frequency-based 
#' representation may be applied for individual data points. 
#' Assuming a categorical dimension j, we calculate the frequency \eqn{f_{jk},} of each category k of dimension j. 
#' The respective axis vector  \eqn{\mathbf{a}_{j},}  is then divided into according blocks, whose size represent the relative frequency (or probability) 
#' \eqn{\frac{f_{jk}}{\sum_{l=1}^m f_{jl}},} of each of the m categories of dimension j.
#' 
#' In summary, given an order for each categorical dimension, the Equation \eqn{g(),} above can be extended  to SC for mixed data by:
#' \deqn{ g_j(\mathbf{x}_i) = F_j(x_{ij}) - \frac{P_j(x_{ij})}{2}   , } if categorical/ordinal
#' \deqn{ g_j(\mathbf{x}_i) = \frac{x_{ij} - min_j}{max_j - min_j}  ,} if numerical
#' 
#' where \eqn{F_j,} is the cumulative density function for (categorical/ordinal) dimension j and \eqn{P_j,} its probability function. 
#' @param df  A dataframe with the data to explore. It should contain only numeric or factor columns.
#' @param color column where labels from the data are extracted.
#' @param approach Standard approach as defined by Kandogan, or Orthographic Star Coordinates (OSC) with a recondition as defined by Lehmann and Thiesel
#' @param numericRepresentation if true attempt to convert all factors to numeric representation, otherwise used mixed representation as defined in Hinted Star Coordinates
#' @param meanCentered center the projection at the mean of the values. May allow for easier value estimation
#' @param projMatrix a pre-defined projection matrix as an initial configuration. Should be defined in the same fashion as the output
#' @param clusterFunc function to define hints, assume increase in value of the function is an increase in quality of the projection. The function will be called with two parameters (points, labels)
#' @return A list with the projection matrix, coordinates of the projected samples and a logical vector with the selected samples
#' @references 
#' Kandogan, E. (2001, August). Visualizing multi-dimensional clusters, trends, and outliers using star coordinates. In Proceedings of the seventh ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 107-116).
#' 
#' Lehmann, D. J., & Theisel, H. (2013). Orthographic star coordinates. IEEE Transactions on Visualization and Computer Graphics, 19(12), 2615-2624.
#' 
#' Rubio-Sánchez, M., & Sanchez, A. (2014). Axis calibration for improving data attribute estimation in star coordinates plots. IEEE transactions on visualization and computer graphics, 20(12), 2013-202
#' 
#' Matute, J., & Linsen, L. (2020, February). Hinted Star Coordinates for Mixed Data. In Computer Graphics Forum (Vol. 39, No. 1, pp. 117-133).
#' 
#' @export
#' @examples
#' if (interactive()) {
#'  library(RadialVisGadgets)
#'  library(datasets)
#'  data(iris)
#'  StarCoordinates(iris, "Species")
#' }
#' 
StarCoordinates <- function(df, color = NULL,  approach="Standard", numericRepresentation=TRUE, meanCentered = TRUE, projMatrix=NULL, clusterFunc = NULL) {

  #######################################################################
  # Simple Pre-processing, validate input to function and create data that
  # won't change throughout the running of the gadget
  colorVar <- color #prefer nomenclature for the parameters as color, but semantically prefer colorVar
  # check for simple errors and stop in case error if bad enough
  inputCheck <- inputValidationSC(df,  colorVar,  approach, projMatrix, clusterFunc )
  if (!is.null(inputCheck)){
      if (inputCheck$stop)  stop(inputCheck$errorMessage)
      else warning(inputCheck$errorMessage)
  }
  odata <- df # Store original data somewhere
  df <- removeNonZeroAndMissing(df)

  # 
  frequencyList <- NULL
  if (!is.null(colorVar)  && (colorVar %in% colnames(df))){
      df[[colorVar]] <- NULL
  }
  else if (!is.null(colorVar)){
       colorVar <- NULL
  }
  
  # if it is a numerical representation, simply try to convert 
  if(numericRepresentation){
       df <-  mutate_if(df, is.factor, as.numeric)
  }
  else {
       frequencyList <- getFrequencyList(df)
  }
  ######################################################################
  # throttle given the type of interaction with the vectors 
  ui <- miniPage(
    useShinyjs(),
    gadgetTitleBar("Star Coordinates"),
    miniContentPanel(
      plotOutput("plot", height = "100%", brush="plotBrush", hover=hoverOpts(id="plot_hover", delayType="throttle", delay= 150),
                                          dblclick = "plotDblClick")
    ),
    miniButtonBlock(
      actionButton("zoomIn","", icon = icon("search-plus")),
      actionButton("zoomOut","", icon = icon("search-minus")),
      actionButton("screenShot","", icon = icon("camera-retro")),
      if( !is.null(colorVar) && !is.null(clusterFunc))
             actionButton("hint","Hint", icon = icon("map-signs")),
    )
  )

  #####################################################################
  server <- function(input, output, session) {
      # Create Status variables, helper data inside the server 
      highlightedIdx <- reactiveVal(-1)
      highlightedCat <- reactiveVal(-1)
      highlightedCatValue <- reactiveVal(value = NULL)
      helperValues <- reactiveValues()
      helperValues$clicked <- FALSE
      helperValues$movingDim <- -1
      helperValues$plotLimit <- 1.5
      helperValues$projectionMatrix <- initProjectionMatrix(df, projMatrix)
    
      helperValues$lbl <- NULL
      helperValues$cumulativeList <-  NULL
      helperValues$rangedData <- NULL
      helperValues$hints <- NULL


      observeEvent( input$plot_hover, {
           helperValues$curX <- input$plot_hover$x
           helperValues$curY <- input$plot_hover$y
           currentClosest <- closestDimension(helperValues$projectionMatrix, helperValues$curX, helperValues$curY )
           if (helperValues$clicked ){
                  # if we are moving the dimensional vectors, deactivate the brushing 
                  session$resetBrush("plotBrush")

                  helperValues$projectionMatrix[helperValues$movingDim,1:2] <- c(helperValues$curX, helperValues$curY)
                  if (approach == "OSC"){ # biggest difference is the recondition criterion when creating the star coordinates 
                       helperValues$projectionMatrix <- OSCRecondition(helperValues$projectionMatrix)
                  }
                  highlightedIdx(helperValues$movingDim)
           }
           else {
                highlightedIdx(currentClosest)
           }

           if (!is.null(highlightedCatValue())){
                helperValues$lbl <- NULL
           }

      })

      observeEvent(input$plotDblClick,{

          if(!numericRepresentation){
            # If categories are present in the representation, then we need to check 
            # for swapping or highlighting the name. 
            if ( highlightedCat() != -1){
              
                nameDim <- colnames(df)[highlightedCat()]
                inside <- insideCategoricalValueList(helperValues$projectionMatrix, nameDim, helperValues$cumulativeList, helperValues$curX, helperValues$curY )
                helperValues$lbl <-  insideCategoricalValueGG(helperValues$projectionMatrix, nameDim, helperValues$cumulativeList, inside)
                currentOrder <- helperValues$cumulativeList[[nameDim]]$order

                if (sum(inside) > 0){
                     selectedCatValue <-  currentOrder[which(inside)]
                      if(is.null( highlightedCatValue() )){
                           highlightedCatValue(selectedCatValue)
                      }
                      else if (!(selectedCatValue %in% highlightedCatValue()) ){
                           newOrder <- swap(currentOrder, selectedCatValue, highlightedCatValue())
                           helperValues$cumulativeList[[nameDim]] <- getCumulativeTable(frequencyList[[nameDim]], newOrder)
                           helperValues$rangedData <- getRangedData(df, numericRepresentation, meanCentered, frequencyList, helperValues$cumulativeList)
                           highlightedCatValue(NULL)
                           helperValues$lbl <-  NULL
                    }
                     else {
                           helperValues$lbl <-  NULL
                           highlightedCatValue(NULL)
                     }

                }
            }

              currentClosest <- closestDimension(helperValues$projectionMatrix, input$plotDblClick$x, input$plotDblClick$y )
              if ( currentClosest != -1){
                       nameDim <- colnames(df)[currentClosest]
                       if( is.factor(df[[nameDim]])){
                             if  ( highlightedCat() != -1){
                                    highlightedCat(-1)
                                    helperValues$lbl <- NULL
                                    highlightedCatValue(NULL)
                             }
                             else {
                                   highlightedCat(currentClosest)
                            }
                       }
              }
          }

      })

      mouseDown <- function(){
         # detect where the mouse was pressed, shiny doesn't divide click in 3 parts, so we have to use the javascript 
         # events 
         currentClosest <- closestDimension(helperValues$projectionMatrix, helperValues$curX, helperValues$curY )
         if ( currentClosest != -1){
               helperValues$movingDim <- currentClosest
               helperValues$clicked <- TRUE
         }
         helperValues$hints = NULL
      }

      mouseUp <- function(){
        helperValues$clicked <- FALSE
        helperValues$movingDim <- -1
      }
      

      initCumulative <- function(){
        # The placement of the points is dependent on the cumulative function in the case of categories
        # it is faster to have smaller tables that transform the frequency tables to cumulative and simply
        # change the cumulative sum afterwards
        if (!is.null( frequencyList) && is.null(helperValues$cumulativeList)){
          helperValues$cumulativeList <- list()
          catVars <- names(frequencyList)
          for( catName in catVars){
            cumulativeTable <- getCumulativeTable(frequencyList[[catName]])
            helperValues$cumulativeList[[catName]] <- cumulativeTable
          }
          helperValues$rangedData <- getRangedData(df, numericRepresentation, meanCentered, frequencyList, helperValues$cumulativeList)
        }
        
        if (approach == "OSC"){
          helperValues$projectionMatrix <- OSCRecondition(helperValues$projectionMatrix)
        }

        # if it's numerical, let re-range the values 
        if( numericRepresentation && is.null( helperValues$rangedData )){
             helperValues$rangedData <- getRangedData(df, numericRepresentation, meanCentered)
        }
      }


      output$plot <- renderPlot({
          # if we don't have the data scaled, then do the first call
          if (is.null( helperValues$rangedData )) initCumulative()

           # background
          curPlot <- drawDimensionVectors(helperValues$projectionMatrix, highlightedIdx(), highlightedCat(), highlightedCatValue(),  helperValues$cumulativeList )
          
          # hovering 
          if (!is.null(helperValues$lbl))  curPlot <- curPlot + helperValues$lbl
          
          # hints
          if (!is.null(helperValues$hints))   curPlot <- curPlot + helperValues$hints

          # actual points ... also add the limit so that the movement doesn't re-scales the image everytime 
          curPlot + getGGProjectedPoints(helperValues$projectionMatrix, helperValues$rangedData , odata, colorVar) +
          coord_cartesian(xlim = c(-helperValues$plotLimit,helperValues$plotLimit ), ylim = c(-helperValues$plotLimit,helperValues$plotLimit ))
      },res=96)

      observeEvent(input$cancel, {
        stopApp(NULL)
      })



      observeEvent(input$hint, {
             actions <- getActionHints(helperValues$projectionMatrix, helperValues$rangedData, odata, colorVar ,approach, clusterFunc)

             updateActionButton(session, "hint", label = paste0("Hint. Max(",formatC(max(actions$diff), format = "e",  digits = 2), ")"))
             if (max(actions$diff) > 0){
                  helperValues$hints <- getGGHints(helperValues$projectionMatrix, actions)

             }
             else {
                  helperValues$hints <- NULL
             }
      })

      # helping events ... zooming and out to an extent and create a screenshot of a current state
      onevent("mouseup", "plot", mouseUp())
      onevent("mousedown","plot", mouseDown())
      observeEvent(input$zoomIn,{ helperValues$plotLimit <-  helperValues$plotLimit*0.9 })
      observeEvent(input$zoomOut,{ helperValues$plotLimit <-  helperValues$plotLimit*1.1 })
      observeEvent(input$screenShot, { screenshot()})
      
      observeEvent(input$done, {
        projMatrix <- getCleanProjectionMatrix(helperValues$projectionMatrix, colorVar, odata, df)
        projectedPoints <- getProjectedPoints(helperValues$projectionMatrix, helperValues$rangedData )
        selected <- brushedPoints(df=projectedPoints, brush=input$plotBrush, xvar="V1", yvar="V2", allRows = TRUE)
        result <- list(projMatrix, selected$selected_, projectedPoints)
        names(result) <- c("Proj.Matrix","Selection","Projected.Points")

        stopApp(result)
      })
  }
   viewer <- paneViewer(300)
   runGadget(ui, server, viewer=viewer)
}



