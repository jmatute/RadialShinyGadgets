source("./R/utils.R")


#' Star Coordinates Gadget
#'
#' Creates a RShiny Gadget for Star Coordinates
#' @param df  A dataframe with the data to explore. It should contain only numeric or factor columns.
#' @param colorVar column where labels from the data are extracted.
#' @param approach Standard approach as defined by Kandogan, or Orthographic Star Coordinates (OSC) with a recondition as defined by Lehmann and Thiesel
#' @param numericRepresentation if true attempt to convert all factors to numeric representation, otherwise used mixed representation as defined in Hinted Star Coordinates
#' @param meanCentered center the projection at the mean of the values. May allow for easier value estimation
#' @param projMatrix a pre-defined projection matrix as an initial configuration. Should be defined in the same fashion as the output
#' @param clusterFunc function to define hints, assume increase in value of the function is an increase in quality of the projection. The function will be called with two parameters (points, labels)
#' @return A list with the projection matrix, coordinates of the projected samples and a logical vector with the selected samples
#' @export
#' @examples
#' library(datasets)
#' data(iris)
#' # StarCoordinates(iris, "Species")

StarCoordinates <- function(df, colorVar = NULL,  approach="Standard", numericRepresentation=TRUE, meanCentered = TRUE, projMatrix=NULL, clusterFunc = NULL) {

  #######################################################################
  # Simple Pre-processing, validate input to function and create data that
  # won't change throughout the running of the gadget
  inputCheck <- inputValidation(df,  colorVar,  approach, projectionMatrix, clusterFunc )
  if (!is.null(inputCheck)){
      if (inputCheck$stop)  stop(inputCheck$errorMessage)
      else warning(inputCheck$errorMessage)
  }
  odata <- df # Store original data somewhere
  df <- removeNonZeroAndMissing(df)

  frequencyList <- NULL
  if (!is.null(colorVar)  && (colorVar %in% colnames(df))){
      df[[colorVar]] <- NULL
  }
  else if (!is.null(colorVar)){
       colorVar <- NULL
  }

  if(numericRepresentation){
       df <-  mutate_if(df, is.factor, as.numeric)
  }
  else {
       frequencyList <- getFrequencyList(df)
  }
  ######################################################################
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
      if( !is.null(colorVar) && !is.null(clusterFunc))
             actionButton("hint","Hint", icon = icon("map-signs")),
    )
  )

  #####################################################################
  server <- function(input, output, session) {
      # Create Status variables
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
                  session$resetBrush("plotBrush")

                  helperValues$projectionMatrix[helperValues$movingDim,1:2] <- c(helperValues$curX, helperValues$curY)
                  if (approach == "OSC"){
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
                                    highlightedCat(-y1)
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
        if (!is.null( frequencyList) && is.null(helperValues$cumulativeList)){
          helperValues$cumulativeList <- list()
          catVars <- names(frequencyList)
          for( catName in catVars){
            cumulativeTable <- getCumulativeTable(frequencyList[[catName]])
            helperValues$cumulativeList[[catName]] <- cumulativeTable
          }
          helperValues$rangedData <- getRangedData(df, numericRepresentation, meanCentered, frequencyList, helperValues$cumulativeList)
        }

        if( numericRepresentation && is.null( helperValues$rangedData )){
             helperValues$rangedData <- getRangedData(df, numericRepresentation, meanCentered)
        }
      }


      output$plot <- renderPlot({
          initCumulative()
          curPlot <- drawDimensionVectors(helperValues$projectionMatrix, highlightedIdx(), highlightedCat(), highlightedCatValue(),  helperValues$cumulativeList )
          if (!is.null(helperValues$lbl))
                curPlot <- curPlot + helperValues$lbl
          if (!is.null(helperValues$hints))
                curPlot <- curPlot + helperValues$hints

          curPlot + getGGProjectedPoints(helperValues$projectionMatrix, helperValues$rangedData , odata, colorVar) +
          coord_cartesian(xlim = c(-helperValues$plotLimit,helperValues$plotLimit ), ylim = c(-helperValues$plotLimit,helperValues$plotLimit ))
      },res=96)

      observeEvent(input$cancel, {
        stopApp(NULL)
      })



      observeEvent(input$hint, {
             actions <- getActionHints(helperValues$projectionMatrix, helperValues$rangedData, odata, colorVar , clusterFunc)

             updateActionButton(session, "hint", label = paste0("Hint. Max(",formatC(max(actions$diff), format = "e",  digits = 2), ")"))
             if (max(actions$diff) > 0){
                  helperValues$hints <- getGGHints(helperValues$projectionMatrix, actions)

             }
             else {
                  helperValues$hints <- NULL
             }
      })

      onevent("mouseup", "plot", mouseUp())
      onevent("mousedown","plot", mouseDown())

      observeEvent(input$zoomIn,{ helperValues$plotLimit <-  helperValues$plotLimit*0.9 })
      observeEvent(input$zoomOut,{ helperValues$plotLimit <-  helperValues$plotLimit*1.1 })

      observeEvent(input$done, {
        projMatrix <- getCleanProjectionMatrix(helperValues$projectionMatrix, colorVar, odata, df)
        projectedPoints <- getProjectedPoints(helperValues$projectionMatrix, helperValues$rangedData )
        selected <- brushedPoints(projectedPoints, input$plotBrush, allRows = TRUE )
        result <- list(projMatrix, selected$selected_, projectedPoints)
        names(result) <- c("Proj.Matrix","Selection","Projected.Points")

        stopApp(result)
      })
   }
   runGadget(ui, server)
}



