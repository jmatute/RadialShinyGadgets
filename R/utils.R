
#' Create Complete Case
#' 
#' How to handle missing data and zero variance attributes is undefined. As such, they are removed from the dataset and a warning is given to the user
#' @name  removeNonZeroAndMissing
#' @param data dataframe to check its state
#' @return A cleaned dataset from missing data and zero-variance attributes
#' @noRd
removeNonZeroAndMissing <- function(data){
  # The approach assumes a complete case scenario
  # Remove any incomplete cases.
  ok <- complete.cases(data)
  if (sum(!ok) > 0){
    warning("Incomplete cases are removed.")
  }
  data <- data[ok,]
  # The values are also scaled to 0..1 per dimension
  # if column contains near zero variance it is useless for the analysis
  nzv <- nearZeroVar(data)
  if (length(nzv) > 0){
    warning("Non-zero variance variables removed.")
    data <- data[, -nzv]
  }
  invisible(data)
}

#' Initialize Projection Matrix
#'
#' Initialize the projection matrix by dividing the unit circle
#' into locations, the size and shape of each point is defined.
#' @param data data to create the original projection atrmix from 
#' @param initProjMatrix optional initial projection matrix 
#' @return initial projection matrix
#' @noRd
initProjectionMatrix <- function(data, initProjMatrix = NULL) {
  # Interaction can change the properties of the points
  ncols <- dim(data)[2]
  initMatrix <- matrix(0L, nrow = 2, ncol = ncols)
  stepSize <- (2.0*pi)/(ncols)
  for(i in 1:ncols) initMatrix[,i] = c(cos((i-1)*stepSize), sin((i-1)*stepSize))
  szs <- rep.int(3, ncols)     # could be rewritten as initMatrix['szs'] <- 3
  shapes <- rep.int(21, ncols) # but this is more explicit
  initMatrix <- as.data.frame(cbind(t(initMatrix),szs, shapes))
  row.names(initMatrix) <- colnames(data)
  # If a projection matrix is given then we initialize the positions accordingly

  if ( !is.null(initProjMatrix)){
    for(name in rownames(initMatrix)){
      initMatrix[name,1:2] <- initProjMatrix[name,]
    }
  }
  invisible(initMatrix)
}


#' Input Validation
#'
#' Check if the input is valid
#' @param data  A dataframe with the data to explore. It should contain only numeric or factor columns.
#' @param colorVar column where labels from the data are extracted.
#' @param approach Standard approach as defined by Kandogan, or Orthographic Star Coordinates (OSC) with a recondition as defined by Lehmann and Thiesel
#' @param projectionMatrix a pre-defined projection matrix as an initial configuration. Should be defined in the same fashion as the output
#' @param clusterFunc function to define hints, assume increase in value of the function is an increase in quality of the projection. The function will be called with two parameters (points, labels)
#' @return list of found errors, NULL if no error was found.
#' @noRd
inputValidation <- function(data, colorVar , approach, projectionMatrix, clusterFunc){
  # Simple input validation for the function. Properties needed are described below
  error <- NULL
  if( !is.data.frame(data)){
    error <- list("Data should be a data frame.", TRUE)
    names(error) <- c("errorMessage", "stop")
  }

  if (is.null(error) && !(approach %in% c("Standard","OSC") )){
    error <- list("Selected approach was not found. Options are: Standard, OSC. ", TRUE)
    names(error) <- c("errorMessage", "stop")
  }

  if (is.null(error) &&  !(sum(sapply(data,  function(x) {is.numeric(x)|| is.factor(x)})) == dim(data)[2]) ){
    error <- list("Data frame columns should be numeric or factor. ", TRUE)
    names(error) <- c("errorMessage", "stop")
  }

  if (is.null(colorVar) && !is.null(clusterFunc) ){
    error <- list("Clustering function needs a labels dimension (colorVar). ", TRUE)
    names(error) <- c("errorMessage", "stop")
  }

  if (is.null(error) && dim(data)[2] <= 1){
    error <- list("Minimum 2 columns.", TRUE)
    names(error) <- c("errorMessage", "stop")
  }

  if ( !is.null(colorVar) && !(colorVar %in% colnames(data))){
    if (is.null(error)){
      error <- list("colorVar: Dimension not found for colouring. Defaulting to NULL. ", FALSE)
      names(error) <- c("errorMessage", "stop")
    }
    else {
      error$errorMessage = paste0( error$errorMessage,"\n","colorVar: Dimension not found for colouring. Defaulting to NULL. " )
    }
  }

  invisible(error)
}


######################################################################
#  Helper functions for drawing and interaction related functions


#' Circle Points
#' 
#' Create the points for the background 
#' @param center center of circle
#' @param radius radius of the background circle
#' @param npoints number of points to be created to define the circle
#' @return points that define the background circle 
#' @noRd
circleFun <- function(center = c(0,0), radius = 1, npoints = 100){
  # A simple generator of the coordinates for a circle that appears as 
  # background
  r = radius
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#' Is Point inside a categorical block
#' 
#' A helper function for checking per categorical block whether the click 
#' was performed inside or out the block
#' @param row a row from categorical value list (internal)
#' @return Boolean w.r.t. point exists within a categorical block
#' @noRd
isPointInside <- function(row){
  # Test whether the point is inside a categorical block
  # 0.03 is the width.
  p1 <- as.numeric(c(row[6], row[7]))
  p2 <- as.numeric(c(row[4], row[5]))
  x0 <-   as.numeric(c(row[8], row[9]))
  n <- (p2 - p1)/sqrt(sum((p2-p1)^2))
  temp <- (p1 -x0)  - (sum((p1-x0)*n))*n
  d <- sqrt(sum(temp^2))
  e1 <- p2 - p1
  recArea <- sum(e1^2)
  e2 <- x0 - p1
  val <- sum(e1*e2)
  invisible(d < 0.03 && val > 0 && val < recArea)
}

#' ggplot2 hints
#' 
#' Creation of the Hints as  ggplot2 segments per action (subsets)
#' @param projectionMatrix the current running projection matrix
#' @param possibleActions the hint increase value of the possible action subset
#' @return ggplot2 segments with different widths based on the increase in quality
#' @noRd
getGGHints <- function(projectionMatrix, possibleActions){
  # Creation of the Hints as  ggplot2 segments
  # per action, 2 points with xy coordinates are created
  # depending on the action different segments are done
  # ref. Hinted Star Coordinates for Mixed Data
  newProjectionMatrix <- projectionMatrix
  newProjectionMatrix["vector"] <- row.names(projectionMatrix)
  actions <- left_join(possibleActions, newProjectionMatrix, by="vector")

  newActions <- actions %>% rowwise() %>% mutate(stx=switch(as.character(.data$action), "inc" = .data$V1*0.5, "dec" = 0.0, "cw" = .data$V1, "ccw" = .data$V1)) %>%
    mutate(sty=switch(as.character(.data$action), "inc" = .data$V2*0.5, "dec" = 0.0, "cw" = .data$V2, "ccw" = .data$V2)) %>%
    mutate(etx=switch(as.character(.data$action), "inc" = .data$V1, "dec" = .data$V1*0.5,  "cw" = .data$x, "ccw" = .data$x)) %>%
    mutate(ety=switch(as.character(.data$action), "inc" = .data$V2, "dec" = .data$V2*0.5,  "cw" = .data$y, "ccw" = .data$y)) %>%
    mutate(width = 0.5 + .data$rangedDiff) # Change the sit

  segments <- geom_segment(aes(x = .data$stx, y = .data$sty, xend = .data$etx, yend = .data$ety, size=.data$width), lineend = "round", data = newActions, show.legend = FALSE)
  invisible(segments)
}


#' ggplot2 label inside categorical value
#' 
#' Create the label for the factor value in the selected factor dimension
#' @param projectionMatrix the current projection matrix
#' @param catName the name of the explored categorical dimension
#' @param cumulativeList a list with the cumulative values,to avoid  unnecessary computation
#' @param isInside for the category catName a logical vector whether is inside for each categorical value
#' @return geom_label is the element is inside a categorical value block
#' @noRd
insideCategoricalValueGG  <- function(projectionMatrix, catName, cumulativeList, isInside){

  dirVector <- projectionMatrix[catName,]
  curCumulativeTable <- cumulativeList[[catName]]$table
  element <- NULL
  if(sum(isInside) > 0){
    curCumulativeTable <- curCumulativeTable[isInside,]
    cumulativeName <- paste0(catName,".Cumulative")
    cumulativeFreq <- paste0(catName,".freq")
    curCumulativeTable  <- curCumulativeTable %>% mutate(bcx= (.data[[cumulativeName]] - .data[[cumulativeFreq]]*0.5) * dirVector$V1 ) %>%
                                                  mutate(bcy= (.data[[cumulativeName]] - .data[[cumulativeFreq]]*0.5) * dirVector$V2 )
    element <- geom_label( data=curCumulativeTable,  aes(.data$bcx, .data$bcy,label=.data[[catName]]))
  }
  invisible(element)
}

#' Closest Dimension to click
#' 
#' Given mouse coordinates, find the closest dimension vector
#' and return if it is under a certain distance
#' @param projectionMatrix current projection matrix
#' @param curX coordinate x of the current location in the plot
#' @param curY coordinate y of the current location in the plot
#' @param threshold distance threshold to be defined as close enough to a dimension
#' @return return index of the closest dimension, -1 if distance is over a threshold
#' @noRd
closestDimension <- function(projectionMatrix, curX, curY,threshold=0.05){

  distances <- apply(projectionMatrix, 1, function(x)  sqrt((x[1]-curX)*(x[1]-curX) + (x[2]-curY)*(x[2]-curY)))
  closest <- which.min(distances)
  if ( distances[closest] < threshold){
    closestDim <- closest
  }
  else  closestDim <- -1
  invisible(closestDim)
}


#' get the ggplot2 points 
#' 
#' Compute the projected points and create the geom_point based on the colouring attribute
#' @param projectionMatrix current projection matrix in order to project the data
#' @param rangedData scaled and mapped data
#' @param oData original data, provides a way to color according to the original values
#' @param colorVar which attribute to be used to color the variables in the dataset
#' @return geom_point of the projected points
#' @noRd
getGGProjectedPoints <- function(projectionMatrix, rangedData, oData,  colorVar=NULL ){
  
  projectedPoints <- getProjectedPoints(projectionMatrix, rangedData)
  drawnPoints <- NULL
  if (!is.null(colorVar) && !(colorVar %in% colnames(oData))){
    colorVar <- NULL
  }

  if (is.null(colorVar))
    drawnPoints <- geom_point(data=projectedPoints, aes(.data$V1,.data$V2),colour = "black", shape=20, size=2, alpha=0.3)
  else {
    drawnPoints <-  geom_point(data=projectedPoints, aes_string(x="V1",y="V2", colour = oData[[colorVar]]), shape=20, size=2, alpha=0.3)
  }
  invisible(drawnPoints)
}



######################################################################
#  Processing Functions


#' Get Value of Clustering Function to Optimize
#' 
#' Given a clustering validation function return the value when applied to the projected points
#' @param projectionMatrix current projection matrix in order to project the data
#' @param rangedData scaled and mapped data
#' @param oData original data, provides a way to color according to the original values
#' @param colorVar which attribute to be used to color the variables in the dataset
#' @param clusterFunc function to define hints, assume increase in value of the function is an increase in quality of the projection.
#'                   The function will be called with two parameters (points, labels)
#' @return value given the cluster validity function 
#' @noRd
getValueAtFunc <- function(projectionMatrix, rangedData, oData, colorVar , clusterFunc){
  projectedPoints <- getProjectedPoints(projectionMatrix, rangedData )
  labels <- as.numeric(oData[[colorVar]])
  value <- clusterFunc(projectedPoints, labels)
  invisible(value)
}


#' Get Value of Clustering Function to all actions
#' 
#' According to Hinted Star Coordinates Paper the hints can be computed with 20% change in angle or size
#' Per dimension, define where the new vector should end based on the four actions: increase(inc), decrease(dec)
#' clockwise rotation (cw) and counter-clockwise rotation(ccw). The total hints to be evaluated are then 4xdimensions
#' @param projectionMatrix current projection matrix in order to project the data
#' @param rangedData scaled and mapped data
#' @param oData original data, provides a way to color according to the original values
#' @param colorVar which attribute to be used to color the variables in the dataset
#' @param approach approach used for star coordinates 
#' @param clusterFunc function to define hints, assume increase in value of the function is an increase in quality of the projection.
#'                   The function will be called with two parameters (points, labels)
#' @return a dataframe containing the change in value due to the execution of the defined actions 
#' @noRd
getActionHints <- function(projectionMatrix, rangedData, oData, colorVar , approach, clusterFunc){
  
  increase <- function(x){x*1.20}
  decrease <- function(x){x*0.80}
  clockwise <- function(x){ x  - 0.20*3.14159*0.5 }
  counterClockwise <- function(x){ x  + 0.20*3.14159*0.5 }

  changesInPts <- projectionMatrix %>% mutate_at(c("V1","V2"), list(inc=increase, dec=decrease)) %>%
    mutate(angle= atan2(.data$V2,.data$V1)) %>% mutate_at( c("angle"), list(cw=clockwise, ccw = counterClockwise)) %>%
    mutate(V1_cw =  cos(.data$cw)*sqrt(.data$V1^2 + .data$V2^2)  )  %>%   mutate(V2_cw =  sin(.data$cw)*sqrt(.data$V1^2 + .data$V2^2)  )  %>%
    mutate(V1_ccw =  cos(.data$ccw)*sqrt(.data$V1^2 + .data$V2^2)  )  %>%   mutate(V2_ccw =  sin(.data$ccw)*sqrt(.data$V1^2 + .data$V2^2)  )

  row.names(changesInPts) <- row.names(projectionMatrix)

  actions <- c("inc","dec","cw","ccw")
  vectors <- row.names(changesInPts)

  origVal <- getValueAtFunc(projectionMatrix, rangedData , oData, colorVar , clusterFunc)
  # Create a data frame where each action per vector is a row
  possibleActions <- data.frame(actions=c(outer(vectors, actions, FUN=paste)) ) %>% tidyr::separate(actions, into=c("vector","action"),sep=" ")
  possibleActions["val"] <- NA
  possibleActions["x"] <- NA
  possibleActions["y"] <- NA


  for ( rowIdx in 1:dim( possibleActions)[1]){
    newProjection <- projectionMatrix
    vec <- possibleActions[rowIdx,]$vector
    action <- possibleActions[rowIdx,]$action
    newProjection[vec,"V1"] <-  changesInPts[vec,paste0("V1_",action)]
    newProjection[vec,"V2"] <-  changesInPts[vec,paste0("V2_",action)]
    
    if ( approach == "OSC"){
        newProjection <-  OSCRecondition(newProjection)
    }
    possibleActions[rowIdx,]$x <-  changesInPts[vec,paste0("V1_",action)]
    possibleActions[rowIdx,]$y <-  changesInPts[vec,paste0("V2_",action)]
    
    # Get the value with the new projection matrix
    possibleActions[rowIdx,]$val <- getValueAtFunc(newProjection, rangedData , oData, colorVar , clusterFunc)
  }
  # The clustering function assumes that a larger value is better.
  possibleActions <- possibleActions %>% mutate(diff= ifelse(.data$val - origVal>0, .data$val - origVal,0)) %>% mutate(rangedDiff= diff/max(diff))
  invisible(possibleActions)
}


#' Check whether is inside categorical value 
#' 
#' Check whether the point selected is inside a categorical value, works together with
#' insideCategoricalValueGG to create the label of the selected categorical value
#' @param projectionMatrix the current projection matrix
#' @param catName the name of the explored categorical dimension
#' @param cumulativeList a list with the cumulative values,to avoid  unnecessary computation
#' @param px mouse click location on x coordinates
#' @param  py mouse click location on y coordinates
#' @return logical vector whether the clicked element is inside a categorical value block
#' @noRd
insideCategoricalValueList <- function(projectionMatrix, catName, cumulativeList, px, py ){
  # Given a factor attribute that is drawn. Check whether the interaction
  # points to a specific categorical value block
  origin <- data.frame(x = 0, y = 0)
  lines <- do.call(rbind, apply(projectionMatrix, 1, function(x) {rbind(x,  origin )}))
  dirVector <- lines[paste0(catName,".1"),]
  cumulativeName <- paste0(catName,".Cumulative")
  cumulativeFreq <- paste0(catName,".freq")

  ndirVector <- dirVector / sqrt(sum(dirVector^2))
  ortho <- c(-ndirVector$y, ndirVector$x)
  curCumulativeTable <- cumulativeList[[catName]]$table
  curCumulativeTable  <- curCumulativeTable %>% mutate(cx= .data[[cumulativeName]] * dirVector$x ) %>% mutate(cy= .data[[cumulativeName]] * dirVector$y )
  curCumulativeTable  <- curCumulativeTable %>% mutate(bcx= (.data[[cumulativeName]] - .data[[cumulativeFreq]]) * dirVector$x ) %>% mutate(bcy= (.data[[cumulativeName]] - .data[[cumulativeFreq]]) * dirVector$y )
  curCumulativeTable  <- curCumulativeTable %>% mutate(x = px) %>% mutate(y = py)
  isInside <- apply(curCumulativeTable, 1, function(row){ isPointInside(row) })

  invisible(isInside)
}


#' Swap elements in a vector
#' 
#' Swap two elements in a vector
#' @param vector where the elements exist
#' @param elem1 first element to swap
#' @param elem2 second element to swap
#' @return vector with swapped elements 
#' @noRd
swap <- function(vector,  elem1, elem2){
  # Simple swap of two elements in a vector
  p1 <- match(elem1, vector)
  p2 <- match(elem2, vector)
  vector[p1] <- elem2
  vector[p2] <- elem1
  invisible(vector)
}

#' Insert a row into an existing dataframe
#' 
#' Insert a row into an existing dataframe
#' @param existingDF existing dataframe to insert row to
#' @param newrow row to be inserted
#' @param r index to insert the new row
#' @return dataframe with the new row at index
#' @noRd
insertRow <- function(existingDF, newrow, r) {
  # Insert a row in a existing data frame
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

#' Orthographic Star Coordinates Recondition
#' 
#' Reconditioning based on the Orthographic Star Coordinates
#' Maintains the orthogonality of the projection.
#' @param projectionMatrix current projection matrix
#' @return new projection matrix where the orthogonality is maintained
#' @noRd
OSCRecondition <- function(projectionMatrix){
  normX <- norm(projectionMatrix$V1, type="2")
  dotP  <- sum((projectionMatrix$V1/normX)*projectionMatrix$V2)
  newY  <- projectionMatrix$V2 - dotP* (projectionMatrix$V1/normX)
  projectionMatrix$V1 <- projectionMatrix$V1/normX
  projectionMatrix$V2 <- newY/norm(newY, type="2")
  invisible(projectionMatrix)
}

#' get projected points
#' 
#' Project points given the projection matrix, simply matrix multiplication
#' @param projectionMatrix projection matrix to multiply with the ranged data
#' @param rangedData scaled and mapped data
#' @return the coordinates of the projected points 
#' @noRd
getProjectedPoints <- function(projectionMatrix, rangedData){
  # 
  A <- t(as.matrix(projectionMatrix[,1:2]))
  projectedPoints <- as.data.frame(t(A %*% rangedData))
  invisible(projectedPoints)
}


#############################################################
# Helper functions for frequency related analysis i.e. numeric representation = FALSE

#' Create a frequency table
#' 
#' Given a category define the frequency i.e. % of occurences
#' @param df dataframe as basis
#' @param varName factor dimension to create the frequency
#' @return frequency table of a categorical dimension
#' @noRd
getFrequencyTable <- function(df, varName){
  freqName <- paste0(varName,".freq")
  freqTable <- df %>% group_by(.data[[varName]]) %>% tally()  %>% mutate(!!freqName := n / sum(n)) %>% select(c(varName,freqName))
  invisible(freqTable)
}

#' Create a cumulative table
#' 
#' Given a frequency table create the cumulative function return a list
#' with the order used for the cumulative function and the cumulative table result
#' @param frequencyTable result from getFrequencyTable 
#' @param order order of the categorical value to be summed over 
#' @return list with the table and order of the cumulative results
#' @noRd
getCumulativeTable <- function(frequencyTable, order=NULL){
  if (!is.null(order)){
    indices <- match(order, frequencyTable[[colnames(frequencyTable)[1] ]])
    frequencyTable <- frequencyTable[indices,]
  }
  else {
    order <- as.character(frequencyTable[[ colnames(frequencyTable)[1] ]])
  }
  cumulative <- cbind(frequencyTable, cumsum(frequencyTable[,2]))
  colnames(cumulative)[3] <- paste0(colnames(cumulative)[1], ".Cumulative")

  tableWithOrder <- list(cumulative, order)
  names(tableWithOrder) <- c("table","order")
  invisible(tableWithOrder)
}

#' Get Frequency List
#' 
#' Per factor dimension create a frequency table and return as list
#' @param data dataset to check for factor attributes 
#' @return list of frequencies for all categories
#' @noRd
getFrequencyList <- function(data){
  l = list()
  for( name in colnames(data)){
    if ( is.factor( data[[name]] )){
      l[[name]] <- getFrequencyTable(data, name)
    }
  }
  invisible(l)
}

#' Draw Dimension Vectors
#' 
#' Draw the background and vectors to interact in the Star Coordinates approach.
#' Modify the vectors according to different interaction
#' @param projectionMatrix current projection matrix
#' @param highlightedIdx if mouse is over a dimension vector, highlight the dimension vector
#' @param highlightedCat a highlighted category displays the frequency blocks and allows interaction between blocks
#' @param highlightedCatValue which element within a categorical dimension has been selected
#' @param cumulativeList cumulative table for all categorical dimensions
#' @return ggplot with the dimension vectors, background and possible interactions with dimensional vectors 
#' @noRd
drawDimensionVectors <- function(projectionMatrix, highlightedIdx, highlightedCat, highlightedCatValue,  cumulativeList){
   # Create the dimensional vectors and any visual modifications to the interaction area
   sizes <- projectionMatrix$szs
   shapes <- projectionMatrix$shapes
   circle <- circleFun()

   if ( highlightedIdx != -1) {
         sizes[highlightedIdx] = 6
         shapes[highlightedIdx] = 20
   }
   if (highlightedCat != -1){
        shapes[highlightedCat] = 20
   }
   # Draw the labels + the point that can be selected to move the vectors
   pts <- geom_point(data=projectionMatrix, aes(.data$V1,.data$V2), colour = "black", fill = "white", shape=shapes, size=sizes)
   labels <-  geom_label( data=projectionMatrix,  aes(.data$V1+.data$V1*0.2,.data$V2+.data$V2*0.2,label=rownames(projectionMatrix)))
   #draw the background circle
   backgroundCircle <- geom_path(data=circle, aes(.data$x, .data$y), alpha=0.2)
   # Create a path for the vectors, by simply adding (0,0) in between the locations
   origin <- data.frame(x = 0, y = 0)
   lines <- do.call(rbind, apply(projectionMatrix, 1, function(x) {rbind(x,  origin )}))
   vectors <- geom_path(data=lines, aes(.data$x,.data$y))
   ## if it is non-numeric representation

   curPlot <- ggplot() + pts+ labels +vectors+ backgroundCircle
   if( length(cumulativeList) > 0){
        catVars <- names(cumulativeList)

        for( catName in catVars){
          dirVector <- lines[paste0(catName,".1"),]
          cumulativeName <- paste0(catName,".Cumulative")
          cumulativeFreq <- paste0(catName,".freq")

          ndirVector <- dirVector / sqrt(sum(dirVector^2))
          ortho <- c(-ndirVector$y, ndirVector$x)
          curCumulativeTable <- cumulativeList[[catName]]$table

          curCumulativeTable  <- curCumulativeTable %>% mutate(cx= .data[[cumulativeName]] * dirVector$x ) %>% mutate(cy= .data[[cumulativeName]] * dirVector$y )
          curCumulativeTable  <- curCumulativeTable %>% mutate(tsx = .data$cx  +  0.02 * ortho[1]) %>% mutate(tex = .data$cx +  -0.02 * ortho[1])
          curCumulativeTable  <- curCumulativeTable %>% mutate(tsy = .data$cy  +  0.02 * ortho[2]) %>% mutate(tey = .data$cy +  -0.02 * ortho[2])

          drawCat <-  highlightedIdx != -1 && catName %in%  rownames(projectionMatrix)[highlightedIdx]
          drawCat <- drawCat || (  highlightedCat != -1 && catName %in%  rownames(projectionMatrix)[highlightedCat] )

          if ( drawCat  ){

            curCumulativeTable  <- curCumulativeTable %>% mutate(bcx= (.data[[cumulativeName]] - .data[[cumulativeFreq]]) * dirVector$x ) %>% mutate(bcy= (.data[[cumulativeName]] - .data[[cumulativeFreq]]) * dirVector$y )
            curCumulativeTable  <- curCumulativeTable %>% mutate(bsx = .data$bcx  +  0.02 * ortho[1]) %>% mutate(bex = .data$bcx +  -0.02 * ortho[1])
            curCumulativeTable  <- curCumulativeTable %>% mutate(bsy = .data$bcy  +  0.02 * ortho[2]) %>% mutate(bey = .data$bcy +  -0.02 * ortho[2])

            curCumulativeTable  <- curCumulativeTable %>%  tidyr::gather( "orig_x", "x" , c(.data$tsx, .data$tex, .data$bex, .data$bsx)) %>%
                                                                  gather( "orig_y", "y" , c(.data$tsy, .data$tey, .data$bey, .data$bsy)) %>% 
                                                           mutate( orig_x=substr(.data$orig_x,1,2)) %>% mutate( orig_y=substr(.data$orig_y,1,2)) %>% filter(.data$orig_x == .data$orig_y)

            

            if ( !is.null(highlightedCatValue)){
                 curCumulativeTable <- curCumulativeTable %>% mutate(fillC = ifelse( .data[[catName]] == highlightedCatValue, "1", "0")  )
                 curPlot <- curPlot + geom_polygon(data=curCumulativeTable, colour="black",  alpha=0.4, show.legend=FALSE,mapping=aes(x=.data$x, y=.data$y,  group=.data[[catName]], fill=.data$fillC))
            }
            else
                curPlot <- curPlot + geom_polygon(data=curCumulativeTable, colour="black",  fill="white", alpha=0.4, mapping=aes(x=.data$x, y=.data$y,  group=.data[[catName]]))
          }
          else {
               curPlot <- curPlot + geom_segment(aes(x = .data$tsx, y = .data$tsy, xend = .data$tex, yend = .data$tey), data = curCumulativeTable)
          }
        } # end of checking the list of categorical
   }
   ##
   invisible( curPlot + theme_void())
}

#' Get Ranged Data
#' 
#' In order to multiply by the projection matrix the values should be scaled
#' if categorical, they should have a numerical value to be used.
#' @param data original dataset to scale the values
#' @param numericRepresentation assume that all attributes are numerical
#' @param meanCentered application of axis calibration for  value estimation 
#' @param frequencyList result from getFrequencyList (not used if numerical representation)
#' @param cumulativeList result from calling getCumulativeTable per category (not used if numerical representation)
#' @return a dataframe of the scaled and mapped dataset
#' @noRd
getRangedData <- function(data, numericRepresentation=TRUE, meanCentered = TRUE, frequencyList = NULL, cumulativeList =NULL){
  # Modify the original data so it can be used for matrix multiplication

  rangedData <- NULL

  if (numericRepresentation){
        if ( meanCentered){
          rangedData <-  as.matrix(t(as.data.frame(apply(data, 2, function(x){
            r <- max( max(x) - mean(x), mean(x) - min(x))
            (x - mean(x))/r
          }))))
        }
        else
          rangedData <-  as.matrix(t(as.data.frame(apply(data, 2, function(x){(x-min(x))/(max(x)-min(x))}))))
  }
  else {
        numericCols <- sapply(data,  function(x) {is.numeric(x)})
        if ( meanCentered){
             data[,numericCols] <- apply(data[,numericCols], 2,   function(x){
               r <- max( max(x) - mean(x), mean(x) - min(x))
               (x - mean(x))/r
             })
             rangedData <- as.matrix(t(as.data.frame(data)))
        }
        else {
          data[,numericCols] <- apply(data[,numericCols], 2,   function(x){
            (x-min(x))/(max(x)-min(x))
          })
        }

          if (!is.null( frequencyList)){
                 catVars <- names(frequencyList)
                 for( catName in catVars){
                       cumulativeTable <- cumulativeList[[catName]]
                       data <- left_join(data, cumulativeTable$table, by=catName)
                       data[[catName]] <- data[[paste0(catName,".Cumulative")]]  -  data[[paste0(catName,".freq")]]*0.5
                       data[[paste0(catName,".Cumulative")]] <- NULL
                       data[[paste0(catName,".freq")]] <- NULL
                 }
          }
          rangedData <-  as.matrix(t(as.data.frame(data)))
  }

  invisible(rangedData)
}

#####################################################################
#   Output Functions


#' Get Clean Projection Matrix
#' 
#' Get the projection matrix as a result without the sizes and shapes of points
#' If a color Variable was selected then insert into the projection matrix at the
#' same location where it should have been based on the dimensions
#' @param projectionMatrix projection matrix to clean
#' @param colorVar attribute used for colouring
#' @param odata original dataset without any modifications
#' @param data data post-preprocessing 
#' @return cleaned projection matrix 
#' @noRd
getCleanProjectionMatrix <- function(projectionMatrix, colorVar, odata, data){

  projMatrix  <- projectionMatrix[,1:2]
  if ( !is.null(colorVar)  && (colorVar %in% colnames(odata))){
    originalIdx <- which( colnames(odata) == colorVar)
    if ( originalIdx == length(colnames(odata)))
      projMatrix  <- rbind( projMatrix , rep.int(0,2))
    else if ( originalIdx == 1)
      projMatrix  <- rbind( rep.int(0,2),  projMatrix )
    else
      projMatrix <- rbind( projMatrix[1:(originalIdx-1),], rep.int(0,2), projMatrix[-(1:(originalIdx-1)),])
  }
  rownames( projMatrix) <- colnames(odata)
  colnames( projMatrix) <- c("x","y")
  invisible(projMatrix)
}


