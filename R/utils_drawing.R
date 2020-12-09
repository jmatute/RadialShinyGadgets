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

#' get the ggplot2 points with computed projected points 
#' 
#' Given the projected points and create the geom_point based on the colouring attribute
#' @param projectedPoints projected points to create the gg drawing from 
#' @param oData original data, provides a way to color according to the original values
#' @param colorVar which attribute to be used to color the variables in the dataset
#' @return geom_point of the projected points
#' @noRd
drawGGProjectedPoints <- function(projectedPoints, oData,  colorVar=NULL ){
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



