

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

#' get projected points Star Coordinates 
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

#' get projected points RadViz
#' 
#' Project points given the projection matrix, simply matrix multiplication
#' @param anchorPoints anchor points to multiply with the ranged data
#' @param rangedData scaled and mapped data
#' @return the coordinates of the projected points 
#' @noRd
getProjectedPointsRadViz <- function(anchors, rangedData){
  
  A <- t(as.matrix(anchors[,1:2]))
  elemSum <- apply(rangedData,2, sum)
  pts <- A %*% rangedData
  pts <- pts / elemSum
  projectedPoints <- as.data.frame(t(pts))
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



