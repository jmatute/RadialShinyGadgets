
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


#' Input Validation for Star Coordinates
#'
#' Check if the input is valid
#' @param data  A dataframe with the data to explore. It should contain only numeric or factor columns.
#' @param colorVar column where labels from the data are extracted.
#' @param approach Standard approach as defined by Kandogan, or Orthographic Star Coordinates (OSC) with a recondition as defined by Lehmann and Thiesel
#' @param projectionMatrix a pre-defined projection matrix as an initial configuration. Should be defined in the same fashion as the output
#' @param clusterFunc function to define hints, assume increase in value of the function is an increase in quality of the projection. The function will be called with two parameters (points, labels)
#' @return list of found errors, NULL if no error was found.
#' @noRd
inputValidationSC <- function(data, colorVar , approach, projectionMatrix, clusterFunc){
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


#' Input Validation for Rad Viz
#'
#' Check if the input is valid
#' @param data  A dataframe with the data to explore. It should contain only numeric  columns.
#' @return list of found errors, NULL if no error was found.
#' @noRd
inputValidationRadViz <- function(data, colorVar){
  # Simple input validation for the function. Properties needed are described below
  error <- NULL
  if( !is.data.frame(data)){
    error <- list("Data should be a data frame.", TRUE)
    names(error) <- c("errorMessage", "stop")
  }
  
  totalNumeric <- sum(sapply(data,  function(x) {is.numeric(x)}))
  if (!is.null(colorVar)){
     if(  is.factor( data[[colorVar]]))
        totalNumeric <- totalNumeric + 1
  }
  
  if (is.null(error) &&  !(totalNumeric == dim(data)[2]) ){
    error <- list("Data frame columns should be numeric only. ", TRUE)
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

  error
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


