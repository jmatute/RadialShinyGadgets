

## Star Coordinates

The following package is creates a RShiny Gadget for Star Coordinates which allows for initial data exploration. 


* **param** *df:*  A dataframe with the data to explore. It should contain only numeric or factor columns.
* **param** *colorVar:* column where labels from the data are extracted.
* **param** *approach:* Standard approach as defined by Kandogan, or Orthographic Star Coordinates (OSC) with a recondition as defined by Lehmann and Thiesel
* **param** *numericRepresentation:* if true attempt to convert all factors to numeric representation, otherwise used mixed representation as defined in Hinted Star Coordinates
* **param** *meanCentered:* center the projection at the mean of the values. May allow for easier value estimation
* **param** *projMatrix:* a pre-defined projection matrix as an initial configuration. Should be defined in the same fashion as the output
* **param** *clusterFunc:* function to define hints, assume increase in value of the function is an increase in quality of the projection. The function will be called with two parameters (points, labels)


* **returns** A list with the projection matrix, coordinates of the projected samples and a logical vector with the selected samples


### Theory 

Traditional Star Coordinates are defined for multidimensional numerical data sets  $X =\{ \mathbf{p}_{1}, \dots, \mathbf{p}_{d}  \}$ for $N$ data points $\mathbf{x}_{i} \in \mathbf{R}^{d}$ of dimensionality $d$. This approach was first proposed by Kandogan [1]. 


Let $A =\{ \mathbf{a}_{1}, \dots, \mathbf{a}_{d}  \}$ be a set of (typically 2D) vectors, each corresponding to one of the $d$ dimensions
.
The projection $\mathbf{p}_i' \in \mathbf{R}^{2}$ of a multidimensional point $\mathbf{p}_i = (p_{i1},\ldots,p_{id}) \in \mathbf{R}^{d}$ in SC is then defined as:
\begin{equation}
  \mathbf{x}_i' = \sum_{j=1}^{d} \mathbf{a}_{j} g_j( \mathbf{p}_i)
\end{equation}
with
\begin{equation} \label{eq:g}
 g_j(\mathbf{p}_i) = \frac{p_{ij} - min_j}{max_j - min_j}
\end{equation}
and $(min_j,max_j)$ denoting the value range of dimension $j$. 

### Star Coordinates for Mixed Data 

The Star Coordinates approach can be extended to handle mixed data i.e. numerical and categorical dimensions [4].  Equally spaced categorical points may not reflect the true nature of the data. Instead, frequency-based representation for individual data points can be used.  Assuming a categorical dimension $j$, we calculate the frequency $f_{jk}$ of each category $k$ of dimension $j$.  The respective axis vector $\mathbf{a}_{j}$  is then divided into according blocks, whose size represent the relative frequency (or probability) $\frac{f_{jk}}{\sum_{l=1}^m f_{jl}}$ of each of the $m$ categories of dimension $j$.

In summary, given an order for each categorical dimension, we can extend the equation above to SC for mixed data by:
\begin{equation} \label{eq:gcat}
 g_j(\mathbf{x}_i) = 
\begin{cases}
F_j(x_{ij}) - \frac{P_j(x_{ij})}{2}  & \text{if j is categorical/ordinal}, \\
\frac{x_{ij} - min_j}{max_j - min_j} & \text{if j is numerical}
\end{cases}
\end{equation}
where $F_j$ is the cumulative density function for (categorical/ordinal) dimension $j$ and $P_j$ its probability function. 


### References

1. Kandogan, E. (2001, August). Visualizing multi-dimensional clusters, trends, and outliers using star coordinates. In Proceedings of the seventh ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 107-116).

2. Lehmann, D. J., & Theisel, H. (2013). Orthographic star coordinates. IEEE Transactions on Visualization and Computer Graphics, 19(12), 2615-2624.

3. Rubio-Sánchez, M., & Sanchez, A. (2014). Axis calibration for improving data attribute estimation in star coordinates plots. IEEE transactions on visualization and computer graphics, 20(12), 2013-2022

4. Matute, J., & Linsen, L. (2020). Hinted Star Coordinates for Mixed Data. In Computer Graphics Forum (Vol. 39, No. 1, pp. 117-133).


### Examples

#### Traditional Star Coordinates 

    library(datasets)
    data(iris)
    StarCoordinates(iris)
    
![*Traditional Star Coordinates*](imgs/standard.gif)


#### Star Coordinates with Class Exploration

    StarCoordinates(iris, colorVar ="Species ")

![*Star Coordinates with a selected class*](imgs/species.gif)


#### Mean Centering

Mean centering is enabled by default. 

    StarCoordinates(iris, colorVar ="Species" , meanCentered = FALSE)
    
![*Star Coordinates Uncentered *](imgs/uncentered.gif)


#### Non-numeric Representation

By setting numeric representation to false, Star Coordinates uses the mixed data approach defined in [4] for the mapping.

    StarCoordinates(iris, numericRepresentation = FALSE )

![*Star Coordinates with Mixed Data *](imgs/mixed.gif)




   
  
