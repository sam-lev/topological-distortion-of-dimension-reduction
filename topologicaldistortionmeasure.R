#
# Optimizing dimension reduction through projection of high dimension data 
# onto low dimensional subspace where quality of dimension reduction is
# determined based of off minimal distortion.
#
# Author: Samuel Leventhal
# Email: samlev@cs.utah.edu

#
# Projection of high dimensional vector b onto subspace of dimension A
#


dir.create(paste(getwd(),"/componentpca",sep=""))
dir.create(paste(getwd(),"/cylinder",sep=""))
dir.create(paste(getwd(),"/pointclouds",sep=""))
dir.create(paste(getwd(),"/randomprojdistortion",sep=""))
dir.create(paste(getwd(),"/RingFigures",sep=""))
dir.create(paste(getwd(),"/ShortLoop",sep=""))
dir.create(paste(getwd(),"/swissroll",sep=""))
dir.create(paste(getwd(),"/tetrahedron",sep=""))
dir.create(paste(getwd(),"/torus",sep=""))
dir.create(paste(getwd(),"/wasserstein",sep=""))
#
# Returns a list of the higher dimensional point cloud pointcloud 
# projected onto the subspace defined by the principle components
# for each 1-cycle component identified by ShortLoop. The list of
# vertices defining each loop component can be passed in the variable
# components if already loaded from the pre-built *_loop.OFF file. If
# components is not provided pcacomp will call componentpoints to
# load the list of components from the pre-built .OFF file. The
# variable plot plots each the projections of pointcloud onto each
# subspace defined by 1-cycle vertices principle components.
# Input: 
#   components: (optional) List of matrices where each point
#     is a vertex defining a loop feature identified by ShortLoop
#   pointcloud: Matrix of points defining a 3-dimensional pointcloud
#   plot: (optional) Boolean, False by default, provides plots of projected
#     point cloud onto principle components of 1-cycles. All figures are
#     saved in componentpca folder
#
#
pcacomp <-
  function(components, pointcloud, plot, plotname){
    
    if(missing(plot)){
      plot <- FALSE
    }
    
    if(missing(plotname)){
      plotname <- "unnamedpointcloud"
    }
    
    
    if (!require(package = "stats")) {
      install.packages("stats")
    }
    library("stats")
    
    #
    # If the pointcloud is the name of a pre-built pointcloud 
    # then load pointcloud.
    #
    if(typeof(pointcloud) == "character"){
      plotname <- pointcloud
      pointcloud <- as.matrix(readpointcloud(paste(pointcloud,".txt",sep="")))
    }
    
    #
    # If list of vertex points defining 1-homological features is missing
    # call componentpoints to load list of vertex points for each 
    # component from pre-made *_loop.OFF file
    #
    if(missing(components)){
      components <- componentpoints(pointcloud, paste(pointcloud,".txt",sep=""), FALSE)
    }
    #
    # Perform Principele component ananlysis
    #
    pcapoints<- list()
    pcarotations <- list()
    
    for( i in 1:length(components)){
      pca <- prcomp(components[[i]], center = TRUE, scale. = FALSE)
      pcapoints[[i]] <- pca
      pcarotations[[i]] <- pca$rotation
    }
    
    #
    # Now project the entire point cloud in the directions optimal for 
    # each component
    #
    proj_pc <- list()
    for( i in 1:length(pcarotations)){
      rotatedpc <- as.data.frame(pointcloud%*%pcarotations[[i]])
      proj_pc[[i]] <- rotatedpc  
    }
    
    if(plot){
      
      for (i in 1:length(proj_pc)){
          pdf(paste("./componentpca/",plotname,"pointcloudprojection",toString(i),"pca", toString(i),".pdf", sep=""),"pdf")
          plot(proj_pc[[i]])
          dev.off()
        
      }

    sink()
    sink(file = NULL)
    }
    
    return(proj_pc)
  }


#
# Calculate the topological distortion between higher and lower
# dimensional point cloud where the lower dimensional point cloud
# is the result of PCA on loop components of point cloud. The 
# dimension reduced data (lower dimensional projection based on PCA)
# is returned from the function pcacomp. The function computes the
# persistence diagrams of the point cloud and all PCA reduced 
# point clouds of interest, stores barcode diagrams based on these
# filtration diagrams, and computes the bottleck and p-wasserstein
# distance which can be written to file 
#
# Input:
#   plot: Boolean (optional) false by default otherwise records data 
#     and figures in filtrationlooppca folder
# 
# Output: 
#

#   Example call to function using functions within code: 
#
#       # load point cloud from file
#       pointcloud_matrix <- as.matrix(readpointcloud(paste(pointcloud,".txt",sep="")))
#
#       # identify loop features found by shortloop. Requires shortloop _loop.off file and 
#       #  stored point cloud matrix in .txt file (can use pointcloudtotext() ) 
#       #  pointcloud is name of .off file such as pointcloud.off
#       components <- componentpoints(pointcloud, paste(pointcloud,".txt",sep=""), FALSE)
#
#       # obtain list of point cloud projects based on PCA from interior loops
#       loopPCA <- pcacomp(components, pointcloud, FALSE)
#
#       # call distortionloopPCA()
#       distortionloopPCA(pointcloud, loopPCA, numfirstprin, numsecondprin)

distortionloopPCA <-
  function(pointcloud, loopPCA, firstprincomp, secondprincomp, plot, plotname)
  {

    if(missing(plot)){
      plot <- FALSE
    }
    if(missing(plotname)){
      plotname <- "unnamedpointcloud"
    }
    
    if (!require(package = "TDA")) {
      install.packages("TDA")
    }
    library("TDA")
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    if(!require(package = "stats")){
      install.packages("stats")
    }
    library("stats")
    

    #
    # If pointcloud already stored, load pointcloud matrix, parse components
    # from pre-compiled  .OFF file into list of matrices whose points are 
    # vertices to each 1-cycle in point cloud, and obtain projections of 
    # point cloud onto principle component subspaces defined by PCA on 
    # those 1-cycles.
    #
    if(typeof(pointcloud) == "character"){
      plotname <- plotname
      pointcloud_matrix <- as.matrix(readpointcloud(paste(pointcloud,".txt",sep="")))
      components <- componentpoints(pointcloud, paste(pointcloud,".txt",sep=""), FALSE)
      loopPCA <- pcacomp(components ,pointcloud,FALSE)
    }else{ 
      pointcloud_matrix <- pointcloud
    }    

    if(missing(loopPCA)){
      loopPCA <-  pcacomp(components , pointcloud_matrix, FALSE)
    }
    
    # 
    #Perform the ripsfiltration of the d dimensional point cloud PC, obtain rips diagram.
    #
    DiagPC <- ripsDiag(X = as.data.frame(pointcloud_matrix), maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram

    
    # 
    #Perform the ripsfiltration of the d dimensional point cloud PC, obtain rips diagram.
    #
    DiagPCproj <- list()
    for( i in 1:length(loopPCA)){
      
    DiagPCproj[[i]] <- ripsDiag(X = as.data.frame(cbind(loopPCA[[i]][,firstprincomp],loopPCA[[i]][,secondprincomp])), maxdimension = 1, maxscale = 1, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    
    }

    #
    # Calculate bottleneck and p-wasserstein distances between 
    # persistence diagrams of entire point cloud and pca projected 
    # point cloud
    #
    for( i in 1:length(DiagPCproj)){
      bDistance = bottleneck(DiagPC, DiagPCproj[[i]], dimension=1)
      w1Distance = wasserstein(DiagPC, DiagPCproj[[i]], p=1, dimension=1)
      w2Distance = wasserstein(DiagPC, DiagPCproj[[i]], p=2, dimension=1)
      w3Distance = wasserstein(DiagPC, DiagPCproj[[i]], p=3, dimension=1)
      w4Distance = wasserstein(DiagPC, DiagPCproj[[i]], p=4, dimension=1)
      w5Distance = wasserstein(DiagPC, DiagPCproj[[i]], p=5, dimension=1)
    }
    
    #
    # Plot results
    #
    if(TRUE){
      sink(file=NULL)
      # 
      # Plot Persistence Diagram of Point Cloud
      #
      pdf(paste("./filtrationloopPCA/",plotname,"pointcloudbarcode.pdf",sep=""),"pdf")
      plot(DiagPC, barcode = TRUE)
      dev.off()
      par(mar=c(1,1,1,1))
      
      # 
      # Plot Persistence Diagram of Point Cloud
      #
      for( i in 1:length(DiagPCproj)){
        pdf(paste("./filtrationloopPCA/",plotname,"projonloopPCA",toString(i),".pdf",sep=""),"pdf")
        plot(DiagPCproj[[i]], barcode = TRUE)
        dev.off()
      }
      
      ProjPCAResults = paste("./filtrationloopPCA/",plotname,"projloopcomponentsPCAresults.txt",sep="")
      sink(ProjPCAResults)
    
      cat(paste("Pointcloud = ",plotname,"\n",sep=""),append = TRUE)
    
      for(i in 1:length(DiagPCproj)){
      cat(paste("Number ", toString(i), " of ", toString(length(DiagPCproj))
                , " distance measures between point cloud and projection based on 1-cylce connected component \n \n"), sep="")
      cat("bottleneck distance:\n", append=TRUE)
      cat(paste(toString(bDistance), "\n", sep=""))
      cat("\n")
    
      cat("1-Wasserstein Distance:\n", append = TRUE)
      cat(paste( "  ",toString(w1Distance), "\n", sep=""))
      cat("\n")
      cat("2-Wasserstein Distance:\n", append = TRUE)
      cat(paste( "  ",toString(w2Distance), "\n", sep=""))
      cat("\n")
      cat("3-Wasserstein Distance:\n", append = TRUE)
      cat(paste( "  ",toString(w3Distance), "\n", sep=""))
      cat("\n")
      cat("4-Wasserstein Distance:\n", append = TRUE)
      cat(paste("  ", toString(w4Distance), "\n", sep=""))
      cat('\n')
      cat("5-Wasserstein Distance:\n", append = TRUE)
      cat(paste( "  ",toString(w5Distance), "\n", sep=""))
      cat("\n")
      cat("\n")
      }
      sink()
      sink()
      sink(file=NULL)
    }
    
    return(list(bDistance, w1Distance,w2Distance,w3Distance,w4Distance, w5Distance))

  }

#
# Helper function for randprojmindistortion. Projects matrix A theta along
# the x component, phi along the y, and psi along the z. Returns the 
# rotated pointcloud.
#
vecnorm <- function(x) {sqrt(sum(x^2))}
proj <-
  function( A , theta, phi, psi ){
    
    A <- as.data.frame(A)
    Arot <- matrix(NA, nrow = nrow(A), ncol= ncol(A))
    rotatexvec_mat <- matrix(c(1, 0, 0, 0, cos(theta), -1*sin(theta)
                               , 0, sin(theta), cos(theta)), nrow=3, ncol=3)
    rotateyvec_mat <- matrix(c(cos(phi), 0, sin(phi), 0, 1, 0
                               , -1*sin(phi), 0, cos(phi)), nrow=3, ncol=3)
    rotatezvec_mat <- matrix(c(cos(psi), -1*sin(psi), 0, sin(psi), cos(psi), 0
                               ,0,0,1), nrow=3, ncol=3)
    for(row in 1:nrow(A)){
      Arot[row,] <- rotatexvec_mat%*%rotateyvec_mat%*%rotatezvec_mat%*%t(as.matrix(A[row,])) 
    }
    
    return(Arot)
  }

#
# Helper Function for randprojmindistortion: Uses the randomly generated angles
# in radians to project the point cloud under question in those directions.
# Input: matrix consisting of points in point cloud, list of 3D directions.
# Output: list of projected point clouds as data frame.
#
#
projmatrix <-
  function(datamatrix , Angles){
    
    #
    # Project point cloud onto m random subspaces of dimension n
    # p = projected point cloud
    Plist <- list(proj(datamatrix, pi/2, 0, 0) )
    #P <- matrix( data = NA, nrow = nrow(as.data.frame(datamatrix)), ncol = ncol(as.data.frame(datamatrix)) )
    
    for(i in 2:nrow(Angles)){
      P<-proj(datamatrix, Angles[i,1], Angles[i,2], Angles[i,3]) # orPCMatrix?
      Plist[[i]] <- P
    }
    
    return(Plist)
  }

#
# Calculate the projection with the minimal bottelneck, 1-Wasserstein, and
# 2-Wasserstein distance between high dimensional point cloud persistence
# and 'numDirections' point clouds projected onto the xy plane after rotated
# by (theta ,phi,psi) from x,y,z axis respectively.
# Input:
#   PC: high dimensional point cloud data
#   numDirections: number directions to rotate point cloud
#   plot: boolean, Default False-no plot or recording min topological distances
#         all plots and data are saved in randomprojdistortion folder
# Output:
#   minB: minimum bottleneck distance, minW1: minimum 1-wasserstein distance, minW2: min 2-wasserstein distance
#   projPC: list of all projections after rotation, ripsProjPC: list of rips filtration between PC and each projection
#   optProjdB: matrix of projected point cloud with minimum bottleneck, optProjW1: matrix of projected point 
#   with minimum 1-wasserstein, optProjW2: matrix of projected pointcloud with minimum 2-wasserstein, 
#   PCAdbottle: bottleneck distance between  between persistence diagrams of PC and PC projected onto its principle components
#   PCAdW1: 1-wasserstein distance between PC and projected point cloud onto principle components
#   PCAdW2: 2-wasserstein between persistence diagrams of PC and PC projected onto its principle components
# 
randprojmindistortion <-
  function(PC, numDirections, plot)
  {
    #if (!require(package = "sp")) {
    #  install.packages("sp")
    #}
    #library("sp")
    
    if (!require(package = "TDA")) {
      install.packages("TDA")
    }
    library("TDA")
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    if(!require(package = "stats")){
      install.packages("stats")
    }
    library("stats")
    
    #
    # If using a stored point cloud matrix, load it
    #
    if(typeof(PC) == "character"){
      
      pcname <- PC
      PC <- as.matrix(readpointcloud(paste(PC,".txt",sep="")))

      
    } 
    
    #
    # Set defaults for max homological feature dimension to 2
    # default number of directions 12
    # default plot boolean to FALSE
    #
    if(missing(numDirections)){
      numDirections <- 12
    }
    if(missing(plot)){
      plot <- FALSE
    }
    
    
    #
    # Convert coordinate point cloud data to matrix
    #
    PCMatrix <- as.matrix(as.data.frame(PC))
    
    
    #
    # Generate m random angles (theta, phi, psi) representing m planes
    # PC will be projected onto.
    #
    Angles <- list()
    Plist <- list()

    for(i in 1:numDirections){
      
      theta <- runif(1,0,1)*pi
      phi <- runif(1,0,1)*pi
      psi <- runif(1,0,1)*pi

      rproj <- proj(PCMatrix, theta, phi, psi)
      Plist[[i]] <- matrix(c(rproj[,1], rproj[,2]),nrow = nrow(rproj), ncol = 2)
      Angles[[i]] <- as.matrix(c(theta,phi,psi))
    }


    #
    # Perform Principele component ananlysis
    #
    pc.pca <- prcomp(PCMatrix, center = TRUE, scale. = FALSE)
    pcarotations <- pc.pca$rotation
    #
    # Projection of matrix onto principle component subplane
    #
 
    pc.pca <- as.data.frame(PCMatrix%*%as.matrix(pcarotations))
    pc.pca <- as.matrix(cbind(pc.pca$PC1, pc.pca$PC2))

    # 
    #Perform the ripsfiltration of the d dimensional point cloud PC, obtain rips diagram.
    #
    DiagPC <- ripsDiag(X=as.data.frame(PC), maxdimension = 1, maxscale=2, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    DiagPCA <- ripsDiag(X=as.data.frame(pc.pca), maxdimension = 1, maxscale=2, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    #
    # Find the persistence diagram for the projected point clouds
    # 
    
    DiagPlist <- list()
    for( i in 1:length(Plist)){
      
      print(paste("Persistence",toString(i),sep=" "))
      
      DiagP <- ripsDiag(X=as.data.frame(Plist[[i]]), maxdimension = 1, maxscale=2, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
      DiagPlist[[i]] <- DiagP
      
    }
    
    
    #
    # Determine all bottleneck, 1-wasserstein, and 2-wasserstein distannce between 
    # persistence diagram of PC and all projected point clouds
    #
    bottleneckDistances <- list()
    w1Distances <- list()
    w2Distances <- list()
    mindB <- Inf
    mindW1 <- Inf
    mindW2 <- Inf
    optProjBottleneck <- NA
    optProjPCdB <- matrix(NA ,nrow = nrow(Plist[[1]]), ncol = 2)
    optProjW1 <- NA
    optProjPCdW1 <- matrix(NA ,nrow = nrow(Plist[[1]]), ncol = 2)
    optProjW2 <- NA
    optProjPCdW2 <- matrix(NA ,nrow = nrow(Plist[[1]]), ncol = 2)
    optBottleneckPersistence <- list()
    opt1WassersteinPersistence <- list()
    opt2WassersteinPersistence <- list()
    indexbottle <- 1
    indexw1 <- 1
    indexw2 <- 1
    
    for(i in 1:length(DiagPlist)){
      
      print(paste("Number ", toString(i), "of distance measures", sep = ""))
      bDistance = bottleneck(DiagPC, DiagPlist[[i]], dimension=1)
      w1Distance = wasserstein(DiagPC, DiagPlist[[i]], p=1, dimension=1)
      w2Distance = wasserstein(DiagPC, DiagPlist[[i]], p=2, dimension=1)
      if(bDistance < mindB){
        mindB <- bDistance
        optProjBottleneck <- Angles[[i]]
        optProjPCdB <- Plist[[i]]
        optBottleneckPersistence <- DiagPlist[[i]]
        indexbottle <- i
      }
      if(w1Distance < mindW1){
        mindW1 <- w1Distance
        optProjW1 <- Angles[[i]]
        optProjPCdW1 <- Plist[[i]]
        opt1WassersteinPersistence <- DiagPlist[[i]]
        indexw1 <- i
      }
      if(w2Distance < mindW2){
        mindW2 <- w2Distance
        optProjW2 <- Angles[[i]]
        optProjPCdW2 <- Plist[[i]]
        opt2WassersteinPersistence <- DiagPlist[[i]]
        indexw2 <- i
      }
      bottleneckDistances[i] <- bDistance
      w1Distances[i] <- w1Distance
      w2Distances[i] <- w2Distance
      
    }
    
    #
    # Collect the bottleneck and wasserstein distances
    # between the the rips filtrations of the original 
    # point cloud and the principle component rips filtration
    #
    
    PCAbDistance = bottleneck(DiagPC, DiagPCA, dimension=1)
    PCAw1Distance = wasserstein(DiagPC, DiagPCA, p=1, dimension=1)
    PCAw2Distance = wasserstein(DiagPC, DiagPCA, p=2, dimension=1)
    
    # Gather figures and data for Paper
    #
    if(plot){
      print("plotting")
      
      # reset figure margins
      par(mar=c(1,1,1,1))
      # file to write
      RandomProjResults = paste("./randomprojdistortion/", pcname, "randomProjResults.txt",sep="")
      sink(RandomProjResults)
      
      cat("Uniform Perpendicular rings. Bottom Ring Radius = 4, z oriented radius = 2","\n",
          "order bottom = 400. order z oriented = 200. Number directions = 12","\n",append = TRUE)
   
      cat("Minimum Bottleneck Distance: ", mindB, " Onto plane with angle (theta, phi, psi)= "
          , optProjBottleneck ," from (x,y,z) Respectively \n", append= TRUE)   
      cat("Minimum  1st-Wasserstein Distance: ", mindW1, " Onto plane with angle (theta, phi, psi)= "
          , optProjW1 ," from (x,y,z) Respectively","\n", append=TRUE)
      cat("Minimum 2nd-Wasserstein Distance: ", mindW2, " Onto plane with angle (theta, phi, psi)= "
          , optProjW2," from (x,y,z) Respectively", "\n", append=TRUE)
      
      cat("PCA Bottleneck Distance: ", PCAbDistance,"\n", append= TRUE)   
      
      cat("PCA 1st-Wasserstein Distance: ", PCAw1Distance,"\n", append=TRUE)
      
      cat("PCA 2nd-Wasserstein Distance: ", PCAw2Distance, "\n",append=TRUE)
      
      #
      # plot the resulting projection from PCA
      #
      pdf(paste("./randomprojdistortion/", pcname, "PointCloudPCAproj.pdf",sep=""),"pdf")
      plot(as.data.frame(pc.pca))
      dev.off()
      
 
      
      #
      # Perpendicular Ring point cloud
      #
      pdf(paste("./randomprojdistortion/", pcname, "pointcloud.pdf",sep=""),"pdf")
      plot(as.data.frame(PCMatrix))
      dev.off()
      
      #
      # High dimensional barcode
      #
      pdf(paste("./randomprojdistortion/", pcname, "PointCloudBarcode.pdf",sep=""),"pdf")
      plot(DiagPC, barcode = TRUE)
      dev.off()
      
      # 
      # PCA barcode
      #
      pdf(paste("./randomprojdistortion/", pcname, "PCAbarcode.pdf",sep=""),"pdf")
      plot(DiagPCA, barcode = TRUE)
      dev.off()
      
      #
      # Barcode resulting from projection with optimal bottleneck distance
      #
      pdf(paste("./randomprojdistortion/", pcname, "optprojBottleneckBarcode.pdf",sep=""),"pdf")
      plot(optBottleneckPersistence, barcode = TRUE)
      dev.off()
      
      #
      # Projection which provided optimal barcode
      #
      pdf(paste("./randomprojdistortion/", pcname, "optprojBottleneckProjection.pdf",sep=""),"pdf")
      plot(as.data.frame(Plist[[indexbottle]]))
      dev.off()
      
      #
      # Barcode for optimal 1th-wasserstein
      #
      pdf(paste("./randomprojdistortion/", pcname, "optproj1WBarcode.pdf",sep=""),"pdf")
      plot(opt1WassersteinPersistence, barcode = TRUE)
      dev.off()
      
      pdf(paste("./randomprojdistortion/", pcname, "optprojdW1PC.pdf",sep=""),"pdf")
      plot(as.data.frame(Plist[[indexw1]]))
      dev.off()
      
      #
      # Barcode for optimal 2nd-wasserstein
      #
      pdf(paste("./randomprojdistortion/", pcname, "optproj2WBarcode.pdf",sep=""),"pdf")
      plot(opt2WassersteinPersistence, barcode = TRUE)
      dev.off()
      
      pdf(paste("./randomprojdistortion/", pcname, "optprojdW2PC.pdf",sep=""),"pdf")
      plot(as.data.frame(Plist[[indexw2]]))
      dev.off()
      
      sink()
      sink(file=NULL)
    }
    result <- list(minB = mindB, minW1 = mindW1, minW2 = mindW2, projPC = Plist, ripsProjPC = DiagPlist
                   , optProjdB = Plist[[indexbottle]], optProjW1 = Plist[[indexw1]], optProjW2 = Plist[[indexw2]]
                   ,PCAdbottle = PCAbDistance, PCAdW1 = PCAw1Distance, PCAdW2 = PCAw2Distance)
    return(result)
  }



##############################################################################################################################
#                                               Begin Example Point Cloud Generators                                         #
#
#
#               Included are the uniform cylinder, the helix, the uniform circle in 2d, the uniform circle in 3d with       #
#                       with rotations, uniform rings in 2d and 3d with interior loops, rings in 3d perpendicular to each   #
#                       other at one point, and the torus. 
#                                                                                                                           #
#
#
##############################################################################################################################

#
# Constructs a point cloud on the surface of a torus
#
uniformtorus <- 
  function(order, innerradius, outerradius){
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    if (!require(package = "TDA")) {
      install.packages("TDA")
    }
    library("TDA")
    
    torus <- torusUnif(n = order, a = innerradius, c = outerradius)
    if (!require(package = "rgl")) {
      install.packages(pkgs = "rgl")
    }
    library(rgl)
    plot3d(torus)
    
    return(torus)
  }

tetrahedron_edges <-
  function(x){
    
  if (!require(package = "stats")) {
    install.packages("stats")
  }
  library("stats")
  
    
  if (!require(package = "rgl")) {
      install.packages("rgl")
  }
  library("rgl")
 
    
  if (!require(package = "scatterplot3d")) {
      install.packages("scatterplot3d")
  }
  library("scatterplot3d") 
  
  names <- lapply(1:(6*x), toString)
  tetrahedron_pc <- as.data.frame(matrix(NA, ncol = 3, nrow = 6*x), paste(names, sep=''))
  tetrahedron_pc_hack <- matrix(NA, ncol = 3, nrow = x)
  
  tetrahedron_pc <- as.data.frame(matrix(NA, ncol = 3, nrow = 6*x), paste(names, sep=''))
  tetrahedron_pc <- as.data.frame(matrix(NA, ncol = 3, nrow = 6*x), paste(names, sep=''))
  tetrahedron_pc <- as.data.frame(matrix(NA, ncol = 3, nrow = 6*x), paste(names, sep=''))
  
  dist <- runif(x, 0,1)
  disttwo <- runif(x,0,1)
  distthree <- runif(x,0,1)
  

  for (i in 1:length(dist)){
    
    # Add on axis points
    tetrahedron_pc[toString(i),] =  c(disttwo[i],0,0)    # Add edge along x axis
    tetrahedron_pc[toString(i+x),] =  c(0, dist[i], 0)      # Add edge along y axis
    tetrahedron_pc[toString(i+2*x),] =  c(0,0, distthree[i])       # Add edge along z axis

    # Now add edges connecting two axes
    tetrahedron_pc[toString(i+3*x),] = c(tan(pi/2)*(1-dist[i]) , dist[i], 0 )     # Add angled xy edge
    tetrahedron_pc[toString(i+4*x),] = c(disttwo[i], 0, tan(pi/2)*(1-disttwo[i]))    # Add angled xz edge
    tetrahedron_pc[toString(i+5*x),] = c(0, dist[i], tan(pi/2)*(1-dist[i]))          # Add angled yz edge

  }


  return()#tetrahedron_pc)
}

uniformcylinder<-
  function(radius, order, height, density){
    
    
    if (!require(package = "stats")) {
      install.packages("stats")
    }
    library("stats")
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    base <- uniformRing3d(radius,order, 0, .8)$ring3d
    cylinder <- base
    
    for( k in 1:height){
      for( i in 1:order){
        for(j in 1:density){
          heightelm <- k+ runif(1)
          cylinder[i,3]  <-  heightelm
        }
      }
    }
    return(cylinder)
  }


HelixExample <-
  function(){
    library("scatterplot3d")
    library("TDA")
    
    
    HelixDistanceFile = "./HelixFigures/persistenceDistance.txt"
    sink(HelixDistanceFile)
    
    z = seq(-5, 5, 0.07)
    x = cos(2*z)
    y = sin(2*z)
    pdf("./HelixFigures/pointCloud3dHelix.pdf")
    helix3d <- scatterplot3d(x, y, z, highlight.3d = TRUE, col.axis = "blue",
                             col.grid = "lightblue", main = "Helix Point Cloud", pch = 20)
    dev.off()
    helix2d.coords = helix3d$xyz.convert(x, y, z)
    
    #
    # Plot our 3D Helix scatter plot and its corresponding 2D projection
    #
    pdf("./HelixFigures/pointCloud2dHelix.pdf")
    plot(helix2d.coords, xlab = "x",ylab = "z", main ="2D Point Cloud After Projection onto xz-plane")
    dev.off()
    
    
    helix3dmatrix = cbind(x, y, z)
    helix2dmatrix = cbind(helix2d.coords$x, helix2d.coords$y)
    
    
    #
    # we calculate te rips diagram with a max limit of filtration of 11 and max dimension
    # of homological features as 3 and 2 respectively. Dimensionality is considered 
    # as 0 for components, 1 for loops, 2 for voids, ect..
    #
    Diag3d= ripsDiag(X=as.matrix(helix3dmatrix), maxdimension=2, maxscale=6, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    Diag2d= ripsDiag(X=helix2dmatrix, maxdimension=2, maxscale=3, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    
    
    #
    # Plot the barcode diagrams of the rips filtration
    #
    pdf("./HelixFigures/barcode3dHelix.pdf")
    plot(Diag3d, barcode = TRUE)
    dev.off()
    #
    pdf("./HelixFigures/barcode2dHelix.pdf")
    plot(Diag2d, barcode = TRUE)
    dev.off()
    
    #
    # Plot the KDE diagram
    #
    pdf("./HelixFigures/KDE3dHelix.pdf")
    plot(x = Diag3d, main = "KDE Diagram")
    dev.off()
    #
    pdf("./HelixFigures/KDE2dHelix.pdf")
    plot(x = Diag2d, main = "KDE Diagram")
    dev.off()
    
    #
    # calculate the topological distance between persistence diagrams
    # with dimension = 1. We consider 1 dimensional features such as loops
    #
    cat("Bottleneck and Wasserstein with Dimension = 0\n")
    
    bDistance = bottleneck(Diag3d, Diag2d, dimension=0)
    wDistance = wasserstein(Diag3d, Diag2d, p=2, dimension=0)
    
    cat(bDistance, wDistance,"\n", append = TRUE)
    
    # Homological features of dim = 1
    cat("Bottleneck and Wasserstein with Dimension = 1\n")
    
    bDistance = bottleneck(Diag3d, Diag2d, dimension=1)
    wDistance = wasserstein(Diag3d, Diag2d, p=2, dimension=1)
    
    cat(bDistance, wDistance,"\n", append = TRUE)
    
    
    
    #                                                                               #
    #                                                                               #
    #########          Now We Preform a Rotation to Compare Projections     #########
    #                                                                               #
    #                                                                               #
    y_ = seq(-5, 5, 0.07)
    x_ = cos(2*z)
    z_ = sin(2*z)
    pdf("./HelixFigures/pointCloud3dHelixRotated.pdf")
    helix3d_ <- scatterplot3d(x_, y_, z_, highlight.3d = TRUE, col.axis = "blue",
                              col.grid = "lightblue", main = "Helix Point Cloud Rotated", pch = 20)
    dev.off()
    helix2d_.coords = helix3d_$xyz.convert(x_, y_, z_)
    
    #
    # Plot our 3D Helix scatter plot and its corresponding 2D projection
    #
    pdf("./HelixFigures/pointCloud2dHelixRotated.pdf")
    plot(helix2d_.coords, xlab = "x",ylab = "z", main ="2D Points Cloud After Projection of Rotated 3D Helix")
    dev.off()
    
    
    helix3dmatrix_ = cbind(x_, y_, z_)
    helix2dmatrix_ = cbind(helix2d_.coords$x, helix2d_.coords$y)
    
    #
    # we calculate te rips diagram with a max limit of filtration of 11 and max dimension
    # of homological features as 3 and 2 respectively. Dimensionality is considered 
    # as 0 for components, 1 for loops, 2 for voids, ect..
    #
    Diag3d_= ripsDiag(X=as.matrix(helix3dmatrix_), maxdimension=2, maxscale=6, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    Diag2d_= ripsDiag(X=helix2dmatrix_, maxdimension=2, maxscale=3, library="GUDHI", location = TRUE, printProgress=FALSE)$diagram
    
    
    #
    # Plot the barcode diagrams of the rips filtration
    #
    pdf("./HelixFigures/barcode3dHelixRotated.pdf")
    plot(Diag3d_, barcode = TRUE)
    dev.off()
    #
    pdf("./HelixFigures/barcode2dHelixRotated.pdf")
    plot(Diag2d_, barcode = TRUE)
    dev.off()
    
    #
    # Plot the KDE diagram
    #
    pdf("./HelixFigures/KDE3dHelixRotated.pdf")
    plot(x = Diag3d_, main = "KDE Diagram")
    dev.off()
    #
    pdf("./HelixFigures/KDE2dHelixRotated.pdf")
    plot(x = Diag2d_, main = "KDE Diagram")
    dev.off()
    
    #
    # calculate the topological distance between persistence diagrams
    # with dimension = 1. We consider 1 dimensional features such as loops
    #
    cat("Bottleneck and Wasserstein After Rotation with Dimension = 0\n")
    
    bDistance_ = bottleneck(Diag3d_, Diag2d_, dimension=0)
    wDistance_ = wasserstein(Diag3d_, Diag2d_, p=2, dimension=0)
    
    cat(bDistance_, wDistance_,"\n", append = TRUE)
    
    # Homological features of dim = 1
    cat("Bottleneck and Wasserstein After Rotation with Dimension = 1\n")
    
    bDistance_ = bottleneck(Diag3d_, Diag2d_, dimension=1)
    wDistance_ = wasserstein(Diag3d_, Diag2d_, p=2, dimension=1)
    
    cat(bDistance_, wDistance_,"\n", append = TRUE)
    
    
    #######
    sink()
  }

uniformPerpRings <-
  function(radiusXY, radiusZ, orderXY, orderZ, shrink_sizeXY, shrink_sizeZ, plot){
    
    
    if (!require(package = "TDA")) {
      install.packages("TDA")
    }
    library("TDA")
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    #
    # Make shrink and option to plot defaults .8 and FALSE
    #
    if(missing(shrink_sizeXY)){
      shrink_sizeXY <- 0.8
    }
    if(missing(shrink_sizeZ)){
      shrink_sizeZ <- 0.8
    }
    if(missing(plot)){
      plot <- FALSE
    }
    
    
    xyRing <- uniformRing3d(radiusXY, orderXY, -pi/2, .8)$ring3d
    zRing <- uniformRing3d(radiusZ, orderZ, pi/2, .8)$ring3d
    
    for(i in 1:orderZ){
      
      zRing[i,3] <- zRing[i,3] #+ 2*radiusXY#+ sqrt(radiusXY*radiusXY + radiusZ*radiusZ)
      zRing[i,1] <- zRing[i,1] - 2*radiusXY #- sqrt(radiusXY*radiusXY + radiusZ*radiusZ)# radiusXY*radiusXY + radiusZ*radiusZ)
      
    }
    if(plot){
      plot3d(rbind(xyRing,zRing))
      rgl.postscript("./RingFigures/doubleParallel3DRing.pdf","pdf") 
      
    }
    result = list(xyparallel = xyRing, zparallel = zRing, perprings = rbind(xyRing,zRing))
    
    return(result)
  }

uniformRing3d <-
  function(radius, order, angle , shrink_size){
    
    
    if (!require(package = "TDA")) {
      install.packages("TDA")
    }
    library("TDA")
    
    #
    # Gather ring point cloud in 2D and seperate to gather x, y coordinates
    #
    ring <- circleUnif(order, r=radius)
    ringMatrix <- as.array(ring)
    ring2DMatrix <- matrix(ringMatrix, ncol = 2)
    
    #
    # Gather x and y coordinates of 2D point cloud
    #
    x <- ring2DMatrix[,1]
    y <- ring2DMatrix[,2]
    z <- c()
    #
    # Manipulate point cloud data to form an interior ring.
    #
    x_interior <- c()
    y_interior <- c()
    z_interior <- c()
    #
    # Form interior ring and z axis component
    #
    for (i in 0:length(x)){
      x_interior[i] = x[i]*shrink_size
      y_interior[i] = y[i]*shrink_size
    }
    #
    # z axis component and interior
    #
    for(i in 0:length(x))
    {
      z[i] = y[i]*angle
      z_interior[i] = z[i]*shrink_size
    }
    
    #
    # Form matrix of possible 3D point clouds
    #
    ring3DMatrix <- as.matrix(cbind(x,y,z))
    ring3DMatrix_interior <- as.matrix(cbind(c(x,x_interior), c(y,y_interior), c(z,z_interior)))
    ring2DMatrix_interior <- as.matrix(cbind(c(x,x_interior), c(y,y_interior)))
    
    result = list(ring3d = ring3DMatrix, ring2d = ring2DMatrix
                  , ring3Dinterior = ring3DMatrix_interior
                  ,ring2Dinterior = ring2DMatrix_interior)
    return(result)
  }


UniformRing <-
  function(){
    
    if (!require(package = "scatterplot3d")) {
      install.packages("scatterplot3d")
    }
    library("scatterplot3d")
    
    if (!require(package = "TDA")) {
      install.packages("TDA")
    }
    library("TDA")
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    #
    # Parameters for angle, density, and radius of point cloud
    #
    radius <- 4
    order <- 300
    angle <- 70#*(pi/180)
    shrink_size <- 0.5
    
    
    ring3DMatrix <- uniformRing3D(radius,order,angle,shrink_size)$ring3D
    ring2DMatrix <- uniformRing3D(radius,order,angle,shrink_size)$ring2D
    rind3DMatrix_interior <- uniformRing3D(radius,order,angle,shrink_size)$ring3Dinterior
    ring2DMatrix_interior <- uniformRing3D(radius,order,angle,shrink_size)$ring2Dinterior
    
  }


#
# To identify both loops use a rips filtration alpha value of 2
#
uniformswissroll <-
  function(order, h, angle, plot){
    
    if (!require(package = "scatterplot3d")) {
      install.packages("scatterplot3d")
    }
    library("scatterplot3d")
    
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    
    if (!require(package = "cems")){
      install.packages("cems")
    }
    library("cems")
    
    if(missing(plot)){
      plot <- FALSE
    }
    
    
    #
    #create 'order' samples with different parameters
    #phi - number of revolutions
    #nstd - std of normal noise added orthogonally
    #
    sr <- swissroll(order, nstd = .3, height = h, phi = angle)$Xn

    srhole <- sr[1,]
    
    for( i in 1:order){
      point <- c()
      # select all points not in square region
      
      if((!(sr[i,1] > -6 & sr[i,1] < 3) & !(h/4 < sr[i,2] & sr[i,2] < 3*h/4) |  !(3 < sr[i,3] & sr[i,3] < 7) )) {
        point <- cbind(point, sr[i,1])
        point <- cbind(point, sr[i,2])
        point <- cbind(point, sr[i,3])
        rownames(point) <- i
      }
      
      
      
      srhole <- rbind(srhole, point)
      
    }
    
    if(plot){
      pdf("./swissroll/swissroll.pdf","pdf")
      plot3d(srhole)
      dev.off()
    }
    return(srhole)
    
  }


##########################################################################################################################
#
#                                        Functionality for ShortLoop
# Note: From ReadMe.txt from shortloop
#           Terminal command: ShortLoop <input.off> -v -t -a 0.01
#           On parameters: Here -a 0.01 defines the alpha value for the Rips complex.
#           Output: In both cases "ShortLoop" will produce three files:
#                   <input_loops.txt> - description of the loops in the computed basis;
#                   <input_loops.off> - visual representation of the loops;
#                   <input_timing.txt> - timing information.
############################################################################################################################

#
# Gather the coordinates of the points identified by shortloop as belonging to a unique component. 
# in the persistence of the pointcloud.
# OFFfile = name of OFFfile_loop.off, OFFfile.off, ect... file
# pointcloudfile = file with pointcoordinate at unique line of file and number of line = number of row of coordinate
# output: coordinates of each point in each component where i is unique component [[1]...[i]...[n]]
componentpoints <-
  function(OFFfile, pointcloudfile, plot){
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")

    #
    # Read point cloud file
    #
    pointcloud <- readpointcloud(pointcloudfile)
    pointcloud[1,] <- pointcloud[3,]
    pointcloud[2,] <- pointcloud[3,]

    #
    # get number of loops, line number of OFF file points located 
    #
    offFile <- read.table(paste(getwd(),"/ShortLoop/",OFFfile,".off",sep=""), header = FALSE, sep = "\n")

    loopFile <- read.table(paste(getwd(),"/ShortLoop/",OFFfile,"_loops.off",sep=""), header = FALSE, sep = "\n")
    loopFileText <- read.table(paste(getwd(),"/ShortLoop/",OFFfile,"_loops.txt",sep=""), header = FALSE, sep = "\n", skip = 1)
    numloops <- strtoi(strsplit(toString(read.table(paste(getwd(),"/ShortLoop/",OFFfile,"_loops.txt",sep=""), header = FALSE, sep = "\n", skip = 0)[1,], split = " "), split = " ")[[1]][1])
    #text <- makeOFF("perpRings", uniformPerpRings(3,2,0,300,200,.8,.8,FALSE)$perprings)
    
    loopOFFcoord <- list() # a list of matrices where the i'th matrix is the coordinates in the OFF file of the points in loop i
    for( i in 1:numloops){
      loopimatrix <- strsplit(toString(loopFileText[i,]), split = " ")
      loopOFFcoord[i] <- loopimatrix
    }
    
    loopeuclidean <- list() # a list of matrices where the i'th matrix the the euclidean coordinates of the i'th loop
    for(j in 1:length(loopOFFcoord)){
      loopi <- c()
      loopicoord <- loopOFFcoord[[j]] # the coordinate matrix for the i'th loop in the loopcoord list
      for (i in 6:length(loopicoord)){
        loopi <- c(loopi, strtoi(loopicoord[i]))  # the list of euclidean coordinates for the vertices in the i'th loop
      }
      loopeuclidean[[j]] <- loopi
    }
    
    loopvertex <- list()
    for( i in 1:length(loopeuclidean)){
      loopjvertex <- c()
      for( j in 1:length(loopeuclidean[[i]])){
        loopjvertex <- c(loopjvertex, loopFile[j,])
      }
      loopvertex[[i]] <- loopjvertex
    }
    
    #
    # Collect the coordinates from the pointcloud file
    #
    componentcoords <- list(NULL)
    for(i in 1:length(loopeuclidean)){
      index <- loopeuclidean[[i]]
      component <- c()
      for( j in 1:(length(index))){
        point <- pointcloud[index[j-1],]
        component <- rbind(component, point)
      }
      component <- as.data.frame(component)
      componentcoords[[i]] <- component
      
    }
    

    
    if(plot==TRUE){
      plot3d(pointcloud, col = "gray", alpha = 0.7)
      for(i in 1:length(componentcoords)){
        print(i*25+47)
        points3d(componentcoords[[i]], col=colors()[i*50+47])
      }
      
      #lapply(componentcoords, function(x) points3d(x, col = 1:length(componentcoords)))
      rgl.postscript(paste(getwd(),"/loopvertexplot/", OFFfile,"plot.pdf", sep=""),"pdf")
    }
    return(componentcoords) 
  }


makeOFF<-
  function(filename, pointcloud){
    file.create(paste(getwd(),"/ShortLoop/",filename,".off",sep=""), overwrite = TRUE)
    OFFfile = paste(getwd(),"/ShortLoop/" , filename , ".off",sep="")
    #"/Documents/phd/TopologicalDistortionProject/DistortionCode/ShortLoop")
    sink(OFFfile)
    
    cat("OFF", "\n",append = TRUE)
    cat(paste(toString(nrow(pointcloud)),"0","0",sep="\t"), "\n",append = TRUE)
    for(i in 1:(nrow(pointcloud)-1)){
      
      cat(paste(toString(pointcloud[i,1]),toString(pointcloud[i,2]),toString(pointcloud[i,3]),sep="\t"),"\n", append = TRUE)
      #write(toString(i),"", append=TRUE)
    }
    cat(paste(toString(pointcloud[nrow(pointcloud),1]),toString(pointcloud[nrow(pointcloud),2])
              ,toString(pointcloud[nrow(pointcloud),3]),sep="\t"), append = TRUE)
    
    sink(file = NULL)
    sink(file = NULL)
    
    return()
  }


pointcloudtotext <-
  function(filename, pointcloud){
    
    if (!require(package = "MASS")) {
      install.packages("MASS")
    }
    library("MASS")
    
    file.create(paste(getwd(),"/pointclouds/",filename,".txt", sep=""),overwrite = TRUE)
    pointcloudfile = paste(getwd(), "/pointclouds/",filename, ".txt",sep="")
    write.matrix(pointcloud, pointcloudfile, sep = " ")
  }

readpointcloud <-
  function(filename){
    
    pointcloudmatrix <- read.table(paste(getwd(),"/pointclouds/", filename, sep=""), header = TRUE, sep="") 
    
    return(pointcloudmatrix)
  }
#
# Write persistence diagrams to .txt file for calculating the pth-wasserstein 
# distance and bottleneck in C++
#
diagramtotext <-
  function(diagram1, diagram2){
    
    diagram1text = paste(getwd(),"/wasserstein/diagram1.txt", sep="")
    sink(diagram1text)
    write("#Persistence Diagram 1", sep="", append = TRUE)
    for(i in 1:nrow(diagram1)){
      cat(paste(toString(diagram1[i,2])," ",toString(diagram1[i,3]), sep=""),"\n",append = TRUE)
    }
    sink()
    sink(file=NULL)
    
    diagram2text = paste(getwd(),"/wasserstein/diagram2.txt", sep="")
    sink(diagram2text)
    write("#Persistence Diagram 2", sep="", append = TRUE)
    for(i in 1:nrow(diagram2)){
      cat(paste(toString(diagram2[i,2])," ",toString(diagram2[i,3]), sep =""),"\n",append = TRUE)
    }
    
    sink(file=NULL)
    
  }

###########################################################################################################
#   
#                           Method Used For Obtaining Plots and Recording Results
#
##########################################################################################################




main <-
  function( which ){
    
    if (!require(package = "rgl")) {
      install.packages("rgl")
    }
    library("rgl")
    
    
    #
    # Output: uniformRing.off, uniformRingPC.txt, uniformPerpRing.off, uniformPerpRingPC.txt 
    #         and point clouf plots
    # 
    if( which == 1){
      uniformRing <- uniformRing3D(5,300,pi/4,.8)
      plot3d(as.data.frame(uniformRing))
      rgl.postscript(paste(getwd(),"\\DistortionCode\\progressreport\\uniformRing.pdf",sep=""),"pdf")
      
      makeOFF("uniformRing", uniformRing$ring3D)
      pointcloudtotext("uniformRingPC.txt", uniformRing$ring3D)
      
      
      uniformRingPerp <- uniformPerpRings(5,3,0,400,300,.8,.8,FALSE)
      plot3d(as.data.frame(uniformRingPerp$perprings))
      rgl.postscript(paste(getwd(),"\\DistortionCode\\progressreport\\uniformPerpRing.pdf",sep=""),"pdf")
      
      makeOFF("uniformPerpRing", uniformRingPerp$perprings)
      pointcloudtotext("uniformPerpRingPC.txt", uniformRingPerp$perprings)
    }
    
    #
    #  Plot principle components from which == 1
    #
    if(which == 2 ){
      
      pointsUniformRing <- componentpoints("uniformRing", "uniformRingPC.txt")
      uniformRing <- readpointcloud("uniformRingPC.txt")
      plotcomponents(pointsUniformRing,uniformRing, "uniformRingLoops.pdf")
      
      sink(file = NULL)
      
      pointsUniformPerpRing <- componentpoints("uniformPerpRing", "uniformPerpRingPC.txt")
      uniformPerpRing <- readpointcloud("uniformPerpRingPC.txt")
      plotcomponents(pointsUniformPerpRing,uniformPerpRing, "uniformPerpRingLoops.pdf")
      
      sink(file = NULL)
      
    }
    
    #
    # Output: .OFF files for ring with interior loop, make matrix file, list of vertex
    #   points associated to each 1-cycle connetcted component, plot each component list.
    #
    if(which == 3){
      sink(file = NULL)
      ringinterior <- uniformRing3D(5,400,pi/4,.7)$ring3Dinterior
      makeOFF("ringinterior",ringinterior)
      pointcloudtotext("ringinteriorPC.txt",ringinterior)
      ringinteriorcomponents <- componentpoints("ringinterior", "ringinteriorPC.txt")
      plotcomponents(ringinteriorcomponents,ringinterior,"ringInterior.pdf")
      
      sink(file = NULL)
    }
    
    if(which == 4){
      ringperpcomp <- componentpoints("uniformPerpRing", "uniformPerpRing.txt", FALSE)
      uniformcylindercomp <- componentpoints("cylinder", "cylinder.txt",FALSE)
      swisscomp <- componentpoints("swissrolltwoholes","swissrolltwoholes.txt",FALSE)
      toruscomp <- componentpoints("torus", "torus.txt",FALSE)
      ringwithinteriorcomp <- componentpoints("ringwithinterior", "ringwithinterior.txt",FALSE)
      
      pcacomp(ringwithinteriorcomp,"ringwithinterior",FALSE)
      pcacomp(toruscomp,"torus",FALSE)
      pcacomp(swisscomp,"swissrolltwoholes",FALSE)
      pcacomp(uniformcylindercomp,"cylinder",FALSE)
      pcacomp(ringperpcomp,"perprings",FALSE)
      
            ringperpcomp <- componentpoints("uniformPerpRing", "uniformPerpRing.txt", FALSE)
            pcacomp(ringperpcomp,"perprings",FALSE)
    }
    
    if(which == 5){
      #
      # Calculate distortion between point cloud and component based PCA projection
      # and plot barcodes for uniform ring with interior loop point cloud. 
      #
          pointcloud = "ringwithinterior"
          numfirstprin = 1
          numsecondprin = 2
      
      #
      # Calculate distortion and plot barcodes for swiss roll with two holes point 
      # cloud point. 
      #
          pointcloud = "swissrolltwoholes"
          numfirstprin = 1
          numsecondprin = 3
      
      #
      # Calculate distortion and plot barcodes for point cloud on torus surface point 
      # cloud point. 
      # overflow error if maxdimension > 1 and maxscale > 3 for ripsfiltration of pointcloud
      # overflow error when computing the rips diagram for projected point cloud
          #pointcloud = "torus"
          #numfirstprin = 1
          #numsecondprin = 2
      
      #
      # Calculate distortion and plot barcodes for uniform perpendicular point 
      # cloud point set. 
      #
          #pointcloud = "perprings"
          #numfirstprin = 1
          #numsecondprin = 3

      # load point cloud from file
      pointcloud_matrix <- as.matrix(readpointcloud(paste(pointcloud,".txt",sep="")))
      
      # identify loop features found by shortloop. Requires shortloop _loop.off file and 
      #  stored point cloud matrix in .txt file (can use pointcloudtotext() )
      #  pointcloud is name of .off file such as pointcloud.off
      components <- componentpoints(pointcloud, paste(pointcloud,".txt",sep=""), FALSE)
      
      # obtain list of point cloud projects based on PCA from interior loops
      loopPCA <- pcacomp(components ,pointcloud, FALSE)
      
      # call distortionloopPCA()
      distortionloopPCA(pointcloud, loopPCA, numfirstprin, numsecondprin)
    }
    return()
  }



