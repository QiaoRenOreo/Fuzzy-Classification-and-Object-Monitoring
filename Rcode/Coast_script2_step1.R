#==================================================================================
# R code for processing one landsat image 10 times with 5 classes
#==================================================================================

#======================================================================================
# Block 1: variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))
graphics.off()

require(rgdal)

# Show intermediate results? (slow)
ShowAll <- TRUE

# Fuzzy parameter
m <- 3 #2.0 #1.5-2.5

# Number of classes
Ncl <- 4

# Set path to directories
Root <- "D:/Study/Module8/IndividualAssignment"  

# Input directory
Path_in <- paste(Root, "/Input",sep="")

# Output directory
Path_out <- paste(Root,"/ObjectMonit_step1_output_m1",sep="")
dir.create(Path_out, recursive = TRUE,showWarnings=FALSE)

# Name of the input file (excluding file extenstion .img)
#imagefn.in <- "lc820130827_s5"
imagefn.in <- "lc820170619_s5"

imageext.in <- "img"
imagefn.out <- paste("FCM_",imagefn.in,sep="")

histstretch<-function(data){
  cur.lim<-quantile(data,c(0.025,0.975),na.rm=TRUE)
  data<-pmax(cur.lim[1],pmin(cur.lim[2],data))
  data<-floor(255*(data-cur.lim[1])/(cur.lim[2]-cur.lim[1]))
  data[is.na(data)]<-0
  data
}

# Combine path, file name and file extension, then read the file
MS <- readGDAL(paste(Path_in,"/",imagefn.in,".",imageext.in,sep="")) 

# Number of bands
Nb <- dim(MS@data)[2]

MS.image <- MS

# Set RGB composition 
#543 for falsecolor 8 band landsat8 (multispectral bands)
#432 for falsecolor of 6 band landsat8 (only the ones comparable to multispectral bands of landsat 7), input images with "r"at the end, or ordinary Landsat 7
nR <- 5
nG <- 4
nB <- 3

MS.image@data$red   <-histstretch(MS.image@data[,nR])
MS.image@data$green <-histstretch(MS.image@data[,nG])
MS.image@data$blue  <-histstretch(MS.image@data[,nB])

windows(7,7)
image(MS.image,red="red",green="green",blue="blue",axes=TRUE)
title(paste(imagefn.in, " RGB= ",nR,":",nG,":",nB,sep=""))

# Image dimensions
d <- MS@grid@cells.dim
psize <- MS@grid@cellsize

M <- d[1]
N <- d[2]

x <- as.matrix(MS@data)
#======================================================================================
# 2. Initialization of membership values. 
# Function sample_frac produces random membership value vector with sum of elements = 1
#======================================================================================
Nruns <- 1
Run_names <- array("0",Nruns)

Run_names[] <- 1:Nruns

for(run in 1:Nruns){
  run_name <- paste("_run_",Run_names[run],sep="")
  
  n <- nrow(MS)
  V <- array(0,c(Ncl,Nb))  
  
  U <- array(0,c(n,Ncl))
  Uref <- U
  
  V[,] <- matrix(rep(colMeans(x),Ncl),byrow=TRUE,nrow=Ncl)
  
  # Add random vector to the mean value of each class 
  amplitude <- array(0,Nb)
  
  for(i in 1:Nb) amplitude[i] <- sd(x[,i])
  
  for(k in 1:Ncl)
    for(i in 1:Nb)
      V[k,i] <- V[k,i] + runif(1,min=-amplitude[i],max=amplitude[i])
  
  Md <- array(0,c(Ncl,n))
  Mdk <- array(0,c(Ncl,Ncl,n))
  
  # Update pixel membership values
  for(k in 1:Ncl){
    #Md[l,,] <- colSums((D-V[l,])^2,1)
    Md[k,] <- rowSums((t(t(x)-V[k,]))^2,1)
  }
  
  for(l in 1:Ncl)
    for (k in 1:Ncl){
      Mdk[k,l,] <- (Md[l,]/Md[k,])^(1/(m-1))
    }
  
  U <- t(1.0/colSums(Mdk,1))
  #======================================================================================
  # 3. Display initial membership values
  #======================================================================================
  FCM <- MS
  FCM@data <- data.frame(U)
  
  if(ShowAll){
    #windows(title="Iinitial membership values")
    Nrow <- round(sqrt(Ncl))
    par(mfrow=c(Nrow,ceiling(Ncl/Nrow)))
    
    for(k in 1:Ncl){
      image(FCM,attr=k,col=gray((0:255)/255),axes=TRUE)
      title( main = paste("Class ",k,sep=""))
    }
    
    "Initial mean values:"
    V
  }
  
  #======================================================================================
  # 4. Fuzzy C means computations
  #======================================================================================
  
  max_val <- max(x)
  
  fcm_objfun <- function(){
    val <- 0
    
    for(k in 1:Ncl){
      
      tmp <- rowSums((t(t(x)-V[k,]))^2,1)
      val <- val + sum(tmp*(U[,k]^m))/(n*(max_val^2))
    }
    
    return(val)
  }
  
  Nit <- 1000
  eps <- 1.0e-03
  
  Md <- array(0,c(Ncl,n))
  Mdk <- array(0,c(Ncl,Ncl,n))
  
  parr <- rep(18,Ncl)
  #colarr <- 2:(1+Ncl)
  
  
  col_arr <- rainbow(Ncl)
  
  if(ShowAll){
    windows(12,6)
    par(mfrow=c(1,2))
  }
  
  Errcount <- 0
  Err <- 1
  
  ind <- sample(n,1000,replace=FALSE)
  
  Err_hist <- array(0,Nit)
  
  for(iter in 1:Nit){
    Uref <- U
    Vref <- V
    
    # Update pixel membership values
    for(k in 1:Ncl){
      #Md[l,,] <- colSums((D-V[l,])^2,1)
      Md[k,] <- rowSums((t(t(x)-V[k,]))^2,1)
    }
    
    for(l in 1:Ncl)
      for (k in 1:Ncl){
        Mdk[k,l,] <- (Md[l,]/Md[k,])^(1/(m-1))
      }
    
    U <- t(1.0/colSums(Mdk,1))
    
    # Update cluster mean values
    for (l in 1:Nb)
      for (k in 1:Ncl){
        V[k,l] <- sum((U[,k]^m)*x[,l])/sum(U[,k]^m)
      }
    
    #Err <- sum(abs(Uref-U))/(n*Ncl)
    #Err <- sum(abs(Vref-V))/(Ncl*Nb)
    Err_new <- max(sqrt(rowMeans((U-Uref)^2)))
    #Err_new <- fcm_objfun()
    Err_hist[iter] <- Err_new
    
    if(abs(Err_new-Err)<eps) Errcount<-Errcount+1 else Errcount=0
    if(Errcount>=3) break
    
    Err <- Err_new
    
    if(ShowAll){
      FCM@data[,1] <- U[,1]
      image(FCM, attr=1, col=gray((0:255)/255))
      title(main = paste("Iteration=",iter, " Mean band ",nR,"=",round(V[1,nR],3),sep=""))
      
      plot(x[ind,nR],x[ind,nG], xlab = paste("Band",nR,sep=""), ylab = paste("Band",nG,sep=""),cex=0.3,pch=16)
      for(i in 1:Ncl){
        points(V[i,nR],V[i,nG],col=col_arr[i],cex=1.1,pch=parr[i])
        points(V[i,nR],V[i,nG],col=col_arr[i],cex=3,pch=3)
      }
    }
  }
  
  # Convergence plot, for testing purpose
  if(ShowAll){
    windows()
    Err_hist[1] <-1
    plot(Err_hist[1:iter])
    lines(Err_hist[1:iter])
  }
  
  #======================================================================================
  # 5. Display results
  #======================================================================================
  FCM <- MS
  FCM@data <- data.frame(U)
  
  cl_names <- "Class1"
  for(j in 2:Ncl)cl_names <- c(cl_names,paste("Class",j,sep=""))
  
  names(FCM) <- cl_names
  
  #windows(title="FCM result: membership values")
  Nrow <- round(sqrt(Ncl))
  par(mfrow=c(Nrow,ceiling(Ncl/Nrow)))
  
  for(k in 1:Ncl){
    image(FCM, attr=k, col=gray((0:255)/255),axes=TRUE)
    title(main = paste("Class ",k,sep=""))
  }
  
  #"Number of iterations:"
  #iter
  
  #"Mean values:"
  #V
  
  #======================================================================================
  # 6. Display feature space
  #======================================================================================
  if(ShowAll){
    k<-nR
    l<-nB
    
    windows(title="Feature space with legend")
    
    plot(x[ind,k],x[ind,l], xlab = paste("Band",k,sep=""), ylab = paste("Band",l,sep=""),cex=0.1)
    
    text_str <- "Class1"
    
    for(i in 1:Ncl){
      points(V[i,k],V[i,l],col=col_arr[i],cex=2,pch=parr[i])
      if(i>1) text_str <- c(text_str,paste("Class",i,sep=""))
    }
    
    legend("right",text_str,fill=col_arr)
  }
  #======================================================================================
  # 7. Display confidence image
  #======================================================================================
  if(ShowAll){
    Win1 <- array(0,n)
    Win2 <- array(0,n)
    
    Win1p <- array(0,n)
    Win2p <- array(0,n)
    
    for(i in 1:n){
      k <- which.max(U[i,])
      Win1[i] <- k
      Win1p[i] <- U[i,k]
      
      l <- which.max(U[i,-k])
      Win2[i] <- l
      Win2p[i] <- (U[i,-k])[l]
    }
    
    m12 <- 1 - (Win2p/Win1p)
    
    Win <- FCM
    Win$conf <- m12
    Win$first <- Win1
    
    col_arr <- rainbow(Ncl) 
    
    #windows()
    par(mfrow=c(1,2))
    image(Win, attr="first",col=col_arr,axes=TRUE)
    title("Winner 1")
    plot.new()
    legend("center",text_str,fill=col_arr)
    
    
    #windows()
    par(mfrow=c(1,2))
    image(Win,attr="conf", col=gray((0:255)/255), axes=TRUE)
    title(main="FCM result: confidence image")
    hist(m12,main="Histogram of confidence image")
  }
  
  #======================================================================================
  # 8. Output the results
  #======================================================================================
  
  write.table(V, file = paste(Path_out,"/",imagefn.in,run_name,"_FCM_Mean.txt",sep=""),append=FALSE,quote=TRUE,sep =" ",eol="\n",na="NA",dec=".",row.names=FALSE,col.names=FALSE,qmethod=c("escape","double"))
  write.table(c(iter,fcm_objfun()), file = paste(Path_out,"/",imagefn.in,run_name,"_FCM_metadata.txt",sep=""),append=FALSE,quote=TRUE,sep =" ",eol="\n",na="NA",dec=".",row.names=FALSE,col.names=FALSE,qmethod=c("escape","double"))
  
  # Merge classes
  #FCM$water <- rowSums(FCM@data[,c(2)])
  ##windows()
  #image(FCM, attr="water", col=gray((0:255)/255),axes=TRUE)
  #title("Merged classes")
  
  
  #Which class is water class?
  cl <-  which.min(rowSums(V[,4:6]))
  FCM$water <- FCM@data[,cl]
  
  # Merge classes (manual input class numbers!!!)
  # cl1 <- 1
  # cl2 <- 2
  # FCM$water <- FCM@data[,cl1]+FCM@data[,cl2]
  
  image(FCM, attr="water", col=gray((0:255)/255),axes=TRUE)
  title(main = "Water class")
  
  OUT <- FCM
  
  #OUT@data <- data.frame(FCM@data[,cl])
  OUT@data <- data.frame(FCM$water)
  
  setwd(Path_out)
  
  OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="Float32")
  saveDataset(OUT.tif,paste(imagefn.out,run_name,".tif",sep=""))
  GDAL.close(OUT.tif)
  
  setwd(Root)
}

# The End
