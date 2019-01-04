#======================================================================================
# Block 1: variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))
graphics.off()

require(rgdal)

Root <-"D:/Study/Module8/IndividualAssignment"

Path_in <- paste(Root, "/ObjectMonit_step1_output",sep="")
Path_out <- paste(Root, "/ObjectMonit_step1_output",sep="")

# Name of the Landsat input file (excluding file extenstion .img)
#image.fn <- "FCM_lc820130827_s5_run_1"
image.fn <- "FCM_lc820170619_s5_run_1"
imageext.in <- "tif"

imagefn.out <- paste("Water_Object_",image.fn,sep="")

histstretch<-function(data) {
  cur.lim<-quantile(data,c(0.025,0.975),na.rm=TRUE)
  data<-pmax(cur.lim[1],pmin(cur.lim[2],data))
  data<-floor(255*(data-cur.lim[1])/(cur.lim[2]-cur.lim[1]))
  data[is.na(data)]<-0
  data
}

FCM <- readGDAL(paste(Path_out,"/",image.fn,".",imageext.in,sep="")) 

# Number of bands
Ncl <- dim(FCM@data)[2]

if(Ncl>1){
  print("ERROR. Input file contains more than one band. Replace by a proper file.")
  FCM <- NULL
}

windows()
image(FCM, col=gray((0:255)/255),axes=TRUE)
title(main = "Membership of water class")

# Image dimensions
d <- FCM@grid@cells.dim
# get pixel size
psize <- FCM@grid@cellsize

M <- d[1]
N <- d[2]
n <- nrow(FCM)

Path_out <- paste(Root,"/ObjectMonit_step2_output_thsh70",sep="")
dir.create(Path_out, recursive = TRUE,showWarnings=FALSE)

wm <- as.vector(FCM@data[1:n,1])

windows()
par(mfrow=c(1,2))
image(FCM,col=gray((0:255)/255),axes=TRUE)
title("Membership for class water")
hist(wm)

# Thresholding

thr <- 0.70

wmt <- wm

wmt[wmt<thr] <- 0
wmt[wmt>0] <- 1

FCM$wmt <- wmt

windows()
image(FCM,attr="wmt",col=gray((0:255)/255),axes=TRUE)
title(main=paste("thresholded membership, threshold=", thr,sep=""))

OUT <- FCM

setwd(Path_out)

OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="Float32")
saveDataset(OUT.tif,paste(imagefn.out,".tif",sep=""))
GDAL.close(OUT.tif)
#==================================================================
# Block 2: Define necessary neighbourhood relationships for region growing
#==================================================================
Nm	<- array(0, c(M, N))
Nsx	<- array(0, c(M, N, 8))
Nsy	<- array(0, c(M, N, 8))

for(i in 2:(M-1))
  for(j in 2:(N-1))
  {
    Nm[i, j] <- 8
    Nsx[i, j, 1:Nm[i, j]] <- c(i-1,i+1,i  ,i  ,i-1,i-1,i+1,i+1)
    Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j  ,j-1,j+1,j-1,j+1,j-1,j+1)
  }

i <- 1
for(j in 2:(N-1))
{
  Nm[i, j] <- 5
  Nsx[i, j, 1:Nm[i, j]] <- c(i+1,i  ,i  ,i+1,i+1)
  Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j-1,j+1,j-1,j+1)
}

i <- M
for(j in 2:(N-1))
{
  Nm[i, j] <- 5
  Nsx[i, j, 1:Nm[i, j]] <- c(i-1,i  ,i  ,i-1,i-1)
  Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j-1,j+1,j-1,j+1)
}

j <- 1
for(i in 2:(M-1))
{
  Nm[i, j] <- 5
  Nsx[i, j, 1:Nm[i, j]] <- c(i-1,i+1,i  ,i-1,i+1)
  Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j  ,j+1,j+1,j+1)
}

j <- N
for(i in 2:(M-1))
{
  Nm[i, j] <- 5
  Nsx[i, j, 1:Nm[i, j]] <- c(i-1,i+1,i  ,i-1,i+1)
  Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j  ,j-1,j-1,j-1)
}

i<-1
j<-1
Nm[i, j] <- 3
Nsx[i, j, 1:Nm[i, j]] <- c(i+1,i  ,i+1)
Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j+1,j+1)

i<-1
j<-N
Nm[i, j] <- 3
Nsx[i, j, 1:Nm[i, j]] <- c(i+1,i  ,i+1)
Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j-1,j-1)

i<-M
j<-1
Nm[i, j] <- 3
Nsx[i, j, 1:Nm[i, j]] <- c(i-1,i  ,i-1)
Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j+1,j+1)

i<-M
j<-N
Nm[i, j] <- 3
Nsx[i, j, 1:Nm[i, j]] <- c(i-1,i  ,i-1)
Nsy[i, j, 1:Nm[i, j]] <- c(j  ,j-1,j-1)
#======================================================================================
# Block 3: Segments an image F(MxN) with Ncl classes, for class Cl
#======================================================================================
Segment <- function(A,Cl)
{
  if(sum(A==Cl)==0)
  {
    report <- list(0, NA)
    names(SReport) <- c("area","pn")
    return(report)
  }
  
  Nleft <- sum(A==Cl)
  
  Lseg_x <- array(0,c(Nleft))
  Lseg_y <- array(0,c(Nleft))
  Lseg   <- array(0,c(Nleft,2))
  
  Nobj <- 0
  Area <- array(0,Nleft)
  Pix_count <- 0
  Offset_seg <- 0
  
  while(Nleft>=1) 
  {
    Nobj <- Nobj+1	
    Pix_count <- Pix_count +1
    t0 <- which(A==Cl,arr.ind=TRUE)
    i <- t0[1,1]
    j <- t0[1,2]
    
    Ncount <- 1
    
    Lseg[Pix_count,] <- c(i,j)
    A[i,j] <- A[i,j] + 1
    
    
    Grow <- TRUE
    while(Grow)
    {
      Found <- FALSE
      
      k<-1
      while(k<=Ncount)
      {
        i <- Lseg[Offset_seg+k,1]
        j <- Lseg[Offset_seg+k,2]
        
        #		  	   	Check neighbours
        for(l in 1:Nm[i,j])
        {
          if(A[Nsx[i,j,l],Nsy[i,j,l]]==Cl)
          {
            Ncount <- Ncount + 1
            Pix_count <- Pix_count +1
            Lseg[Pix_count,] <- c(Nsx[i,j,l],Nsy[i,j,l]) 
            A[Nsx[i,j,l],Nsy[i,j,l]] <- A[Nsx[i,j,l],Nsy[i,j,l]] +1
            Found <- TRUE
          }
        }
        k<-k+1
      }
      if(!Found) Grow <- FALSE
    }
    
    Area[Nobj] <- Ncount
    Offset_seg <- Pix_count
    
    A[A==(Cl+1)] <- 0
    
    Nleft <- sum(A==Cl)
  }
  
  Area_max <- max(Area)
  RArea <- array(0,Nobj)
  RArea <- Area[1:Nobj]
  
  # Sort the segments by area (descending)
  
  R_sorted <- array(0,c(sum(RArea),2))	# Sorted array of pixel positions per segment
  Sort_count <- 0
  SPix_count <- 0
  SArea<- array(0,Nobj)		# Sorted array of areas
  TArea<- RArea			# Temporary array
  
  while(Sort_count<Nobj)
  {
    w <- which(TArea==max(TArea))[1]
    Sort_count <- Sort_count +1
    n <- TArea[w]
    SArea[Sort_count] <- n
    Rstart <- 0
    if(w>1)Rstart<-sum(RArea[1:(w-1)])
    
    R_sorted[(SPix_count+1):(SPix_count+n),] <- Lseg[(Rstart+1):(Rstart+n),]
    SPix_count <- SPix_count+n
    TArea[w] <- 0
  }
  
  SReport <- list(Nobj, SArea, R_sorted)
  names(SReport) <- c("nobj","area","xy")
  
  return(SReport)
}

Dm <- matrix(wmt,nrow=M,ncol=N,byrow=FALSE)

seg <- Segment(Dm,1)

# Get the largest segment
xy <- seg$xy[1:(seg$area[1]),]

#display the largest segment

FCM$seg <- array(0,n)

ind <- xy[,1]+(xy[,2]-1)*M
FCM$seg[ind] <- 1

windows()
image(FCM,attr="seg",col=gray((0:255)/255),axes=TRUE)
title("largest segment")

# Compute area in m^2
A <- seg$area[1]* psize[1]*psize[2]
print(paste("Area of the largest segment is",A, "m^2",sep=" "))

#==============================================================================
# Get the second largest segment

xy <- seg$xy[(seg$area[1]+1):(sum(seg$area[1:2])),]

#display the 2nd largest segment

FCM2 <-FCM
FCM2$seg <- array(0,n)

ind <- xy[,1]+(xy[,2]-1)*M
FCM2$seg[ind] <- 1

windows()
image(FCM2,attr="seg",col=gray((0:255)/255),axes=TRUE)
title("2nd largest segment")

# Compute area in m^2
A <- seg$area[2]* psize[1]*psize[2]
print(paste("Area of the 2nd largest segment is",A, "m^2",sep=" "))

#======================================================================================
# Block 4: Output the results
#======================================================================================

OUT <- FCM

OUT@data <- data.frame(FCM$seg)

setwd(Path_out)

OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="Float32")
saveDataset(OUT.tif,paste(imagefn.out,".tif",sep=""))
GDAL.close(OUT.tif)

setwd(Root)

# The End
