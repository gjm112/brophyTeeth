############################################################################################################
#How can we (or can we at all?) we use both the MTurk images and Juliet's dissertation 
############################################################################################################
#install_github("vbonhomme/Momocs")
library(geomorph)
library(EBImage)
library(shapes)
library(devtools)
library(Momocs)
library(plyr)
library(caret)
library(randomForest)
library(e1071)
library(caret)
numHarm<-20

################################################
#Mechanical Turk data.  Black and white images.
################################################
imgsList <- list()
#vector of image names
dscn <- list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/resultsClean")

for (d in dscn){
  print(d)
  
  files <- list.files(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/resultsClean/",d))
  files <- files[grep("jpg",files)] #Pull out only .jpg
  
  #Import all tracings of a tooth.  Up to three copies.  
  imgsList[[d]] <- import_jpg(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/resultsClean/",d,"/",files,sep=""))
  
  #Manual removal of this MTurk image
  imgsList[[d]][which(names(imgsList[[d]])=="3A1COHJ8NJU7FB7Z1T98O6918IG8H5")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3C8HJ7UOP7T8RL9X1GPYTVE1ODMZMA")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3R0T90IZ1SBVRI21YZ7V5STJJ9PGCL")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3NJM2BJS4W514VV01IXIZ17BK0BCPK")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="37C0GNLMHF23ZHJ9MITKD7YCASHD6T")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3N1FSUEFL5ZPKUFV3U05G9EYE454DZ")]<-NULL
  
  #Remove images with too few observations.  These are clearly mistakes.  
  remove <- names(which(unlist(lapply(imgsList[[d]], nrow)) <= 50))
  if (length(remove) > 0){
    for (q in remove){
      imgsList[[d]][[q]] <- NULL
    }
  }
  
  
}

#List of image names that have been loaded into R.  
dscn <- names(imgsList)

#Do procrustes and then find the average tooth.  
imgsListMean <- list()
for (d in dscn){print(d)
  
  test <- efourier(Out(imgsList[[d]]), nb.h = numHarm)
  numPoints <- 200
  ptsList <- list()
  arr <- array(NA, dim=c(numPoints, 2, dim(test$coe)[1]))
  for (j in 1:dim(test$coe)[1]){
    #arr[,,j]<-ptsList[[j]]<-efourier_shape(test$coe[j,grep("A",colnames(test$coe))],test$coe[j,grep("B",colnames(test$coe))],test$coe[j,grep("C",colnames(test$coe))],test$coe[j,grep("D",colnames(test$coe))],plot=TRUE,nb.pts=numPoints)
    efa <- list(an = test$coe[j,grep("A",colnames(test$coe))],
              bn = test$coe[j, grep("B", colnames(test$coe))],
              cn = test$coe[j, grep("C", colnames(test$coe))],
              dn = test$coe[j, grep("D", colnames(test$coe))],
              a0 = 0,c0= 0)
    arr[, , j] <- ptsList[[j]] <- efourier_i(efa, nb.pts = numPoints)
  }
  
  #set of shapes to be predicted.  
  imgsListMean[[d]] <- mshape <- as.matrix(procGPA(arr)$mshape)
  
}





################################################
#Training Data
################################################
ptsTrainList<-list()
refFile<-data.frame()

setwd("/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/")

tribeVec <- list.files()
for (t in tribeVec){print(t)
  speciesVec<-list.files(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/",t))
  for (s in speciesVec){print(s)
    teethVec<-list.files(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/",t,"/",s))
    for (tooth in teethVec){
      filesVec<-list.files(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/",t,"/",s,"/",tooth))
      for (f in filesVec){
        file<-(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/",t,"/",s,"/",tooth,"/",f))
        
        if (file != "/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/Antilopini/Antidorcas marsupialis/LM2 A marsupialis/DSCN3190"){
          temp<-matrix(unlist(read.table(file))[-c(1:20)],ncol=13,byrow=TRUE)[,-13]
          
          pts<-data.frame(x=c(t(temp[,seq(1,11,2)])),y=c(t(temp[,seq(1,11,2)+1])))
        }
        
        if (file == "/Users/gregorymatthews/Dropbox/brophyTeeth/Raw data/Point files/Antilopini/Antidorcas marsupialis/LM2 A marsupialis/DSCN3190"){
          pts<-read.table(file)
          names(pts)<-c("x","y")
        }
        
        #Check if the tooth is in the test data set.  If it is don't add it to the training data set.
        if (!f%in%toupper(names(imgsListMean))){
          ptsTrainList[[substring(tooth,1,3)]][[f]]<-as.matrix(pts)
          refFile<-rbind(refFile,data.frame(ref=f,tooth=substring(tooth,1,3),tribe=t,species=s))
        }
        
      }
    }
  }
}

ptsTrainListNorm <- list()
for (toothname in names(ptsTrainList)){
ptsTrainListNorm[[toothname]] <- list()

#List of image names that have been loaded into R.  
dscn <- names(ptsTrainList[[toothname]])

for (d in dscn){print(d)
  test <- efourier(Out(ptsTrainList[[toothname]][[d]]), nb.h = numHarm, start = TRUE)
  numPoints <- 200
  ptsList <- list()
    efa <- list(an = test$coe[j,grep("A",colnames(test$coe))],
                bn = test$coe[j, grep("B", colnames(test$coe))],
                cn = test$coe[j, grep("C", colnames(test$coe))],
                dn = test$coe[j, grep("D", colnames(test$coe))],
                a0 = 0, c0= 0)
    #this is a TEST!!!
    if (efa$dn[1] < 0){
      efa$dn <- -1*efa$dn
      efa$cn <- -1*efa$cn
    }
    ptsTrainListNorm[[toothname]][[d]] <- efourier_i(efa, nb.pts = numPoints)
}

}


#################################################
#Testing out some stuff down here: 
#################################################

toothname <- "LM1"
for (d in names(ptsTrainListNorm[[toothname]])){print(d)
plot(ptsTrainListNorm[[toothname]][[d]],type="l", xlab = "x", ylab = "y", main = d)
}

for (d in names(imgsListMean)){print(d)
  plot(imgsListMean[[d]],type="l", xlab = "x", ylab = "y", main = d)
}

plot(imgsListMean[[d]])
temp<-imgsListMean[[d]]
test <- efourier(Out(temp), nb.h = 8)
test2 <- efourier(Out(cbind(temp[,1],-temp[,2])), nb.h = 8)

plot(Out(cbind(temp[,1],-temp[,2]))[[1]][[1]])
plot(Out(temp)[[1]][[1]])

test <- efourier(Out(temp), nb.h = 8)
test <- efourier(Out(cbind(temp[,1],-temp[,2])), nb.h = 8)

numPoints <- 200
efa <- list(an = test$coe[j,grep("A",colnames(test$coe))],
            bn = test$coe[j, grep("B", colnames(test$coe))],
            cn = test$coe[j, grep("C", colnames(test$coe))],
            dn = test$coe[j, grep("D", colnames(test$coe))],
            a0 = 0, c0= 0)
if (efa$dn[1] < 0){
  efa$dn <- -1*efa$dn
  efa$cn <- -1*efa$cn
}
print(efa)
greg <- efourier_i(efa, nb.pts = numPoints)
plot(greg,type="l")

################################################################################################
#Juliet dissertation data: ptsTrainListNorm
#Mech Turk data: imgsListMean
################################################################################################
save(imgsListMean, file = "/Users/gregorymatthews/Dropbox/brophyTeeth/imgsListMean.RData")
save(ptsTrainListNorm, file = "/Users/gregorymatthews/Dropbox/brophyTeeth/ptsTrainListNorm.RData")

plot(imgsListMean[[1]])
plot(ptsTrainListNorm[["LM1"]][[1]])
plot(ptsTrainList[["LM1"]][[1]])



