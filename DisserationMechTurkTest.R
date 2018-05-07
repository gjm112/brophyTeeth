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
plot(imgsListMean[[3]])
plot(ptsTrainListNorm[["LM1"]][[1]])
plot(ptsTrainList[["LM1"]][[1]])




#Note: 
##Before doing efourier and training the models and doing prediction.
##I need to rotate all shapes in training and test using procrustes so that they all have the same alignment.  



#This has been run
# 
# #Import MTurk results
# folders<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/results")
# 
# for (j in folders){print(j)
# #Pick out one folder
# folderTemp <- paste("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/results/",j,sep="")
# fff<-list.files(folderTemp,full=TRUE)
# 
# #Shell script
# #ImageMagick
# system(paste("cd ",folderTemp,'; for i in *.png ; do convert "$i" -background white -alpha remove "$i" ; done',sep=" "))
# #Convert to .jpg
# system(paste("cd ",folderTemp,'; for i in *.png ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
# system(paste("cd ",folderTemp,'; for i in *.jpeg ; do convert "$i" "${i%.*}.jpg" ; done',sep=" "))
# #Convert to gray scale
# system(paste("cd ",folderTemp,'; for i in *.jpg; do convert "$i" -colorspace Gray "$i"; done',sep=" "))
# #trim
# system(paste("cd ",folderTemp,'; for i in *.jpg; do convert -trim "$i" "$i"; done',sep=" "))
# }

################################################
#Links DSCN numbers to Taxa
################################################
require(gdata)
keys<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/Extant teeth Images and Key/Key to bovid teeth/",full=TRUE)
keys<-keys[grep(".csv",keys)]
keysList<-list()
for (k in keys){
  temp<-read.csv(k)
  if(names(temp)[3]=="sub"){
    temp<-subset(temp,select=-sub)
  }
  names(temp)[1:10]<-c("Tribe","Genus","Species","Institution","No","Origin","Element","M3","M2","M1")
  temp<-apply(temp,2,as.character)
  temp<-temp[temp[,"Tribe"]!="",]
  keysList[[k]]<-temp[,1:10]
}

keyFile<-as.data.frame(do.call(rbind,keysList))
keyFile$M3<-paste0("dscn",as.character(keyFile$M3))
keyFile$M2<-paste0("dscn",as.character(keyFile$M2))
keyFile$M1<-paste0("dscn",as.character(keyFile$M1))

keyFile$M1[keyFile$M1=="dscn"]<-NA
keyFile$M2[keyFile$M2=="dscn"]<-NA
keyFile$M3[keyFile$M3=="dscn"]<-NA

keyFile$M1<-toupper(keyFile$M1)
keyFile$M2<-toupper(keyFile$M2)
keyFile$M3<-toupper(keyFile$M3)
################################################
#test Data from mechanical turk.  
#This is the data to be classified
################################################

imgsList<-list()
dscn<-list.files("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/resultsClean")

for (d in dscn){print(d)
  
  files<-list.files(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/resultsClean/",d))
  files <- files[grep("jpg",files)]
  
  imgsList[[d]]<-import_jpg(paste0("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/resultsClean/",d,"/",files,sep=""))
  
  #Manual removal of this MTurk image
  imgsList[[d]][which(names(imgsList[[d]])=="3A1COHJ8NJU7FB7Z1T98O6918IG8H5")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3C8HJ7UOP7T8RL9X1GPYTVE1ODMZMA")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3R0T90IZ1SBVRI21YZ7V5STJJ9PGCL")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3NJM2BJS4W514VV01IXIZ17BK0BCPK")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="37C0GNLMHF23ZHJ9MITKD7YCASHD6T")]<-NULL
  imgsList[[d]][which(names(imgsList[[d]])=="3N1FSUEFL5ZPKUFV3U05G9EYE454DZ")]<-NULL
  
  #Remove images with too few observations.  These are clearly mistakes.  
  
  remove<-names(which(unlist(lapply(imgsList[[d]],nrow))<=50))
  if (length(remove)>0){
    for (q in remove){
      imgsList[[d]][[q]]<-NULL
    }
  }
  
  
}

dscn<-names(imgsList)
#Do procrustes and then find the average tooth.  
imgsListMean<-list()
for (d in dscn){print(d)
  
  test <- efourier(Out(imgsList[[d]]),nb.h=numHarm)
  numPoints<-200
  ptsList<-list()
  arr<-array(NA, dim=c(numPoints, 2, dim(test$coe)[1]))
  for (j in 1:dim(test$coe)[1]){
    #arr[,,j]<-ptsList[[j]]<-efourier_shape(test$coe[j,grep("A",colnames(test$coe))],test$coe[j,grep("B",colnames(test$coe))],test$coe[j,grep("C",colnames(test$coe))],test$coe[j,grep("D",colnames(test$coe))],plot=TRUE,nb.pts=numPoints)
    efa<-list(an=test$coe[j,grep("A",colnames(test$coe))],
              bn=test$coe[j,grep("B",colnames(test$coe))],
              cn=test$coe[j,grep("C",colnames(test$coe))],
              dn=test$coe[j,grep("D",colnames(test$coe))],
              a0=0,c0=0)
    arr[,,j]<-ptsList[[j]]<-efourier_i(efa,nb.pts=numPoints)
  }
  
  #set of shapes to be predicted.  
  imgsListMean[[d]]<-mshape<-as.matrix(procGPA(arr)$mshape)
  
}



#save(imgsList,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsList.RData")
#load("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsList.RData")
#save(imgsListMean,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsListMean.RData")
#load("/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/imgsListMean.RData")
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


set.seed(1234)
outList <- list()
ptsList <- mod <- modSpecies <- list()
predTempProbsSpecies <- list()
rm(tooth)
#Remove observatiobs from the training set that are in the test set.  
#No UM2 in the test data set.  
for (tooth in c("LM1","LM2","LM3","UM1","UM3")){print(tooth)
  #Align the training and test sets using procrustes so that their Fourier coefficients all have the same interpretation
  numPoints<-200
  modSpecies[[tooth]] <- list()
  
  
  #Converting all of the training data set to have 60 points.  
  #Note that the efourier_i function is modified to remove rotation.  
  test <- efourier(Out(ptsTrainList[[tooth]]),nb.h=20)
  for (j in 1:length(ptsTrainList[[tooth]])){print(j)
    efa<-list(an=test$coe[j,grep("A",colnames(test$coe))],
              bn=test$coe[j,grep("B",colnames(test$coe))],
              cn=test$coe[j,grep("C",colnames(test$coe))],
              dn=test$coe[j,grep("D",colnames(test$coe))],
              a0=0,c0=0)
    ptsTrainList[[tooth]][[j]]<-efourier_i(efa,nb.pts = 200,nb.h=20)
  }
  
  #pull out only the worn images for a specific target tooth
  imgsListMeanTemp<-imgsListMean[which(toupper(names(imgsListMean))%in%keyFile[[substring(tooth,2,3)]])]
  
  
  # for (j in 1:length(imgsListMeanTemp)){print(j)
  #   efa <- efourier(imgsListMeanTemp[[j]],nb.h=20)
  #   imgsListMeanTemp[[j]]<-efourier_i(efa,nb.pts = 200,nb.h=20)
  # }
  
  
  
  arr<-array(NA, dim=c(numPoints, 2, length(ptsTrainList[[tooth]])+length(imgsListMeanTemp)))
  for (j in 1:length(ptsTrainList[[tooth]])){
    arr[,,j]<-ptsTrainList[[tooth]][[j]]
  }
  
  for (j in 1:length(imgsListMeanTemp)){
    arr[,,j+length(ptsTrainList[[tooth]])]<-imgsListMeanTemp[[j]]
  }
  
  
  # ptsList[[tooth]]<-list()
  # test<-efourier(Out(ptsTrainList[[tooth]]),nb.h=20)
  # for (j in 1:length(ptsTrainList[[tooth]])){
  #   efa<-list(an=test$coe[j,grep("A",colnames(test$coe))],
  #        bn=test$coe[j,grep("B",colnames(test$coe))],
  #        cn=test$coe[j,grep("C",colnames(test$coe))],
  #        dn=test$coe[j,grep("D",colnames(test$coe))],
  #        a0=0,c0=0)
  # 
  # arr[,,j]<-ptsList[[tooth]][[j]]<-efourier_i(efa,nb.pts=numPoints,nb.h=20)
  # }
  
  
  procrusted<-procGPA(arr)
  #plot(procrusted$rotated[,,1])
  #for (greg in 1:dim(procrusted$rotated)[3]){
  #  points(procrusted$rotated[,,greg],type="l")
  #}
  
  outCooTrain<-Out(procrusted$rotated[,,1:length(ptsTrainList[[tooth]])])
  outCooTest<-Out(procrusted$rotated[,,-c(1:length(ptsTrainList[[tooth]]))])
  
  efTrain<-as.data.frame(efourier(outCooTrain,nb.h=numHarm)$coe)
  efTest<-as.data.frame(efourier(outCooTest,nb.h=numHarm)$coe)
  
  
  pc<-princomp(rbind(efTrain,efTest))
  efTrain<-cbind(efTrain,pc$scores[1:nrow(efTrain),1:20])
  efTest<-cbind(efTest,pc$scores[-c(1:nrow(efTrain)),1:20])
  
  efTrain$ref<-names(ptsTrainList[[tooth]])
  efTest$ref<-toupper(names(imgsListMeanTemp))
  
  efTrain<-merge(efTrain,refFile,by.x="ref",by.y="ref",all.x=TRUE)
  
  #####################
  #Train the model
  #####################
  form<-as.formula(paste0("tribe~",paste(c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20)),collapse="+")))
  #Omit PCS
  #form<-as.formula(paste0("tribe~",paste(c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm)),collapse="+")))
  
  #tune.control <- tune.control(cross=6)
  #TribeInnerTune.svm<-tune.svm(y=efTrain$tribe,x=efTrain[,c(paste0("A",1:numHarm),paste0("B",2:numHarm),paste0("C",2:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],gamma=2^seq(-10,-3,length=15),cost=10^seq(-2,1.5,length=15),tune.control=tune.control)
  #print(c(TribeInnerTune.svm$best.parameters$gamma,TribeInnerTune.svm$best.parameters$cost))
  
  mod[[tooth]]<-svmFitted<-randomForest(form,data=efTrain,ntree=5000)
  
  
  formula<-as.formula(paste0("species~",paste(c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20)),collapse="+")))
  tribes <- as.character(sort(unique(efTrain$tribe)))
  for (j in c(1,4:7)){
    temp <- subset(efTrain, efTrain$tribe == tribes[j])
    temp$species<-as.factor(as.character(temp$species))
    
    #sampsize <- round(dim(temp)[1]*sampSizeVecSpecies[m])
    # #################################
    # #Tuning for Species
    # #################################
    # SpecInnerTune.rf <- tune.randomForest(formula,data = temp,ntree = 2000,mtry = seq(10,50,10))
    # print(SpecInnerTune.rf$best.parameters$mtry)
    # tuningParmList[[tooth[m]]][["Species"]][[j]][[i]]<-SpecInnerTune.rf$best.parameters
    
    modSpecies[[tooth]][[j]] <- cart <- randomForest(formula,data=temp,ntree=5000)
    
  }
  
  #################################
  #Fit with best tunes parameters
  #################################
  #svmFitted <- svm(y=efTrain$tribe,x=efTrain[,c(paste0("A",1:numHarm),paste0("B",2:numHarm),paste0("C",2:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],probability=TRUE,gamma = TribeInnerTune.svm$best.parameters$gamma,cost = TribeInnerTune.svm$best.parameters$cost)
  
  #Which observations in the test set are a specific tooth.  
  #ids<-which(toupper(names(imgsListMean))%in%toupper(keyFile$M1) |
  #             toupper(names(imgsListMean))%in%toupper(keyFile$M2) |
  #             toupper(names(imgsListMean))%in%toupper(keyFile$M3)) 
  
  if (substring(tooth,1,1)=="U"){
    ids<-which(toupper(names(imgsListMean))%in%toupper(subset(keyFile,Element=="maxilla")[[substring(tooth,2,3)]]))
  }
  
  if (substring(tooth,1,1)=="L"){
    ids<-which(toupper(names(imgsListMean))%in%toupper(subset(keyFile,Element=="mandible")[[substring(tooth,2,3)]]))
  }
  
  #predict the tribes
  predTempProbs<-predict(svmFitted,efTest[,c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],type="prob")
  predTempProbs<-t(apply(predTempProbs,1,function(x){
    x[x==1]<-0.999;x[x==0]<-0.001
    #Now normalize
    x <- x/sum(x)
    return(x)
  }))
  
  
  #predict the species
  predTempProbsSpecies[[tooth]]<-list()
  for (i in c(1:7)){
    
    if (i %in% c(2:3)){
      predTempProbsSpecies[[tooth]][[i]] <- as.matrix(predTempProbs[,tribes[i]])
    }
    
    if (i %in% c(1,4:7)){
      predTempProbsSpecies[[tooth]][[i]]<-predict(modSpecies[[tooth]][[i]] ,efTest[,c(paste0("A",1:numHarm),paste0("B",1:numHarm),paste0("C",1:numHarm),paste0("D",1:numHarm),paste0("Comp.",1:20))],type="prob")
      for (q in 1:nrow(predTempProbsSpecies[[tooth]][[i]])){  
        predTempProbsSpecies[[tooth]][[i]][q,] <- predTempProbsSpecies[[tooth]][[i]][q,] * predTempProbs[q,tribes[i]]
      }
    }
  }
  
  
  predSpecies <- do.call(cbind,predTempProbsSpecies[[tooth]])
  
  outList[[tooth]]<-data.frame(ref=toupper(names(imgsListMeanTemp)),predTempProbs,predSpecies)
  
}

######################################################
#Pull in ID's on the chipps vs worn teeth.  
######################################################
fils<-toupper(list.files("/Users/gregorymatthews/Dropbox/Wear Test/chipped/"))
fils<-fils[grep("JPG",fils)]
filsChipped<-substring(fils,1,nchar(fils)-4)

fils<-c(toupper(list.files("/Users/gregorymatthews/Dropbox/Wear Test/notchipped/")),list.files("/Users/gregorymatthews/Dropbox/Wear Test/notchipped/Incorrect/"))
fils<-fils[grep("JPG",fils)]
filsNotChipped<-substring(fils,1,nchar(fils)-4)

######################################################################
#Now pull out only the teeth we need from each
#Prediction are made on all test teeth but we only need the ones for that specific model
outList[["LM1"]]<-subset(merge(outList[["LM1"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M1")],by.x="ref",by.y="M1",all.x=TRUE),Element=="mandible")
outList[["LM2"]]<-subset(merge(outList[["LM2"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M2")],by.x="ref",by.y="M2",all.x=TRUE),Element=="mandible")
outList[["LM3"]]<-subset(merge(outList[["LM3"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M3")],by.x="ref",by.y="M3",all.x=TRUE),Element=="mandible")

outList[["UM1"]]<-subset(merge(outList[["UM1"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M1")],by.x="ref",by.y="M1",all.x=TRUE),Element=="maxilla")
#outList[["UM2"]]<-subset(merge(outList[["UM2"]][,1:8],keyFile[,c("Tribe","Genus","Species","Element","M2")],by.x="ref",by.y="M2",all.x=TRUE),Element=="mandible")
outList[["UM3"]]<-subset(merge(outList[["UM3"]][,1:28],keyFile[,c("Tribe","Genus","Species","Element","M3")],by.x="ref",by.y="M3",all.x=TRUE),Element=="maxilla")


outList[["LM1"]]$tooth <- "LM1"
outList[["LM2"]]$tooth <- "LM2"
outList[["LM3"]]$tooth <- "LM3"

outList[["UM1"]]$tooth <- "UM1"
#outList[["UM2"]]$tooth <- "UM2"
outList[["UM3"]]$tooth <- "UM3"

preds<-do.call(rbind,outList)

names(preds)[13:14] <- c("Antidorcas.marsupialis","Syncerus.cafer")

names(preds)[9:28] <- do.call(rbind,strsplit(names(preds)[9:28],"[.]"))[,2]

preds$Species <- gsub(" ","",preds$Species)

################################################
#Merge on indicator for chipped vs not chipped.
################################################
preds$Chipped<-preds$ref%in%filsChipped+0
preds$Chipped[preds$ref%in%filsChipped]<-"Chipped"
preds$Chipped[!preds$ref%in%filsChipped]<-"Not Chipped"


preds$probTribe <- NA
for (tr in as.character(unique(preds$Tribe))){
  preds$probTribe[preds$Tribe==tr] <- preds[[tr]][preds$Tribe==tr]
}

preds$probSpecies <- NA
for (sp in as.character(unique(preds$Species))){
  preds$probSpecies[preds$Species==sp] <- preds[[sp]][preds$Species==sp]
}


# write.csv(preds,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/results20180208_RandomForest_withSpecies.csv")
# save(preds,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/preds_results20180208_RandomForest_withSpecies.RData")
# save(mod,file="/Users/gregorymatthews/Dropbox/brophyTeeth/MTurk-WearTest/mod_results20180208_RandomForest_withSpecies.RData")



temp<-subset(preds,Tribe=="Alcelaphini")
boxplot(temp$Alcelaphini~temp$Chipped)
library(dplyr)
group_by(preds,tooth,Tribe,Chipped)  %>% summarize(mn=mean(prob),n=n())
group_by(preds,tooth,Chipped)  %>% summarize(mn=mean(prob),n=n())

boxplot(preds$prob~preds$tooth)
boxplot(preds$prob~preds$Tribe)

preds$correct<-(preds$predictedTribe==preds$Tribe)+0



ggplot(data=preds) + geom_boxplot(aes(x=Tribe,y=prob)) + geom_point(aes(y=prob,x=Tribe,colour=as.factor(correct))) 
ggplot(data=preds) + geom_boxplot(aes(x=Tribe,y=prob)) + geom_point(aes(y=prob,x=Tribe,colour=as.factor(correct))) + facet_wrap(~Chipped)
ggplot(data=preds) + geom_boxplot(aes(x=predictedTribe,y=prob)) + geom_point(aes(y=prob,x=predictedTribe,colour=as.factor(correct))) 
ggplot(data=preds) + geom_boxplot(aes(x=predictedTribe,y=prob)) + geom_point(aes(y=prob,x=predictedTribe,colour=as.factor(correct))) + facet_wrap(~Chipped)




##############################
#Overall
##############################
#which.max(preds[1,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")])
preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]
table(preds$predictedTribe,preds$Tribe)
confusionMatrix(preds$predictedTribe,preds$Tribe)


table(preds$predictedTribe[preds$tooth=="LM1"],preds$Tribe[preds$tooth=="LM1"])
table(preds$predictedTribe[preds$tooth=="LM2"],preds$Tribe[preds$tooth=="LM2"])
table(preds$predictedTribe[preds$tooth=="LM3"],preds$Tribe[preds$tooth=="LM3"])
table(preds$predictedTribe[preds$tooth=="UM1"],preds$Tribe[preds$tooth=="UM1"])
table(preds$predictedTribe[preds$tooth=="UM3"],preds$Tribe[preds$tooth=="UM3"])

##############################
#Chipped
##############################
#which.max(preds[1,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")])
preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]

table(subset(preds,Chipped==1)$predictedTribe,subset(preds,Chipped==1)$Tribe)
confusionMatrix(subset(preds,Chipped==1)$predictedTribe,subset(preds,Chipped==1)$Tribe)


##############################
#Not Chipped
##############################
#which.max(preds[1,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")])
preds$predictedTribe<-c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")[apply(preds[,c("Alcelaphini","Antilopini","Bovini","Hippotragini","Neotragini","Reduncini","Tragelaphini")],1,which.max)]

table(subset(preds,Chipped==0)$predictedTribe,subset(preds,Chipped==0)$Tribe)
confusionMatrix(subset(preds,Chipped==0)$predictedTribe,subset(preds,Chipped==0)$Tribe)


