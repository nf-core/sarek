# ASCAT 2.5
# author: Peter Van Loo
# PCF and ASPCF: Gro Nilsen
# GC correction: Jiqiu Cheng
# ASmultiPCF: Edith Ross

#' @title ascat.loadData
#' @description Function to read in SNP array data
#' @details germline data files can be NULL - in that case these are not read in
#' @param Tumor_LogR_file file containing logR of tumour sample(s)
#' @param Tumor_BAF_file file containing BAF of tumour sample(s)
#' @param Germline_LogR_file file containing logR of germline sample(s), NULL
#' @param Germline_BAF_file file containing BAF of germline sample(s), NULL
#' @param chrs a vector containing the names for the chromosomes (e.g. c(1:22,"X"))
#' @param gender a vector of gender for each cases ("XX" or "XY"). Default = all female ("XX")
#' @param sexchromosomes a vector containing the names for the sex chromosomes
#'
#' @return ascat data structure containing:\cr
#' 1. Tumor_LogR data matrix\cr
#' 2. Tumor_BAF data matrix\cr
#' 3. Tumor_LogR_segmented: placeholder, NULL\cr
#' 4. Tumor_BAF_segmented: placeholder, NULL\cr
#' 5. Germline_LogR data matrix\cr
#' 6. Germline_BAF data matrix\cr
#' 7. SNPpos: position of all SNPs\cr
#' 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]],] will output the Tumor_LogR data of chromosome 13\cr
#' 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)\cr
#' 10. gender: a vector of gender for each cases ("XX" or "XY"). Default = NULL: all female ("XX")\cr
#'
#' @export
#'
ascat.loadData = function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, Germline_BAF_file = NULL, chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y")) {

  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)

  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA

  Germline_LogR = NULL
  Germline_BAF = NULL
  if(!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t", check.names=F)

    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }

  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X,Y (or whatever is given in the input value of chrs)
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]

  # if some chromosomes have no data, just remove them
  chrs = intersect(chrs,unique(SNPpos[,1]))

  Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  # make sure it is all converted to numerical values
  for (cc in 1:dim(Tumor_LogR)[2]) {
    Tumor_LogR[,cc]=as.numeric(as.vector(Tumor_LogR[,cc]))
    Tumor_BAF[,cc]=as.numeric(as.vector(Tumor_BAF[,cc]))
  }

  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos),c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[rownames(SNPpos),c(-1,-2),drop=F]
    for (cc in 1:dim(Germline_LogR)[2]) {
      Germline_LogR[,cc]=as.numeric(as.vector(Germline_LogR[,cc]))
      Germline_BAF[,cc]=as.numeric(as.vector(Germline_BAF[,cc]))
    }
  }

  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]

  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  }

  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chunks...")
  chr = split_genome(SNPpos)

  if (is.null(gender)) {
    gender = rep("XX",dim(Tumor_LogR)[2])
  }
  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
              samples = colnames(Tumor_LogR), gender = gender,
              sexchromosomes = sexchromosomes,
              failedarrays = NULL))
}



#' @title ascat.plotRawData
#' @description Plots SNP array data
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#'
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#'
#' @export
ascat.plotRawData = function(ASCATobj) {
  print.noquote("Plotting tumor data")
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    png(filename = paste(ASCATobj$samples[i],".tumour.png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
    plot(c(1,dim(ASCATobj$Tumor_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[,i],col="red")
    #points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,dim(ASCATobj$Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[,i],col="red")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }

  if(!is.null(ASCATobj$Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(ASCATobj$Germline_LogR)[2]) {
      png(filename = paste(ASCATobj$samples[i],".germline.png",sep=""), width = 2000, height = 1000, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
      plot(c(1,dim(ASCATobj$Germline_LogR)[1]), c(-1,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(ASCATobj$Germline_BAF)[1]), c(0,1), type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_BAF[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
}


#' @title ascat.GCcorrect
#' @description Corrects logR of the tumour sample(s) with genomic GC content
#' @param ASCATobj an ASCAT object
#' @param GCcontentfile File containing the GC content around every SNP for increasing window sizes
#' @details Note that probes not present in the GCcontentfile will be lost from the results
#' @return ASCAT object with corrected tumour logR
#'
#' @export
ascat.GCcorrect = function(ASCATobj, GCcontentfile = NULL) {
  if(is.null(GCcontentfile)) {
    print.noquote("Error: no GC content file given!")
  }
  else {
    GC_newlist<-read.table(file=GCcontentfile,header=TRUE,as.is=TRUE)
    colnames(GC_newlist)[c(1,2)] = c("Chr","Position")
    GC_newlist$Chr<-as.character(GC_newlist$Chr)
    GC_newlist$Position<-as.numeric(as.character(GC_newlist$Position))

    ovl = intersect(row.names(ASCATobj$Tumor_LogR),row.names(GC_newlist))

    GC_newlist<-GC_newlist[ovl,]

    SNPpos = ASCATobj$SNPpos[ovl,]
    Tumor_LogR = ASCATobj$Tumor_LogR[ovl,,drop=F]
    Tumor_BAF = ASCATobj$Tumor_BAF[ovl,,drop=F]

    chrs = intersect(ASCATobj$chrs,unique(SNPpos[,1]))

    Germline_LogR = NULL
    Germline_BAF = NULL
    if(!is.null(ASCATobj$Germline_LogR)) {
      Germline_LogR = ASCATobj$Germline_LogR[ovl,,drop=F]
      Germline_BAF = ASCATobj$Germline_BAF[ovl,,drop=F]
    }

    last = 0;
    ch = list();
    for (i in 1:length(ASCATobj$chrs)) {
      chrke = SNPpos[SNPpos[,1]==ASCATobj$chrs[i],]
      chrpos = chrke[,2]
      names(chrpos) = rownames(chrke)
      chrpos = sort(chrpos)
      ch[[i]] = (last+1):(last+length(chrpos))
      last = last+length(chrpos)
    }

    for (s in 1:length(ASCATobj$samples)) {
      print.noquote(paste("Sample ", ASCATobj$samples[s], " (",s,"/",length(ASCATobj$samples),")",sep=""))
      Tumordata = Tumor_LogR[,s]
      names(Tumordata) = rownames(Tumor_LogR)

      # Calculate weighted correlation
      length_tot<-NULL
      corr_tot<-NULL
      for(chrindex in unique(SNPpos[,1])) {
        GC_newlist_chr<-GC_newlist[GC_newlist$Chr==chrindex,]
        td_chr<-Tumordata[GC_newlist$Chr==chrindex]

        flag_nona<-(complete.cases(td_chr) & complete.cases(GC_newlist_chr))

        #only work with chromosomes that have variance
        if(length(td_chr[flag_nona])>0 & var(td_chr[flag_nona])>0){
          corr<-cor(GC_newlist_chr[flag_nona,3:ncol(GC_newlist_chr)],td_chr[flag_nona])
          corr_tot<-cbind(corr_tot,corr)
          length_tot<-c(length_tot,length(td_chr))
        }
      }
      corr<-apply(corr_tot,1,function(x) sum(abs(x*length_tot))/sum(length_tot))
      index_1M<-c(which(names(corr)=="X1M"),which(names(corr)=="X1Mb"))
      maxGCcol_short<-which(corr[1:(index_1M-1)]==max(corr[1:(index_1M-1)]))
      maxGCcol_long<-which(corr[index_1M:length(corr)]==max(corr[index_1M:length(corr)]))
      maxGCcol_long<-(maxGCcol_long+(index_1M-1))

      cat("weighted correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")
      cat("Short window size: ",names(GC_newlist)[maxGCcol_short+2],"\n")
      cat("Long window size: ",names(GC_newlist)[maxGCcol_long+2],"\n")

      # Multiple regression
      flag_NA<-(is.na(Tumordata))|(is.na(GC_newlist[,2+maxGCcol_short]))|(is.na(GC_newlist[,2+maxGCcol_long]))
      td_select<-Tumordata[!flag_NA]
      GC_newlist_select <- GC_newlist[!flag_NA,]
      x1<-GC_newlist_select[,2+maxGCcol_short]
      x2<-GC_newlist_select[,2+maxGCcol_long]
      x3<-(x1)^2
      x4<-(x2)^2
      model<-lm(td_select~x1+x2+x3+x4,y=TRUE)

      GCcorrected<-Tumordata
      GCcorrected[]<-NA
      GCcorrected[!flag_NA] <- model$residuals

      Tumor_LogR[,s] = GCcorrected

      chr = split_genome(SNPpos)
    }

    # add some plotting code for each sample while it is generated!!!!

    return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
                Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
                Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
                SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
                samples = colnames(Tumor_LogR), gender = ASCATobj$gender,
                sexchromosomes = ASCATobj$sexchromosomes))
  }
}


#' @title ascat.aspcf
#' @description run ASPCF segmentation
#' @details This function can be easily parallelised by controlling the selectsamples parameter\cr
#' it saves the results in LogR_PCFed[sample]_[segment].txt and BAF_PCFed[sample]_[segment].txt
#' @param ASCATobj an ASCAT object
#' @param selectsamples a vector containing the sample number(s) to PCF. Default = all
#' @param ascat.gg germline genotypes (NULL if germline data is available)
#' @param penalty penalty of introducing an additional ASPCF breakpoint (expert parameter, don't adapt unless you know what you're doing)
#'
#' @return output: ascat data structure containing:\cr
#' 1. Tumor_LogR data matrix\cr
#' 2. Tumor_BAF data matrix\cr
#' 3. Tumor_LogR_segmented: matrix of LogR segmented values\cr
#' 4. Tumor_BAF_segmented: list of BAF segmented values; each element in the list is a matrix containing the segmented values for one sample (only for probes that are germline homozygous)\cr
#' 5. Germline_LogR data matrix\cr
#' 6. Germline_BAF data matrix\cr
#' 7. SNPpos: position of all SNPs\cr
#' 8. ch: a list containing vectors with the indices for each chromosome (e.g. Tumor_LogR[ch[[13]],] will output the Tumor_LogR data of chromosome 13\cr
#' 9. chr: a list containing vectors with the indices for each distinct part that can be segmented separately (e.g. chromosome arm, stretch of DNA between gaps in the array design)\cr
#'
#' @export
#'
ascat.aspcf = function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg = NULL, penalty = 25) {
  #first, set germline genotypes
  gg = NULL
  if(!is.null(ascat.gg)) {
    gg = ascat.gg$germlinegenotypes
  }
  else {
    gg = ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 0.7
  }
  # calculate germline homozygous stretches for later resegmentation
  ghs = predictGermlineHomozygousStretches(ASCATobj$chr, gg)

  segmentlengths = unique(c(penalty,25,50,100,200,400,800))
  segmentlengths = segmentlengths[segmentlengths>=penalty]

  Tumor_LogR_segmented = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  rownames(Tumor_LogR_segmented) = rownames(ASCATobj$Tumor_LogR)
  colnames(Tumor_LogR_segmented) = colnames(ASCATobj$Tumor_LogR)
  Tumor_BAF_segmented = list();
  for (sample in selectsamples) {
    print.noquote(paste("Sample ", ASCATobj$samples[sample], " (",sample,"/",length(ASCATobj$samples),")",sep=""))
    logrfilename = paste(ASCATobj$samples[sample],".LogR.PCFed.txt",sep="")
    baffilename = paste(ASCATobj$samples[sample],".BAF.PCFed.txt",sep="")
    logRPCFed = numeric(0)
    bafPCFed = numeric(0)
    for (segmentlength in segmentlengths) {
      logRPCFed = numeric(0)
      bafPCFed = numeric(0)
      tbsam = ASCATobj$Tumor_BAF[,sample]
      names(tbsam) = rownames(ASCATobj$Tumor_BAF)
      homosam = gg[,sample]
      for (chrke in 1:length(ASCATobj$chr)) {
        lr = ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]],sample]
        #winsorize to remove outliers
        #this has a problem with NAs
        lrwins = vector(mode="numeric",length=length(lr))
        lrwins[is.na(lr)] = NA
        lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)],2.5,25)$ywin
        baf = tbsam[ASCATobj$chr[[chrke]]]
        homo = homosam[ASCATobj$chr[[chrke]]]
        Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)
        bafsel = baf[Select_het]
        # winsorize BAF as well (as a safeguard)
        bafselwinsmirrored = madWins(ifelse(bafsel>0.5,bafsel,1-bafsel),2.5,25)$ywin
        bafselwins = ifelse(bafsel>0.5,bafselwinsmirrored,1-bafselwinsmirrored)
        indices = which(Select_het)
        logRaveraged = NULL;
        if(length(indices)!=0) {
          averageIndices = c(1,(indices[1:(length(indices)-1)]+indices[2:length(indices)])/2,length(lr)+0.01)
          startindices = ceiling(averageIndices[1:(length(averageIndices)-1)])
          endindices = floor(averageIndices[2:length(averageIndices)]-0.01)
          if(length(indices)==1) {
            startindices = 1
            endindices = length(lr)
          }
          nrIndices = endindices - startindices + 1
          logRaveraged = vector(mode="numeric",length=length(indices))
          for(i in 1:length(indices)) {
            if(is.na(endindices[i])) {
              endindices[i]=startindices[i]
            }
            logRaveraged[i]=mean(lrwins[startindices[i]:endindices[i]], na.rm=T)
          }
        }
        # if there are no probes in the segment (after germline homozygous removal), don't do anything, except add a LogR segment
        if(length(logRaveraged)>0) {
          logRASPCF = NULL
          bafASPCF = NULL
          if(length(logRaveraged)<6) {
            logRASPCF = rep(mean(logRaveraged),length(logRaveraged))
            bafASPCF = rep(mean(bafselwins),length(logRaveraged))
          }
          else {
            PCFed = fastAspcf(logRaveraged,bafselwins,6,segmentlength)
            logRASPCF = PCFed$yhat1
            bafASPCF = PCFed$yhat2
          }
          names(bafASPCF)=names(indices)
          logRc = numeric(0)
          for(probe in 1:length(logRASPCF)) {
            if(probe == 1) {
              logRc = rep(logRASPCF[probe],indices[probe])
            }
            # if probe is 1, set the beginning, and let the loop go:
            if(probe == length(logRASPCF)) {
              logRc = c(logRc,rep(logRASPCF[probe],length(lr)-indices[probe]))
            }
            else if(logRASPCF[probe]==logRASPCF[probe+1]) {
              logRc = c(logRc, rep(logRASPCF[probe],indices[probe+1]-indices[probe]))
            }
            else {
              #find best breakpoint
              d = numeric(0)
              totall = indices[probe+1]-indices[probe]
              for (bp in 0:(totall-1)) {
                dis = sum(abs(lr[(1:bp)+indices[probe]]-logRASPCF[probe]), na.rm=T)
                if(bp!=totall) {
                  dis = sum(dis, sum(abs(lr[((bp+1):totall)+indices[probe]]-logRASPCF[probe+1]), na.rm=T), na.rm=T)
                }
                d = c(d,dis)
              }
              breakpoint = which.min(d)-1
              logRc = c(logRc,rep(logRASPCF[probe],breakpoint),rep(logRASPCF[probe+1],totall-breakpoint))
            }
          }
          #2nd step: adapt levels!
          logRd = numeric(0)
          seg = rle(logRc)$lengths
          startprobe = 1
          endprobe = 0
          for (i in 1:length(seg)) {
            endprobe = endprobe+seg[i]
            level = mean(lr[startprobe:endprobe], na.rm=T)
            logRd = c(logRd, rep(level,seg[i]))
            startprobe = startprobe + seg[i]
          }
          logRPCFed = c(logRPCFed,logRd)
          bafPCFed = c(bafPCFed,bafASPCF)
        }
        # add a LogR segment
        else {
          level = mean(lr,na.rm=T)
          reps = length(lr)
          logRPCFed = c(logRPCFed,rep(level,reps))
        }
        # correct wrong segments in germline homozygous stretches:
        homsegs = ghs[[sample]][ghs[[sample]][,1]==chrke,]
        startchr = min(ASCATobj$chr[[chrke]])
        endchr = max(ASCATobj$chr[[chrke]])
        # to solve an annoying error when homsegs has length 1:
        if(length(homsegs)==3) {
          homsegs=t(as.matrix(homsegs))
        }
        if(!is.null(homsegs)&&!is.na(homsegs)&&dim(homsegs)[1]!=0) {
          for (i in 1:dim(homsegs)[1]) {
            # note that only the germline homozygous segment is resegmented, plus a bit extra (but this is NOT replaced)
            startpos = max(homsegs[i,2],startchr)
            endpos = min(homsegs[i,3],endchr)
            # PCF over a larger fragment
            startpos2 = max(homsegs[i,2]-100,startchr)
            endpos2 = min(homsegs[i,3]+100,endchr)
            # take into account a little extra (difference between startpos2 and startpos3 is not changed)
            startpos3 = max(homsegs[i,2]-5,startchr)
            endpos3 = min(homsegs[i,3]+5,endchr)
            # note that the parameters are arbitrary, but <100 seems to work on the ERBB2 example!
            # segmentlength is lower here, as in the full data, noise on LogR is higher!
            # do this on Winsorized data too!
            towins = ASCATobj$Tumor_LogR[startpos2:endpos2,sample]
            winsed = madWins(towins[!is.na(towins)],2.5,25)$ywin
            pcfed = vector(mode="numeric",length=length(towins))
            pcfed[!is.na(towins)] = exactPcf(winsed,6,floor(segmentlength/4))
            pcfed2 = pcfed[(startpos3-startpos2+1):(endpos3-startpos2+1)]
            dif = abs(pcfed2-logRPCFed[startpos3:endpos3])
            #0.3 is hardcoded here, in order not to have too many segments!
            #only replace if enough probes differ (in order not to get singular probes with different breakpoints)
            if(!is.na(dif)&&sum(dif>0.3)>5) {
              #replace a bit more to make sure no 'lone' probes are left (startpos3 instead of startpos)
              logRPCFed[startpos3:endpos3]=ifelse(dif>0.3,pcfed2,logRPCFed[startpos3:endpos3])
            }
          }
        }
      }
      #fill in NAs (otherwise they cause problems):
      #some NA probes are filled in with zero, replace those too:
      logRPCFed = fillNA(logRPCFed, zeroIsNA=TRUE)

      #adapt levels again
      seg = rle(logRPCFed)$lengths
      logRPCFed = numeric(0)
      startprobe = 1
      endprobe = 0
      prevlevel = 0
      for (i in 1:length(seg)) {
        endprobe = endprobe+seg[i]
        level = mean(ASCATobj$Tumor_LogR[startprobe:endprobe,sample], na.rm=T)
        #making sure no NA's get filled in...
        if(is.nan(level)) {
          level=prevlevel
        }
        else {
          prevlevel=level
        }
        logRPCFed = c(logRPCFed, rep(level,seg[i]))
        startprobe = startprobe + seg[i]
      }
      #put in names and write results to files
      names(logRPCFed) = rownames(ASCATobj$Tumor_LogR)

      # if less than 800 segments: this segmentlength is ok, otherwise, rerun with higher segmentlength
      if(length(unique(logRPCFed))<800) {
        break
      }
    }

    write.table(logRPCFed,logrfilename,sep="\t",col.names=F)
    write.table(bafPCFed,baffilename,sep="\t",col.names=F)
    bafPCFed = as.matrix(bafPCFed)
    Tumor_LogR_segmented[,sample] = logRPCFed
    Tumor_BAF_segmented[[sample]] = 1-bafPCFed
  }

  ASCATobj = list(Tumor_LogR = ASCATobj$Tumor_LogR,
                  Tumor_BAF = ASCATobj$Tumor_BAF,
                  Tumor_LogR_segmented = Tumor_LogR_segmented,
                  Tumor_BAF_segmented = Tumor_BAF_segmented,
                  Germline_LogR = ASCATobj$Germline_LogR,
                  Germline_BAF = ASCATobj$Germline_BAF,
                  SNPpos = ASCATobj$SNPpos,
                  ch = ASCATobj$ch,
                  chr = ASCATobj$chr,
                  chrs = ASCATobj$chrs,
                  samples = colnames(ASCATobj$Tumor_LogR), gender = ASCATobj$gender,
                  sexchromosomes = ASCATobj$sexchromosomes, failedarrays = ascat.gg$failedarrays)
  return(ASCATobj)
}


#' @title ascat.plotSegmentedData
#' @description plots the SNP array data before and after segmentation
#'
#' @param ASCATobj an ASCAT object (e.g. from ascat.aspcf)
#'
#' @return png files showing raw and segmented tumour logR and BAF
#'
#' @export
#'
ascat.plotSegmentedData = function(ASCATobj) {
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    Select_nonNAs = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) = rownames(ASCATobj$Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    png(filename = paste(ASCATobj$samples[arraynr],".ASPCF.png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2)
    r = ASCATobj$Tumor_LogR_segmented[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]),arraynr]
    beta = ASCATobj$Tumor_BAF_segmented[[arraynr]][,,drop=FALSE]
    plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
    points(r,col="green")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]],HetIDs);
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(ASCATobj$Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
    points(beta, col = "green")
    points(1-beta, col = "green")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]],HetIDs);
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }
}



#' @title ascat.runAscat
#' @description ASCAT main function, calculating the allele-specific copy numbers
#' @param ASCATobj an ASCAT object from ascat.aspcf
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100\% aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
#' @param pdfPlot Optional flag if nonrounded plots and ASCAT profile in pdf format are desired. Default=F
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#  @param textFlag Optional flag to add the positions of fragments located outside of the plotting area to the plots. Default=F
#' @param circos Optional file to output the non-rounded values in Circos track format. Default=NA
#' @param rho_manual optional argument to override ASCAT optimization and supply rho parameter (not recommended)
#' @param psi_manual optional argument to override ASCAT optimization and supply psi parameter (not recommended)
#' @details Note: for copy number only probes, nA contains the copy number value and nB = 0.
#' @return an ASCAT output object, containing:\cr
#' 1. nA: copy number of the A allele\cr
#' 2. nB: copy number of the B allele\cr
#' 3. aberrantcellfraction: the aberrant cell fraction of all arrays\cr
#' 4. ploidy: the ploidy of all arrays\cr
#' 5. failedarrays: arrays on which ASCAT analysis failed\cr
#' 6. nonaberrantarrays: arrays on which ASCAT analysis indicates that they show virtually no aberrations\cr
#' 7. segments: an array containing the copy number segments of each sample (not including failed arrays)\cr
#' 8. segments_raw: an array containing the copy number segments of each sample without any rounding applied\cr
#' 9. distance_matrix: distances for a range of ploidy and tumor percentage values
#'
#' @export
#'
ascat.runAscat = function(ASCATobj, gamma = 0.55, pdfPlot = F, y_limit = 5, circos=NA, rho_manual = NA, psi_manual = NA) {
  goodarrays=NULL
  res = vector("list",dim(ASCATobj$Tumor_LogR)[2])
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    print.noquote(paste("Sample ", ASCATobj$samples[arraynr], " (",arraynr,"/",length(ASCATobj$samples),")",sep=""))
    lrr=ASCATobj$Tumor_LogR[,arraynr]
    names(lrr)=rownames(ASCATobj$Tumor_LogR)
    baf=ASCATobj$Tumor_BAF[,arraynr]
    names(baf)=rownames(ASCATobj$Tumor_BAF)
    lrrsegm = ASCATobj$Tumor_LogR_segmented[,arraynr]
    names(lrrsegm) = rownames(ASCATobj$Tumor_LogR_segmented)
    bafsegm = ASCATobj$Tumor_BAF_segmented[[arraynr]][,,drop=FALSE]
    names(bafsegm) = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    failedqualitycheck = F
    if(ASCATobj$samples[arraynr]%in%ASCATobj$failedarrays) {
      failedqualitycheck = T
    }
    ending = ifelse(pdfPlot, "pdf", "png")
    circosName=NA
    if(!is.na(circos)){
      circosName=paste(circos,"_",ASCATobj$samples[arraynr],sep="")
    }
    if(is.na(rho_manual)) {
      res[[arraynr]] = runASCAT(lrr,baf,lrrsegm,bafsegm,ASCATobj$gender[arraynr],ASCATobj$SNPpos,ASCATobj$ch,ASCATobj$chrs,ASCATobj$sexchromosomes, failedqualitycheck,
                                paste(ASCATobj$samples[arraynr],".sunrise.png",sep=""),paste(ASCATobj$samples[arraynr],".ASCATprofile.", ending ,sep=""),
                                paste(ASCATobj$samples[arraynr],".rawprofile.", ending ,sep=""),NA,
                                gamma,NA,NA,pdfPlot, y_limit, circosName)
    } else {
      res[[arraynr]] = runASCAT(lrr,baf,lrrsegm,bafsegm,ASCATobj$gender[arraynr],ASCATobj$SNPpos,ASCATobj$ch,ASCATobj$chrs,ASCATobj$sexchromosomes, failedqualitycheck,
                                paste(ASCATobj$samples[arraynr],".sunrise.png",sep=""),paste(ASCATobj$samples[arraynr],".ASCATprofile.", ending,sep=""),
                                paste(ASCATobj$samples[arraynr],".rawprofile.", ending,sep=""),NA,
                                gamma,rho_manual[arraynr],psi_manual[arraynr], pdfPlot, y_limit, circosName)
    }
    if(!is.na(res[[arraynr]]$rho)) {
      goodarrays[length(goodarrays)+1] = arraynr
    }
  }

  if(length(goodarrays)>0) {
    n1 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    n2 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
    rownames(n1) = rownames(ASCATobj$Tumor_LogR)
    rownames(n2) = rownames(ASCATobj$Tumor_LogR)
    colnames(n1) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
    colnames(n2) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      n1[,i] = res[[goodarrays[i]]]$nA
      n2[,i] = res[[goodarrays[i]]]$nB
    }

    distance_matrix = vector("list",length(goodarrays))
    names(distance_matrix) <- colnames(ASCATobj$Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      distance_matrix[[i]] = res[[goodarrays[i]]]$distance_matrix
    }

    tp = vector(length=length(goodarrays))
    psi = vector(length=length(goodarrays))
    ploidy = vector(length=length(goodarrays))
    goodnessOfFit = vector(length=length(goodarrays))
    naarrays = NULL
    for (i in 1:length(goodarrays)) {
      tp[i] = res[[goodarrays[i]]]$rho
      psi[i] = res[[goodarrays[i]]]$psi
      ploidy[i] = mean(res[[goodarrays[i]]]$nA+res[[goodarrays[i]]]$nB,na.rm=T)
      goodnessOfFit[i] = res[[goodarrays[i]]]$goodnessOfFit
      if(res[[goodarrays[i]]]$nonaberrant) {
        naarrays = c(naarrays,ASCATobj$samples[goodarrays[i]])
      }
    }
    fa = colnames(ASCATobj$Tumor_LogR)[-goodarrays]
    names(tp) = colnames(n1)
    names(ploidy) = colnames(n1)
    names(psi) = colnames(n1)
    names(goodnessOfFit) = colnames(n1)

    seg = NULL
    for (i in 1:length(goodarrays)) {
      segje = res[[goodarrays[i]]]$seg
      seg = rbind(seg,cbind(ASCATobj$samples[goodarrays[i]],as.vector(ASCATobj$SNPpos[segje[,1],1]),
                            ASCATobj$SNPpos[segje[,1],2],
                            ASCATobj$SNPpos[segje[,2],2],segje[,3],segje[,4]))
    }
    colnames(seg) = c("sample","chr","startpos","endpos","nMajor","nMinor")
    seg = data.frame(seg,stringsAsFactors=F)
    seg[,3]=as.numeric(seg[,3])
    seg[,4]=as.numeric(seg[,4])
    seg[,5]=as.numeric(seg[,5])
    seg[,6]=as.numeric(seg[,6])

    seg_raw = NULL
    for (i in 1:length(goodarrays)) {
      segje = res[[goodarrays[i]]]$seg_raw
      seg_raw = rbind(seg_raw,cbind(ASCATobj$samples[goodarrays[i]],as.vector(ASCATobj$SNPpos[segje[,1],1]),
                                    ASCATobj$SNPpos[segje[,1],2],
                                    ASCATobj$SNPpos[segje[,2],2],segje[,3],segje[,4:ncol(segje)]))

    }
    colnames(seg_raw) = c("sample","chr","startpos","endpos","nMajor","nMinor","nAraw","nBraw")

    seg_raw = data.frame(seg_raw,stringsAsFactors=F)
    seg_raw[,3]=as.numeric(seg_raw[,3])
    seg_raw[,4]=as.numeric(seg_raw[,4])
    seg_raw[,5]=as.numeric(seg_raw[,5])
    seg_raw[,6]=as.numeric(seg_raw[,6])
    seg_raw[,7]=as.numeric(seg_raw[,7])
    seg_raw[,8]=as.numeric(seg_raw[,8])

  }
  else {
    n1 = NULL
    n2 = NULL
    tp = NULL
    ploidy = NULL
    psi = NULL
    goodnessOfFit = NULL
    fa = colnames(ASCATobj$Tumor_LogR)
    naarrays = NULL
    seg = NULL
    seg_raw = NULL
    distance_matrix = NULL
  }

  return(list(nA = n1, nB = n2, aberrantcellfraction = tp, ploidy = ploidy, psi = psi, goodnessOfFit = goodnessOfFit,
              failedarrays = fa, nonaberrantarrays = naarrays, segments = seg, segments_raw = seg_raw, distance_matrix = distance_matrix))
}

# helper function to split the genome into parts
split_genome = function(SNPpos) {

  # look for gaps of more than 5Mb (arbitrary treshold to account for big centremeres or other gaps) and chromosome borders
  bigHoles = which(diff(SNPpos[,2])>=5000000)+1
  chrBorders = which(SNPpos[1:(dim(SNPpos)[1]-1),1]!=SNPpos[2:(dim(SNPpos)[1]),1])+1

  holes = unique(sort(c(bigHoles,chrBorders)))

  # find which segments are too small
  #joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  #while (1 %in% joincandidates) {
  #  holes=holes[-1]
  #  joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  #}
  #while ((length(holes)+1) %in% joincandidates) {
  #  holes=holes[-length(holes)]
  #  joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  #}

  #while(length(joincandidates)!=0) {
  # the while loop is because after joining, segments may still be too small..

  #startseg = c(1,holes)
  #endseg = c(holes-1,dim(SNPpos)[1])

  # for each segment that is too short, see if it has the same chromosome as the segments before and after
  # the next always works because neither the first or the last segment is in joincandidates now
  #previoussamechr = SNPpos[endseg[joincandidates-1],1]==SNPpos[startseg[joincandidates],1]
  #nextsamechr = SNPpos[endseg[joincandidates],1]==SNPpos[startseg[joincandidates+1],1]

  #distanceprevious = SNPpos[startseg[joincandidates],2]-SNPpos[endseg[joincandidates-1],2]
  #distancenext = SNPpos[startseg[joincandidates+1],2]-SNPpos[endseg[joincandidates],2]

  # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
  #joins = ifelse(previoussamechr&nextsamechr,
  #               ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
  #               ifelse(nextsamechr, joincandidates, joincandidates-1))

  #holes=holes[-joins]

  #joincandidates=which(diff(c(0,holes,dim(SNPpos)[1]))<200)
  #}
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary, but I'm sure this is no problem.

  startseg = c(1,holes)
  endseg = c(holes-1,dim(SNPpos)[1])

  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }

  return(chr)
}



#' @title predictGermlineHomozygousStretches
#' @description helper function to predict germline homozyguous segments for later re-segmentation
#' @param chr a list containing vectors with the indices for each distinct part that can be segmented separately
#' @param hom germline genotypes
#'
#' @keywords internal
#'
#' @return germline homozyguous segments
#'
#'
predictGermlineHomozygousStretches = function(chr, hom) {

  # contains the result: a list of vectors of probe numbers in homozygous stretches for each sample
  HomoStretches = list()

  for (sam in 1:dim(hom)[2]) {
    homsam = hom[,sam]

    perchom = sum(homsam,na.rm=T)/sum(!is.na(homsam))

    # NOTE THAT A P-VALUE THRESHOLD OF 0.001 IS HARDCODED HERE
    homthres = ceiling(log(0.001,perchom))

    allhprobes = NULL
    for (chrke in 1:length(chr)) {
      hschr = homsam[chr[[chrke]]]

      hprobes = vector(length=0)
      for(probe in 1:length(hschr)) {
        if(!is.na(hschr[probe])) {
          if(hschr[probe]) {
            hprobes = c(hprobes,probe)
          }
          else {
            if(length(hprobes)>=homthres) {
              allhprobes = rbind(allhprobes,c(chrke,chr[[chrke]][min(hprobes)],chr[[chrke]][max(hprobes)]))
            }
            hprobes = vector(length=0)
          }
        }
      }
      # if the last probe is homozygous, this is not yet accounted for
      if(!is.na(hschr[probe]) & hschr[probe]) {
        if(length(hprobes)>=homthres) {
          allhprobes = rbind(allhprobes,c(chrke,chr[[chrke]][min(hprobes)],chr[[chrke]][max(hprobes)]))
        }
      }

    }

    HomoStretches[[sam]]=allhprobes

  }

  return(HomoStretches)
}


#' @title make_segments
#' @description Function to make segments of constant LRR and BAF.\cr
#' This function is more general and does not depend on specific ASPCF output, it can also handle segmention performed on LRR and BAF separately
#' @param r segmented logR
#' @param b segmented BAF
#' @return segments of constant logR and BAF including their lengths
#' @keywords internal
#' @export
make_segments = function(r,b) {
  m = matrix(ncol = 2, nrow = length(b))
  m[,1] = r
  m[,2] = b
  m = as.matrix(na.omit(m))
  pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
  colnames(pcf_segments) = c("r","b","length");
  index = 0;
  previousb = -1;
  previousr = 1E10;
  for (i in 1:dim(m)[1]) {
    if (m[i,2] != previousb || m[i,1] != previousr) {
      index=index+1;
      count=1;
      pcf_segments[index, "r"] = m[i,1];
      pcf_segments[index, "b"] = m[i,2];
    }
    else {
      count = count + 1;
    }
    pcf_segments[index, "length"] = count;
    previousb = m[i,2];
    previousr = m[i,1];
  }
  pcf_segments = as.matrix(na.omit(pcf_segments))[,,drop=FALSE]
  return(pcf_segments);
}



# function to create the distance matrix (distance for a range of ploidy and tumor percentage values)
# input: segmented LRR and BAF and the value for gamma
create_distance_matrix = function(segments, gamma) {
  s = segments
  psi_pos = seq(1,6,0.05)
  rho_pos = seq(0.1,1.05,0.01)
  d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) = psi_pos
  colnames(d) = rho_pos
  dmin = 1E20;
  for(i in 1:length(psi_pos)) {
    psi = psi_pos[i]
    for(j in 1:length(rho_pos)) {
      rho = rho_pos[j]
      nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
      # choose the minor allele
      nMinor = NULL
      if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
        nMinor = nA
      }
      else {
        nMinor = nB
      }
      d[i,j] = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T)
    }
  }
  return(d)
}




#' @title runASCAT
#' @description the ASCAT main function
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param gender a vector of gender for each cases ("XX" or "XY"). Default = NULL: all female ("XX")
#' @param SNPpos position of all SNPs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param chrnames a vector containing the names for the chromosomes (e.g. c(1:22,"X"))
#' @param sexchromosomes a vector containing the names for the sex chromosomes
#' @param failedqualitycheck did the sample fail any previous quality check or not?
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file
#' @param aberrationreliabilitypng aberration reliability score is plotted if filename is given
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100\% aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays)
#' @param rho_manual optional argument to override ASCAT optimization and supply rho parameter (not recommended)
#' @param psi_manual optional argument to override ASCAT optimization and supply psi parameter (not recommended)
#' @param pdfPlot Optional flag if nonrounded plots and ASCAT profile in pdf format are desired. Default=F
#' @param circos Optional file to output the non-rounded values in Circos track format. Default=NA
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#'
#' @keywords internal
#'
#' @return list containing optimal purity and ploidy
#'
#' @import RColorBrewer
#'
#' @export
runASCAT = function(lrr, baf, lrrsegmented, bafsegmented, gender, SNPpos, chromosomes, chrnames, sexchromosomes, failedqualitycheck = F,
                    distancepng = NA, copynumberprofilespng = NA, nonroundedprofilepng = NA, aberrationreliabilitypng = NA, gamma = 0.55,
                    rho_manual = NA, psi_manual = NA, pdfPlot = F, y_limit = 5, circos=NA) {
  ch = chromosomes
  chrs = chrnames
  b = bafsegmented
  r = lrrsegmented[names(bafsegmented)]

  SNPposhet = SNPpos[names(bafsegmented),]
  autoprobes = !(SNPposhet[,1]%in%sexchromosomes)

  b2 = b[autoprobes]
  r2 = r[autoprobes]

  s = make_segments(r2,b2)
  d = create_distance_matrix(s, gamma)

  TheoretMaxdist = sum(rep(0.25,dim(s)[1]) * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1),na.rm=T)

  # flag the sample as non-aberrant if necessary
  nonaberrant = F
  MINABB = 0.03
  MINABBREGION = 0.005

  percentAbb = sum(ifelse(s[,"b"]==0.5,0,1)*s[,"length"])/sum(s[,"length"])
  maxsegAbb = max(ifelse(s[,"b"]==0.5,0,s[,"length"]))/sum(s[,"length"])
  if(percentAbb <= MINABB & maxsegAbb <= MINABBREGION) {
    nonaberrant = T
  }


  MINPLOIDY = 1.5
  MAXPLOIDY = 5.5
  MINRHO = 0.2
  MINGOODNESSOFFIT = 80
  MINPERCZERO = 0.02
  MINPERCZEROABB = 0.1
  MINPERCODDEVEN = 0.05
  MINPLOIDYSTRICT = 1.7
  MAXPLOIDYSTRICT = 2.3

  nropt = 0
  localmin = NULL
  optima = list()

  if(!failedqualitycheck && is.na(rho_manual)) {

    # first, try with all filters
    for (i in 4:(dim(d)[1]-3)) {
      for (j in 4:(dim(d)[2]-3)) {
        m = d[i,j]
        seld = d[(i-3):(i+3),(j-3):(j+3)]
        seld[4,4] = max(seld)
        if(min(seld) > m) {
          psi = as.numeric(rownames(d)[i])
          rho = as.numeric(colnames(d)[j])
          nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
          nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho

          # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
          ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);

          percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])

          goodnessOfFit = (1-m/TheoretMaxdist) * 100

          if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & percentzero > MINPERCZERO) {
            nropt = nropt + 1
            optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
            localmin[nropt] = m
          }
        }
      }
    }

    # if no solution, drop the percentzero > MINPERCZERO filter (allow non-aberrant solutions - but limit the ploidy options)
    if (nropt == 0) {
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i,j]
          seld = d[(i-3):(i+3),(j-3):(j+3)]
          seld[4,4] = max(seld)
          if(min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);

            perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }

            goodnessOfFit = (1-m/TheoretMaxdist) * 100

            if (ploidy > MINPLOIDYSTRICT & ploidy < MAXPLOIDYSTRICT & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT & perczeroAbb > MINPERCZEROABB) {
              nropt = nropt + 1
              optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }

    # if still no solution, allow solutions with 100% aberrant cells (include the borders with rho = 1), but in first instance, keep the percentzero > 0.01 filter
    if (nropt == 0) {
      #first, include borders
      cold = which(as.numeric(colnames(d))>1)
      d[,cold]=1E20
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i,j]
          seld = d[(i-3):(i+3),(j-3):(j+3)]
          seld[4,4] = max(seld)
          if(min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);

            percentzero = (sum((round(nA)==0)*s[,"length"])+sum((round(nB)==0)*s[,"length"]))/sum(s[,"length"])
            percOddEven = sum((round(nA)%%2==0&round(nB)%%2==1|round(nA)%%2==1&round(nB)%%2==0)*s[,"length"])/sum(s[,"length"])
            perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }

            goodnessOfFit = (1-m/TheoretMaxdist) * 100

            if (!nonaberrant & ploidy > MINPLOIDY & ploidy < MAXPLOIDY & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT &
                (perczeroAbb > MINPERCZEROABB | percentzero > MINPERCZERO | percOddEven > MINPERCODDEVEN)) {
              nropt = nropt + 1
              optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }

    # if still no solution, drop the percentzero > MINPERCENTZERO filter, but strict ploidy borders
    if (nropt == 0) {
      for (i in 4:(dim(d)[1]-3)) {
        for (j in 4:(dim(d)[2]-3)) {
          m = d[i,j]
          seld = d[(i-3):(i+3),(j-3):(j+3)]
          seld[4,4] = max(seld)
          if(min(seld) > m) {
            psi = as.numeric(rownames(d)[i])
            rho = as.numeric(colnames(d)[j])
            nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
            nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho

            # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
            ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);

            perczeroAbb = (sum((round(nA)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1))+sum((round(nB)==0)*s[,"length"]*ifelse(s[,"b"]==0.5,0,1)))/sum(s[,"length"]*ifelse(s[,"b"]==0.5,0,1))
            # the next can happen if BAF is a flat line at 0.5
            if (is.na(perczeroAbb)) {
              perczeroAbb = 0
            }

            goodnessOfFit = (1-m/TheoretMaxdist) * 100

            if (ploidy > MINPLOIDYSTRICT & ploidy < MAXPLOIDYSTRICT & rho >= MINRHO & goodnessOfFit > MINGOODNESSOFFIT) {
              nropt = nropt + 1
              optima[[nropt]] = c(m,i,j,ploidy,goodnessOfFit)
              localmin[nropt] = m
            }
          }
        }
      }
    }
  }

  if (!is.na(rho_manual)) {

    rho = rho_manual
    psi = psi_manual

    nA = (rho-1-(s[,"b"]-1)*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho
    nB = (rho-1+s[,"b"]*2^(s[,"r"]/gamma)*((1-rho)*2+rho*psi))/rho

    # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
    ploidy = sum((nA+nB) * s[,"length"]) / sum(s[,"length"]);

    nMinor = NULL
    if (sum(nA,na.rm=T) < sum(nB,na.rm=T)) {
      nMinor = nA
    }
    else {
      nMinor = nB
    }
    m = sum(abs(nMinor - pmax(round(nMinor),0))^2 * s[,"length"] * ifelse(s[,"b"]==0.5,0.05,1), na.rm=T)

    goodnessOfFit = (1-m/TheoretMaxdist) * 100

    nropt = 1
    optima[[1]] = c(m,rho,psi,ploidy,goodnessOfFit)
    localmin[1] = m

  }


  if (nropt>0) {
    if (is.na(rho_manual)) {
      optlim = sort(localmin)[1]
      for (i in 1:length(optima)) {
        if(optima[[i]][1] == optlim) {
          psi_opt1 = as.numeric(rownames(d)[optima[[i]][2]])
          rho_opt1 = as.numeric(colnames(d)[optima[[i]][3]])
          if(rho_opt1 > 1) {
            rho_opt1 = 1
          }
          ploidy_opt1 = optima[[i]][4]
          goodnessOfFit_opt1 = optima[[i]][5]
        }
      }
    } else {
      rho_opt1 = optima[[1]][2]
      psi_opt1 = optima[[1]][3]
      ploidy_opt1 = optima[[1]][4]
      goodnessOfFit_opt1 = optima[[1]][5]
    }
  }

  if(nropt>0) {
    #plot Sunrise
    if (!is.na(distancepng)) {
      png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
    }
    ascat.plotSunrise(d,psi_opt1,rho_opt1)
    if (!is.na(distancepng)) {
      dev.off()
    }

    rho = rho_opt1
    psi = psi_opt1
    SNPposhet = SNPpos[names(bafsegmented),]
    haploidchrs = unique(c(substring(gender,1,1),substring(gender,2,2)))
    if(substring(gender,1,1)==substring(gender,2,2)) {
      haploidchrs = setdiff(haploidchrs,substring(gender,1,1))
    }
    diploidprobes = !(SNPposhet[,1]%in%haploidchrs)
    nullchrs = setdiff(sexchromosomes,unique(c(substring(gender,1,1),substring(gender,2,2))))
    nullprobes = SNPposhet[,1]%in%nullchrs

    nAfull = ifelse(diploidprobes,
                    (rho-1-(b-1)*2^(r/gamma)*((1-rho)*2+rho*psi))/rho,
                    ifelse(nullprobes,0,
                           ifelse(b<0.5,(rho-1+((1-rho)*2+rho*psi)*2^(r/gamma))/rho,0)))
    nBfull = ifelse(diploidprobes,
                    (rho-1+b*2^(r/gamma)*((1-rho)*2+rho*psi))/rho,
                    ifelse(nullprobes,0,
                           ifelse(b<0.5,0,(rho-1+((1-rho)*2+rho*psi)*2^(r/gamma))/rho)))
    nA = pmax(round(nAfull),0)
    nB = pmax(round(nBfull),0)

    if(!is.na(circos)){
      frame<-cbind(SNPposhet,nAfull,nBfull)
      chrSegmA<-rle(frame$nAfull)
      chrSegmB<-rle(frame$nBfull)
      if(all(chrSegmA$lengths==chrSegmB$lengths)){
        start=1
        for(i in 1:length(chrSegmA$values)){
          valA<-chrSegmA$values[i]
          valB<-chrSegmB$values[i]
          size<-chrSegmA$lengths[i]
          write(c(paste("hs",frame[start,]$Chr,sep=""),frame[start,]$Position,frame[(start+size-1),]$Position,valA), file = paste(circos,"_major",sep=""), ncolumns = 4, append = TRUE, sep = "\t")
          write(c(paste("hs",frame[start,]$Chr,sep=""),frame[start,]$Position,frame[(start+size-1),]$Position,valB), file = paste(circos,"_minor",sep=""), ncolumns = 4, append = TRUE, sep = "\t")
          start=start+size
        }
      }
      else{
        print("Major and minor allele copy numbers are segmented differently.")
      }
    }

    if (is.na(nonroundedprofilepng)) {
      dev.new(10,5)
    }
    else {
      if(pdfPlot){
        pdf(file = nonroundedprofilepng, width = 20, height = y_limit, pointsize=20)
      }
      else{
        png(filename = nonroundedprofilepng, width = 2000, height = (y_limit*100), res = 200)
      }
    }

    ascat.plotNonRounded(ploidy_opt1, rho_opt1, goodnessOfFit_opt1, nonaberrant, nAfull, nBfull, y_limit, bafsegmented, ch,lrr, chrnames)

    if (!is.na(nonroundedprofilepng)) {
      dev.off()
    }

    rho = rho_opt1
    psi = psi_opt1

    diploidprobes = !(SNPpos[,1]%in%haploidchrs)
    nullprobes = SNPpos[,1]%in%nullchrs

    #this replaces an occurrence of unique that caused problems
    #introduces segment spanning over chr ends, when two consecutive probes from diff chr have same logR!
    # build helping vector
    chrhelp = vector(length=length(lrrsegmented))
    for (chrnr in 1:length(ch)) {
      chrke = ch[[chrnr]]
      chrhelp[chrke] = chrnr
    }

    tlr2 = rle(lrrsegmented)
    tlr.chr= rle(chrhelp)

    tlrstart = c(1,cumsum(tlr2$lengths)+1)
    tlrstart = tlrstart[1:(length(tlrstart)-1)]
    tlrend = cumsum(tlr2$lengths)

    tlrstart.chr= c(1,cumsum(tlr.chr$lengths)+1)
    tlrstart.chr = tlrstart.chr[1:(length(tlrstart.chr)-1)]
    tlrend.chr = cumsum(tlr.chr$lengths)

    tlrend<-sort(union(tlrend, tlrend.chr))
    tlrstart<-sort(union(tlrstart, tlrstart.chr))

    tlr=NULL
    for(ind in tlrstart){
      val<-lrrsegmented[ind]
      tlr<-c(tlr, val)
    }

    # For each LRR probe, find the matching BAF probe
    # and its position in bafsegmented
    probeLookup = data.frame(
      lrrprobe = names(lrrsegmented),
      bafpos = match(names(lrrsegmented), names(bafsegmented)),
      stringsAsFactors=F
    )

    seg = NULL
    for (i in 1:length(tlr)) {
      logR = tlr[i]
      #pr = which(lrrsegmented==logR) # this was a problem
      pr = tlrstart[i]:tlrend[i]
      start = min(pr)
      end = max(pr)

      bafpos = probeLookup$bafpos[pr]
      bafpos = bafpos[!is.na(bafpos)]
      bafke  = bafsegmented[bafpos][1]

      #if bafke is NA, this means that we are dealing with a germline homozygous stretch with a copy number change within it.
      #in this case, nA and nB are irrelevant, just their sum matters
      if(is.na(bafke)) {
        bafke=0
      }

      nAraw = ifelse(diploidprobes[start],
                     (rho-1-(bafke-1)*2^(logR/gamma)*((1-rho)*2+rho*psi))/rho,
                     ifelse(nullprobes[start],0,
                            (rho-1+((1-rho)*2+rho*psi)*2^(logR/gamma))/rho))
      nBraw = ifelse(diploidprobes[start],(rho-1+bafke*2^(logR/gamma)*((1-rho)*2+rho*psi))/rho,0)
      # correct for negative values:
      if (nAraw+nBraw<0) {
        nAraw = 0
        nBraw = 0
      }
      else if (nAraw<0) {
        nBraw = nAraw+nBraw
        nAraw = 0
      }
      else if (nBraw<0) {
        nAraw = nAraw+nBraw
        nBraw = 0
      }
      # when evidence for odd copy number in segments of BAF = 0.5, assume a deviation..
      limitround = 0.5
      nA = ifelse(bafke==0.5,
                  ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
                         round(nAraw)+1,
                         ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround,
                                round(nAraw),
                                round(nAraw))),
                  round(nAraw))
      nB = ifelse(bafke==0.5,
                  ifelse(nAraw+nBraw>round(nAraw)+round(nBraw)+limitround,
                         round(nBraw),
                         ifelse(nAraw+nBraw<round(nAraw)+round(nBraw)-limitround,
                                round(nBraw)-1,
                                round(nBraw))),
                  round(nBraw))
      if (is.null(seg)) {
        seg = t(as.matrix(c(start,end,nA,nB)))
        seg_raw = t(as.matrix(c(start,end,nA,nB,nAraw,nBraw)))
      }
      else {
        seg = rbind(seg,c(start,end,nA,nB))
        seg_raw = rbind(seg_raw,c(start,end,nA,nB,nAraw,nBraw))
      }
    }
    colnames(seg)=c("start","end","nA","nB")
    colnames(seg_raw)=c("start","end","nA","nB","nAraw","nBraw")

    # every repeat joins 2 ends. 20 repeats will join about 1 million ends..
    for (rep in 1:20) {
      seg2=seg
      seg = NULL
      skipnext = F
      for(i in 1:dim(seg2)[1]) {
        if(!skipnext) {
          if(i != dim(seg2)[1] && seg2[i,"nA"]==seg2[i+1,"nA"] && seg2[i,"nB"]==seg2[i+1,"nB"] &&
             chrhelp[seg2[i,"end"]]==chrhelp[seg2[i+1,"start"]]) {
            segline = c(seg2[i,"start"],seg2[i+1,"end"],seg2[i,3:4])
            skipnext = T
          }
          else {
            segline = seg2[i,]
          }

          if (is.null(seg)) {
            seg = t(as.matrix(segline))
          }
          else {
            seg = rbind(seg,segline)
          }
        }
        else {
          skipnext = F
        }
      }
      colnames(seg)=colnames(seg2)
    }
    rownames(seg)=NULL

    nMajor = vector(length = length(lrrsegmented))
    names(nMajor) = names(lrrsegmented)
    nMinor = vector(length = length(lrrsegmented))
    names(nMinor) = names(lrrsegmented)

    for (i in 1:dim(seg)[1]) {
      nMajor[seg[i,"start"]:seg[i,"end"]] = seg[i,"nA"]
      nMinor[seg[i,"start"]:seg[i,"end"]] = seg[i,"nB"]
    }

    n1all = vector(length = length(lrrsegmented))
    names(n1all) = names(lrrsegmented)
    n2all = vector(length = length(lrrsegmented))
    names(n2all) = names(lrrsegmented)

    # note: any of these can have length 0
    NAprobes = which(is.na(lrr))
    heteroprobes = setdiff(which(names(lrrsegmented)%in%names(bafsegmented)),NAprobes)
    homoprobes = setdiff(setdiff(which(!is.na(baf)),heteroprobes),NAprobes)
    CNprobes = setdiff(which(is.na(baf)),NAprobes)

    n1all[NAprobes] = NA
    n2all[NAprobes] = NA
    n1all[CNprobes] = nMajor[CNprobes]+nMinor[CNprobes]
    n2all[CNprobes] = 0
    heteroprobes2 = names(lrrsegmented)[heteroprobes]
    n1all[heteroprobes] = ifelse(baf[heteroprobes2]<=0.5,nMajor[heteroprobes], nMinor[heteroprobes])
    n2all[heteroprobes] = ifelse(baf[heteroprobes2]>0.5,nMajor[heteroprobes], nMinor[heteroprobes])
    n1all[homoprobes] = ifelse(baf[homoprobes]<=0.5,nMajor[homoprobes]+nMinor[homoprobes],0)
    n2all[homoprobes] = ifelse(baf[homoprobes]>0.5,nMajor[homoprobes]+nMinor[homoprobes],0)


    # plot ASCAT profile
    if (is.na(copynumberprofilespng)) {
      dev.new(10,2.5)
    }
    else {
      if(pdfPlot){
        pdf(file = copynumberprofilespng, width = 20, height = y_limit, pointsize=20)
      }
      else{
        png(filename = copynumberprofilespng, width = 2000, height = (y_limit*100), res = 200)
      }
    }
    #plot ascat profile
    ascat.plotAscatProfile(n1all, n2all, heteroprobes, ploidy_opt1, rho_opt1, goodnessOfFit_opt1, nonaberrant,y_limit, ch, lrr, bafsegmented, chrnames)


    if (!is.na(copynumberprofilespng)) {
      dev.off()
    }


    if (!is.na(aberrationreliabilitypng)) {
      png(filename = aberrationreliabilitypng, width = 2000, height = 500, res = 200)
      par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)

      diploidprobes = !(SNPposhet[,1]%in%haploidchrs)
      nullprobes = SNPposhet[,1]%in%nullchrs

      rBacktransform = ifelse(diploidprobes,
                              gamma*log((rho*(nA+nB)+(1-rho)*2)/((1-rho)*2+rho*psi),2),
                              # the value for nullprobes is arbitrary (but doesn't matter, as these are not plotted anyway because BAF=0.5)
                              ifelse(nullprobes,-10,gamma*log((rho*(nA+nB)+(1-rho))/((1-rho)*2+rho*psi),2)))

      bBacktransform = ifelse(diploidprobes,
                              (1-rho+rho*nB)/(2-2*rho+rho*(nA+nB)),
                              ifelse(nullprobes,0.5,0))

      rConf = ifelse(abs(rBacktransform)>0.15,pmin(100,pmax(0,100*(1-abs(rBacktransform-r)/abs(r)))),NA)
      bConf = ifelse(diploidprobes & bBacktransform!=0.5, pmin(100,pmax(0,ifelse(b==0.5,100,100*(1-abs(bBacktransform-b)/abs(b-0.5))))), NA)
      confidence = ifelse(is.na(rConf),bConf,ifelse(is.na(bConf),rConf,(rConf+bConf)/2))
      maintitle = paste("Aberration reliability score (%), average: ", sprintf("%2.0f",mean(confidence,na.rm=T)),"%",sep="")
      plot(c(1,length(nAfull)), c(0,100), type = "n", xaxt = "n", main = maintitle, xlab = "", ylab = "")
      points(confidence,col="blue",pch = "|")
      abline(v=0,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (i in 1:length(ch)) {
        chrk = ch[[i]];
        chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk_hetero)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,5,chrs[i], pos = 1, cex = 2)
        abline(v=vpos,lty=1,col="lightgrey")
      }
      dev.off()
    }

    return(list(rho = rho_opt1, psi = psi_opt1, goodnessOfFit = goodnessOfFit_opt1, nonaberrant = nonaberrant,
                nA = n1all, nB = n2all, seg = seg, seg_raw = seg_raw, distance_matrix = d))

  }

  else {

    name=gsub(".sunrise.png","",basename(distancepng))

    png(filename = distancepng, width = 1000, height = 1000, res = 1000/7)
    ascat.plotSunrise(d,0,0)
    dev.off()

    warning(paste("ASCAT could not find an optimal ploidy and cellularity value for sample ", name, ".\n", sep=""))
    return(list(rho = NA, psi = NA, goodnessOfFit = NA, nonaberrant = F, nA = NA, nB = NA, seg = NA, seg_raw = NA, distance_matrix = NA))
  }

}

#' ascat.plotSunrise
#'
#' @param d distance matrix for a range of ploidy and tumour percentage values
#' @param psi_opt1 optimal ploidy
#' @param rho_opt1 optimal aberrant cell fraction
#' @param minim when set to true, optimal regions in the sunrise plot are depicted in blue; if set to false, colours are inverted and red corresponds to optimal values (default: TRUE)
#'
#' @return plot visualising range of ploidy and tumour percentage values
#' @export
ascat.plotSunrise<-function(d, psi_opt1, rho_opt1, minim=T){

  par(mar = c(5,5,0.5,0.5), cex=0.75, cex.lab=2, cex.axis=2)

  if(minim){
    hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256))
  } else {
    hmcol = colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256)
  }
  image(log(d), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")

  ploidy_min<-as.numeric(rownames(d)[1])
  ploidy_max<-as.numeric(rownames(d)[nrow(d)])
  purity_min<-as.numeric(colnames(d)[1])
  purity_max<-as.numeric(colnames(d)[ncol(d)])

  axis(1, at = seq(0, 1, by = 1/(ploidy_max-1)), labels = seq(ploidy_min, ploidy_max, by = 1))
  axis(2, at = seq(0, 1/purity_max, by = 1/3/purity_max), labels = seq(purity_min, purity_max, by = 3/10))

  if(psi_opt1>0 && rho_opt1>0){
    points((psi_opt1-ploidy_min)/(ploidy_max-1),(rho_opt1-purity_min)/(1/purity_max),col="green",pch="X", cex = 2)
  }
}


#' @title ascat.plotNonRounded
#'
#' @description Function plotting the unrounded ASCAT copy number over all chromosomes
#' @param ploidy ploidy of the sample
#' @param rho purity of the sample
#' @param goodnessOfFit estimated goodness of fit
#' @param nonaberrant boolean flag denoting non-aberrated samples
#' @param nAfull copy number major allele
#' @param nBfull copy number minor allele
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#  @param textFlag Optional flag to add the positions of fragments located outside of the plotting area to the plots. Default=F
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param chrs a vector containing the names for the chromosomes (e.g. c(1:22,"X"))
#' @param ch a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#'
#' @return plot showing the nonrounded copy number profile, using base plotting function
#'
#' @export
ascat.plotNonRounded <- function(ploidy, rho, goodnessOfFit, nonaberrant, nAfull, nBfull,y_limit=5,bafsegmented,ch,lrr, chrs){
  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit),"%", ifelse(nonaberrant,", non-aberrant",""),sep="")
  nBfullPlot<-ifelse(nBfull<y_limit, nBfull, y_limit+0.1)
  nAfullPlot<-ifelse((nAfull+nBfull)<y_limit, nAfull+nBfull, y_limit+0.1)
  colourTotal = "purple"
  colourMinor = "blue"
  base.gw.plot(bafsegmented,nAfullPlot,nBfullPlot,colourTotal,colourMinor,maintitle,ch,lrr,chrs,y_limit,twoColours=TRUE)
}


#  @title create.bb.plot.average
#
#  @param BAFvals B Allele Frequency for every probe with position
#  @param subclones Segments with chromosomal locations
#  @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#  @param ploidy ploidy of the sample
#  @param rho purity of the sample
#  @param goodnessOfFit estimated goodness of fit
#  @param pos_min
#  @param pos_max
#  @param segment_states_min Vector containing copy number per segment minor allele
#  @param segment_states_tot Vector containing copy number per segment total copy number
#  @param chr.segs Vector containing chromosome segments
#
#  @return plot showing Battenberg average copy number profile using base plotting function
#
# create.bb.plot.average = function(BAFvals, subclones, bafsegmented, ploidy, rho, goodnessOfFit, pos_min, pos_max, segment_states_min, segment_states_tot, chr.segs){
#   maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit*100),"%",sep="")
#   nTotal = array(NA, nrow(BAFvals))
#   nMinor = array(NA, nrow(BAFvals))
#   for (i in 1:nrow(subclones)) {
#     segm_chr = subclones$chr[i] == BAFvals$Chromosome & subclones$startpos[i] < BAFvals$Position & subclones$endpos[i] >= BAFvals$Position
#     pos_min = min(which(segm_chr))
#     pos_max = max(which(segm_chr))
#     nTotal = c(nTotal, segment_states_tot[pos_min[i]:pos_max[i]])
#     nMinor = c(nMinor, segment_states_min[pos_min[i]:pos_max[i]])
#   }
#   colourTotal = "#E69F00"
#   colourMinor = "#2f4f4f"
#   base.gw.plot(bafsegmented,nTotal,nMinor,colourTotal,colourMinor,maintitle,chr.segs,y_limit,textFlag)
# }

#' @title base.gw.plot
#' @description Basis for the genome-wide plots
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param nAfullPlot Total segment copy number
#' @param nBfullPlot Segment copy number minor allele
#' @param colourTotal Colour to plot total copy number
#' @param colourMinor Colour to plot minor allele
#' @param maintitle Title comprising ploidy, rho, goodness of fit
#' @param chr.segs Vector comprising chromosome segments
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param chr.names Vector giving the names of the chromosomes as displayed on the figure
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#' @param twoColours Optional flag to specify colours, if TRUE colour is paler for CN values > y_limit
#' @keywords internal
#' @return basic plot containing chromosome positions and names, plots copy number for either ASCAT non rounded or BB average
#' @export

base.gw.plot = function(bafsegmented,nAfullPlot,nBfullPlot,colourTotal,colourMinor,maintitle,chr.segs,lrr,chr.names,y_limit,twoColours=FALSE){
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2.5)
  ticks=seq(0, y_limit, 1)
  plot(c(1,length(nAfullPlot)), c(0,y_limit), type = "n", xaxt = "n", yaxt="n", main = maintitle, xlab = "", ylab = "")
  axis(side = 2, at = ticks)
  abline(h=ticks, col="lightgrey", lty=1)

  A_rle<-rle(nAfullPlot)
  start=0
  #plot total copy number
  for(i in 1:length(A_rle$values)){
    val<-A_rle$values[i]
    size<-A_rle$lengths[i]
    rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal,red.f=0.75,green.f=0.75,blue.f=0.75), colourTotal), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourTotal,red.f=0.75,green.f=0.75,blue.f=0.75), colourTotal))
    start=start+size
  }

  B_rle<-rle(nBfullPlot)
  start=0
  #plot minor allele copy number
  for(i in 1:length(B_rle$values)){
    val<-B_rle$values[i]
    size<-B_rle$lengths[i]
    rect(start, (val-0.07), (start+size-1), (val+0.07), col=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75,green.f=0.75,blue.f=0.75), colourMinor), border=ifelse((twoColours & val>=y_limit), adjustcolor(colourMinor, red.f=0.75,green.f=0.75,blue.f=0.75), colourMinor))
    start=start+size
  }

  chrk_tot_len = 0
  abline(v=0,lty=1,col="lightgrey")
  for (i in 1:length(chr.segs)) {
    chrk = chr.segs[[i]];
    chrk_hetero = intersect(names(lrr)[chrk],names(bafsegmented))
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk_hetero)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,y_limit,chr.names[i], pos = 1, cex = 2)
    abline(v=vpos,lty=1,col="lightgrey")
  }

  #   #add text to too high fragments
  #   if(textFlag){
  #     #      rleB<-rle(nBfullPlot>y_limit)
  #     #      pos<-0
  #     #      for(i in 1:length(rleB$values)){
  #     #        if(rleB$values[i]){
  #     #          xpos=pos+(rleB$lengths[i]/2)
  #     #          text(xpos,y_limit+0.1,sprintf("%1.2f",nBfull[pos+1]), pos = 1, cex = 0.7)
  #     #        }
  #     #        pos=pos+rleB$lengths[i]
  #     #      }
  #
  #     rleA<-rle(nAfullPlot>y_limit)
  #     pos<-0
  #     for(i in 1:length(rleA$values)){
  #       if(rleA$values[i]){
  #         xpos=pos+(rleA$lengths[i]/2)
  #         text(xpos,y_limit+0.1,sprintf("%1.2f",(nAfull[pos+1]+nBfull[pos+1])), pos = 1, cex = 0.7)
  #       }
  #       pos=pos+rleA$lengths[i]
  #     }
  #   }
}

#' @title ascat.plotAscatProfile
#'
#' @description Function plotting the rounded ASCAT profiles over all chromosomes
#' @param n1all copy number major allele
#' @param n2all copy number minor allele
#' @param heteroprobes probes with heterozygous germline
#' @param ploidy ploidy of the sample
#' @param rho purity of the sample
#' @param goodnessOfFit estimated goodness of fit
#' @param nonaberrant boolean flag denoting non-aberrated samples
#' @param y_limit Optional parameter determining the size of the y axis in the nonrounded plot and ASCAT profile. Default=5
#' @param ch a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chrs a vector containing the names for the chromosomes (e.g. c(1:22,"X"))
#'
#' @return plot showing the ASCAT profile of the sample
#'
#' @export
#'
ascat.plotAscatProfile<-function(n1all, n2all, heteroprobes, ploidy, rho, goodnessOfFit, nonaberrant, y_limit=5, ch, lrr, bafsegmented, chrs){
  nA2 = n1all[heteroprobes]
  nB2 = n2all[heteroprobes]
  nA = ifelse(nA2>nB2,nA2,nB2)
  nB = ifelse(nA2>nB2,nB2,nA2)

  nBPlot<-ifelse(nB<=y_limit, nB+0.1, y_limit+0.1)
  nAPlot<-ifelse(nA<=y_limit, nA-0.1, y_limit+0.1)

  colourTotal="red"
  colourMinor="green"

  maintitle = paste("Ploidy: ",sprintf("%1.2f",ploidy),", aberrant cell fraction: ",sprintf("%2.0f",rho*100),"%, goodness of fit: ",sprintf("%2.1f",goodnessOfFit),"%", ifelse(nonaberrant,", non-aberrant",""),sep="")
  base.gw.plot(bafsegmented,nAPlot,nBPlot,colourTotal,colourMinor,maintitle,ch,lrr,chrs,y_limit,twoColours=TRUE)
}



#
# Enhanced bivariate PCF filter for aCGH data (v. 08.02.2010)
# Whole chromosomes/chromosome arms wrapper function
#

fastAspcf <- function(logR, allB, kmin, gamma){

  N <- length(logR)
  w <- 1000 #w: windowsize
  d <- 100

  startw = -d
  stopw = w-d

  nseg = 0
  var2 = 0
  breakpts = 0
  larger = TRUE
  repeat{
    from <- max(c(1,startw))
    to  <-  min(c(stopw,N))
    logRpart <- logR[from:to]
    allBpart <- allB[from:to]
    allBflip <- allBpart
    allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]

    sd1 <- getMad(logRpart)
    sd2 <- getMad(allBflip)

    #Must check that sd1 and sd2 are defined and != 0:
    sd.valid <- c(!is.na(sd1),!is.na(sd2),sd1!=0,sd2!=0)
    if(all(sd.valid)){
      #run aspcfpart:
      #part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      part.res <- aspcfpart(logRpart=logRpart, allBflip=allBflip, a=startw, b=stopw, d=d, sd1=sd1, sd2=sd2, N=N, kmin=kmin, gamma=gamma)
      breakptspart <- part.res$breakpts
      # the 'larger' is (occasionally) necessary in the last window of the segmentation!
      larger = breakptspart>breakpts[length(breakpts)]
      breakpts <- c(breakpts, breakptspart[larger])
      var2 <- var2 + sd2^2
      nseg = nseg+1
    }

    if(stopw < N+d){
      startw <- min(stopw-2*d + 1,N-2*d)
      stopw <- startw + w
    }else{
      break
    }

  }#end repeat
  breakpts <- unique(c(breakpts, N))
  if(nseg==0){nseg=1}  #just in case the sd-test never passes.
  sd2 <- sqrt(var2/nseg)

  # On each segment calculate mean of unflipped B allele data
  frst <- breakpts[1:length(breakpts)-1] + 1
  last <- breakpts[2:length(breakpts)]
  nseg <- length(frst)

  yhat1 <- rep(NA,N)
  yhat2 <- rep(NA,N)

  for(i in 1:nseg){
    yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), last[i]-frst[i]+1)
    yi2 <- allB[frst[i]:last[i]]
    # Center data around zero (by subtracting 0.5) and estimate mean
    if(length(yi2)== 0){
      mu <- 0
    }else{
      mu <- mean(abs(yi2-0.5))
    }

    # Make a (slightly arbitrary) decision concerning branches
    # This may be improved by a test of equal variances
    if(sqrt(sd2^2+mu^2) < 2*sd2){
      mu <- 0
    }
    yhat2[frst[i]:last[i]] <- rep(mu+0.5,last[i]-frst[i]+1)
  }

  return(list(yhat1=yhat1,yhat2=yhat2))

}#end fastAspcf



aspcfpart <- function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma){

  from <- max(c(1,a))
  usefrom <- max(c(1,a+d))
  useto <- min(c(N,b-d))

  N <- length(logRpart)
  y1 <- logRpart
  y2 <- allBflip

  #Check that vectors are long enough to run algorithm:
  if(N < 2*kmin){
    breakpts <- 0
    return(list(breakpts=breakpts))
  }

  # Find initSum, initKvad, initAve for segment y[1..kmin]
  initSum1 <- sum(y1[1:kmin])
  initKvad1 <- sum(y1[1:kmin]^2)
  initAve1 <- initSum1/kmin
  initSum2 <- sum(y2[1:kmin])
  initKvad2 <- sum(y2[1:kmin]^2)
  initAve2 <- initSum2/kmin

  # Define vector of best costs
  bestCost <- rep(0,N)
  cost1 <- (initKvad1 - initSum1*initAve1)/sd1^2
  cost2 <- (initKvad2 - initSum2*initAve2)/sd2^2
  bestCost[kmin] <- cost1 + cost2

  # Define vector of best splits
  bestSplit <- rep(0,N)

  # Define vector of best averages
  bestAver1 <- rep(0,N)
  bestAver2 <- rep(0,N)
  bestAver1[kmin] <- initAve1
  bestAver2[kmin] <- initAve2


  #Initialize
  Sum1 <- rep(0,N)
  Sum2 <- rep(0,N)
  Kvad1 <- rep(0,N)
  Kvad2 <- rep(0,N)
  Aver1 <- rep(0,N)
  Aver2 <- rep(0,N)
  Cost <- rep(0,N)

  # We have to treat the region y(1..2*kmin-1) separately, as it
  # cannot be split into two full segments
  kminP1 <- kmin+1
  for (k in (kminP1):(2*kmin-1)) {
    Sum1[kminP1:k] <- Sum1[kminP1:k]+y1[k]
    Aver1[kminP1:k] <- Sum1[kminP1:k]/((k-kmin):1)
    Kvad1[kminP1:k] <- Kvad1[kminP1:k]+y1[k]^2
    Sum2[kminP1:k] <- Sum2[kminP1:k]+y2[k]
    Aver2[kminP1:k] <- Sum2[kminP1:k]/((k-kmin):1)
    Kvad2[kminP1:k] <- Kvad2[kminP1:k]+y2[k]^2


    bestAver1[k] <- (initSum1+Sum1[kminP1])/k
    bestAver2[k] <- (initSum2+Sum2[kminP1])/k
    cost1 <- ((initKvad1+Kvad1[kminP1])-k*bestAver1[k]^2)/sd1^2
    cost2 <- ((initKvad2+Kvad2[kminP1])-k*bestAver2[k]^2)/sd2^2

    bestCost[k] <- cost1 + cost2

  }


  for (n in (2*kmin):N) {

    nMkminP1=n-kmin+1

    Sum1[kminP1:n] <- Sum1[kminP1:n]+ y1[n]
    Aver1[kminP1:n] <- Sum1[kminP1:n]/((n-kmin):1)
    Kvad1[kminP1:n] <- Kvad1[kminP1:n]+ (y1[n])^2

    cost1 <- (Kvad1[kminP1:nMkminP1]-Sum1[kminP1:nMkminP1]*Aver1[kminP1:nMkminP1])/sd1^2

    Sum2[kminP1:n] <- Sum2[kminP1:n]+ y2[n]
    Aver2[kminP1:n] <- Sum2[kminP1:n]/((n-kmin):1)
    Kvad2[kminP1:n] <- Kvad2[kminP1:n]+ (y2[n])^2
    cost2 <- (Kvad2[kminP1:nMkminP1]-Sum2[kminP1:nMkminP1]*Aver2[kminP1:nMkminP1])/sd2^2

    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)] + cost1 + cost2

    Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
    cost <- Cost[Pos] + gamma

    aver1 <- Aver1[Pos]
    aver2 <- Aver2[Pos]
    totAver1 <- (Sum1[kminP1]+initSum1)/n
    totCost1 <- ((Kvad1[kminP1]+initKvad1) - n*totAver1*totAver1)/sd1^2
    totAver2 <- (Sum2[kminP1]+initSum2)/n
    totCost2 <- ((Kvad2[kminP1]+initKvad2) - n*totAver2*totAver2)/sd2^2
    totCost <- totCost1 + totCost2

    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver1 <- totAver1
      aver2 <- totAver2
    }
    bestCost[n] <- cost
    bestAver1[n] <- aver1
    bestAver2[n] <- aver2
    bestSplit[n] <- Pos-1


  }#endfor


  # Trace back
  n <- N
  breakpts <- n
  while(n > 0){
    breakpts <- c(bestSplit[n], breakpts)
    n <- bestSplit[n]
  }#endwhile

  breakpts <- breakpts + from -1
  breakpts <- breakpts[breakpts>=usefrom & breakpts<=useto]

  return(list(breakpts=breakpts))

}#end aspcfpart






div <- function(a, b, c){

  if(nargs() < 3){
    c <- 0
  }#endif

  if(b > 0){
    v <-  a/b
  }else{
    v <- c
  }#endif

  return(v)
}#endfunction


#function v = div(a, b, c)
# if nargin < 3
#     c = 0;
# end
# if b > 0
#     v = a/b;
# else
#     v = c;
# end
#end



#Get mad SD (based on KL code)
getMad <- function(x,k=25){

  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]

  #Calculate runMedian
  runMedian <- medianFilter(x,k)

  dif <- x-runMedian
  SD <- mad(dif)

  return(SD)
}





exactPcf <- function(y, kmin, gamma) {
  ## Implementaion of exact PCF by Potts-filtering
  ## x: input array of (log2) copy numbers
  ## kmin: Mininal length of plateaus
  ## gamma: penalty for each discontinuity
  N <- length(y)
  yhat <- rep(0,N);
  if (N < 2*kmin) {
    yhat <- rep(mean(y),N)
    return(yhat)
  }
  initSum <- sum(y[1:kmin])
  initKvad <- sum(y[1:kmin]^2)
  initAve <- initSum/kmin;
  bestCost <- rep(0,N)
  bestCost[kmin] <- initKvad - initSum*initAve
  bestSplit <- rep(0,N)
  bestAver <- rep(0,N)
  bestAver[kmin] <- initAve
  Sum <- rep(0,N)
  Kvad <- rep(0,N)
  Aver <- rep(0,N)
  Cost <- rep(0,N)
  kminP1=kmin+1
  for (k in (kminP1):(2*kmin-1)) {
    Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
    Aver[kminP1:k] <- Sum[kminP1:k]/((k-kmin):1)
    Kvad[kminP1:k] <- Kvad[kminP1:k]+y[k]^2
    bestAver[k] <- (initSum+Sum[kminP1])/k
    bestCost[k] <- (initKvad+Kvad[kminP1])-k*bestAver[k]^2
  }
  for (n in (2*kmin):N) {
    yn <- y[n]
    yn2 <- yn^2
    Sum[kminP1:n] <- Sum[kminP1:n]+yn
    Aver[kminP1:n] <- Sum[kminP1:n]/((n-kmin):1)
    Kvad[kminP1:n] <- Kvad[kminP1:n]+yn2
    nMkminP1=n-kmin+1
    Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]+Kvad[kminP1:nMkminP1]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
    Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
    cost <- Cost[Pos]
    aver <- Aver[Pos]
    totAver <- (Sum[kminP1]+initSum)/n
    totCost <- (Kvad[kminP1]+initKvad) - n*totAver*totAver
    if (totCost < cost) {
      Pos <- 1
      cost <- totCost
      aver <- totAver
    }
    bestCost[n] <- cost
    bestAver[n] <- aver
    bestSplit[n] <- Pos-1
  }
  n <- N
  while (n > 0) {
    yhat[(bestSplit[n]+1):n] <- bestAver[n]
    n <- bestSplit[n]
  }
  return(yhat)
}


#Perform MAD winsorization:
madWins <- function(x,tau,k){
  xhat <- medianFilter(x,k)
  d <- x-xhat
  SD <- mad(d)
  z <- tau*SD
  xwin <- xhat + psi(d, z)
  outliers <- rep(0, length(x))
  outliers[x > xwin] <- 1
  outliers[x < xwin] <- -1
  return(list(ywin=xwin,sdev=SD,outliers=outliers))
}


#########################################################################
# Function to calculate running median for a given a window size
#########################################################################

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1

  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }

  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)

}


fillNA = function(vec, zeroIsNA=TRUE) {
  if (zeroIsNA) {vec[vec==0] <- NA}
  nas = which(is.na(vec))

  if(length(nas) == 0) {
    return(vec)
  }

  # Find stretches of contiguous NAs
  starts = c(1, which(diff(nas)>1)+1)
  ends = c(starts[-1] - 1, length(nas))

  starts = nas[starts]
  ends = nas[ends]

  # Special-case: vec[1] is NA
  startAt = 1
  if(starts[1]==1) {
    vec[1:ends[1]] = vec[ends[1]+1]
    startAt = 2
  }

  if (startAt > length(starts)) {
    return(vec)
  }

  # Special-case: last element in vec is NA
  endAt = length(starts)
  if(is.na(vec[length(vec)])) {
    vec[starts[endAt]:ends[endAt]] = vec[starts[endAt]-1]
    endAt = endAt-1
  }

  if (endAt < startAt) {
    return(vec)
  }

  # For each stretch of NAs, set start:midpoint to the value before,
  # and midpoint+1:end to the value after.
  for(i in startAt:endAt) {
    start = starts[i]
    end = ends[i]
    N = 1 + end-start
    if (N==1) {
      vec[start] = vec[start-1]
    } else {
      midpoint = start+ceiling(N/2)
      vec[start:midpoint] = vec[start-1]
      vec[(midpoint+1):end] = vec[end+1]
    }
  }

  return(vec)
}


psi <- function(x,z){
  xwin <- x
  xwin[x < -z] <- -z
  xwin[x > z] <- z
  return(xwin)
}







#' @title ascat.predictGermlineGenotypes
#' @description predicts the germline genotypes of samples for which no matched germline sample
#' is available
#' @param ASCATobj an ASCAT object
#' @param platform used array platform
#' @details Currently possible values for platform:\cr
#' AffySNP6 (default)\cr
#' Custom10k\cr
#' Illumina109k\cr
#' IlluminaCytoSNP\cr
#' Illumina610k\cr
#' Illumina660k\cr
#' Illumina700k\cr
#' Illumina1M\cr
#' Illumina2.5M\cr
#' IlluminaOmni5\cr
#' Affy10k\cr
#' Affy100k\cr
#' Affy250k_sty\cr
#' Affy250k_nsp\cr
#' AffyOncoScan\cr
#' AffyCytoScanHD\cr
#' HumanCNV370quad\cr
#' HumanCore12\cr
#' HumanCoreExome24\cr
#' HumanOmniExpress12\cr
#' IlluminaOmniExpressExome\cr
#'
#' @return predicted germline genotypes
#'
#' @export
ascat.predictGermlineGenotypes = function(ASCATobj, platform = "AffySNP6") {
  Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
  rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)

  if (platform=="Custom10k") {
    maxHomozygous = 0.05
    proportionHetero = 0.59
    proportionHomo = 0.38
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform=="Illumina109k") {
    maxHomozygous = 0.05
    proportionHetero = 0.35
    proportionHomo = 0.60
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform=="IlluminaCytoSNP") {
    maxHomozygous = 0.05
    proportionHetero = 0.28
    proportionHomo = 0.62
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="Illumina610k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina660k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina700k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina1M") {
    maxHomozygous = 0.05
    proportionHetero = 0.22
    proportionHomo = 0.74
    proportionOpen = 0.02
    segmentLength = 100
    #previousvalues:
    #proportionHetero = 0.24
    #proportionOpen = 0.01
    #segmentLength = 60
  }
  else if (platform=="Illumina2.5M") {
    maxHomozygous = 0.05
    proportionHetero = 0.21
    proportionHomo = 0.745
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="IlluminaOmni5") {
    maxHomozygous = 0.05
    proportionHetero = 0.13
    proportionHomo = 0.855
    proportionOpen = 0.01
    segmentLength = 100
  }
  else if (platform=="Affy10k") {
    maxHomozygous = 0.04
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 20
  }
  else if (platform=="Affy100k") {
    maxHomozygous = 0.05
    proportionHetero = 0.27
    proportionHomo = 0.62
    proportionOpen = 0.09
    segmentLength = 30
  }
  else if (platform=="Affy250k_sty") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="Affy250k_nsp") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="AffySNP6") {
    maxHomozygous = 0.05
    proportionHetero = 0.25
    proportionHomo = 0.67
    proportionOpen = 0.04
    segmentLength = 100
  }
  else if (platform=="AffyOncoScan") {
    maxHomozygous = 0.04
    proportionHetero = 0.355
    proportionHomo = 0.605
    proportionOpen = 0.025
    segmentLength = 30
    # maxHomozygous = 0.05
    # proportionHetero = 0.24
    # proportionHomo = 0.69
    # proportionOpen = 0.04
    # segmentLength = 30
  }
  else if (platform=="AffyCytoScanHD") {
    # maxHomozygous = 0.05
    # proportionHetero = 0.26
    # proportionHomo = 0.69
    # proportionOpen = 0.03
    # segmentLength = 30
    maxHomozygous = 0.04
    proportionHetero = 0.32
    proportionHomo = 0.60
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="HumanCNV370quad") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 20
  }
  else if (platform=="HumanCore12") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 20
  }
  else if (platform=="HumanCoreExome24") {
    maxHomozygous = 0.05
    proportionHetero = 0.175
    proportionHomo = 0.79
    proportionOpen = 0.02
    segmentLength = 100
  }
  else if (platform=="HumanOmniExpress12") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 100
  }
  else if (platform=="IlluminaOmniExpressExome") {
    maxHomozygous = 0.05
    proportionHetero = 0.35
    proportionHomo = 0.60
    proportionOpen = 0.03
    segmentLength = 100
  }
  #   else if (platform=="OmniZhonghua8") {
  #     maxHomozygous = 0.05
  #     proportionHetero = 0.295
  #     proportionHomo = 0.67
  #     proportionOpen = 0.015
  #     segmentLength = 100
  #   }
  else {
    print("Error: platform unknown")
  }

  failedarrays = NULL

  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {

    Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[,i]),i]
    names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[,i])]
    Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA),i]
    names(Tumor_LogR_noNA) = names(Tumor_BAF_noNA)

    chr_noNA = list()
    prev = 0
    for(j in 1:length(ASCATobj$chr)) {
      chrke = ASCATobj$chr[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,i]))
      chr_noNA[[j]] = (prev+1):next2
      prev = next2
    }

    ch_noNA = list()
    prev = 0
    for(j in 1:length(ASCATobj$ch)) {
      chrke = ASCATobj$ch[[j]]
      next2 = prev + sum(!is.na(ASCATobj$Tumor_BAF[chrke,i]))
      ch_noNA[[j]] = (prev+1):next2
      prev = next2
    }

    tbsam = Tumor_BAF_noNA
    # sample, mirrored
    bsm = ifelse(tbsam<0.5, tbsam, 1-tbsam)

    homoLimit = max(sort(bsm)[round(length(bsm)*proportionHomo)],maxHomozygous)

    if(homoLimit>0.25) {
      failedarrays = c(failedarrays,ASCATobj$samples[i])
    }

    Hom = ifelse(bsm<homoLimit,T,NA)

    Homo = sum(Hom==T, na.rm=T)
    Undecided = sum(is.na(Hom))

    extraHetero = round(min(proportionHetero * length(Tumor_BAF_noNA),Undecided-proportionOpen*length(Tumor_BAF_noNA)))

	Hetero = 0

    if(extraHetero>0) {

      allProbes=1:length(Tumor_BAF_noNA)
      nonHomoProbes = allProbes[is.na(Hom)|Hom==F]

      lowestDist = NULL

      # bsm, with homozygous replaced by NA
      bsmHNA=bsm
      bsmHNA[!is.na(Hom)&Hom]=NA

      for (chrke in chr_noNA) {

        chrNonHomoProbes = intersect(nonHomoProbes,chrke)

        # there must be a minimum number of probes on the chromosome, otherwise these are called homozygous anyway
        if (length(chrNonHomoProbes)>5) {

          #make sure we're not going over any borders..
          segmentLength2 = min(length(chrNonHomoProbes)-1,segmentLength)

          chrNonHomoProbesStartWindowLeft = c(rep(NA,segmentLength2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2)])
          chrNonHomoProbesEndWindowLeft = c(NA,chrNonHomoProbes[1:(length(chrNonHomoProbes)-1)])
          chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)],NA)
          chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2))
          chrNonHomoProbesStartWindowMiddle = c(rep(NA,segmentLength2/2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2/2)])
          chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2/2))

          chrLowestDist = NULL

          for (probeNr in 1:length(chrNonHomoProbes)) {
            probe = chrNonHomoProbes[probeNr]
            if(!is.na(chrNonHomoProbesStartWindowLeft[probeNr])&!is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
              medianLeft = median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], na.rm=T)
            }
            else {
              medianLeft = NA
            }
            if(!is.na(chrNonHomoProbesStartWindowRight[probeNr])&!is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
              medianRight = median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]], na.rm=T)
            }
            else {
              medianRight = NA
            }

            if(!is.na(chrNonHomoProbesStartWindowMiddle[probeNr])&!is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
              medianMiddle = median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]],
                                      bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]), na.rm=T)
            }
            else {
              medianMiddle = NA
            }

            chrLowestDist[probeNr] = min(abs(medianLeft-bsm[probe]),abs(medianRight-bsm[probe]),abs(medianMiddle-bsm[probe]),Inf,na.rm=T)
          }
        }

        # if too few probes on the chromosome
        else {
          chrLowestDist = NULL
          if (length(chrNonHomoProbes)>0) {
            # 1 is higher than any practical distance
            chrLowestDist[1:length(chrNonHomoProbes)] = 1
          }
        }

        lowestDist = c(lowestDist,chrLowestDist)
      }

      lowestDistUndecided = lowestDist[is.na(Hom[nonHomoProbes])]
      names(lowestDistUndecided)=names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]

      sorted = sort(lowestDistUndecided)
      Hom[names(sorted[1:min(length(sorted),extraHetero)])] = F

      Hetero = sum(Hom==F, na.rm=T)
      Homo = sum(Hom==T, na.rm=T)
      Undecided = sum(is.na(Hom))

    }

    png(filename = paste("tumorSep",colnames(ASCATobj$Tumor_LogR)[i],".png",sep=""), width = 2000, height = 500, res = 200)
    title = paste(paste(colnames(ASCATobj$Tumor_BAF)[i], Hetero), Homo)
    ascat.plotGenotypes(ASCATobj,title,Tumor_BAF_noNA, Hom, ch_noNA)
    dev.off()

    # set all Undecided to homozygous
    Hom[is.na(Hom)] = T
    Homozygous[names(Hom),i] = Hom
  }

  return(list(germlinegenotypes = Homozygous, failedarrays = failedarrays))

}


#' ascat.plotGenotypes
#'
#' @param ASCATobj an ASCAT object
#' @param title main title of the plot
#' @param Tumor_BAF_noNA B-allele frequencies of the tumour sample with removed NA values
#' @param Hom Boolean vector denoting homozygous SNPs
#' @param ch_noNA vector of probes per chromosome (NA values excluded)
#'
#' @return plot showing classified BAF per sample, with unused SNPs in green, germline homozygous SNPs in blue and all others in red
#' @export
#'
ascat.plotGenotypes<-function(ASCATobj, title, Tumor_BAF_noNA, Hom, ch_noNA){
  par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
  plot(c(1,length(Tumor_BAF_noNA)), c(0,1), type = "n", xaxt = "n", main = title, xlab = "", ylab = "")
  points(Tumor_BAF_noNA,col=ifelse(is.na(Hom),"green",ifelse(Hom,"blue","red")))

  abline(v=0.5,lty=1,col="lightgrey")
  chrk_tot_len = 0
  for (j in 1:length(ch_noNA)) {
    chrk = ch_noNA[[j]];
    chrk_tot_len_prev = chrk_tot_len
    chrk_tot_len = chrk_tot_len + length(chrk)
    vpos = chrk_tot_len;
    tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
    abline(v=vpos+0.5,lty=1,col="lightgrey")
  }
}
