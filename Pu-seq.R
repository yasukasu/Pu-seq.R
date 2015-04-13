version="Pu-seq analysis v1.0"

##functions
moving.ave <- function(data, n){
  r <- c()
  for(i in 1:length(data)){ 
    d.data = data[max(1,i-n): min(length(data), i+n)]
    d.data <- d.data[!is.infinite(d.data)] # exluding "Inf"
    d.data <- d.data[!is.nan(d.data)] # exluding "NaN"
    r[i]=mean(d.data)
  }
  (r)
}

Closeori<-function(pos, bin){
  
  x<-c(1:length(pos)) 
  remove<-c()
  replace<-c()
  for(i in 1:(length(pos)-1)){
    if (abs(pos[i]-pos[i+1])==bin)  
    {remove<-c(remove,x[i])
     next
    }else if
    (abs(pos[i]-pos[i+1])==2*bin) 
    {remove<-c(remove,x[i],x[i+1])
     replace<-c(replace,(pos[x[i]]+bin))
     next
    }else if
    (abs(pos[i]-pos[i+1])==3*bin) 
    {remove<-c(remove,x[i],x[i+1])
     replace<-c(replace,(pos[x[i]]+bin))
     next
    }else if
    (abs(pos[i]-pos[i+1])==4*bin) 
    {remove<-c(remove,x[i],x[i+1])
     replace<-c(replace,(pos[x[i]]+2*bin))
    }
    
  }
  posremove<-c(pos[-remove])
  posreplace<-sort(unique(c(posremove,replace)))
  return(posreplace)
  
  
}

orieff<-function(close, ratio, pos){
  #tableclose=peak positions from diff, ratio=pol usage ratio, pos=all position)
  maxpos<-c(match(close,pos))
  
  ratiomin<-c(tail((which(diff(c(FALSE,diff(ratio[1:maxpos[1]])>0,TRUE))>0)),n=1))
  
  for (i in 1:(length(maxpos)-1)){
    
    ratiomin<-c(ratiomin,(-1+maxpos[i]+tail((which(diff(c(FALSE,diff(ratio[maxpos[i]:maxpos[i+1]])>0,TRUE))>0)),n=1)))
  }
  
  ratiomax<-c()
  
  for (i in 1:(length(maxpos)-1)){
    ratiomax<-c(ratiomax,(-1+maxpos[i]+head((which(diff(c(TRUE,diff(ratio[maxpos[i]:maxpos[i+1]])>=0,FALSE))<0)),n=1)))
  }
  ratiomax<-c(ratiomax, (-1+maxpos[length(maxpos)]+(head((which(diff(c(TRUE,diff(ratio[maxpos[i]:maxpos[i+1]])>=0,FALSE))<0)),n=1))))
  #which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
  
  orieff<-c((ratio[ratiomax]-ratio[ratiomin])*100)
  oriefftable<-cbind(close,orieff)
  
  return(oriefftable)
}

Findlocalmax<-function(diffdata,position,percentile){
  
  max<-c(which(diff(c(TRUE,diff(diffdata)>=0,FALSE))<0 & diffdata>0) )
  tableall<-cbind(position[max], diffdata[max])
  perc<-quantile(tableall[,2],percentile)
  per<-which(tableall[,2]<=perc)
  table<-tableall[-per,]
  return(table)
}

orieff_merge<-function(orieff_ef, orieff_ef_pos, orieff_dr, orieff_dr_pos,chro,bin){
  
  value<-c()
  valuepos<-c()
  drpaired<-c()
  efpaired<-c()
  for (i in 1:length(orieff_ef_pos)){
    if 
    (orieff_ef_pos[i] %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[match(orieff_ef_pos[i],orieff_dr_pos)])))
      valuepos<-c(valuepos, orieff_ef_pos[i])
      drpaired<-c(drpaired, orieff_dr_pos[match(orieff_ef_pos[i],orieff_dr_pos)] )
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
    }else if
    ((orieff_ef_pos[i]+bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, orieff_ef_pos[i])
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])   
      next
      
    }else if
    ((orieff_ef_pos[i]-bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, orieff_ef_pos[i])
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])   
      next
    }else if
    ((orieff_ef_pos[i]+2*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+2*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]+bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+2*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i]) 
      next
    }else if
    ((orieff_ef_pos[i]-2*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-2*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]-bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-2*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i]) 
      next
    }else if
    ((orieff_ef_pos[i]+3*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+3*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]+2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+3*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
      
    }else if
    ((orieff_ef_pos[i]-3*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-3*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]-2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-3*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
    }else if
    ((orieff_ef_pos[i]+4*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]+4*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]+2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]+4*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      next
      
    }else if
    ((orieff_ef_pos[i]-4*bin) %in% orieff_dr_pos){
      value<-c(value, mean(c(orieff_ef[i], orieff_dr[(match((orieff_ef_pos[i]-4*bin),orieff_dr_pos))])))
      valuepos<-c(valuepos, (orieff_ef_pos[i]-2*bin))
      drpaired<-c(drpaired, (orieff_dr_pos[match((orieff_ef_pos[i]-4*bin),orieff_dr_pos)] ))
      efpaired<-c(efpaired, orieff_ef_pos[i])
      
    }
  }
  
  efpairedno<-match(efpaired,orieff_ef_pos)
  drpairedno<-match(drpaired,orieff_dr_pos)
  efunpaired<-orieff_ef_pos[-efpairedno]
  drunpaired<-orieff_dr_pos[-drpairedno]
  
  if((length(efunpaired) != 0) & (length(drunpaired) != 0)){
    
    
    
    for (i in 1:length(drunpaired)){ 
      if 
      (drunpaired[i] %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match(drunpaired[i],orieff_ef_pos)])))
        valuepos<-c(valuepos, drunpaired[i])
        next
      }else if
      ((drunpaired[i]+bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, drunpaired[i])
        next
        
      }else if
      ((drunpaired[i]-bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, drunpaired[i])
        next
        
      }else if
      ((drunpaired[i]+2*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+2*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]+bin))
        next
        
      }else if
      ((drunpaired[i]-2*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-2*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]-bin))
        next
        
      }else if
      ((drunpaired[i]+3*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+3*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]+2*bin))
        next
        
      }else if
      ((drunpaired[i]-3*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-3*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]-2*bin))
        next
      }else if
      ((drunpaired[i]+4*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]+4*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]+2*bin))
        next
        
      }else if
      ((drunpaired[i]-4*bin) %in% efunpaired){
        value<-c(value, mean(c(orieff_dr[match(drunpaired[i],orieff_dr_pos)], orieff_ef[match((drunpaired[i]-4*bin),orieff_ef_pos)])))
        valuepos<-c(valuepos, (drunpaired[i]-2*bin))
        
      } 
      chromosome<-rep(chro,length(valuepos))
      orilist<-cbind(chromosome,valuepos,round(value,3))
      colnames(orilist)<-c("chromosome","maxpos","efficiency")
    }
  }else if(length(efunpaired)==0 | length(drunpaired)==0){
    
    chromosome<-rep(chro,length(valuepos))
    orilist<-cbind(chromosome,valuepos,round(value,3))
    colnames(orilist)<-c("chromosome","maxpos","efficiency")}
  
  return(orilist) 
}


### extacting parmeters from the command-line

args <- commandArgs(trailingOnly = T)
#cat("args:", args, "\n")

n.d.f <- grep("-df", args)
n.d.r <- grep("-dr", args)
n.e.f <- grep("-ef", args)
n.e.r <- grep("-er", args)
if(length(n.d.f)!=1 || length(n.d.r)!=1 || length(n.e.f)!=1 || length(n.e.r)!=1){stop("Error in inputting files!")}

path.d.f = args[n.d.f+1]
path.d.r = args[n.d.r+1]
path.e.f = args[n.e.f+1]
path.e.r = args[n.e.r+1]
paths <- c(path.d.f, path.d.r, path.e.f, path.e.r)
if(length(grep("\\.csv$", paths))!=4){stop("Please indiates csv files")}

# the parameter for moving ave (2N+1) to smooth genome-wide data
NN = 3 
if(length(grep("-ma1", args))){
  n.ma1 = grep("-ma1", args)
  NN = as.integer(args[n.ma1+1])
  if(NN < 0){stop("The value MA must be positive.")}
}

# the parameter for moving ave (2N+1) to smooth differencials (ori detection)
NN2 = NN 
if(length(grep("-ma2", args))){
  n.ma2 = grep("-ma2", args)
  NN2 = as.integer(args[n.ma2+1])
  if(NN2 < 0){stop("The value MA for ori detection must be positive.")}
}

# prefix will be used in output files.
if(length(grep("-prefix", args))){
  n.prefix = grep("-prefix", args)
  prefix =  args[n.prefix+1]
  if(!is.character(prefix)){stop("The prefix includes characters.")}
}

# To enable Origin detection
orimode=F
if(length(grep("-ori", args))){
  orimode=T
}

cat("Usage of Pol-delta on the watson str.  =>", path.d.f, "\n")
cat("Usage of Pol-delta on the crick str.   =>", path.d.r, "\n")
cat("Usage of Pol-epsilon on the watson str.=>", path.e.f, "\n")
cat("Usage of Pol-epsilon on the crick str. =>", path.e.r, "\n")
cat("Moving ave.(Pol-usage ratio) => ", 2*NN+1, " (2x", NN, "+1) bins\n", sep="")
if(orimode){cat("Moving ave.(diff. of Pol-usage ratio) => ", 2*NN2+1, " (2x", NN2, "+1) bins\n", sep="")}
cat("Prefix => ", prefix, "\n")
cat("\n")

### Preparing data to be used 

Pol.d.f.df <- read.csv(path.d.f)
Pol.d.r.df <- read.csv(path.d.r)
Pol.e.f.df <- read.csv(path.e.f)
Pol.e.r.df <- read.csv(path.e.r)

Pol.d.f.count.sum = sum(Pol.d.f.df[,c("count")])
Pol.d.r.count.sum = sum(Pol.d.r.df[,c("count")])
Pol.e.f.count.sum = sum(Pol.e.f.df[,c("count")])
Pol.e.r.count.sum = sum(Pol.e.r.df[,c("count")])

if(all(nrow(Pol.d.f.df) == nrow(Pol.d.r.df), 
       nrow(Pol.d.f.df) == nrow(Pol.e.f.df),
       nrow(Pol.d.f.df) == nrow(Pol.e.r.df))){cat("Lengths of 4 datasets from files are identical.\n") 
} else { stop("Lengths of 4 datasets from files are NOT identical\n", sep=" ")}

chromos1 = as.character(levels(factor(Pol.d.f.df$chro)))
chromos2 = as.character(levels(factor(Pol.d.r.df$chro)))
chromos3 = as.character(levels(factor(Pol.e.f.df$chro)))
chromos4 = as.character(levels(factor(Pol.e.r.df$chro)))

if(all(identical(chromos1,chromos2), 
       identical(chromos1,chromos3), 
       identical(chromos1,chromos4))){ cat("Detected chromosome name:", chromos1, "\n\n", sep=" ")
} else { stop("Chromosome names in datasets did not match.\n", sep=" ")}

Pol.d.f.cal.all = c()
Pol.d.r.cal.all = c()
Pol.e.f.cal.all = c()
Pol.e.r.cal.all = c()

# header of wig file

header1 = function(name, h.line, color){
  paste('track type=wiggle_0 name="', name,
    '" description="', version, ';', date(),' MA:', 2*NN+1, '(2x', NN, '+1)', 
    '" visibility=full autoScale=off color=', color,
    ' yLineOnOff=on yLineMark=', h.line, 
    ' priority=10', sep="")
}

out.path.d.f=paste(prefix, "_Pol-delta.watson.wig", sep="")
out.path.d.r=paste(prefix, "_Pol-delta.crick.wig", sep="")
out.path.e.f=paste(prefix, "_Pol-epsilon.watson.wig", sep="")
out.path.e.r=paste(prefix, "_Pol-epsilon.crick.wig", sep="")

write(header1(paste(prefix, "_Pol-delta_watson-str", sep=""), 0.5, "0,0,255"), file=out.path.d.f, append=F)
write(header1(paste(prefix, "_Pol-delta_crick-str", sep=""), 0.5, "0,0,255"), file=out.path.d.r, append=F)
write(header1(paste(prefix, "_Pol-epsilon_watson-str", sep=""), 0.5, "255,0,0"), file=out.path.e.f, append=F)
write(header1(paste(prefix, "_Pol-epsilon_crick-str", sep=""), 0.5, "255,0,0"), file=out.path.e.r, append=F)

### Calculation 
if(orimode){orilist=c()}

for (chromo in chromos1){
  cat("Processing ", chromo, "...\n" ,sep="")
  
  Pol.d.f.df.chr = Pol.d.f.df[Pol.d.f.df$chro==chromo,]
  Pol.d.r.df.chr = Pol.d.r.df[Pol.d.r.df$chro==chromo,]
  Pol.e.f.df.chr = Pol.e.f.df[Pol.e.f.df$chro==chromo,]
  Pol.e.r.df.chr = Pol.e.r.df[Pol.e.r.df$chro==chromo,]

  if(all(nrow(Pol.d.f.df.chr) == nrow(Pol.d.r.df.chr), 
         nrow(Pol.d.f.df.chr) == nrow(Pol.e.f.df.chr),
         nrow(Pol.d.f.df.chr) == nrow(Pol.e.r.df.chr))){cat("Lengths of", chromo, "in 4 datasets are identical.\n") 
  } else { stop("Lengths of", chromo, "in 4 datasets are NOT identical\n", sep=" ")}
  
  Pol.d.f.count = Pol.d.f.df.chr[,c("count")]
  Pol.d.r.count = Pol.d.r.df.chr[,c("count")]
  Pol.e.f.count = Pol.e.f.df.chr[,c("count")]
  Pol.e.r.count = Pol.e.r.df.chr[,c("count")]
  
  Pol.d.f.pos   = Pol.d.f.df.chr[,c("pos")]
  bin.size = Pol.d.f.pos[2]-Pol.d.f.pos[1]
  ### normalization uisng Pol.d and epsilon both
  
  Pol.d.f.ratio <- Pol.d.f.count/Pol.d.f.count.sum
  Pol.d.r.ratio <- Pol.d.r.count/Pol.d.r.count.sum
  Pol.e.f.ratio <- Pol.e.f.count/Pol.e.f.count.sum
  Pol.e.r.ratio <- Pol.e.r.count/Pol.e.r.count.sum
  
  Pol.d.f.cal <- Pol.d.f.ratio/(Pol.d.f.ratio+Pol.e.f.ratio)
  Pol.e.f.cal <- Pol.e.f.ratio/(Pol.d.f.ratio+Pol.e.f.ratio)
  Pol.d.r.cal <- Pol.d.r.ratio/(Pol.d.r.ratio+Pol.e.r.ratio)
  Pol.e.r.cal <- Pol.e.r.ratio/(Pol.d.r.ratio+Pol.e.r.ratio)
  
  Pol.d.f.cal.ma <- moving.ave(Pol.d.f.cal, NN)
  Pol.e.f.cal.ma <- moving.ave(Pol.e.f.cal, NN)
  Pol.d.r.cal.ma <- moving.ave(Pol.d.r.cal, NN)
  Pol.e.r.cal.ma <- moving.ave(Pol.e.r.cal, NN) 

  # outputting data of chormosome 
  header2 = paste('fixedStep chrom=', chromo,
                  ' start=0 step=', bin.size,
                  ' span=', bin.size, sep="")
  
  write(header2, file=out.path.d.f, append=T) 
  write(header2, file=out.path.d.r, append=T) 
  write(header2, file=out.path.e.f, append=T) 
  write(header2, file=out.path.e.r, append=T)

  write.table(round(Pol.d.f.cal.ma,3), file=out.path.d.f, row.names=F, col.names=F,append=T)
  write.table(round(Pol.d.r.cal.ma,3), file=out.path.d.r, row.names=F, col.names=F,append=T)
  write.table(round(Pol.e.f.cal.ma,3), file=out.path.e.f, row.names=F, col.names=F,append=T)
  write.table(round(Pol.e.r.cal.ma,3), file=out.path.e.r, row.names=F, col.names=F,append=T)

  cat(chromo, " => ", length(Pol.d.f.cal.ma), " bins (1-", bin.size*length(Pol.d.f.cal.ma)," bp)\n", sep="")

  ### origin detection
  if(!orimode){cat("\n");next}

  dat_ef.chr.diff<-moving.ave(diff(Pol.e.f.cal.ma),NN2)
  dat_dr.chr.diff<-moving.ave(diff(Pol.d.r.cal.ma),NN2)
  
  #find local maxima - p= percentile treshold
  p=0.3
  
  dat_ef.chr.table<-Findlocalmax(dat_ef.chr.diff,Pol.d.f.pos, p)
  dat_dr.chr.table<-Findlocalmax(dat_dr.chr.diff,Pol.d.f.pos, p)
  
  
  #merge peaks within 4 bins
  
  dat_ef.chr.table.merged<-Closeori(as.numeric(as.vector(dat_ef.chr.table[,1])),bin.size)
  dat_dr.chr.table.merged<-Closeori(as.numeric(as.vector(dat_dr.chr.table[,1])),bin.size)
  
  #calculate efficiency
  
  dat_ef.chr.orieff<-orieff(dat_ef.chr.table.merged, Pol.e.f.cal.ma, Pol.d.f.pos)
  dat_dr.chr.orieff<-orieff(dat_dr.chr.table.merged, Pol.d.r.cal.ma, Pol.d.f.pos)
  
  #find peaks that are in both within plusminus 4 bins
  
  orilist.chr<-orieff_merge(dat_ef.chr.orieff[,2],dat_ef.chr.orieff[,1],dat_dr.chr.orieff[,2],dat_dr.chr.orieff[,1],chromo,bin.size)
  ex=is.na(as.numeric(orilist.chr[,3])) | as.numeric(orilist.chr[,3])==0 
  orilist<-rbind(orilist,orilist.chr[!ex,])

  cat(nrow(orilist), "origins are detected.\n\n")
  
}

if(orimode){
  bedgraph.df <- cbind(orilist[,1], orilist[,2], orilist[,2], as.numeric(orilist[,3])/100)
  
  out.path.ori=paste(prefix, "_origin.bedgraph", sep="")
  header = paste('track name=', prefix, '_Origins description="', version, ';', date(), ' MA(Pol-ratio):', 2*NN+1, '(2x', NN, '+1) MA(diff):', 2*NN2+1, '(2x', NN2, '+1)"', sep="")
  write(header, file=out.path.ori, append = F)
  write.table(bedgraph.df, file=out.path.ori, append = T, quote=F, col.names=F, row.names=F)
}

