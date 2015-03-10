####Sliding window for downsampling high-res EV data
####KS 02/2015
####Specific goal is to downsample Gacu PHASTCONs data from UCSC genome
####input file format:- csv: lg,pos1,pos2,phastcons (or any ev)

#great R library or greatest R library?
library(data.table)

#WINDOWS
home.dir<-"E:/Genome Meta Analysis/evs/nonwindow/conservation"

#MAC OS
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/nonwindow/conservation"
setwd(home.dir)

#find all the csv files in the home dir
file.list<-grep(".csv",value=TRUE,list.files())

#define window size (in base pairs)

window.size<-75000
ev.name<-"phastcons"

#calc windowed stats for each coordinate run

ev.chr.window.list<-list()
for (i in 1: length(file.list)){
  
  ev.chr<-fread(file.list[i])
  num.windows<-as.integer(max(ev.chr$pos1)/window.size)
  nubbin<-max(ev.chr$pos2)-(window.size*num.windows)
  if(nubbin>0){
    num.windows<-num.windows+1
  }
  
  print (paste("processing",num.windows,"windows","for lg",i))
  lg<-rep(NA,num.windows)
  ev.chr.window<-data.frame(lg=lg,pos1=NA,pos2=NA,ev=NA)
  for (j in 1:(num.windows)){
    
    if (j==1){
      window.range <- ((0*window.size+1):((0*window.size)+window.size))
    } else if ((j==num.windows)&&(nubbin>0)){
      window.range<-((max(ev.chr$pos2)-nubbin):(max(ev.chr$pos2)))
      print (paste("lg",i,"has a nubbin of",nubbin,"base pairs"))
    } else{
      window.range <- (((j-1)*window.size):((j*window.size)))+(j-1)
    }
    #print(paste("window",j,min(window.range),":",max(window.range)))  
    
    ev.chr.window$lg[j]<-i
    ev.chr.window$ev[j]<-mean(ev.chr[window.range]$phastcons)
    ev.chr.window$pos1[j]<-min(window.range)
    ev.chr.window$pos2[j]<-max(window.range)
    
  }
  names(ev.chr.window)[4]<-ev.name
  ev.chr.window.list[[i]]<-ev.chr.window
}


ev.chr.out<-data.frame(lg=i,pos1=ev.chr.wind.pos1,pos2=ev.chr.wind.pos2,ev=ev.chr.means)
names(ev.chr.out)<-ev.name
write.table(ev.chr.out,file=paste("lg",i,"window",sep="_"))



#for (i in 1:length(file.list)){
  
}
