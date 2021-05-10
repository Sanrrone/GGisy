from __future__ import with_statement 

# ==============================================================================
# 						GGisy (python v2.7)
#
# Author: Sandro Valenzuela (sandrolvalenzuead@gmail.com)
# Bugs and errors: https://github.com/Sanrrone/GGisy/issues
#
# Please type "python GGisy.py -h" for usage help
#
# ==============================================================================

__author__ = 'Sandro Valenzuela (sandrolvalenzuead@gmail.com)'
__version__ = '1.0'

import sys, os, subprocess, glob, csv, collections
from optparse import OptionParser
from operator import itemgetter
from Bio import SeqIO


def main():

	parser = OptionParser(usage = "Usage: python GGisy.py -r genome1.fna -q genome2.fna")
	parser.add_option("-r","--reference",dest="genome1",help="First genome to be used as reference", default=None)
	parser.add_option("-q","--query",dest="genome2",help="Second genome to be used as query against the first genome (-r)", default=None)
	parser.add_option("-l","--alignmentLength",dest="alignL",help="Aligment length cutoff in blast output [default: 1000]",default=1000)
	parser.add_option("-e","--evalue",dest="evalue",help="E-value cutoff for blastn search [default: 1e-3]",default=1e-3)
	parser.add_option("-i","--identity",dest="Identity",help="Identity cutoff on the blastn alignment to consider the region [default: 50]",default=50)
	parser.add_option("-t","--threads",dest="Threads",help="Number of threads to be used for blast [default: 4]",default=4)
	parser.add_option("-b","--blastout",dest="Blastout",help="Blast output file to be used instead doing it [default: none]",default=None)
	parser.add_option("-c","--clean",dest="clean",help="clean files after execution [default: True]",default=True)

	(options,args) = parser.parse_args()

	genome1 = str(options.genome1)
	genome2 = str(options.genome2)
	alignL= int(options.alignL)
	evalue= str(options.evalue)
	Identity= int(options.Identity)
	threads= str(options.Threads) #for subcallproccess must be str()
	blastout= options.Blastout #dont cast to str
	cleanf=options.clean

	#check variables
	if not genome1 or genome1 is None:
		print("* No genome was provided (-g1), use -h for help")
		sys.exit()
	else:
		if os.path.isfile(genome1) == False:
			print("*",genome1," doesn't exist")
			sys.exit()

	if not genome2 or genome2 is None:
		print("* its mandatory provide 2 genomes (-g2), use -h for help")
		sys.exit()
	else:
		if os.path.isfile(genome2) == False:
			print("* ",genome2," doesn't exist")
			sys.exit()

	if blastout != None:
		if os.path.isfile(blastout) == False:
			print("* ", blastout, "not found, check if file exist or let the program do the blast omiting this option (-b)")
			sys.exit()

	blastBIN=which("blastn")
	if blastBIN == None:
		print("No blastn was found, install it before continue (make sure is in your $PATH)")
		sys.exit()

	makeblastBIN=which("makeblastdb")
	if makeblastBIN == None:
		print("No makeblastdb was found, install it from blast+ (make sure is in your $PATH)")
		sys.exit()

	rscriptBIN=which("Rscript")
	if rscriptBIN == None:
		print("No Rscript was found, make sure is in your $PATH")
		sys.exit()

	Inputs = collections.namedtuple('Inputs', ['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8'])
	I = Inputs(genome1, genome2, alignL, evalue, Identity, threads, blastout, cleanf)
	return I

def which(program): #function to check if some program exists 
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

def blasting(genome1, genome2, evalue, threads):
	#searching for blast binaries

	subprocess.call(["makeblastdb", "-in", genome1, "-input_type", "fasta", "-dbtype", "nucl", "-out", "ref"])
	subprocess.call(["blastn", "-query", genome2, "-db", "ref", 
		"-evalue", evalue, "-outfmt", "6", "-strand", "both", 
		"-num_threads", threads, "-out", "tmp.tsv"])

	return str("tmp.tsv")

def filterBlastOutput(blastout,alignL,evalue,identity):

	PARSED=open("synteny.tsv",'w') #overwrite if exist
	with open(blastout) as tsvfile:
		tsvreader = csv.reader(tsvfile, delimiter="\t")
		for line in tsvreader:

			#formula line [n-1:n]
			toint = [int(i) for i in line[3:4]]
			if toint[0] >= alignL:
				toint = [float(i) for i in line[2:3]]
				if toint[0] >= float(identity):
					PARSED.write("\t".join(map(str, line[0:3]+line[6:10]))+"\n")

	PARSED.close()


def parsingGenomes(genome):

	gname = genome.split('/')[-1]
	PARSED=open(str(gname+"_info.tsv"),'w') #overwrite if exist
	fasta_sequences = SeqIO.parse(open(genome),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		lengthSeq= len(sequence)

		PARSED.write("%s\t1\t%s\n" % (name, lengthSeq))

	PARSED.close
	return str(gname+"_info.tsv")

def handleR(conn, reference, query, alignL):

	plotstep=open("handle.R", 'w')
	plotstep.write("""rm(list=ls());
library(OmicCircos)
library(RColorBrewer)
library(varhandle)
args<-commandArgs()
handlefile<-as.character(args[6])
refname<-as.character(args[7])
queryname<-as.character(args[8])
filterl<-as.numeric(args[9])
handle<-read.table(handlefile,sep = "\\t",stringsAsFactors = F,check.names = F)
ref<-read.table(refname,sep = "\\t",stringsAsFactors = F,check.names = F)
query<-read.table(queryname,sep = "\\t", stringsAsFactors = F,check.names = F)
rownames(ref)<-ref$V1
rownames(query)<-query$V1
qryUniq<-unique(sort(handle$V1))
refUniq<-unique(sort(handle$V2))
ref<-ref[refUniq,]
ref<-ref[with(ref, order(-V3, V1)), ]
query<-query[qryUniq,]
query<-query[with(query, order(+V3, V1)), ]
data<-rbind(ref,query)
refname<-unlist(strsplit(refname,"_info.tsv"))[1]
queryname<-unlist(strsplit(queryname,"_info.tsv"))[1]
lowId<-min(handle$V3)
fhand<-handle[handle$V6<handle$V7,]
rhand<-handle[handle$V6>handle$V7,]
linkf<-data.frame(seg1=fhand$V1, start1=fhand$V4, end1=fhand$V5, seg2=fhand$V2, start2=fhand$V6, end2=fhand$V7, stringsAsFactors = F)
linkr<-data.frame(seg1=rhand$V1, start1=rhand$V4, end1=rhand$V5, seg2=rhand$V2, start2=rhand$V6, end2=rhand$V7, stringsAsFactors = F)
#fix reverse positions
for(i in 1:nrow(linkr)){
  contign<-linkr[i,4]
  contigl<-ref[contign,3]
  linkr[i,5]<- contigl-linkr[i,5]+1
  linkr[i,6]<- contigl-linkr[i,6]+1
}
data["V5"]<-data["V4"]<-1
colnames(data)<- c("chr", "start", "end","V4","V5")
tocir <- segAnglePo(data, seg=data$chr)
gl<-sum(data$end)+nrow(data)
maxangr<-270+(350/gl)*sum(ref$V3)
spacer<-maxangr/(maxangr-270)/nrow(ref)

for(i in 1:nrow(ref)){
  #358 is the total angles (aviable) for all
  tocir[i,"angle.end"]<-as.character(as.numeric(tocir[i,"angle.start"]) + (350/gl)*as.numeric(tocir[i,7]))
  tocir[i+1,"angle.start"]<-as.character(as.numeric(tocir[i,"angle.end"])+spacer)
}
tocir[i+1,"angle.start"]<-as.character(as.numeric(tocir[i+1,"angle.start"])+2.5)
tocir[i+1,"angle.end"]<-as.character(as.numeric(tocir[i+1,"angle.start"]) + (350/gl)*as.numeric(tocir[i+1,7]))

maxangq<-628-maxangr
spacer<-628/maxangq/nrow(query)
if(nrow(ref)+2>=nrow(tocir)){
  i<-nrow(tocir)
  tocir[i,"angle.start"]<-as.character(as.numeric(tocir[i-1,"angle.end"])+spacer)
  tocir[i,"angle.end"]<-as.character(628)
}else{
  for(i in (nrow(ref)+2):nrow(tocir)-1){
    #358 is the total angles (aviable) for all
    tocir[i,"angle.end"]<-as.character(as.numeric(tocir[i,"angle.start"]) + (350/gl)*as.numeric(tocir[i,7]))
    tocir[i+1,"angle.start"]<-as.character(as.numeric(tocir[i,"angle.end"])+spacer)
  } 
}
refang<-as.numeric(tocir[1:nrow(ref),2])
qryang<-as.numeric(tocir[(nrow(ref)+1):(nrow(ref)+nrow(query)),2])
maxangr<-max(refang)
maxangq<-max(qryang)
faketocir <- tocir
faketocir[,1]<-""
maxangr<-max(refang)
for(i in 1:nrow(tocir)){
  if(270+(maxangr-270)/2<as.numeric(tocir[i,2])){
    break
  }
}
faketocir[i,1]<-refname
maxangq<-max(qryang)
for(i in 1:nrow(tocir)){
  if(maxangr+(maxangq-maxangr)/2<as.numeric(tocir[i,2])){
    break
  }
}
faketocir[i,1]<-queryname
colors<-rev(colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(20))
delta<-(100-lowId)/20
scaleColors<- function(x){
  cArray<-c()
  for(id in x){
    for(i in 1:20){
      if(id>=100-(delta*i)){
        break
      }
    }
    cArray<-c(cArray,colors[i])
  }
  return(cArray)
}
addalpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))
}
black<-addalpha("#000000",0.7)
colors<-addalpha(colors,1)
try({
  linkf[,"colors"]<-addalpha(scaleColors(fhand$V3),1)
},silent = T)
try({
  linkr[,"colors"]<-addalpha(scaleColors(rhand$V3),1)
},silent = T)


pdf(file="synteny.pdf", width = 10, height =10)

if(nrow(data)<=20){
  par(mar=c(2,2,2,2))
  xorigin=700
  yorigin=1000
  plot(c(0,2000), c(0,2000), type="n", axes=FALSE, xlab="", ylab="", main="")
  
  circos(R=450, cir=tocir, W=10,type="chr", print.chr.lab=T, scale=F,xc = xorigin,yc = yorigin,
           col = c(rep("dark blue",nrow(ref)),rep("#FEE496",nrow(query))),cex = 5)

  if(nrow(linkf)>0){
    circos(R=440, cir=tocir, mapping=linkf , type="link.pg", lwd=0.5, col=linkf$colors,xc = xorigin,yc = yorigin)
  }  
  if(nrow(linkr)>0){
    circos(R=440, cir=tocir, mapping=linkr , type="link.pg", lwd=0.5, col=linkr$colors,xc = xorigin,yc = yorigin)
    newlinkr<-linkr
    newlinkr$start1<-newlinkr$start1+as.integer((newlinkr$end1-newlinkr$start1)/2)+1
    newlinkr$start2<-newlinkr$start2+as.integer((newlinkr$end2-newlinkr$start2)/2)-1
    circos(R=440, cir=tocir, W=10, mapping=newlinkr , type="link", lwd=0.6, col=black,xc = xorigin,yc = yorigin)  
  }

  
  legend(x = 1500, y=1700, legend = c(refname,queryname),
         ncol = 1, cex = 0.8,  bty="n",
         fill=c("dark blue","#FEE496"),
         border = c("dark blue","#FEE496"),text.width=c(0.5,0.5),
         title="Sequences")
  legend(x = 1430, y=1500, legend = c(paste("Reference: ", nrow(ref), " (", sum(ref$V3), " bp)", sep = ""), paste("Query: ",nrow(query), " (", sum(query$V3), " bp)", sep="")),
         ncol = 1, cex = 0.8,  bty="n",
         fill=c("dark blue","#FEE496"),
         border = c("dark blue","#FEE496"),text.width=c(0.5,0.5),
         title=paste("Contigs align >= ", filterl, " bp", sep=""))
  legend(x = 1520, y=1300, legend = c("Forward","Reverse"),lty = c(0,1),merge=T,seg.len = 0.6,
         ncol = 1, cex = 0.8,  bty="n",
         fill="white",
         border = "black",text.width=c(0.5,0.5),
         title="Strand Match\n(on reference)")
  legend(x = 1505, y=1100, legend = c("100","","","","","","","","","",(100-lowId)/2 + lowId,"","","","","","","","",lowId),
         ncol = 1, cex = 0.8,  bty="n",
         fill=colors,
         border = colors,
         y.intersp = 0.5,
         x.intersp = 0.5,text.width=c(0.5,0.5),
         title="Identity percent\n")
}else{
  par(mar=c(2,2,2,2))
  xorigin=750
  yorigin=550
  plot(c(0,1500), c(0,1500), type="n", axes=FALSE, xlab="", ylab="", main="")

  circos(R=450, cir=faketocir, W=10,type="chr", print.chr.lab=T, scale=F,xc = xorigin,yc = yorigin,
         col = "white")
  circos(R=410, cir=tocir, W=10,type="chr", print.chr.lab=F, scale=F,xc = xorigin,yc = yorigin,
         col = c(rep("dark blue",nrow(ref)),rep("#FEE496",nrow(query))),cex = 5)
  
  if(nrow(linkf)>0){
    highlightr <- c(420, 450, tocir[1,1], 1, tocir[nrow(ref),1], tocir[nrow(ref),7], "dark blue", NA)
    circos(cir=tocir, mapping=highlightr, type="hl",xc = xorigin,yc = yorigin)
    circos(R=400, cir=tocir, mapping=linkf , type="link.pg", lwd=0.5, col=linkf$colors,xc = xorigin,yc = yorigin)
  }
  if(nrow(linkr)>0){
    highlightq <- c(420, 450, query[1,1], 1, query[nrow(query),1], query[nrow(query),3], "#FEE496", NA)
    circos(cir=tocir, mapping=highlightq, type="hl",xc = xorigin,yc = yorigin)
    circos(R=400, cir=tocir, mapping=linkr , type="link.pg", lwd=0.5, col=linkr$colors,xc = xorigin,yc = yorigin)
    newlinkr<-linkr
    newlinkr$start1<-newlinkr$start1+as.integer((newlinkr$end1-newlinkr$start1)/2)+1
    newlinkr$start2<-newlinkr$start2+as.integer((newlinkr$end2-newlinkr$start2)/2)-1
    circos(R=400, cir=tocir, W=10, mapping=newlinkr , type="link", lwd=0.3, col=black,xc = xorigin,yc = yorigin)
    
  }
  legend(x = 210, y=1500, legend = c(paste("Reference: ", nrow(ref), " (", sum(ref$V3), " bp)", sep = ""), paste("Query: ",nrow(query), " (", sum(query$V3), " bp)", sep="")),
         ncol = 1, cex = 0.8,  bty="n",
         fill=c("dark blue","#FEE496"),
         border = c("dark blue","#FEE496"),text.width=c(0.5,0.5),
         title=paste("Contigs align >= ", filterl, " bp", sep=""))
  legend(x = 270, y=1300, legend = c("Forward","Reverse"),lty = c(0,1),merge=T,seg.len = 0.6,
         ncol = 1, cex = 0.8,  bty="n",
         fill="white",
         border = "black",text.width=c(0.5,0.5),
         title="Strand Match\\n(on reference)")
  legend(x = 990, y=1500, legend = c("100","","","","","","","","","",(100-lowId)/2 + lowId,"","","","","","","","",lowId),
         ncol = 1, cex = 0.8,  bty="n",
         fill=colors,
         border = colors,
         y.intersp = 0.5,
         x.intersp = 0.5,text.width=c(0.5,0.5),
         title="Identity percent\\n")
}
dev.off()""")
	plotstep.close()
	subprocess.call(["Rscript", "handle.R", conn, reference, query, str(alignL), "--vanilla"])


def cleanfiles(ginfo1, ginfo2):

	if os.path.isfile("tmp.tsv"):
		os.remove("tmp.tsv")
	if os.path.isfile("ref.nin"):
		os.remove("ref.nin")
	if os.path.isfile("ref.nsq"):
		os.remove("ref.nsq")
	if os.path.isfile("ref.nhr"):
		os.remove("ref.nhr")
	if os.path.isfile("handle.R"):
		os.remove("handle.R")
	if os.path.isfile(ginfo1):
		os.remove(ginfo1)
	if os.path.isfile(ginfo2):
		os.remove(ginfo2)

if __name__ == '__main__':
	mainV=main()
	blastout=mainV.v7
	if blastout is None:
		blastout=blasting(genome1=mainV.v1, genome2=mainV.v2, evalue=mainV.v4, threads=mainV.v6)

	filterBlastOutput(blastout=blastout, alignL=mainV.v3, evalue=mainV.v4, identity=mainV.v5)
	ref=parsingGenomes(genome=mainV.v1)
	que=parsingGenomes(genome=mainV.v2)
	handleR(conn="synteny.tsv",reference=ref, query=que, alignL=mainV.v3)
	
	if mainV.v8 == True:
		cleanfiles(ref,que)

	sys.exit()





