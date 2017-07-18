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

	parser = OptionParser(usage = "Usage: python GGisy.py -g1 genome1.fna -g2 genome2.fna")
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
		print "* No genome was provided (-g1), use -h for help"
		sys.exit()
	else:
		if os.path.isfile(genome1) == False:
			print "*",genome1," doesn't exist"
			sys.exit()

	if not genome2 or genome2 is None:
		print "* its mandatory provide 2 genomes (-g2), use -h for help"
		sys.exit()
	else:
		if os.path.isfile(genome2) == False:
			print "* ",genome2," doesn't exist"
			sys.exit()

	if blastout != None:
		if os.path.isfile(blastout) == False:
			print "* ", blastout, "not found, check if file exist or let the program do the blast omiting this option (-b)"
			sys.exit()

	blastBIN=which("blastn")
	if blastBIN == None:
		print "No blastn was found, install it before continue (make sure is in your $PATH)"
		sys.exit()

	makeblastBIN=which("makeblastdb")
	if makeblastBIN == None:
		print "No makeblastdb was found, install it from blast+ (make sure is in your $PATH)"
		sys.exit()

	rscriptBIN=which("Rscript")
	if rscriptBIN == None:
		print "No Rscript was found, make sure is in your $PATH"
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
	subprocess.call(["blastn", "-query", genome2, "-db", "ref", "-evalue", evalue, "-outfmt", "6", "-strand", "both", "-num_threads", threads, "-out", "tmp.tsv"])

	return str("tmp.tsv")

def filterBlastOutput(blastout,alignL,evalue,identity):

	PARSED=open("parsed.tsv",'w') #overwrite if exist
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

def handleR(conn, reference, query):

	plotstep=open("handle.R", 'w')
	plotstep.write("""rm(list=ls());
library(OmicCircos)
library(RColorBrewer)
args<-commandArgs()
handlefile<-as.character(args[6])
refname<-as.character(args[7])
queryname<-as.character(args[8])
handle<-read.table(handlefile,sep = "\\t",stringsAsFactors = F,check.names = F)
ref<-read.table(refname,sep = "\\t",stringsAsFactors = F,check.names = F)
query<-read.table(queryname,sep = "\\t", stringsAsFactors = F,check.names = F)
rownames(ref)<-ref$V1
rownames(query)<-query$V1
qryUniq<-unique(sort(handle$V1))
refUniq<-unique(sort(handle$V2))
gl<-sum(query[qryUniq,3])
gl<-gl+sum(ref[refUniq,3])
ref<-ref[refUniq,]
query<-query[qryUniq,]
data<-rbind(ref,query)
lowId<-min(handle$V3)
link.pg.v<-data.frame(seg1=handle$V1, start1=handle$V4, end1=handle$V5, seg2=handle$V2, start2=handle$V6, end2=handle$V7)
measure<-data.frame(seg.name=seq(1:8), seg.start=seq(1:8), seg.end=seq(1:8)+1, seg.value=1)
measure2<-data.frame(seg.name=1, seg.start=1, seg.end=signif(2*gl/12/1000), seg.value=round(2*gl/12/1000,0))
measure4<-data.frame(seg.name=1, seg.start=1, seg.end=signif(4*gl/12/1000), seg.value=round(4*gl/12/1000,0))
measure8<-data.frame(seg.name=1, seg.start=1, seg.end=signif(8*gl/12/1000), seg.value=round(8*gl/12/1000,0))
measure10<-data.frame(seg.name=1, seg.start=1, seg.end=signif(10*gl/12/1000), seg.value=round(10*gl/12/1000,0))
measure["V5"]<-measure2["V5"]<-measure4["V5"]<-measure8["V5"]<-measure10["V5"]<-1
data["V5"]<-data["V4"]<-1
colnames(data)<- c("chr", "start", "end","V4","V5")
tocirmeasure<-segAnglePo(measure, seg=c(as.matrix(measure["seg.name"][,])))
tocirmeasure[2,2]<-300;tocirmeasure[2,3]<-620
tocirmeasure[3,2]<-0;tocirmeasure[3,3]<-620
tocirmeasure[4,2]<-420;tocirmeasure[4,3]<-620
tocirmeasure[5,2]<-90;tocirmeasure[5,3]<-620
tocirmeasure[6,2]<-120;tocirmeasure[6,3]<-620
tocirmeasure[7,2]<-180;tocirmeasure[7,3]<-620
tocirmeasure[8,2]<-240;tocirmeasure[8,3]<-620
tocirmeasure2<-segAnglePo(measure2, seg=c(as.matrix(measure2["seg.name"][,])))
tocirmeasure2[1,2]<-330
tocirmeasure4<-segAnglePo(measure4, seg=c(as.matrix(measure4["seg.name"][,])))
tocirmeasure4[1,2]<-395
tocirmeasure8<-segAnglePo(measure8, seg=c(as.matrix(measure8["seg.name"][,])))
tocirmeasure8[1,2]<-146
tocirmeasure10<-segAnglePo(measure10, seg=c(as.matrix(measure10["seg.name"][,])))
tocirmeasure10[1,2]<-210
tocir <- segAnglePo(data, seg=c(as.matrix(data["chr"][,])))
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
colors<-addalpha(colors,0.8)
link.pg.v[,"colors"]<-addalpha(scaleColors(handle$V3),0.8)
pdf(file="synteny.pdf", width = 10, height =10)
par(mar=c(2,2,2,2))
xorigin=600
yorigin=740
plot(c(0,1500), c(0,1500), type="n", axes=FALSE, xlab="", ylab="", main="")
circos(R=450, cir=tocir, W=10,type="chr", print.chr.lab=T, scale=F,xc = xorigin,yc = yorigin,
       col = c(rep("dark blue",nrow(ref)),rep("#FEE496",nrow(query))),cex = 10)
circos(R=445, cir=tocir, W=10, mapping=link.pg.v , type="link.pg", lwd=1, col=link.pg.v$colors,xc = xorigin,yc = yorigin)
legend(x = 1230, y=1300, legend = c(unlist(strsplit(refname,"_info\\\.tsv"))[1],unlist(strsplit(queryname,"_info\\\.tsv"))[1]),
       ncol = 1, cex = 0.8,  bty="n",
       fill=c("dark blue","#FEE496"),
       border = c("dark blue","#FEE496"),text.width=c(0.5,0.5),
       title=expression(bold("Sequences")))
legend(x = 1200, y=1100, legend = c("100","","","","","","","","","",(100-lowId)/2 + lowId,"","","","","","","","",lowId),
       ncol = 1, cex = 0.8,  bty="n",
       fill=colors,
       border = colors,
       y.intersp = 0.5,
       x.intersp = 0.5,text.width=c(0.5,0.5),
       title=expression(bold("Identity percent\\n")))

dev.off()""")
	plotstep.close()
	subprocess.call(["Rscript", "handle.R", conn, reference, query, "--vanilla"])


def cleanfiles(ginfo1, ginfo2):

	if os.path.isfile("parsed.tsv"):
		os.remove("parsed.tsv")
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
	handleR(conn="parsed.tsv",reference=ref, query=que)
	
	if mainV.v8 == True:
		cleanfiles(ref,que)

	sys.exit()




