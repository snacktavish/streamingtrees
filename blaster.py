import sys
from subprocess import call
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

#WTFFFFFF
#infi=sys.argv[1] #query info should be in non interleaved nexus format
#ancseqs=sys.argv[2] #test sequences should be in fasta format

allseq=open("rbcL.nex").readlines()
treseq=open("rbcL_sub2.nex").readlines()


nams=set()
for lin in treseq[6:-3]:
   nams.add(lin.split('\t')[0])


testseqs=[]
for lin in allseq[6:-3]:
   if lin.split('\t')[0] not in nams:
     testseqs.append(lin)


seqnams=set()
identdict={}

#for i,lin in enumerate(open(infi)):
#   if i>=6:
for item in testseqs:
      nam=lin.split()[0]
      if nam not in identdict:
            identdict[nam]={}
      seqnams.add(nam)
      fi=open("query.txt",'w')
      fi.write(lin.split()[1]) 
      fi.close()
      for i in range(1,2):
        blastn_run = NcbiblastnCommandline(query="query.txt", subject="ancestorT{}.fasta".format(i), evalue=0.001, outfmt=5, out="test.xml")
        stdout, stderr = blastn_run()
        result_handle = open("test.xml")
        blast_records = NCBIXML.parse(result_handle)
        for rec in blast_records:
           for alignment in rec.alignments:
               node=alignment.title.split()[-1]
               if node not in identdict[nam]:
                   identdict[nam][node]=[]
               for hsp in alignment.hsps:
                  identdict[nam][node].append(hsp.bits)


for nam in identdict:
    sumdict=[]
    for node in identdict[nam]:
         sumdict.append((sum(identdict[nam][node]),node))
    sumdict.sort()
    print(nam,sumdict[0])
  
