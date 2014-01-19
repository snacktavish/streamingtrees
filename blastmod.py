import sys
from subprocess import call
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import dendropy

#infi=sys.argv[1] #query info should be in non interleaved nexus format
#ancseqs=sys.argv[2] #test sequences should be in fasta format

def blaster(allseq,queryseq): #both nexus files?

     allseq=open("rbcL_noEPI.nex").readlines()
treseq=open("rbcL_noEPI_1.nex").readlines()


nams=set()
for lin in treseq[6:-3]:
   nams.add(lin.split('\t')[0])


testseqs=[]
for lin in allseq[6:-3]:
   if lin.split('\t')[0] not in nams:
     testseqs.append(lin)


seqnams=set()
identdict={}

for item in testseqs:
      nam=item.split()[0]
      if nam not in identdict:
            identdict[nam]={}
      seqnams.add(nam)
      fi=open("query.txt",'w')
      fi.write(item.split()[1]) 
      fi.close()
      for i in range(1,2):
        call(["python","convertformat.py","-i","nexus","fasta","ancestorT{:d}.nex".format(i)])
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
    print(nam,sumdict[-1])


t1=dendropy.Tree.get_from_path("rep_1.tre", "newick")  
nad=t1.find_node_with_label(sumdict[-1][-1])
#.new_child(taxon=dendropy.Taxon(label=nam),edge_length=0.15)
nad.edge.length


t1.write(open("addition.tre",'w'),schema="nexus")

call(["Garli","score.conf"])


