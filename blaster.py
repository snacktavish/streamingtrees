import sys
from subprocess import call
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import dendropy

#infi=sys.argv[1] #query info should be in non interleaved nexus format
#ancseqs=sys.argv[2] #test sequences should be in fasta format

allseq=open("rbcL_noEPI.nex").readlines()
treseq=open("rbcL_noEPI_subA.nex").readlines()


nams=set()
for lin in treseq[6:-2]:
   nams.add(lin.split()[0])


identdict={}
testseqs=[]
for lin in allseq[6:-3]:
   if lin.split('\t')[0] not in nams:
     testseqs.append(lin.split('\t'))
     identdict[lin.split('\t')[0]]={}


#for i,lin in enumerate(open(infi)):
#   if i>=6:
for lin in testseqs:
      nam=lin[0]
      assert(nam in identdict)
      fi=open("query.txt",'w')
      fi.write(lin[1]) 
      fi.close()
      for i in range(1,11):
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
    tlst = dendropy.Tree.get_from_path('test.tre', "newick")
    tlst.find_node_with_label(sumdict[-1]).new_child(taxon=dendropy.dataobject.taxon.Taxon(taxon=None, label=nam, oid=None))


