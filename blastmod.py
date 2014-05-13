import sys
from subprocess import call
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import dendropy



#call(["python","convertformat.py","-i","nexus","fasta","ancestorT{:d}.nex".format(i)])

def blaster(ancseqs,queryseq,nancestors=5): #both fasta files. Ancseqs should be list of fastsfiles
    dna = dendropy.DnaCharacterMatrix.get_from_path(queryseq, 'fasta')
    qseq = dna.as_string('fasta').split()[3]
    fi=open("query.txt",'w')
    fi.write(qseq) 
    fi.close()
    identdict={}
    for anc in ancseqs:
        blastn_run = NcbiblastnCommandline(query="query.txt", subject=anc, evalue=0.001, outfmt=5, out="tmp.xml")
        stdout, stderr = blastn_run()
        result_handle = open("tmp.xml")
        blast_records = NCBIXML.parse(result_handle)
        for rec in blast_records:
           for alignment in rec.alignments:
               node=alignment.title.split()[-1]
               if node not in identdict:
                   identdict[node]=[]
               for hsp in alignment.hsps:
                  identdict[node].append(hsp.bits)
    sumdict=[]
    for node in identdict:
         sumdict.append((sum(identdict[node]),node))
    sumdict.sort()
    return(sumdict[-1])


#t1=dendropy.Tree.get_from_path("rep_1.tre", "newick")  
#nad=t1.find_node_with_label(sumdict[-1][-1])
#.new_child(taxon=dendropy.Taxon(label=nam),edge_length=0.15)
#nad.edge.length


#t1.write(open("addition.tre",'w'),schema="nexus")

#call(["Garli","score.conf"])


