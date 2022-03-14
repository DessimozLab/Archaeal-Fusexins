import glob
import subprocess
import shlex

def runHmmScan( proteome , pfam , ncores = 24 ):
    outfile = proteome + 'hmmerscan.csv'
    args = 'hmmscan --tblout '+ outfile +' --cpu '+ str(ncores) +' ' + pfam + ' ' + proteome
    print(args)
    p = subprocess.run( shlex.split(args) )
    return p , outfile

pfampath = '/home/cactuskid13/pfam/Pfam-A.hmm'
proteomes = glob.glob('./*proteo.fasta')
for proteo in proteomes:
    if proteo+'.hmmerscan.csv' not in glob.glob('./*'):
        p, out = runHmmScan(proteo , pfampath)
    else:
        print( 'done' , proteo+'.hmmerscan.csv' )
