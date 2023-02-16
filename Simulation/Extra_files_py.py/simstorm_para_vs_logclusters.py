import argparse
import random
import time
import math

def current_milli_time():
    return round(time.time() * 1000)

#--max-ins 10000
def get_basis(invade):
    return " {0} -no-x-cluins --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --x 0.01 --rr 4,4,4,4,4 --rep 1 --u 0.2 --basepop 100 --steps 5000 --silent".format(invade)

def get_filter():
    return """|grep -v "^Invade"|grep -v "^#" """

def get_rand_para():
    r=random.randint(1,500)
    a=""
    for x in range(0,r):
        a+=f"{x},"
    a=a[:-1]
    a = a + " "
    return f"500:{a}"

def get_rand_clusters():
    r=math.floor(10**random.uniform(3.69899,5.69899))
    return f"{r},{r},{r},{r},{r}"


def run_cluster_negsel(invade,count,output):
    """
    TE invasion that is stopped by cluster insertions and neg selection against TEs
    """
    basis =get_basis(invade) 
    commandlist=[]
    for i in range(0,count):
        x=get_rand_clusters()
        u=get_rand_para()
        tr=current_milli_time()+i
        #command=basis+" --x {0} --paramutation {1} --replicate-offset {2} --file-tally tally{2}.txt --file-sfs sfs{2}.txt --file-mhp mhp{2}.txt --seed {3} ".format(x,u,i,tr)
        command=basis+" --cluster bp:{0} --paramutation {1} --replicate-offset {2} --seed {3} ".format(x,u,i,tr)
        #ri=random.random()
        arr = u.split(',')
        lastWord_para = arr[- 1]
        arr_2 = x.split(',')
        lastWord_cluster = arr_2[- 1]
        command+= "--sampleid \"{0}\t{1}\"  {2} > {3}{4}".format(lastWord_cluster,lastWord_para,get_filter(),output,i)
        commandlist.append(command)
    return commandlist;


parser = argparse.ArgumentParser(description="""           
Description
-----------
    Simulation storm""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

Authors
-------
    Robert Kofler
    Filip Wierzbicki
    Almorò Scarpa
""")


parser.add_argument("--number", type=int, required=True, dest="count", default=None, help="the number of simulations")
parser.add_argument("--threads", type=int, required=True, dest="threads", default=None, help="the threads of simulations")
parser.add_argument("--output", type=str, required=True, dest="output", default=None, help="the outputfile of simulations")
parser.add_argument("--invade", type=str, required=True, dest="invade", default=None, help="the invade.go")
# parser.add_argument("--type", type=str, required=True, dest="type", default=None, help="nostop|ns|nst|cluster|nscluster")
parser.add_argument("--silent",action="store_true", dest="silent", help="be quit; default=False")



args = parser.parse_args()
#0   1    2    3     4     5    6       7    8    9           10     11   1     13   14        15   16   17   18   19   20       21   22        23 
#1	0	|	0.01	1.00	0.01	0.0005	0	|	0.0000	0.00	0.0000	0	|	0.0000	10	0	0.10	0.00	1000	|	u=0.01	ci=1		1
#1	100	|	0.11	1.00	0.12	0.0074	0	|	0.0000	0.00	0.0000	0	|	0.0010	8	0	0.35	0.00	1000	|	u=0.01	ci=1		1
#1	200	|	0.40	1.00	0.50	0.0075	0	|	0.0270	0.03	0.0070	0	|	0.0060	33	2	0.70	0.17	1000	|	u=0.01	ci=1		1
#1	300	|	0.84	1.00	1.76	0.0092	0	|	0.2700	0.30	0.0370	0	|	0.0090	96	4	1.33	0.51	1000	|	u=0.01	ci=1		1

commandlist=run_cluster_negsel(args.invade,args.count,args.output)

"""
if(args.type=="nostop"):
    commandlist=run_nostop(args.invade,args.count,args.output)
elif(args.type=="ns"):
    commandlist=run_negsel(args.invade,args.count,args.output)
elif(args.type=="nst"):
    commandlist=run_negsel_t(args.invade,args.count,args.output)
elif(args.type=="cluster"):
    commandlist=run_cluster(args.invade,args.count,args.output)
elif(args.type=="nscluster"):
    commandlist=run_cluster_negsel(args.invade,args.count,args.output)
else:
    raise Exception("Unknown simulation type")
"""
    
    



 
def submit_job_max_len(commandlist, max_processes):
    import subprocess
    import time
    sleep_time = 10.0
    processes = list()
    for command in commandlist:
        if(not args.silent):
            print ('running {n} processes. Submitting {proc} '.format(n=len(processes),proc=str(command)))
        processes.append(subprocess.Popen(command, shell=True, stdout=None))
        while len(processes) >= max_processes:
            time.sleep(sleep_time)
            processes = [proc for proc in processes if proc.poll() is None]
    while len(processes) > 0:
        time.sleep(sleep_time)
        processes = [proc for proc in processes if proc.poll() is None]
 
submit_job_max_len(commandlist, max_processes=args.threads)
print ("Done")
    


    








