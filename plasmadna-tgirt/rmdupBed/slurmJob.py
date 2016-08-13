#!/bin/env python

import fileinput
import argparse

def writeJob(commandlist, jobnamem, commandRank, numberOfJob, numberOfNode, allocation, queue, time):
    commandFile = 'command_%i.bash' %commandRank
    options = \
	"#!/bin/bash \n" +\
	"#SBATCH -J  %s_%i # Job name \n"                             %(jobname, commandRank) +\
	"#SBATCH -N  %i   # Total number of nodes \n"                 %(numberOfNode)+\
	"#SBATCH -n  %i   # Total number of tasks\n"                  %(numberOfJob)+\
	"#SBATCH -p %s    # Queue name \n"                            %(queue)+\
	"#SBATCH -o %s_%i.o%s # Name of stdout output file \n"        %(jobname,commandRank,'%j')+ \
	"#SBATCH -t  %s           # Run time (hh:mm:ss) \n"           %time +\
	"#SBATCH -A %s \n"                                            %(allocation)  +\
        "export LAUNCHER_JOB_FILE=%s\n" %commandFile +\
        "export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins\n" +\
        "export LAUNCHER_RMI=SLURM\n" +\
        "$LAUNCHER_DIR/paramrun\n"
    with open('launcher_%i.slurm' %(commandRank), 'w') as slurmFile:
	slurmFile.write(options)
    with open(commandFile,'w') as commands:
        commands.write('\n'.join(commandlist) + '\n')
    return 0

def main(jobname, commandFile, numberOfJob, numberOfNode, allocation, queue, time):
	with open(commandFile,'ru') as f:
            commands = f.readlines()
            commandlist = []
            i = 0
            commandRank = 0
            for command in commands: 
                commandlist.append(command.strip())
                i += 1
                if i % numberOfJob == 0: 
    		    writeJob(commandlist, jobname, commandRank, numberOfJob, numberOfNode, allocation, queue, time)
                    commandRank += 1
		    i = 0
                    commandlist=[]
                elif i == len(commands)-1:
    		    writeJob(commandlist, jobname, commandRank, i, numberOfNode, allocation, queue, time)
		    i = 0
                    commandlist=[]
	print 'Written %i scripts' %commandRank

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A script to create slurm scripts from list of commands')
	parser.add_argument('-c', '--cmdlst', help='A list of command, each line is a command', required=True)
	parser.add_argument('-j', '--jobname', default='job',help='Jobname (default: job)')
	parser.add_argument('-N', '--numberOfNode', default=1, type=int, help='Number of node for each command (default: 1)')
	parser.add_argument('-n', '--numberOfJob', default=1, type=int, help='Number of job per node (default: 1)')
	parser.add_argument('-A', '--allocation', default = '2013lambowitz', help= 'Account (default: 2013lambowitz)', choices = {'tRNA-profiling-and-b', '2013lambowitz', 'Exosome-RNA-seq'})
	parser.add_argument('-t', '--time', default='01:00:00', help='Run time (hh:mm:ss) default: 1:00:00')
	parser.add_argument('-q','--queue', default='normal',help='Queue (default: normal)')
	args = parser.parse_args()
	cmdlst = args.cmdlst
	jobname = args.jobname
	numberOfJob = args.numberOfJob
	numberOfNode = args.numberOfNode
	allocation = args.allocation
	queue = args.queue	
	time = args.time
	main(jobname, cmdlst, numberOfJob, numberOfNode, allocation, queue, time)
