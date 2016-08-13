#!/bin/evn python

from multiprocessing import Pool
import subprocess
import sys

def runCommand(command):
    command = command.strip()
    subprocess.call(command, shell=True)
    print 'finished', command
    return 0

def main():
    if len(sys.argv) != 2:
        print 'usage: %s <job file>' %(sys.argv[0])
        sys.exit()
    commandFile = sys.argv[1]
    with open(commandFile,'r') as jobs:
        commands = jobs.readlines()
        result = Pool(24).map(runCommand, commands)
    return 0

if __name__ == '__main__':
    main()
