import platform
import subprocess as sp
import time
import re
from os import system
from os import walk
import numpy as np
from os import listdir
from subprocess import call
from functools import partial
from multiprocessing.dummy import Pool
import numpy as np
import time
import random
import sys
from random import randint
from time import sleep

computer_name = platform.node()

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

################################################################################
# Defining commands to run

#sleep(randint(1, 5))

nb_run = 1

NB_PHOTON_PER_RUN = 10000000
ALTITUDE_LIST = [6,7,8,10,12,14,15]

# defining the commands to be run in parallel
commands = []
executable = './TGF_Propa'

for _ in range(nb_run):
    for alt in ALTITUDE_LIST:
	    commands.append(
			executable
			+ ' ' + str(NB_PHOTON_PER_RUN)
			+ ' ' + str(alt)
			)

#####################################

random.shuffle(commands)

# LOCAL RUN (uses python multiprocessing library)
nb_thread = 7  # number of threads (cpu) to run

# Making an array where each element is the list of command for a given thread

command_number = len(commands)

print('Number of commands required ' + str(command_number))

pool = Pool(nb_thread)  # to be always set to 1 for this case
for i, returncode in enumerate(pool.imap(partial(call, shell=True), commands)):
    if returncode != 0:
        print("%d command failed: %d" % (i, returncode))
