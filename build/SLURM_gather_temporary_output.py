import subprocess
import re
import os
import getpass

##
def execute_bash_command(bashCommand):
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, _ = process.communicate()
    return output

def execute_bash_command_alternative(bashCommand):
    output = os.system(bashCommand)
    print(output)

## Find list of ID of currently running slurm jobs
username = getpass.getuser()
bashCommand = "squeue -u " + username
print(bashCommand)
output = execute_bash_command(bashCommand)
splited_lines = output.split('\n')
# print(splited_lines)
#
re1='(\\d+)'	# Integer Number 1
re2='(\\s+)'	# White Space 1
re3='((?:[a-z][a-z]+))'	# Word 1
rg = re.compile(re1+re2+re3,re.IGNORECASE|re.DOTALL)

list_jobID=[]

for line in splited_lines:
    # if the status is anything except " R "
    if (" PD   " in line) or (" CG   " in line) or (" CD   " in line) or (" F   " in line) or (" PD   " in line) or (" PR   " in line) or (" S   " in line) or (" ST   " in line):
        continue
    if "TGF-TEB-" in line:
        continue

    txt=line
    m = rg.search(txt)
    if m:
        int1=m.group(1)
        ws1=m.group(2)
        word1=m.group(3)
        print(int1)
        list_jobID.append(int1)

## For each job, gather particles and send
for jobID in list_jobID:

    bashCommand = " cat /cluster/work/jobs/{0}/output_ascii/*.out > /cluster/work/jobs/{0}/fused_{0}.out".format(str(jobID))
    print(bashCommand)
    execute_bash_command_alternative(bashCommand)

    dirname = os.path.dirname(os.path.abspath(__file__))
    print(dirname)
    bashCommand = " cp /cluster/work/jobs/{0}/fused_{0}.out {1}/output_ascii ".format(str(jobID),dirname)
    print(bashCommand)
    execute_bash_command_alternative(bashCommand)

    
