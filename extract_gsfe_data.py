'''
Extracts data from calculation to obtain gsfe and store it in data.csv
'''

import numpy as np
import subprocess
import re


# Get area
bash_command = "grep -A 3 'CELL_PARAMETERS' pw_0.in | tail -3"
out = subprocess.run(bash_command, capture_output=True, shell=True)

l_param = np.array(out.stdout.decode().split(), dtype='f')
vec1 = l_param[0:3]
vec2 = l_param[3:6]
area = np.linalg.norm(np.cross(vec1, vec2))/(1E10)**2

f = open("data.csv", "w")
for i in range(0, 21):

    bash_command = "grep ! pw_"+str(i)+".out | tail -1"
    out = subprocess.run(bash_command, capture_output=True, shell=True)
    if(out.stdout.decode()):
        energy = float(re.findall(
            "-\d+\.\d+", out.stdout.decode())[0])*2.179872E-15

        if(i == 0):
            energy_0 = energy

        sfe = (energy-energy_0)/area
        f.write(f'{i/10:.2f},{sfe}\n')

f.close()
