'''
This program reads pw.in file, generates stacking sequence and introduces defects in them.
The output file will be pw_0.in, pw_1.in .... pw_20.in for 0, 0.1 ... 2 fractional displacement respectively.
Run using "python stacking_single_elem.py  <element> <lattice_parameter>" ex. "python stacking_single_elem.py  Al 4.04"
'''

import numpy as np
import math
import subprocess
import sys

element = sys.argv[1]

lattice_parameter = float(sys.argv[2])

# Total layers
layers = 9

# Starting layer for shift
start_l = 6

# Upto this layer the atoms are shifted
end_l = layers



a = b = lattice_parameter/math.sqrt(2)
c = math.sqrt((a*math.sqrt(3)/2)**2-(a*1/(math.sqrt(3)*2))**2)*layers

atom_pos = np.zeros((layers, 3))
layer = 0
while(layer < layers):
    atom_pos[layer, 0] = 0
    atom_pos[layer, 1] = 0
    atom_pos[layer, 2] = layer/9*c

    atom_pos[layer+1, 0] = 0.5*a
    atom_pos[layer+1, 1] = 1/math.sqrt(3)*a/2
    atom_pos[layer+1, 2] = (layer+1)/9*c

    atom_pos[layer+2, 0] = a
    atom_pos[layer+2, 1] = 1/math.sqrt(3)*a
    atom_pos[layer+2, 2] = (layer+2)/9*c

    layer += 3


b_twin = lattice_parameter/math.sqrt(6)

# minimum  shift
shift = 0.1*b_twin

# Minimum shift in x and y direction
shiftX = shift*math.cos(30*math.pi/180)
shiftY = shift*math.sin(30*math.pi/180)


lattice_vectors = np.array([
    [a, 0.0, 0.0],
    [b/2, b*math.sqrt(3)/2, 0.0],
    [0.0, 0.0, c]
])


# returns angle between vectors in radian
def angle_between_vectors(vec1, vec2):
    hat_vec1 = vec1 / np.linalg.norm(vec1)
    hat_vec2 = vec2 / np.linalg.norm(vec2)
    dot_product = np.dot(hat_vec1, hat_vec2)
    return np.arccos(dot_product)


# For converting co-ordinates from cartesian to fractional
def convert_to_fractional(lattice_vectors, atom_pos):

    A, B, C = a, b, np.linalg.norm(lattice_vectors[2])
    alpha = angle_between_vectors(lattice_vectors[1], lattice_vectors[2])
    beta = angle_between_vectors(lattice_vectors[0], lattice_vectors[2])
    gamma = angle_between_vectors(lattice_vectors[0], lattice_vectors[1])

    V = A*B*C*math.sqrt(1-math.cos(alpha)**2-math.cos(beta)**2 -
                        math.cos(gamma)**2+2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))

    transformation_matrix = np.array([
        [1/A, -math.cos(gamma)/(A*math.sin(gamma)),
         B*C*(math.cos(alpha) * math.cos(gamma)-math.cos(beta))/(V*math.sin(gamma))],

        [0, 1/(B*math.sin(gamma)),
         A*C*(math.cos(beta) * math.cos(gamma)-math.cos(alpha))/(V*math.sin(gamma))],

        [0, 0, A*B*math.sin(gamma)/V]
    ])

    atom_pos_fractional = np.zeros((layers, 3))
    for layer in range(0, layers):
        l_v = np.array(atom_pos[layer]).T
        atom_pos_fractional[layer] = np.matmul(transformation_matrix, l_v)

    return atom_pos_fractional


def out_bash_command(bash_command):
    out = subprocess.run(bash_command, capture_output=True, shell=True)
    return out.stdout.decode()


frac_dis_end = 20
for i in range(0, frac_dis_end+1):

    str_lattice_vectors = np.array2string(lattice_vectors, formatter={
        'float_kind': lambda lattice_vectors: "%.8f" % lattice_vectors}).replace("[", " ").replace("]", "")

    atom_pos_fractional = convert_to_fractional(lattice_vectors, atom_pos)
    str_atom_pos = np.array2string(atom_pos_fractional, formatter={
        'float_kind': lambda atom_pos_fractional: "%.8f" % atom_pos_fractional}).replace("[[", " [").replace("[", element+"\t").replace("]]", "]").replace("]", " 0 0 1")

    # Extracting from pw.in
    bash_command = "sed '/CELL_PARAMETERS/q' pw.in"
    s = out_bash_command(bash_command)
    s += (str_lattice_vectors+"\n\n")
    s += ("ATOMIC_POSITIONS {crystal}\n")
    s += (str_atom_pos+"\n\n")
    bash_command = "grep -A 1 'K_POINTS' pw.in"
    s += out_bash_command(bash_command)

    # Writing string to file
    with open(f'pw_{i}.in', 'w') as f:
        f.write(s)

    # Condition for start layer shift
    if(i > 0 and i % 10 == 0):
        start_l += 1

    # Shifting cell
    lattice_vectors[2, 0] += shiftX
    lattice_vectors[2, 1] += shiftY

    # Shifting atoms
    atom_pos[start_l-1:end_l, 0] += shiftX
    atom_pos[start_l-1:end_l, 1] += shiftY
