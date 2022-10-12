#!/bin/bash

QE_BIN=/home/albert20/softwares/qe-6.7/bin
POT_DIR=/home/albert20/pseudopotentials/pslibrary
calc_command="mpirun -np 96 $QE_BIN/pw.x -nk 4"
pot="Cu.pbe-dn-kjpaw_psl.1.0.0.UPF"
element="Cu"
atomic_mass=63.546

pressure=-100
while [ $pressure -le 100 ]
do
    mkdir "pressure_"$pressure
    cd "pressure_"$pressure
    mkdir "cell_relaxation" "gsfe"
    cd "cell_relaxation"
     
    cat > pw.in <<EOF
&control
   calculation = 'vc-relax' ,
   prefix = 'pw' ,
   outdir = './temp/' ,
   pseudo_dir = '$POT_DIR' ,
   disk_io = 'none' ,
   etot_conv_thr = 1.0D-4 ,
   forc_conv_thr = 1.0d-3 ,
/
&system
    ibrav = 0
    ecutwfc = 100 ,
    ecutrho = 400 ,
    occupations = 'smearing' ,
    smearing = 'gaussian' ,
    degauss = 0.01 ,
    nat = 9 ,
    ntyp = 1
/
&electrons
   conv_thr = 1.0d-6,
   diagonalization = 'david' ,
   mixing_mode = 'plain' ,
   startingpot = 'atomic' ,
   startingwfc = 'atomic+random' ,
   mixing_beta = 0.7 ,
/
&ions
   ion_dynamics = 'bfgs' ,
/
&cell
   cell_dynamics = 'bfgs' ,
   cell_factor = 2.0 ,
   press = $pressure ,
   press_conv_thr = 1.0D-2 ,
/

ATOMIC_SPECIES
$element  $atomic_mass  $pot
 
CELL_PARAMETERS angstrom
  2.56538340 0.00000000 0.00000000
  1.28269170 2.22168720 0.00000000
  0.00000000 0.00000000 18.85164099

ATOMIC_POSITIONS {crystal}
 $element	0.00000000 0.00000000 0.00000000
 $element	0.33333333 0.33333333 0.11111111
 $element	0.66666667 0.66666667 0.22222222
 $element	-0.00000000 -0.00000000 0.33333333
 $element	0.33333333 0.33333333 0.44444444
 $element	0.66666667 0.66666667 0.55555556
 $element	-0.00000000 -0.00000000 0.66666667
 $element	0.33333333 0.33333333 0.77777778
 $element	0.66666667 0.66666667 0.88888889

K_POINTS automatic
   27 27 3 0 0 0

EOF
   # Perform relaxation
   $calc_command < pw.in > pw.out

   # Extract volume
   volume=$(grep -oE "[0-9]+.[0-9]+\s+Ang\^3" pw.out | tail -1 | grep -oP "[0-9]+.[0-9]+")

   lattice_parameter=$(echo "e( l($volume/9*4)/3 )" | bc -l)
   echo $lattice_parameter

   cp pw.in ../gsfe
   cd ../gsfe
   cp ../../*.py $PWD

   # Generated defect structures
   python3 stacking_single_elem.py $element $lattice_parameter
   
   # Relaxation for defect structure
   for((i=0;i<=20;i++))
   do
       $calc_command < pw_$i.in > pw_$i.out
   done

   # Extracts energy to obtain GSFE data
   python3 extract_gsfe_data.py 
   
   cd ../../
   pressure=$(( $pressure + 50 ))
done
