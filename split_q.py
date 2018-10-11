"""
This scripts prepares the input files for a dvscf/gkkp calculation with QE.
The calculation is managed by 'manual' parallelisation over q-points (both regular grid or randomly chosen).
The folder structures and the underlying scf/nscf calculations must be already present.
"""
from __future__ import print_function
from yambopy     import *
from qepy        import *
import numpy as np
import os
import argparse
import subprocess

def clean_shit(folder,qpts,pwd,r_or_mp):
	""" Removes files that are not needed in between calculations
	"""
	os.system('rm -rf %s/q%d/*wfc*'%(folder,qpts))
	if r_or_mp == 'mp':
		if qpts==0: os.system('cd %s/q%d/_ph0 ; mv *dvscf* ../ ; rm -rf * ; mv ../*dvscf* ./ ; cd %s'%(folder,qpts,pwd))
		else:	
			os.system('rm -rf %s/q%d/_ph0/*.phsave'%(folder,qpts))
			os.system('cd %s/q%d/_ph0/*.q_%d/ ; mv *dvscf* ../ ; rm -rf * ; mv ../*dvscf* ./ ; cd %s'%(folder,qpts,qpts+1,pwd))
	if r_or_mp == 'r': os.system('cd %s/q%d/_ph0 ; mv *dvscf* ../ ; rm -rf * ; mv ../*dvscf* ./ ; cd %s'%(folder,qpts,pwd))
#
def check_calc(folder,qpts,calc_type):
	""" Check if qpts-th job has finished properly. If it failed, the script will proceed to the next job without exiting.
	"""
	status = True
	check = subprocess.check_output("grep JOB %s/q%d/slurm*"%(folder,qpts), shell=True) # !! Works for SLURM standard output !!
        check = check.strip().split()[-1]
        if check != 'DONE.':
        	print('[WARNING] %s calculation for qpoint %s/q%d FAILED.'%(calc_type,work_dir,iq))
                status=False
	return status
#
def locus_2dhex(vec): #A 2-d hexagon oriented in the "ibrav=4" way is the locus of points that have the same distance according to the norm defined here
	return max(abs(vec[1]),(abs(vec[1])+sqrt(3)*abs(vec[0]))*0.5)
def generate_random_pts(nqpts,scf_inp,work_dir):
	"""Generate nqpts random points inside the Wigner-Seitz cell
	   IMPLEMENTATION IS ONLY FOR HEXAGONAL CELLS ORIENTED IN THE "ibrav=4" WAY!!!!!
	"""
	out_nm = '%s/random_xy_qpts.dat'%work_dir
	if os.path.isfile('%s'%out_nm): rand3d = np.loadtxt('%s'%out_nm)
	else:
		qe = PwIn(scf_inp)
		alat = float(qe.system['celldm(1)'])
		lat_vec = np.array(qe.cell_parameters)/alat #Lattice vectors in alat units
		rlat_vec= rec_lat(lat_vec) #Reciprocal lattice vectors in 2pi/alat units
		norm1, norm2 = [np.linalg.norm(rlat_vec[0]), np.linalg.norm(rlat_vec[1])]
		if abs(norm1-norm2)>1e-5: #This is a check that rec lat vecs actually have the same norm (hexagonal cell)
    			print("ERROR in generate_random_pts: the two lattice vectors have different lengths")
    			exit()
		R_hex = norm1/2. #Radius of the hexagonal "ball"
		i=0
		randxy = []
		print("Generating random q-points...")
		while i < nqpts:
    			r_qpt = np.random.rand(1,2)[0]*2.*norm1-norm1 #Random 2d point in an area that contains the Wigner-Seitz cell
    			if locus_2dhex(r_qpt)<=R_hex: #If the point is inside the WS cell we keep it, otherwise we generate another
        			randxy.append(r_qpt)
        			i+=1
		randxy = np.array(randxy)
		randxy[0]=np.zeros(2) #First point must always be Gamma!
		rand3d = np.c_[randxy, np.zeros(nqpts), np.ones(nqpts)] #Adding z-components and fake weight column
        	if len(rand3d)!=nqpts:
                	print("ERROR in generate_random_pts: number of selected random points < nqpts")
                	exit()
		f = open('%s'%out_nm,'w')
                for q in rand3d: f.write('%12.8f '*4%tuple(q) + '\n')
                f.close()
		print("...Done.")		
	return rand3d
#
parser = argparse.ArgumentParser(description='Script to run electron-phonon calculation with manual qpts parallelisation')
parser.add_argument('-r' ,'--random',      action="store_true", help='Random qpoint grid')
parser.add_argument('-mp' ,'--regular',  action="store_true", help='Monkhorst-Pack qpoint grid')
parser.add_argument('-inp','--inputs', action="store_true", help='dvscf and gkkp input files')
parser.add_argument('-ph','--run_phonons',action="store_true", help='Run the calculations: either dvscf or gkkp')
args = parser.parse_args()

##################### USER INPUT AREA ##############################
pwd=subprocess.check_output("pwd", shell=True).strip()
prefix = 'hbn'
nbnd   = '20'
scf_dir   = os.path.expandvars('%s/scf'%pwd)
scf_nm    = '%s.scf'%prefix
nscf_dir  = os.path.expandvars('%s/nscf/%s'%(pwd,nbnd))
dvscf_inp = 'dvscf.in'
gkkp_inp  = 'gkkp.in'
nq1, nq2, nq3 = [12,12,1] #Regular grid
nq  = 576                 #Random grid: total number of qpoints / Regular grid: number of qpoints in the IBZ
qfirst,qlast = [0,nq]     #Which qpoints will be run (default: all of them)
#work_dir  = os.path.expandvars('%s/%dx%dx%d'%(pwd,nq1,nq2,nq3))
work_dir  = os.path.expandvars('%s/%dr'%(pwd,nq)) #Directory where the q-dependent jobs will be run
""" Temporary note, nq in IBZ for ibrav=4
12x12x1 --> nq=19
12x12x3 --> nq=38
12x12x6 --> nq=76
12x12x9 --> nq=95
12x12x12 --> nq=133
"""
############# END USER INPUT AREA #################################

#check if submission script is present in the same folder as this script
if not os.path.isfile('%s/run_ph.sh'%pwd):
	print('run_ph.sh script absent from main folder')
        exit()
os.system('cp run_ph.sh %s/'%work_dir)

#check if the scf cycle is present
if not os.path.isdir('%s/%s.save' % (scf_dir,prefix)):
    print('scf calculation not found!')
    exit()

#check if the nscf cycle is present
if not os.path.isdir('%s/%s.save' % (nscf_dir,prefix)):
    print('nscf calculation not found!')
    exit()

#check if parser instruction was given
if not args.random and not args.regular:
        print('Select a qpoint generation mode')
        exit()

# Input files for the electron-phonon calculation
phin = PhIn()
phin['prefix']          = "'%s'"%prefix
phin['tr2_ph']          = 1.0e-14
phin['fildyn']          = "'%s.dyn'"%prefix
phin['qplot']           = '.true.'
phin['fildvscf']        = "'dvscf'"
phin['iverbosity']      = 1
phin['ldisp']           = '.false.'
phin['trans']           = '.true.'
phin['electron_phonon'] = "'dvscf'" 
if args.random:
	q_list = generate_random_pts(nq,'%s/%s'%(scf_dir,scf_nm),work_dir)
	whichclean = 'r'
	phin['qplot']           = '.true.'
	phin['ldisp']           = '.false.'
	phin.qpoints 		= q_list
if args.regular:
	whichclean = 'mp'
	phin['qplot']           = '.false.'
	phin['ldisp']           = '.true.'
	phin['nq1']             = nq1
	phin['nq2']             = nq2
	phin['nq3'] 		= nq3
# Writing the full input files
phin.write('%s/%s'%(work_dir,dvscf_inp))
phin['trans']           = '.false.'
phin['electron_phonon'] = "'yambo'"
phin.write('%s/%s'%(work_dir,gkkp_inp))

if args.inputs:
	#Preparing folder structure and q-dependent input files for calculations
	for iq in range(nq):
                phin['start_q'] = iq+1
                phin['last_q']  = iq+1
                if not os.path.isdir('%s/q%d'%(work_dir,iq)):
                        os.system('mkdir -p %s/q%d'%(work_dir,iq))
                        phin['trans']           = '.true.'
                        phin['electron_phonon'] = "'dvscf'"
                        phin.write('%s/q%d/%s'%(work_dir,iq,dvscf_inp)) #dvscf input
                        phin['trans']           = '.false.'
                        phin['electron_phonon'] = "'yambo'"
                        phin.write('%s/q%d/%s'%(work_dir,iq,gkkp_inp)) #gkkp input

if args.run_phonons:
	# Start of manually parallelised calculations
	for iq in range(qfirst,qlast):
		# Copy elphon files to generate yambo DB
		if os.path.isdir('%s/q%d/elph_dir'%(work_dir,iq)):
			if check_calc(work_dir,iq,'gkkp') == False: continue
			clean_shit(work_dir,iq,pwd,r_or_mp=whichclean)
			os.system('mv %s/q%d/elph_dir/s.dbph_%06d %s/'%(work_dir,iq,iq+1,work_dir))
			continue
		# Run dvscf calculations
		if not os.path.isdir('%s/q%d/_ph0'%(work_dir,iq)):
			os.system('ln -s %s/%s.save %s/q%d'%(scf_dir,prefix,work_dir,iq)) #scf save is linked
			os.system('cd %s ; sh run_ph.sh %d %s %s ; cd ../'%(work_dir,iq,'dvscf',dvscf_inp))
		# Run gkkp calculations
		else:
			if check_calc(work_dir,iq,'dvscf') == False: continue
			clean_shit(work_dir,iq,pwd,r_or_mp=whichclean)
			os.system('rm %s/q%d/%s.save'%(work_dir,iq,prefix))
			os.system('ln -s %s/%s.save %s/q%d'%(nscf_dir,prefix,work_dir,iq)) #nscf save is linked
			os.system('cd %s ; sh run_ph.sh %d %s %s ; cd ../'%(work_dir,iq,'gkkp',gkkp_inp))

