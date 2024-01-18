"""
Create all files required to run TAT simulation(s) with set seeds and other geometric arguments. 
simulation & clustering run separately so that multithreading can be used on simulation
"""
import os
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--spacing", nargs="+",  type=float ,default = [1.0], 
                    help="spacing in micrometer for the dna voxels")

parser.add_argument("--n_div_R", type=int,
                    help="n of divisions in R", default = 80)
parser.add_argument("--n_div_Z", type=int,
                    help="n of divisions in R", default = 40)
parser.add_argument("--n_div_theta", type=int,
                    help="n of divisions in R", default = 100)

parser.add_argument("--startR", type=float, default = 10.5)
parser.add_argument("--vessel_halflength", type = float, default = 20.0)
parser.add_argument("--shell_halflength", type = float, default = 1.5)
parser.add_argument("--gun_length", type = float, default = 40.0)
parser.add_argument("--gpsHalfZ", type = float, default = 10.0)
parser.add_argument("--n", type=int, default = 100)
parser.add_argument("--slurm", type=bool, required=False, default = False)
parser.add_argument("--seed", type=int, required=False)
parser.add_argument("--nthreads", type=int, default=4)
parser.add_argument("--mem", type=int, default=20)
parser.add_argument("--continuous", type=bool, default=False)
parser.add_argument("--folder", type=str, required = False)

args = parser.parse_args()



# parameters to set
if args.nthreads:
    numThread = args.nthreads
else:
    numThread = 4

if args.slurm:
    slurm = args.slurm
else:
    slurm = False
if args.seed:
    seed = args.seed
else:
    seed = random.randint(1, 10000)
if args.mem:
    mem = args.mem
else:
    mem = 20
simulation_parent = os.path.abspath(os.curdir)
if args.continuous:
    continuous = args.continuous
else:
    continuous = False

if continuous:
    histoneFile = os.path.join(
        simulation_parent, "geometryFiles", "histonePos_output.bin"
    )
    sugarFile = os.path.join(
        simulation_parent, "geometryFiles", "sugarPos_output.bin"
    )

else:
    histoneFile = os.path.join(
        simulation_parent, "geometryFiles", "histonePositions_4x4_300nm.bin"
    )
    sugarFile = os.path.join(
        simulation_parent, "geometryFiles", "sugarPos_4x4_300nm.bin"
    )

time_decay = "0-10:00:00"
time_DNA = "2-24:00:00"
time_clustering = "0-6:00:00"


targetSize = 300  # nanometers
gpsRadius = 10 # micrometers

if args.gpsHalfZ:
    gpsHalfZ = args.gpsHalfZ  # micrometers

if args.n_div_R:
    n_div_R = args.n_div_R

if args.n_div_theta:
    n_div_theta = args.n_div_theta
    
if args.n_div_Z:
    n_div_Z = args.n_div_Z

if args.spacing:
    spacing = args.spacing

if args.startR:
    startR = args.startR

if args.vessel_halflength:
    vessel_halflength = args.vessel_halflength
if args.shell_halflength:
    shell_halflength = args.shell_halflength
if args.gun_length:
    gun_length = args.gun_length

if args.n:
    numIons = args.n

printProgress = int(numIons/10)

n_string = f"{int(numIons/1000)}k" if numIons/1000 >=1 else f"{numIons}"

if args.folder:
    folder = args.folder
else:
    if continuous:
        folder = f"test_At{n_string}_shell_ps_{startR}um_{n_div_R}R_continuous"
    else:
        folder = f"test_At{n_string}_shell_ps_{startR}um_{n_div_R}R"
# print(f"Creating scripts in folder: {folder} with seed {seed}.")
####

decay_sim_folder = os.path.join(simulation_parent, "decay_simulation")
dna_sim_folder = os.path.join(simulation_parent, "simulation")
clustering_folder = os.path.join(simulation_parent, "Clustering")

makerundir = lambda d: os.path.join(d, "build")
projectcode="phys004501"

#### create folder for current test
test_dir = os.path.join(simulation_parent, f"output/{folder}")

os.makedirs(test_dir, exist_ok=True)

#### create one file per spacing
for s in spacing:
    s_string = f"{int(s)}"

    #  write mac file

    filename_mac = os.path.join(
        test_dir, f"input_Atdecay_{n_string}_length_{gun_length}_{seed}.in"
    )  # to be run from folder created
    with open(filename_mac, "w") as f:
        f.write("/run/verbose 2\n")
        f.write("/control/verbose 2\n")
        f.write("\n")
        f.write(f"/det/set_ndiv_R {n_div_R}\n")
        # f.write(f"/det/set_ndiv_theta {n_div_theta}\n")
        f.write(f"/det/set_ndiv_Z {n_div_Z}\n")
        f.write("\n")
        f.write(f"/det/set_spacing {s} um\n")
        f.write(f"/det/set_startR {startR} um\n")
        f.write(f"/det/set_vessel_halflength {vessel_halflength} um\n")
        f.write(f"/det/set_shell_halflength {shell_halflength} um\n")
        f.write(f"/shell/gun/gun_length {gun_length} um\n")
        f.write("\n")
        f.write("/run/initialize\n")
        f.write("\n")
    
        f.write("/gun/particle ion\n")
        f.write("/gun/ion 85 211\n")
        f.write(f"/gun/energy 0 MeV\n")
        f.write("\n")
        f.write(f"/run/printProgress {printProgress}\n")
        f.write(f"/run/beamOn {numIons}\n")
    
    filename_tat_mac = os.path.join(test_dir, f"tat_{numThread}thread.in")
    with open(filename_tat_mac, "w") as f:
        f.write("/run/verbose 2\n")
        f.write("/control/verbose 2\n")
        f.write("\n")       
        f.write(f"/run/numberOfThreads {numThread}\n")
        f.write("/run/initialize \n")
        f.write("\n")
        f.write("/gun/particle gamma\n")
        f.write("/gun/energy 0\n")
        f.write("/gun/momentumAmp 0 \n")
        f.write("\n")
        f.write("/scheduler/endTime 5 nanosecond\n")
        f.write("\n")
        f.write("/run/printProgress 1\n")
        f.write("\n")
        f.write("/run/beamOn 0\n")

    filename_decay = os.path.join(
    test_dir, f"run_Atdecay_{n_string}_spacing_{s_string}um_length_{gun_length}_{seed}.sh")  # to be run from folder created


    # if slurm:
    #     print(f"sbatch run_Atdecay_{n_string}_spacing_{s_string}um_{seed}.sh")
    with open(filename_decay, "w") as f:
        f.write("#!/bin/bash --login\n")
        if slurm:
            f.write(f"#SBATCH --account={projectcode}\n")
            f.write(f"#SBATCH --job-name=Atdecay_spacing_{s_string}um_{seed}\n")
            f.write(f"#SBATCH --output={test_dir}/Atdecay_spacing_{s_string}um_{seed}.out.%J\n")
            f.write(f"#SBATCH --error={test_dir}/Atdecay_spacing_{s_string}um_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={time_decay}\n")
            f.write("#SBATCH --nodes=1\n")
            # f.write("#SBATCH -p short\n")
            f.write(f"#SBATCH --ntasks-per-node={numThread}\n")
            f.write(f"#SBATCH --mem 2GB\n")
            f.write("\n")
            f.write("module load apps/geant/4.11.1\n")
            f.write("module load apps/root/6.26.00\n")
		
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate dart\n")

        else:
            f.write("conda activate rootpy\n")
            f.write("source /opt/geant4-v11.1.0-install/bin/geant4.sh\n")
        # run decay simulation
        f.write(
            f"time {makerundir(decay_sim_folder)}/decaySim -mac {filename_mac} -out {test_dir}/out_Atdecay_{n_string}_spacing_{s_string}um_{seed} -seed {seed}\n"
        )
        

        filename_DNA = os.path.join(
        test_dir, f"run_AtDNA_{n_string}_spacing_{s_string}um_length_{gun_length}_{seed}.sh"
    )  # to be run from folder created
    # if slurm:
    #     print(f"sbatch run_AtDNA_{n_string}_spacing_{s_string}um_{seed}.sh")
    with open(filename_DNA, "w") as f:
        f.write("#!/bin/bash --login\n")
        if slurm:
            f.write(f"#SBATCH --account={projectcode}\n")
            f.write(f"#SBATCH --job-name=AtDNA_spacing_{s_string}um_{seed}\n")
            f.write(f"#SBATCH --output={test_dir}/AtDNA_spacing_{s_string}um_{seed}.out.%J\n")
            f.write(f"#SBATCH --error={test_dir}/AtDNA_spacing_{s_string}um_{seed}.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={time_DNA}\n")
            f.write("#SBATCH --nodes=1\n")
            # f.write("#SBATCH -p short\n")
            f.write(f"#SBATCH --ntasks-per-node={numThread}\n")
            f.write(f"#SBATCH --mem {mem}GB\n")
            f.write("\n")
            f.write("module load apps/geant/4.11.1\n")
            f.write("module load apps/root/6.26.00\n")
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate dart\n")

        else:
            f.write("conda activate rootpy\n")
            f.write("source /opt/geant4-v11.1.0-install/bin/geant4.sh\n")
        # run DNA simulation
        f.write(
            f"time {makerundir(dna_sim_folder)}/tat -mac {filename_tat_mac} -PS {test_dir}/out_Atdecay_{n_string}_spacing_{s_string}um_{seed}.bin -out {test_dir}/out_AtDNA_{n_string}_spacing_{s_string}um_{seed}.root -sugar {sugarFile} -histone {histoneFile} -seed {seed} \n"
        )

        

        # combinedFilename = os.path.join(test_dir, f"alpha_{angle}deg_{seed}_18.5.root")

        # output = "hadd " + combinedFilename

        # for i in range(0, numThread):
        #     toCombine = f" alpha_{angle}deg_{seed}_18.5_t{i}.root"
        #     output += toCombine

        # f.write(output + "\n")

        # f.write(f"if test -f {combinedFilename}; \n")
        # f.write("\t then \n")
        # fileToDelete = os.path.join(test_dir, f"alpha_{angle}deg_{seed}_18.5_t*.root")
        # f.write("rm {}; fi\n".format(fileToDelete))

    filename_clustering = os.path.join(
        test_dir, f"run_clustering_At_{n_string}_spacing_{s_string}um_length_{gun_length}_{seed}.sh"
    )  # to be run from folder created

    with open(filename_clustering, "w") as f:
        f.write("#!/bin/bash --login\n")

        if slurm:
            f.write(f"#SBATCH --job-name=clustering_At_{n_string}_spacing_{s_string}um\n")
            f.write(f"#SBATCH --account={projectcode}\n")
            f.write(f"#SBATCH --output={test_dir}/clustering_At_{n_string}_spacing_{s_string}um.out.%J\n")
            f.write(f"#SBATCH --error={test_dir}/clustering_At_{n_string}_spacing_{s_string}um.err.%J\n")
            # maximum job time in D-HH:MM
            f.write(f"#SBATCH --time={time_clustering}\n")
            f.write("#SBATCH --nodes=1\n")
            # f.write("#SBATCH -p short\n")
            f.write("#SBATCH --ntasks-per-node=1\n")
            f.write(f"#SBATCH --mem {mem}GB\n")
            f.write("\n")
            f.write("module load apps/root/6.26.00\n")

        f.write("\n")

        # clustering
        clustering_dir = os.path.join(simulation_parent, "Clustering")
        clustering_outdir = os.path.join(test_dir, "clustering_out")
        if not os.path.exists(clustering_outdir):
            os.makedirs(clustering_outdir)

        if slurm:
            f.write("source /user/home/yw18581/.bash_profile\n")
            f.write("source activate clustering\n")
        else:
            f.write("conda activate clustering\n")
        f.write(
            f"python {clustering_dir}/run_up_part.py --filename {test_dir}/out_AtDNA_{n_string}_spacing_{s_string}um_{seed}.root --output {clustering_outdir}/out_AtDNA_{n_string}_spacing_{s_string}um_{seed}_part.csv --sugar {sugarFile} --ndiv_R {n_div_R}  --continuous {continuous}\n"
        )
    filename_runscript = os.path.join(test_dir, f"run_script_At_{n_string}_spacing_{s_string}um_length_{gun_length}_{seed}.sh")
    with open(filename_runscript, "w") as f:
        f.write("#!/bin/bash --login\n")

        if slurm:
            f.write(f"#SBATCH --job-name=run_At{n_string}_sp_{s_string}um_{seed}\n")
            f.write(f"#SBATCH --account={projectcode}\n")
            f.write(f"#SBATCH --output={test_dir}/run_At_{n_string}_spacing_{s_string}um_{seed}.out.%J\n")
            f.write(f"#SBATCH --error={test_dir}/run_At_{n_string}_spacing_{s_string}um_{seed}.err.%J\n")
            f.write(f"#SBATCH --time=0-10:00\n")
            f.write("#SBATCH --nodes=1\n")
            # f.write("#SBATCH -p short\n")
            f.write("#SBATCH --ntasks-per-node=1\n")
            f.write(f"#SBATCH --mem 1GB\n")
            f.write("\n")
            f.write("source /user/home/yw18581/.bash_profile\n")
            # f.write("module load apps/root/6.26.00\n")

        f.write("\n")

        
        if not slurm:
            f.write(f"source {filename_decay}\n")
            f.write(f"source {filename_DNA}\n")
            f.write(f"source {filename_clustering}\n")
        if  slurm:
            f.write(f"job_decay=$(sbatch --parsable {filename_decay})\n")
            f.write(f"job_DNA=$(sbatch --parsable --dependency=afterok:$job_decay {filename_DNA})\n")
            f.write(f"job_clustering=$(sbatch --dependency=afterany:$job_DNA {filename_clustering})\n")
        
if not slurm:     
    # print(f"{seed}", end=' ') 
    print(f"tmux new-session -d -s tat_{n_string}_{seed}\n 'source {filename_runscript}' ")
else:
    # print(f"sbatch {filename_runscript}")
    print(f"{seed}", end=' ')
