import os
import filecmp
from subprocess import Popen, PIPE

# CONFIG
DAS = True         # change whether running on DAS or not
programs = ["matmul_basic", "matmul_opt", "matmul_omp"]
reference = programs[0]
matrices = ["arc130.mtx", "bcspwr06.mtx", "bcsstk13.mtx", "G6.mtx", "tols1090.mtx"]
input_dir = "input/"
output_dir = "output/"

# LOGIC
template_command = ["prun", "-np", "1"] if DAS else []

try: 
    # if output dir does not exist, create it..
    os.mkdir(output_dir)
    print(f"Created directory '{output_dir}'")
except FileExistsError: 
    pass

for program in programs:
    print("Running:", program)
    for matrix in matrices:
        print("\tFile:", matrix)
        in_matrix_path = input_dir + matrix
        if program == reference:
            out_matrix_path = output_dir + "ref_out_" + matrix
        else:
            out_matrix_path = output_dir + "tmp.txt"

        command = template_command + ["./" + program, in_matrix_path, in_matrix_path, out_matrix_path]
        process = Popen(command, stdout=PIPE)
        process.wait()
        output_lines = process.stdout.readlines()

        last_line_idx = len(output_lines) - 1
        gflops = float(output_lines[last_line_idx - 1].split()[-1])
        exec_time = float(output_lines[last_line_idx].split()[-1])

        print("\t\tGFLOP/s: ", gflops)
        print("\t\tSeconds: ", exec_time)
        if program != reference:
            if filecmp.cmp(out_matrix_path, output_dir + "ref_out_" + matrix, shallow=False):
                print("\t\tOutput: Correct")
            else:
                print("\t\tOutput: False")