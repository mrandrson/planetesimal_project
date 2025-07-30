import os

def write_script(v0, N_exp, t, saveN, script_name):
    lines = [
        '#!/bin/bash\n',
        '#SBATCH -n 1\n',
        '#SBATCH -o test.out\n',
        '#SBATCH -e test.err\n',
        '#SBATCH -p phys-k\n',
        f'python testparticle.py {v0} {N_exp} {t} {saveN}\n'
    ]
    with open(script_name, 'w') as arc:
        arc.writelines(lines)

def create_output_directories(base_dirs, sub_dirs):
      for base_dir in base_dirs:
        for sub_dir in sub_dirs:
              full_path = os.path.join(base_dir, sub_dir)
              os.makedirs(full_path, exist_ok=True)
base_dirs = []
v0 = {1: 5, 2: 10, 3: 15}
for i in [1, 2, 3]:
    base_dirs.append(f'./data{i}')
#base_dirs = ['./data1', './data2', './data3']
sub_dirs = ['testacceleration', 'testvelocity', 'testposition', 'testtime', 'testenergy', 'testke', 'testpe']

create_output_directories(base_dirs, sub_dirs)

write_script(5, 31186, 350000, './data1', 'test1.sh')
os.system('sbatch test1.sh')

write_script(10, 62358, 350000, './data2', 'test2.sh')
os.system('sbatch test2.sh')

write_script(15, 93533, 350000, './data3', 'test3.sh')
os.system('sbatch test3.sh')

