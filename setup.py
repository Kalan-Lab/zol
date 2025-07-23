#!/usr/bin/env python
from setuptools import setup
import os
import traceback

if __name__ == "__main__":
    setup()
    conda_dir = os.path.abspath(os.environ['CONDA_PREFIX']) + '/'
    try:
        zol_exec_path = conda_dir + 'bin'
        zol_data_path = conda_dir + 'share/zol/db'
        conda_act_dir = conda_dir + 'etc/conda/activate.d/'
        conda_dea_dir = conda_dir + 'etc/conda/deactivate.d/'
    
        os.system('ZOL_EXEC_PATH=' + zol_exec_path)
        os.system('ZOL_DATA_PATH=' + zol_data_path)
        os.environ['ZOL_DATA_PATH'] = zol_data_path
        os.environ['ZOL_EXEC_PATH'] = zol_exec_path

        os.system("g++ -std=c++14 -o " + zol_exec_path + "/runRBH src/zol/orthologs/runRBH.cpp")
        os.system("g++ -std=c++14 -o " + zol_exec_path + "/splitDiamondResults src/zol/orthologs/splitDiamondResults.cpp")
        os.system("g++ -std=c++14 -o " + zol_exec_path + "/splitDiamondResultsForFai src/zol/splitDiamondResultsForFai.cpp")

        assert(os.path.isfile(zol_exec_path + '/runRBH') and os.path.isfile(zol_exec_path + '/splitDiamondResults') and os.path.isfile(zol_exec_path + '/splitDiamondResultsForFai'))
            
        os.system('mkdir -p ' + zol_data_path)
        os.system('mkdir -p ' + conda_act_dir + ' ' + conda_dea_dir) 
        os.system("echo 'Default conda space for downloading annotation databases.\n' > " + zol_data_path + "/README.txt")
        
        act_outf = open(conda_act_dir + 'zol.sh', 'w')
        act_outf.write('#!/usr/bin/env bash\n\n')
        act_outf.write('export ZOL_DATA_PATH=' + zol_data_path + '\n')
        act_outf.write('export ZOL_EXEC_PATH=' + zol_exec_path + '\n')
        act_outf.close()

        dea_outf = open(conda_dea_dir + 'zol.sh', 'w')
        dea_outf.write('#!/usr/bin/env bash\n\n')
        dea_outf.write('unset ZOL_DATA_PATH\n')
        dea_outf.write('unset ZOL_EXEC_PATH\n')
        dea_outf.close()
    except:
        outf = open('installation_error.txt', 'w')
        outf.write('conda dir is: ' + conda_dir + '\n')
        outf.write('C++ compilation or setup of executables failed or had issues setting up space in conda environment for databases!\n')
        outf.write('traceback.format_exc()\n')
        outf.close()
