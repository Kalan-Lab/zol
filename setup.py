from setuptools import setup
import os

setup(name='zol',
      version='1.2.1',
      description='',
      url='http://github.com/Kalan-Lab/zol/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['zol'],
      scripts=['bin/zol',
               'bin/fai',
               'bin/prepTG',
               'scripts/runProdigalAndMakeProperGenbank.py',
               'scripts/listAllGenomesInDirectory.py',
               'scripts/setup_annotation_dbs.py',
               'scripts/processNCBIGenBank.py',
               'scripts/extractBiG-SCAPEclusters.py',
               'zol/orthologs/findOrthologs.py',
               'scripts/generateSyntenicVisual.py',
               'scripts/convertMiniprotGffToGbkAndProt.py',
               'scripts/expandDereplicatedAlignment.py'],
      zip_safe=False)

# compile RBH/InParanoid-esque programs written in C++
os.system("g++ -std=c++11 -o zol/orthologs/runRBH zol/orthologs/runRBH.cpp")
os.system("g++ -std=c++11 -o zol/orthologs/splitDiamondResults zol/orthologs/splitDiamondResults.cpp")
os.system("g++ -std=c++11 -o zol/splitDiamondResultsForFai zol/splitDiamondResultsForFai.cpp")
