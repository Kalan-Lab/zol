from setuptools import setup
import os

setup(name='zol',
      version='1.00',
      description='',
      url='http://github.com/Kalan-Lab/zol/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['zol'],
      scripts=['bin/zol',
               'bin/fai',
               'scripts/readifyAdditionalGenomes.py',
               'scripts/runProdigalAndMakeProperGenbank.py',
               'scripts/setup_annotation_dbs.py',
               'zol/orthologs/findOrthologs.py',
               'external_tools/Treemmer_v0.3.py'],
      zip_safe=False)

# compile RBH/InParanoid-esque programs written in C++
os.system("g++ -o zol/orthologs/runRBH zol/orthologs/runRBH.cpp")
os.system("g++ -o zol/orthologs/splitDiamondResults zol/orthologs/splitDiamondResults.cpp")
# set up STAG by Emms and Kelly
if not os.path.isdir('STAG_1.0.0/'):
      os.system("axel https://github.com/davidemms/STAG/releases/download/1.0.0/STAG_1.0.0.tar.gz")
      os.system("tar -zxvf STAG_1.0.0.tar.gz")
      os.system("rm -f STAG_1.0.0.tar.gz")