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
               'scripts/setup_annotation_dbs.py',
               'zol/orthologs/findOrthologs.py'],
      zip_safe=False)

os.system("g++ -o zol/orthologs/runRBH zol/orthologs/runRBH.cpp")
os.system("g++ -o zol/orthologs/splitDiamondResults zol/orthologs/splitDiamondResults.cpp")
