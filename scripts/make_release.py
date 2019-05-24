from __future__ import print_function

import os
import shutil

curdir = os.path.dirname(__file__)

basepath = r'D:\downloads\olex2-gui\util\pyUtil\FragmentDB'

# subprocess.call(args=['C:/tools/Python2.7.15/python.exe', r'-c ../DSR/sql_export.py'], cwd=os.path.split(curdir)[0])


shutil.copy('../DSR/fragment-database.sqlite', '/fragment-database.sqlite')
print('database copied from DSR git to FragmentDB.')

files = ['fragmentdb.pyc', 'fragmentdb_handler.pyc', 'helper_functions.pyc', 'refine_model_tasks.pyc']

for file in files:
  print('copy file:', file)
  shutil.copy(os.path.join(basepath, file), './')
