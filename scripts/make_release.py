import os
import shutil
import subprocess

curdir = os.path.dirname(__file__)

basepath = r'D:\downloads\olex2-gui\util\pyUtil\FragmentDB'

#subprocess.call(args=['C:/tools/Python2.7.15/python.exe', r'-c ../DSR/sql_export.py'], cwd=os.path.split(curdir)[0])

shutil.copy('../DSR/fragment-database.sqlite', '/fragment-database.sqlite')

#shutil.copy(os.path.join(basepath, 'drawstyle.cds'), './')
#shutil.copy(os.path.join(basepath, 'fragmentdb.md'), './')
#shutil.copy(os.path.join(basepath, 'fragmentdb.phil'), './')
#shutil.copy(os.path.join(basepath, 'fragmentdb.htm'), './')
shutil.copy(os.path.join(basepath, 'fragmentdb.pyc'), './')
#shutil.copy(os.path.join(basepath, 'fragmentdb.py'), './')
shutil.copy(os.path.join(basepath, 'fragmentdb_handler.pyc'), './')
#shutil.copy(os.path.join(basepath, 'fragmentdb_handler.py'), './')
shutil.copy(os.path.join(basepath, 'helper_functions.pyc'), './')
#shutil.copy(os.path.join(basepath, 'helper_functions.py'), './')
#shutil.copy(os.path.join(basepath, 'inputfrag.htm'), './')
#shutil.copy(os.path.join(basepath, 'plugins.xld'), './')
#shutil.copy(os.path.join(basepath, 'README.md'), './')
shutil.copy(os.path.join(basepath, 'refine_model_tasks.pyc'), './')
#shutil.copy(os.path.join(basepath, 'refine_model_tasks.py'), './')

