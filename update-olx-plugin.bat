@ECHO OFF

rem This file is for development only.

echo %FDBDIR%
rem set FDBDIR="C:\Program Files\Olex2-1.2\util\pyUtil\PluginLib\plugin-Fragment_DB"
set FDBDIR2="D:\downloads\olex2-gui\util\pyUtil\FragmentDB"
set GIT="D:\GitHub\FragmentDB"


rem xcopy /Y %GIT%\Fragment_DB.htm %FDBDIR%
rem xcopy /Y %GIT%\helptext.htm %FDBDIR%
rem xcopy /Y %GIT%\inputfrag.htm %FDBDIR%
rem xcopy /Y %GIT%\Fragment_DB.phil %FDBDIR%

rem xcopy /Y %GIT%\Fragment_DB.py %FDBDIR%
rem xcopy /Y %GIT%\helper_functions.py %FDBDIR%
rem xcopy /Y %GIT%\FragmentDB_handler.py %FDBDIR%


rem xcopy /Y %GIT%\fragment-database.sqlite %FDBDIR%
rem del %FDBDIR%\*.pyc

rem ##########################################

xcopy /Y %GIT%\*.py  %FDBDIR2%
xcopy /Y %GIT%\*.htm  %FDBDIR2%
xcopy /Y %GIT%\*.md  %FDBDIR2%
xcopy /Y %GIT%\fragment-database.sqlite %FDBDIR2%
xcopy /Y %GIT%\drawstyle.cds %FDBDIR2%
xcopy /Y %GIT%\fragmentdb.phil %FDBDIR2%

rem xcopy /Y %GIT%\__init__.py %FDBDIR2%
rem xcopy /Y %GIT%\fragmentdb.htm %FDBDIR2%
rem xcopy /Y %GIT%\fragmentdb.py %FDBDIR2%
rem xcopy /Y %GIT%\fragmentdb_handler.py %FDBDIR2%
rem xcopy /Y %GIT%\helper_functions.py %FDBDIR2%
rem xcopy /Y %GIT%\inputfrag.htm %FDBDIR2%
rem xcopy /Y %GIT%\plugins.xld %FDBDIR2%
rem xcopy /Y %GIT%\README.md %FDBDIR2%
rem xcopy /Y %GIT%\refine_model_tasks.py %FDBDIR2%


pause