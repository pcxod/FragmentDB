@ECHO OFF


echo %FDBDIR%
rem set FDBDIR="C:\Program Files\Olex2-1.2\util\pyUtil\PluginLib\plugin-Fragment_DB"
set FDBDIR2="D:\downloads\olex2-gui\util\pyUtil\PluginLib\plugin-Fragment_DB"
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

xcopy /Y %GIT%\Fragment_DB.htm %FDBDIR2%
xcopy /Y %GIT%\helptext.htm %FDBDIR2%
xcopy /Y %GIT%\inputfrag.htm %FDBDIR2% 
xcopy /Y %GIT%\Fragment_DB.phil %FDBDIR2%
xcopy /Y %GIT%\Fragment_DB.py %FDBDIR2%
xcopy /Y %GIT%\helper_functions.py %FDBDIR2% 
xcopy /Y %GIT%\FragmentDB_handler.py %FDBDIR2%
xcopy /Y %GIT%\refine_model_tasks.py %FDBDIR2%
xcopy /Y %GIT%\natsort.py %FDBDIR2%

xcopy /Y %GIT%\fragment-database.sqlite %FDBDIR2%


xcopy /Y %FDBDIR2%\FragmentDB_handler.pyc %GIT%
xcopy /Y %FDBDIR2%\Fragment_DB.pyc %GIT%
xcopy /Y %FDBDIR2%\helper_functions.pyc %GIT%
xcopy /Y %FDBDIR2%\refine_model_tasks.pyc %GIT%
xcopy /Y %FDBDIR2%\natsort.pyc %GIT%

