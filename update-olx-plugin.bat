@ECHO OFF


echo %FDBDIR%
set FDBDIR="C:\Program Files\Olex2-1.2\util\pyUtil\PluginLib\plugin-Fragment_DB"
set FDBDIR2="D:\downloads\olex2-gui\util\pyUtil\PluginLib\plugin-Fragment_DB"
set GIT="D:\GitHub\DSR-db\plugin-Fragment_DB"


xcopy /Y %GIT%\Fragment_DB.htm %FDBDIR%
xcopy /Y %GIT%\Fragment_DB.phil %FDBDIR%
xcopy /Y %GIT%\Fragment_DB.py %FDBDIR%
xcopy /Y %GIT%\FragmentDB_handler.py %FDBDIR%
xcopy /Y %GIT%\inputfrag.htm %FDBDIR% 

xcopy /Y %GIT%\fragment-database.sqlite %FDBDIR%

rem ##########################################

xcopy /Y %GIT%\Fragment_DB.htm %FDBDIR2%
xcopy /Y %GIT%\Fragment_DB.phil %FDBDIR2%
xcopy /Y %GIT%\Fragment_DB.py %FDBDIR2%
xcopy /Y %GIT%\FragmentDB_handler.py %FDBDIR2%
xcopy /Y %GIT%\inputfrag.htm %FDBDIR2% 

xcopy /Y %GIT%\fragment-database.sqlite %FDBDIR2%



del %FDBDIR%\*.pyc