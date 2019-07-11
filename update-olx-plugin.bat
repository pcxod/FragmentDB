@ECHO OFF

rem This file is for development only.

echo %FDBDIR%
set FDBDIR2="D:\downloads\olex2-gui\util\pyUtil\FragmentDB"
set GIT="D:\GitHub\FragmentDB"


rem ##########################################

xcopy /Y %GIT%\*.py  %FDBDIR2%
xcopy /Y %GIT%\*.htm  %FDBDIR2%
xcopy /Y %GIT%\*.md  %FDBDIR2%
xcopy /Y %GIT%\fragment-database.sqlite %FDBDIR2%
xcopy /Y %GIT%\drawstyle.cds %FDBDIR2%
xcopy /Y %GIT%\fragmentdb.phil %FDBDIR2%



pause