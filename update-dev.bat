
set GIT="F:\GitHub\DSR-db"
set OLEX="F:\Programme\Olex2-1.2-dev\etc\scripts\"

set GIT2="C:\Users\daniel\Documents\GitHub\DSR-db"
set OLEX2="C:\Program Files\Olex2-1.2-dev\etc\scripts\"

xcopy /Y %GIT%\example.py %OLEX%
xcopy /Y %GIT%\fragmentdb.py %OLEX%

xcopy /Y %GIT2%\example.py %OLEX2%
xcopy /Y %GIT2%\fragmentdb.py %OLEX2%

rem pause