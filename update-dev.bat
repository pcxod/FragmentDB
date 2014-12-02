
set GIT="F:\GitHub\DSR-db"
set GIT2="C:\Users\daniel\Documents\GitHub\DSR-db"

set OLEX="D:\Programme\Olex2-1.2-dev\etc\scripts\"
set OLEX2="C:\Program Files\Olex2-1.2-dev\etc\scripts\"
set OLEX3="C:\Program Files\Olex2-1.2-dev\Olex2-1.2-alpha\etc\scripts\"

xcopy /Y %GIT%\example.py %OLEX%
xcopy /Y %GIT%\fragmentdb.py %OLEX%
xcopy /Y %GIT%\dk-database.sqlite %OLEX%

xcopy /Y %GIT2%\example.py %OLEX2%
xcopy /Y %GIT2%\fragmentdb.py %OLEX2%
xcopy /Y %GIT2%\dk-database.sqlite %OLEX2%

xcopy /Y %GIT2%\example.py %OLEX3%
xcopy /Y %GIT2%\fragmentdb.py %OLEX3%
xcopy /Y %GIT2%\dk-database.sqlite %OLEX3%

pause