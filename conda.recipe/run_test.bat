REM nosetests -s -vv --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS
"%PYTHON%" -m nose
if errorlevel 1 exit 1
