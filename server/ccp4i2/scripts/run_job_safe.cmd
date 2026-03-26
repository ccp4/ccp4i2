@echo off
REM Crash-safe job runner wrapper (Windows version).
REM
REM Catches process crashes (e.g., segfault in a C extension library) and marks
REM the job as FAILED in the database, preventing stuck "running" jobs.
REM
REM See run_job_safe.sh for the Unix equivalent.
REM
REM Usage: run_job_safe.cmd <python> <job_uuid>

set PYTHON=%~1
set JOB_UUID=%~2

if "%PYTHON%"=="" goto :usage
if "%JOB_UUID%"=="" goto :usage

"%PYTHON%" -m django run_job -ju "%JOB_UUID%" -y
set EXIT_CODE=%ERRORLEVEL%

if %EXIT_CODE% neq 0 (
    echo [run_job_safe] Job %JOB_UUID%: exited with code %EXIT_CODE% 1>&2
    echo [run_job_safe] Marking job %JOB_UUID% as FAILED 1>&2
    "%PYTHON%" -m django set_job_status -ju "%JOB_UUID%" -s FAILED 2>nul
    if %ERRORLEVEL% neq 0 (
        echo [run_job_safe] WARNING: could not update job status 1>&2
    )
)

exit /b %EXIT_CODE%

:usage
echo [run_job_safe] Usage: run_job_safe.cmd ^<python^> ^<job_uuid^> 1>&2
exit /b 1
