# CCP4i2 Windows Test VM Setup Script
# Run this in an elevated PowerShell session after RDP-ing into the VM.
#
# Usage:
#   1. Download CCP4 development tarball and extract with:
#        tar -xzf ccp4-XXXXXXXX.tar.gz -C C:\Users\ccp4admin\Downloads
#   2. Open PowerShell as Administrator
#   3. Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#   4. .\setup-windows-test-vm.ps1 -CCP4Dir C:\Users\ccp4admin\Downloads\ccp4-XXXXXXXX
#
# This script will:
#   - Install Visual C++ runtime (required for gemmi and other CCP4 DLLs)
#   - Install Git (winget is not available on Windows Server)
#   - Move aside the bundled ccp4i2 to avoid conflicts
#   - Clone the ccp4i2 repo (django branch)
#   - Install Python dependencies into CCP4's Python
#   - Run Django migrations
#   - Download the latest Electron app from CI artifacts

param(
    [Parameter(Mandatory=$true)]
    [string]$CCP4Dir,
    [string]$WorkDir = "C:\ccp4i2",
    [string]$GitHubRepo = "ccp4/ccp4i2",
    [string]$Branch = "django"
)

$ErrorActionPreference = "Stop"

Write-Host "=== CCP4i2 Windows Test VM Setup ===" -ForegroundColor Cyan

# --- Step 1: Validate CCP4 installation ---
Write-Host "`n--- Checking CCP4 installation ---" -ForegroundColor Cyan

# Find ccp4-python - check both bin/ and Python311/ locations
$ccp4Python = $null
$pythonDir = $null

# Check for ccp4-python.bat in bin/
$batPath = Join-Path $CCP4Dir "bin\ccp4-python.bat"
if (Test-Path $batPath) {
    $ccp4Python = $batPath
}

# Also find the raw Python interpreter (needed for some operations)
foreach ($pyDir in @("Python311", "Python312", "Python39")) {
    $candidate = Join-Path $CCP4Dir "$pyDir\python.exe"
    if (Test-Path $candidate) {
        $pythonDir = Join-Path $CCP4Dir $pyDir
        if (-not $ccp4Python) {
            $ccp4Python = $candidate
        }
        break
    }
}

if (-not $ccp4Python) {
    Write-Host @"

CCP4 Python not found in $CCP4Dir.
Expected one of:
  $CCP4Dir\bin\ccp4-python.bat
  $CCP4Dir\Python311\python.exe

Please check the CCP4Dir path and try again.
"@ -ForegroundColor Red
    exit 1
}
Write-Host "[OK] CCP4 found at: $CCP4Dir" -ForegroundColor Green
Write-Host "     Python: $ccp4Python" -ForegroundColor Green

# --- Step 2: Install Visual C++ Redistributable ---
# Required for gemmi_ext and other CCP4 native DLLs on fresh Windows Server
Write-Host "`n--- Installing Visual C++ Redistributable ---" -ForegroundColor Cyan

$vcTestDll = "$env:SystemRoot\System32\msvcp140.dll"
if (Test-Path $vcTestDll) {
    Write-Host "[OK] Visual C++ runtime already installed" -ForegroundColor Green
} else {
    Write-Host "Downloading Visual C++ Redistributable..." -ForegroundColor Yellow
    $vcRedistUrl = "https://aka.ms/vs/17/release/vc_redist.x64.exe"
    $vcRedistPath = "$env:TEMP\vc_redist.x64.exe"
    Invoke-WebRequest -Uri $vcRedistUrl -OutFile $vcRedistPath
    Write-Host "Installing (silent)..."
    Start-Process -FilePath $vcRedistPath -ArgumentList "/install", "/quiet", "/norestart" -Wait
    Write-Host "[OK] Visual C++ runtime installed" -ForegroundColor Green
}

# Verify gemmi loads with CCP4 bin on PATH
Write-Host "Verifying gemmi loads correctly..."
$env:Path = "$(Join-Path $CCP4Dir 'bin');$env:Path"
try {
    & $ccp4Python -c "import gemmi; print('gemmi OK')" 2>&1
    if ($LASTEXITCODE -ne 0) { throw "gemmi import failed" }
    Write-Host "[OK] gemmi loads successfully" -ForegroundColor Green
} catch {
    Write-Host "[WARN] gemmi import failed - CCP4 native libraries may have issues" -ForegroundColor Yellow
    Write-Host "       You may need to check DLL dependencies manually" -ForegroundColor Yellow
}

# --- Step 3: Install Git ---
Write-Host "`n--- Installing Git ---" -ForegroundColor Cyan

$gitPath = Get-Command git -ErrorAction SilentlyContinue
if ($gitPath) {
    Write-Host "[OK] Git already installed: $($gitPath.Source)" -ForegroundColor Green
} else {
    Write-Host "Downloading Git for Windows..." -ForegroundColor Yellow
    $gitInstallerUrl = "https://github.com/git-for-windows/git/releases/download/v2.47.1.windows.2/Git-2.47.1.2-64-bit.exe"
    $gitInstallerPath = "$env:TEMP\git-installer.exe"
    Invoke-WebRequest -Uri $gitInstallerUrl -OutFile $gitInstallerPath
    Write-Host "Installing Git (silent)..."
    Start-Process -FilePath $gitInstallerPath -ArgumentList "/VERYSILENT", "/NORESTART" -Wait
    # Refresh PATH
    $env:Path = [System.Environment]::GetEnvironmentVariable("Path", "Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path", "User")
    Write-Host "[OK] Git installed" -ForegroundColor Green
}

# --- Step 4: Move aside bundled ccp4i2 ---
Write-Host "`n--- Preparing CCP4 Python environment ---" -ForegroundColor Cyan

if ($pythonDir) {
    $bundledCcp4i2 = Join-Path $pythonDir "Lib\site-packages\ccp4i2"
    $bundledCcp4i2Old = Join-Path $pythonDir "Lib\site-packages\ccp4i2_bundled"
    if ((Test-Path $bundledCcp4i2) -and -not (Test-Path $bundledCcp4i2Old)) {
        Write-Host "Moving aside bundled ccp4i2 to avoid conflicts..." -ForegroundColor Yellow
        Rename-Item -Path $bundledCcp4i2 -NewName "ccp4i2_bundled"
        Write-Host "[OK] Bundled ccp4i2 moved to ccp4i2_bundled" -ForegroundColor Green
    } elseif (Test-Path $bundledCcp4i2Old) {
        Write-Host "[OK] Bundled ccp4i2 already moved aside" -ForegroundColor Green
    } else {
        Write-Host "[OK] No bundled ccp4i2 found (clean environment)" -ForegroundColor Green
    }
}

# --- Step 5: Clone the repository ---
Write-Host "`n--- Cloning ccp4i2 repository ---" -ForegroundColor Cyan

if (Test-Path $WorkDir) {
    Write-Host "[OK] $WorkDir already exists, pulling latest..." -ForegroundColor Green
    Push-Location $WorkDir
    git fetch origin
    git checkout $Branch
    git pull origin $Branch
    Pop-Location
} else {
    git clone --branch $Branch "https://github.com/$GitHubRepo.git" $WorkDir
}

# --- Step 6: Install Python dependencies into CCP4 Python ---
Write-Host "`n--- Installing Python dependencies ---" -ForegroundColor Cyan

# Set up CCP4 environment
$env:CCP4 = $CCP4Dir
$env:CCP4I2_TOP = $WorkDir
$env:CCP4I2_BACKEND = "django"
$env:DJANGO_SETTINGS_MODULE = "ccp4i2.config.settings"

Push-Location "$WorkDir\server"

Write-Host "Installing Django and dependencies..."
& $ccp4Python -m pip install `
    "django>=5.0" `
    "djangorestframework>=3.15" `
    "django-cors-headers>=4.0" `
    "dj-database-url>=2.0" `
    "python-dotenv>=1.0" `
    "openpyxl>=3.1.0" `
    "uvicorn[standard]" `
    "pytz"

# Install in editable mode if pyproject.toml exists
if (Test-Path "pyproject.toml") {
    Write-Host "Installing ccp4i2 in editable mode..."
    & $ccp4Python -m pip install --no-cache-dir -e .
}

# --- Step 7: Run Django migrations ---
Write-Host "`n--- Running Django migrations ---" -ForegroundColor Cyan

& $ccp4Python -m django migrate
if ($LASTEXITCODE -eq 0) {
    Write-Host "[OK] Database migrations complete" -ForegroundColor Green
} else {
    Write-Host "[WARN] Migrations failed - server may not start correctly" -ForegroundColor Yellow
}

Pop-Location

# --- Step 8: Download latest Electron app from CI ---
Write-Host "`n--- Downloading Electron app from CI ---" -ForegroundColor Cyan

$electronDir = "$WorkDir\electron-app"
if (-not (Test-Path $electronDir)) {
    New-Item -ItemType Directory -Path $electronDir | Out-Null
}

# Check if gh CLI is available
$ghPath = Get-Command gh -ErrorAction SilentlyContinue
if (-not $ghPath) {
    Write-Host @"

GitHub CLI (gh) is not installed. To download the Electron app:
  1. Download manually from:
     https://github.com/$GitHubRepo/actions/workflows/electron-multiplatform-build.yml
  2. Save the .exe to: $electronDir

"@ -ForegroundColor Yellow
} else {
    $ghAuth = gh auth status 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Host @"

GitHub CLI is not authenticated. To download the Electron app from CI:
  1. Run: gh auth login
  2. Then run: gh run download --repo $GitHubRepo --name windows-installers --dir $electronDir

Or download manually from:
  https://github.com/$GitHubRepo/actions/workflows/electron-multiplatform-build.yml

"@ -ForegroundColor Yellow
    } else {
        Write-Host "Downloading latest Windows Electron build..."
        try {
            gh run download --repo $GitHubRepo --name windows-installers --dir $electronDir
            Write-Host "[OK] Electron app downloaded to: $electronDir" -ForegroundColor Green
        } catch {
            Write-Host "Could not download automatically. Download manually from GitHub Actions." -ForegroundColor Yellow
        }
    }
}

# --- Step 9: Create convenience shortcuts ---
Write-Host "`n--- Creating desktop shortcuts ---" -ForegroundColor Cyan

$startServerScript = @"
@echo off
echo Starting CCP4i2 Django dev server...
set CCP4=$CCP4Dir
set PATH=$CCP4Dir\bin;%PATH%
set CCP4I2_TOP=$WorkDir
set CCP4I2_BACKEND=django
set DJANGO_SETTINGS_MODULE=ccp4i2.config.settings
cd /d $WorkDir\server
"$ccp4Python" manage.py runserver 0.0.0.0:8000
pause
"@

$startServerPath = "$WorkDir\start-server.bat"
Set-Content -Path $startServerPath -Value $startServerScript

# Desktop shortcut for server
$desktop = [Environment]::GetFolderPath("Desktop")
$shortcutPath = Join-Path $desktop "CCP4i2 Dev Server.lnk"
$shell = New-Object -ComObject WScript.Shell
$shortcut = $shell.CreateShortcut($shortcutPath)
$shortcut.TargetPath = $startServerPath
$shortcut.WorkingDirectory = "$WorkDir\server"
$shortcut.Save()

# --- Summary ---
Write-Host @"

=== Setup Complete ===

CCP4:           $CCP4Dir
Repository:     $WorkDir (branch: $Branch)
Electron app:   $electronDir

To start the Django dev server:
  Double-click "CCP4i2 Dev Server" on the desktop
  OR run: $startServerPath

To run the Electron app:
  Run the .exe from $electronDir

To run tests:
  cd $WorkDir\server
  $ccp4Python -m pytest ccp4i2\tests\ -v

"@ -ForegroundColor Green
