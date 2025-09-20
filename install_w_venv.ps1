# Create and activate the Corvus environment
python -m venv "$HOME\.venv\Corvus"
& "$HOME\.venv\Corvus\Scripts\Activate.ps1"

# Install corvus
pip install .

# Create the corvus_shell shortcut on the desktop.
$TargetFile = (Get-Command powershell).Path
$Arguments = "-NoExit -File `"$HOME\.venv\Corvus\Scripts\Activate.ps1`""
$ShortcutFile = "$HOME\Desktop\corvus_shell.lnk"

$WScriptShell = New-Object -ComObject WScript.Shell
$Shortcut = $WScriptShell.CreateShortcut($ShortcutFile)
$Shortcut.TargetPath = $TargetFile
$Shortcut.Arguments = $Arguments
$Shortcut.Save()

# Now download and install scigui
$response = Invoke-WebRequest -Uri "https://github.com/times-software/SciGUI/archive/refs/heads/main.zip" -OutFile SciGUI.zip
# Check if SciGUI folder exists, ask to clobber if it does
if (Test-Path -Path ".\SciGUI" -PathType Container) {
    $response = Read-Host "Directory SciGUI already exists. Overwrite? (Y/N)[Y]"

    if ($response -eq 'N' -or $response -eq 'n') {
        Write-Host "Aborting installation."
        exit
    }
}

Expand-Archive -Force -Path "SciGUI.zip" -DestinationPath "SciGUI"
cd SciGUI/SciGUI-main
pip install .

# Create a startup script that activates the Corvus environment and launches corvus
echo "~\.venv\Corvus\Scripts\Activate.ps1" > ~\.venv\Corvus\Scripts\corvus_launcher.ps1
echo "corvus.exe" >> ~\.venv\Corvus\Scripts\corvus_launcher.ps1

# Finally, create a shortcut for corvus (GUI).
$TargetFile = (Get-Command powershell).Path
$Arguments = "-File `"$HOME\.venv\Corvus\Scripts\corvus_launcher.ps1`"" 
$ShortcutFile = "$HOME\Desktop\corvus.lnk"
$WScriptShell = New-Object -ComObject WScript.Shell
$Shortcut = $WScriptShell.CreateShortcut($ShortcutFile)
$Shortcut.TargetPath = $TargetFile
$Shortcut.Arguments = $Arguments
$Shortcut.Save()

cd "..\.."

# Install jmol if user wants it.
$response = Read-Host "Corvus can use jmol to visualize structures. Download and install? (Y/N)[Y]"

if ($response -eq 'N' -or $response -eq 'n') {
    Write-Host "Not installing jmol. Installation finishing."
    exit
}

Invoke-WebRequest -UserAgent "Wget" -Uri http://sourceforge.net/projects/jmol/files/latest/download -OutFile jmol.zip
Expand-Archive -F jmol.zip
cd jmol
$jmoldir = (ls | Select-Object -ExpandProperty Name)

$source = $jmoldir
$destination = "$HOME\.venv\Corvus\share\"

# Check if the destination file exists
if (Test-Path $destination\$jmoldir) {
    # If jmol directory already found.
    Write-Host "Previous installation of Jmol found. If you want to reinstall"
    echo "please remove the directory: $destination\$jmoldir"
    echo "and re-run the installer."
} else {
    Move-Item -Path "$source" -Destination "$destination\$jmoldir"
}

cd ..

# Create scigui.ini file to tell it where jmol is.
$OutputEncoding = [System.Text.Encoding]::UTF8
echo "[visualization]" > "$HOME\.Corvus\scigui.ini"
echo "jmol_path = `"$HOME\.venv\Corvus\share\$jmoldir\jmol.bat`"" >> "$HOME\.Corvus\scigui.ini" 
echo "Finished installing corvus."
