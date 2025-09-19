# Create and activate the Corvus environment
python -m venv "$HOME\.venv\Corvus"
& "$HOME\.venv\Corvus\Scripts\Activate.ps1"

# Install corvus
pip install .

# Create the corvus_shell shortcut on the desktop.
$TargetFile = (Get-Command powershell).Path
$Arguments = "-NoExit -File `"$HOME\.venv\Corvus\Scripts\Activate.ps1`""
$ShortcutFile = "$HOME\Desktop\corvus.lnk"

$WScriptShell = New-Object -ComObject WScript.Shell
$Shortcut = $WScriptShell.CreateShortcut($ShortcutFile)
$Shortcut.TargetPath = $TargetFile
$Shortcut.Arguments = $Arguments
$Shortcut.Save()