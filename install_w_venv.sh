#!/bin/bash -x
# Create a venv environment
python3 -m venv "$HOME/.venv/Corvus"
source "$HOME/.venv/Corvus/bin/activate"
pip3 install .

# Create a shell script that will start an interactive session when double clicked.
shell=`basename $SHELL`
echo "User shell: $shell"
echo "#!/bin/$shell" > "$HOME/Desktop/corvus_shell"
if [ $shell=="zsh" ]
then
	echo "export ZDOTDIR="$HOME/.venv/Corvus/bin/"" >> "$HOME/Desktop/corvus_shell"
	# Copy activation script to a new startup shell file
	cp "$HOME/.venv/Corvus/bin/activate" "$HOME/.venv/Corvus/bin/.zshrc"
	echo "exec $SHELL -i" >> "$HOME/Desktop/corvus_shell"
elif [ $shell == "bash" ]
then
	echo "exec $SHELL --rcfile "$HOME/.venv/Corvus/bin"" >> "$HOME/Desktop/corvus_shell"
elif [ $shell == "tcsh" ]
then
	echo "exec $shell -c 'source "$HOME/.venv/Corvus/bin/activate"'; tcsh" >> "$HOME/Desktop/corvus_shell"
fi

# Install SciGUI if user want to.
echo 'Would you like to install the corvus graphical user interface? (Y/N)[Y]' 
read a
echo "Response: $a"
if [ ! "X$a"=="XN" ]
then
	echo "Will not install the corvus GUI. Finished installation of corvus."
else

	# Download and install SciGUI
	curl -L "https://github.com/times-software/SciGUI/archive/refs/heads/main.zip" -o SciGUI.zip
	unzip SciGUI.zip
	cd SciGUI-main
	pip3 install .
	
	# Now make the corvus shortcut on the desktop.
	echo "#!/bin/$shell" > "$HOME/Desktop/corvus"
	echo "source \"$HOME/.venv/Corvus/bin/activate\"" >> "$HOME/Desktop/corvus"
	echo "corvus" >> "$HOME/Desktop/corvus"
	cd ..
fi

