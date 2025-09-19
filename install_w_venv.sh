#!/bin/bash -x
python -m venv ~/.venv/Corvus
# Windows powershell script
# ~/.venv/Corvus/Scripts/Activate.ps1
source ~/.venv/Corvus/bin/activate
python -m pip install .

# Create a shell script that will start an interactive session when double clicked.
# Copy activation script to a new startup shell file
cp ~/.venv/Corvus/bin/activate ~/.venv/Corvus/bin/.zshrc
shell=`basename $SHELL`
echo "User shell: $shell"
echo "#!/bin/$shell" > ~/Desktop/corvus_shell
if [ $shell=="zsh" ]
then
	echo "export ZDOTDIR=~/.venv/Corvus/bin/" >> ~/Desktop/corvus_shell
	cp ~/.venv/Corvus/bin/activate ~/.venv/Corvus/bin/.zshrc
	echo "exec $SHELL -i" >> ~/Desktop/corvus_shell
elif [ $shell == "bash" ]
then
	echo "exec $SHELL --rcfile ~/.venv/Corvus/bin" >> ~/Desktop/corvus_shell
elif [ $shell == "tcsh" ]
then
	echo "exec $shell -c 'source ~/.venv/Corvus/bin/activate'; tcsh" >> ~/Desktop/corvus_shell
fi
