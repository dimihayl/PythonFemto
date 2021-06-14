#!/bin/bash

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

cd bin

if ! cmake ../CMake; then
	echo "Configuration of the script ${red}failed${reset}"
	return 3
fi

echo "The configuration of the script was ${green}successful${reset}"
#echo "  Configured using "$1
echo "  To proceed type: make"

cd ..

return 0
