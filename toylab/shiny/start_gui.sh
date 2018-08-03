#!/bin/bash

##
##  Shiny semantic similarity GUI
##

## Press Ctrl-C to stop GUI server

if [[ $# -lt 1 ]]
then
  echo "Usage:  ./start_gui.sh model.rda [model2.rda ...]"
  exit 1
fi

R --no-save -e "shiny::runApp('.', port=8200L, launch.browser=TRUE)" --args "$@"
