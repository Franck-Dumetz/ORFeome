#!/bin/bash
if [[ "$(uname)" == "Linux" ]]; then
  export LD_LIBRARY_PATH="$OLD_LD_LIBRARY_PATH"
  unset OLD_LD_LIBRARY_PATH
fi
