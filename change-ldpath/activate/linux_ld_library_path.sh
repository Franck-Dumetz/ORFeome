#!/bin/bash
if [[ "$(uname)" == "Linux" ]]; then
  export OLD_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
fi
