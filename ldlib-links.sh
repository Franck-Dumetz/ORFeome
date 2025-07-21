mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
mkdir -p "$CONDA_PREFIX/etc/conda/deactivate.d"

cp change-ldpath/activate/linux_ld_library_path.sh "$CONDA_PREFIX/etc/conda/activate.d/"
cp change-ldpath/deactivate/linux_ld_library_path.sh "$CONDA_PREFIX/etc/conda/deactivate.d/"
