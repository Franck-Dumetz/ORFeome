#!/bin/bash

# ORFeome – Analyzing ORFeome screening data
# Copyright (C) 2025 Anushka Shome and Franck Dumetz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see https://www.gnu.org/licenses/.

# Downloads SRA files listed in an input text file and converts them to FASTQ format using prefetch and fastq-dump.

mkdir -p fastqs

while read srr; do
    echo "Processing $srr..."

    prefetch "$srr"

    # Convert to FASTQ format
    fastq-dump --split-3 --outdir fastqs "$srr/$srr.sra"

    rm -r "$srr"

done < $1
