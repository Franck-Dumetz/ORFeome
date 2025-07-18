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
# along with this program.  If not, see https://www.gnu.org/licenses/.

# Converts a two-column CSV of sample names and conditions into a treatment design file for ORFeome analysis.
# Depending on the flag (-u, -m, or both), it appends _m1, _m10, or both to generate condition labels.

import csv
import sys

if sys.argv[1] == "-u":
    with open(sys.argv[2], mode='r', newline='') as i:
        reader = csv.reader(i)
        with open("treatments.csv", mode="w", newline='') as o:
            writer = csv.writer(o)
            writer.writerow(["Sample", "Condition"])
            for row in reader:
                name = row[0].lstrip('\ufeff')
                treat = row[1].lstrip('\ufeff')
                mod_row = [f"{name}_m1", f"{treat}_m1"]
                writer.writerow(mod_row)

elif sys.argv[1] == "-m":
    with open(sys.argv[2], mode='r', newline='') as i:
        reader = csv.reader(i)
        with open("treatments.csv", mode="w", newline='') as o:
            writer = csv.writer(o)
            writer.writerow(["Sample", "Condition"])
            for row in reader:
                name = row[0].lstrip('\ufeff')
                treat = row[1].lstrip('\ufeff')
                mod_row = [f"{name}_m10", f"{treat}_m10"]
                writer.writerow(mod_row)

else:
    with open(sys.argv[2], mode='r', newline='') as i:
        reader = csv.reader(i)
        with open("treatments.csv", mode='w', newline='') as o:
            writer = csv.writer(o)
            writer.writerow(["Sample", "Condition"])
            for row in reader:
                name = row[0].lstrip('\ufeff')
                treat = row[1].lstrip('\ufeff')
                mod_row = [f"{name}_m1", f"{treat}_m1"]
                writer.writerow(mod_row)
                mod_row = [f"{name}_m10", f"{treat}_m10"]
                writer.writerow(mod_row)







