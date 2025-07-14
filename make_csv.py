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







