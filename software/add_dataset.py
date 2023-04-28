import csv
import sys
with open(sys.argv[1], "r", encoding="UTF-8") as handler:
    datasets_list = list(csv.reader(handler))

datasets = {}
for row in datasets_list:
    datasets[row[0]] = row[1]


with open(sys.argv[2], "r", encoding="UTF-8") as handler:
    entries = list(csv.reader(handler))

with open(sys.argv[2], "w", encoding="UTF-8") as handler:
    writer = csv.writer(handler)
    entries[0].append("dataset")
    writer.writerow(entries[0])

    for entry in entries[1:]:
        entry.append(datasets[entry[0]])
        writer.writerow(entry)
