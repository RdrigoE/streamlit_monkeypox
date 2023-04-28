import csv
import sys
from typing import Iterable


def read_csv(file_name: str) -> list[list[str]]:
    with open(file_name) as handler:
        return list(csv.reader(handler))


def transform_into_dict(to_iter: Iterable) -> dict[str, str]:
    new_dict = {}

    for k, v in to_iter:
        new_dict[k] = v
    print(new_dict)
    return new_dict


def main():
    depth_file = sys.argv[1]
    csv_file = sys.argv[2]

    depth_info = read_csv(depth_file)
    csv_info = read_csv(csv_file)

    depth_dict = transform_into_dict(depth_info)

    csv_info[0].append("sample_depth")
    for row in csv_info[1:]:
        row.append(depth_dict["sample_" + row[0]])

    with open(csv_file, "w") as handler:
        writer = csv.writer(handler)
        for row in csv_info:
            writer.writerow(row)


if __name__ == "__main__":
    main()
