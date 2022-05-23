import sys
import os
import math

files = os.listdir(sys.argv[1])

for file in files:
    matrix_file = open(sys.argv[1] + file, 'r')
    lines = matrix_file.readlines()

    idx = 0
    for line in lines:
        if not line.startswith("%"):
            break

        idx += 1

    matrix_properties = lines[idx].split()
    rows_size = int(matrix_properties[0])

    rows = [0] * rows_size

    for line in lines[idx + 1:]:
        rows[int(line.split()[0]) - 1] += 1

    max_rows = max(rows)

    print(f"{file} rows: {max_rows} seconds: "
          f"{math.ceil(rows_size/3072.0) * (max_rows * 0.000006 + max_rows * 0.000006 + rows_size * 0.00000032)}")
