import re
import csv

# Input and output file names
input_file = 'test-stats.txt'
output_file = 'output.csv'

# Regular expression patterns for variables
pattern_header = re.compile(r"\(n, α, _iter\) = \((\d+), ([\d\.]+), (\d+)\)")
pattern_data = re.compile(r"(.*?):\s+([\d\.]+(?:e[\+-]?\d+)?|\d+|OPTIMAL)")

# List to hold all rows of data
rows = []

# Parse the file
with open(input_file, 'r') as file:
    current_row = {}
    for line in file:
        # Check for the header line with independent variables
        header_match = pattern_header.match(line)
        if header_match:
            # Save the current row if it has data
            if current_row:
                rows.append(current_row)
            current_row = {
                'n': int(header_match.group(1)),
                'α': float(header_match.group(2)),
                '_iter': int(header_match.group(3))
            }

        # Check for dependent variable lines
        data_match = pattern_data.match(line)
        if data_match:
            key = data_match.group(1).strip()
            value = data_match.group(2).strip()
            # Convert value to float if possible
            try:
                value = float(value)
            except ValueError:
                pass  # Keep as string if not convertible
            current_row[key] = value

    # Add the last row
    if current_row:
        rows.append(current_row)

# Extract all unique keys (column names)
columns = sorted({key for row in rows for key in row.keys()})

# Write to CSV
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)

print(f"CSV table has been written to {output_file}")
