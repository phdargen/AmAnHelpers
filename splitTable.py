def split_latex_table(input_file, output_file1, output_file2, split_column):
    with open(input_file, 'r') as file:
        lines = [line for line in file.readlines() if line.strip()]  # Ignore empty lines

    # Find the beginning and end of the tabular environment
    start_idx = -1
    end_idx = -1
    for i, line in enumerate(lines):
        if '\\begin{tabular}' in line:
            start_idx = i
        if '\\end{tabular}' in line:
            end_idx = i + 1  # Include the end tabular line

    # If the tabular environment is not found, exit the function
    if start_idx == -1 or end_idx == -1:
        print("No tabular environment found in the input file.")
        return

    # Extract the tabular section
    tabular_content = lines[start_idx:end_idx]

    # Split the header to determine where to divide the columns
    header_parts = tabular_content[2].split('&')
    header1 = '&'.join(header_parts[:split_column + 1]) + ' \\\\\n'
    header2 = header_parts[0].strip() + ' & ' + '&'.join(header_parts[split_column + 1:]).strip() + ' \\\\\n'

    # Create the two parts of the table
    table1 = tabular_content[:3] + [header1]
    table2 = [tabular_content[0]] + [tabular_content[1]] + [header2]

    for line in tabular_content[3:-3]:  # Ignore the last hline
        parts = line.split('&')
        part1 = '&'.join(parts[:split_column + 1]) + ' \\\\\n'
        part2 = parts[0].strip() + ' & ' + '&'.join(parts[split_column + 1:]).strip() + ' \\\\\n'
        table1.append(part1)
        table2.append(part2)

    table1 += tabular_content[-3:]
    table2 += tabular_content[-3:]

    # Save the tables to new files
    with open(output_file1, 'w') as file1, open(output_file2, 'w') as file2:
        file1.writelines(table1)
        file2.writelines(table2)

    print(f"Tables saved as {output_file1} and {output_file2}")


# Split the table after the 6th column (Model 4)
split_latex_table('altModels_table.tex', 'tab1.tex', 'tab2.tex', 6)
