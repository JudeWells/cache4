import re

def read_sdf_file(sdf_file):
    with open(sdf_file, 'r') as f:
        sdf_string = f.read()
    return sdf_string

def reformat_sdf(sdf_string):
    molecule_sections = sdf_string.split('$$$$')
    new_sections = []

    for section in molecule_sections:
        lines = section.split('\n')
        new_lines = []

        for line in lines:
            if '<idnumber>' in line:
                new_lines[0] = line.strip()
                new_lines.append(line)
            else:
                new_lines.append(line)

        new_sections.append('\n'.join(new_lines))

    return '$$$$\n'.join(new_sections)

def write_sdf_file(sdf_string, sdf_file):
    with open(sdf_file, 'w') as f:
        f.write(sdf_string)

sdf_string = read_sdf_file('/Users/judewells/Documents/dataScienceProgramming/enamine_collections/Enamine_legacy_collection_202305.sdf')
new_sdf_string = reformat_sdf(sdf_string)
write_sdf_file(new_sdf_string, '/Users/judewells/Documents/dataScienceProgramming/enamine_collections/Enamine_legacy_collection_202305_reformatted.sdf')
