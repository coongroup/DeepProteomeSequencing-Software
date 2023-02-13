def get_stat_evidence(file_name):
    rcnt = 0
    fcnt = 0
    with open(file_name) as file_stream:
        header = file_stream.readline().rstrip().split('\t')
        for i in range(len(header)):
            print(i, header[i])
        reverse_index = header.index("Reverse")
        contaminant_index = header.index("Potential contaminant")
        for line in file_stream:
            spl = line.rstrip().split("\t")
            if spl[contaminant_index] == "+":
                continue
            if spl[reverse_index] == "+":
                rcnt += 1
            else:
                fcnt += 1
    return rcnt, fcnt


def get_stat_peptides(file_name):
    rcnt = 0
    fcnt = 0
    with open(file_name) as file_stream:
        header = file_stream.readline().rstrip().split('\t')
        for i in range(len(header)):
            print(i, header[i])
        reverse_index = header.index("Reverse")
        contaminant_index = header.index("Potential contaminant")
        # print(header[reverse_index])
        for line in file_stream:
            spl = line.rstrip().split("\t")
            # if spl[reverse_index] != "":
            #     spl[reverse_index]
            if spl[contaminant_index] == "+":
                continue
            if spl[reverse_index] == "+":
                rcnt += 1
            else:
                fcnt += 1
    return rcnt, fcnt


def get_stat_proteins(file_name):
    rcnt = 0
    fcnt = 0
    with open(file_name) as file_stream:
        header = file_stream.readline().rstrip().split('\t')
        for i in range(len(header)):
            print(i, header[i])
        reverse_index = header.index("Reverse")
        contaminant_index = header.index("Potential contaminant")
        mod_group_index = header.index("Only identifcation by site")
        # print(header[reverse_index])
        for line in file_stream:
            spl = line.rstrip().split("\t")
            # if spl[reverse_index] != "":
            #     spl[reverse_index]
            if spl[contaminant_index] == "+" or spl[mod_group_index] == "+":
                continue
            if spl[reverse_index] == "+":
                rcnt += 1
            else:
                fcnt += 1
    return rcnt, fcnt


if __name__ == "__main__":
    evidence_file = "D:\\data\\coon\\standard\\proteomics\\evidence.txt"
    peptides_file = "D:\\data\\coon\\standard\\proteomics\\peptides.txt"
    protein_groups_file = "D:\\data\\coon\\standard\\proteomics\\proteinGroups.txt"

    print("Evidence:", get_stat_evidence(evidence_file))
    # print("Peptides:", get_stat_peptides(peptides_file))
    # print("Proteins:", get_stat_proteins(protein_groups_file))
