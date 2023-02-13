from typing import Dict, List
from scipy.stats import wilcoxon

import statistics
import matplotlib.pyplot as plt

N = 72

experiments = {
    "WT_Tryp_rep1": [(i, i) for i in range(24)],
    "WT_Tryp_rep2": [(i, i) for i in range(24)],
    "delOCT1_Tryp_rep1": [(i, i) for i in range(24)],
}

identifications = [
    int("100"*8+"010"*8+"001"*8, 2),
    int("111"*8+"000"*8+"000"*8, 2),
    int("000"*8+"111"*8+"000"*8, 2),
    int("000"*8+"000"*8+"111"*8, 2),
    int("111"*8+"111"*8+"111"*8, 2)
]

# N = 8
# experiments = {
#     "WT_Tryp_rep1": [(1, 8), (9, 16), (17, 24)],
#     "WT_Tryp_rep2": [(1, 8), (9, 16), (17, 24)],
#     "delOCT1_Tryp_rep1": [(1, 8), (9, 16), (17, 24)],
# }

# identifications = [
#     int("100010001", 2),
#     int("111000000", 2),
#     int("000111000", 2),
#     int("000000111", 2),
#     int("100100100", 2),
#     int("010010010", 2),
#     int("001001001", 2),
#     int("111111111", 2)
# ]


def read_fasta(fasta_file):
    fasta_records = {}
    with open(fasta_file) as file_stream:
        name = ""
        for line in file_stream:
            if line.startswith(">"):
                name = line.rstrip().split(" ")[0].split("|")[1]
                fasta_records[name] = ""
            else:
                fasta_records[name] += line.rstrip()
    return fasta_records


class Peptide:
    def __init__(
            self,
            sequence: str,
            protein2index: Dict[str, List[int]],
            identification: int):
        self.sequence = sequence
        self.protein2index = protein2index
        self.identification = identification


def read_evidences(evidence_file, fasta_records):
    peptides = {}
    with open(evidence_file) as file_stream:
        header = file_stream.readline().rstrip().split('\t')
        sequence_index = header.index("Sequence")
        proteins_index = header.index("Proteins")
        experiment_index = header.index("Experiment")
        fraction_index = header.index("Fraction")
        reverse_index = header.index("Reverse")
        contaminant_index = header.index("Potential contaminant")
        for line in file_stream:
            spl = line.rstrip().split("\t")
            if spl[reverse_index] == "+" or spl[contaminant_index] == "+":
                continue
            sequence = spl[sequence_index]
            experiment = spl[experiment_index]
            fraction = int(spl[fraction_index])
            if experiment not in experiments:
                continue
            index = 0
            index0 = 0
            for e, rs in experiments.items():
                for r in rs:
                    if e == experiment and r[0] <= fraction <= r[1]:
                        index = index0
                    index0 += 1
            if sequence not in peptides:
                protein2index = {}
                for protein in spl[proteins_index].split(';'):
                    if protein not in fasta_records:
                        continue
                    pos = fasta_records[protein].find(sequence)
                    while pos != -1:
                        if protein not in protein2index:
                            protein2index[protein] = []
                        protein2index[protein].append(pos)
                        pos = fasta_records[protein].find(sequence, pos + 1)
                if len(protein2index) == 0:
                    continue
                peptides[sequence] = Peptide(
                        sequence,
                        protein2index,
                        0
                    )
            peptides[sequence].identification |= (1 << ((N - 1) - index))
    return [peptide for seq, peptide in peptides.items()]


def get_coverage(fasta_records: Dict[str, str], peptides: List[Peptide], identification: int):
    coverage = {name: [False for _ in range(len(sequence))] for name, sequence in fasta_records.items()}
    for peptide in peptides:
        if (identification & peptide.identification) == 0:
            continue
        for protein, indexes in peptide.protein2index.items():
            for index in indexes:
                for i in range(index, index + len(peptide.sequence)):
                    coverage[protein][i] = True
    results = {}
    for name, cov in coverage.items():
        if sum(cov) == 0:
            continue
        results[name] = (sum(cov) / len(cov))
    return results


def get_test(data0, data1):
    data = []
    for name, value in data0.items():
        if name not in data1:
            continue
        data.append(data0[name] - data1[name])
    res = wilcoxon(data)
    return res.statistic, res.pvalue


def plot(peptides, fasta_records):
    data = []
    datas = []
    for identification in identifications:
        data0 = get_coverage(fasta_records, peptides, identification)
        data.append([value for name, value in data0.items()])
        datas.append(data0)
        print(statistics.median(data[-1]))
    # print(get_test(datas[1], datas[3]))
    # print(get_test(datas[1], datas[4]))
    # print(get_test(datas[1], datas[5]))
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(7, 5), squeeze=True)
    axs.set_xticklabels(identifications)
    axs.boxplot(data, showfliers=False)
    plt.show()


if __name__ == "__main__":
    evidence_file = "D:\\data\\coon\\yeast_replicates_revision\\evidence.txt"
    fasta_file = "D:\\data\\coon\\yeast_replicates_revision\\559292_SaccharomycesCerevisiae_Uniprot_ISOFORMS_TARGET_20180228.fasta"

    fasta_records = read_fasta(fasta_file)
    peptides = read_evidences(evidence_file, fasta_records)
    plot(peptides, fasta_records)
