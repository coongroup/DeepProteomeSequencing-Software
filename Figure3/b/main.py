import json

from typing import List, Dict, Set, Tuple

import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


class GMutation:
    def __init__(
            self,
            identity: str,
            chr: str,
            nt_before: str,
            nt_pos: int,
            nt_after: str):
        self.identity: str = identity
        self.chr = chr
        self.nt_before: str = nt_before
        self.nt_pos: int = nt_pos
        self.nt_after: str = nt_after

    def __str__(self):
        return f"{self.chr}:{self.nt_before}:{self.nt_pos}:{self.nt_after}"


class PMutation:
    def __init__(
            self,
            identity: str,
            aa_before: str,
            aa_pos: int,
            aa_after: str):
        self.identity: str = identity
        self.aa_before: str = aa_before
        self.aa_pos: int = aa_pos
        self.aa_after: str = aa_after

    def __str__(self):
        return f"{self.aa_before}{self.aa_pos}{self.aa_after}"


def PMutation_convert(line: str) -> PMutation:
    m = line.split(':')
    if len(m) != 2:
        print(line)
    return PMutation(
        m[0],
        m[1][0],
        int(m[1][1:-1]) - 1,  # TODO zero-based
        m[1][-1]
    )


class Mutation:
    def __init__(
            self,
            gmutation: GMutation,
            pmutation: PMutation):
        self.gmutation = gmutation
        self.pmutation = pmutation


class FastaRecord:
    def __init__(
            self,
            name: str,
            mutations: List[PMutation],
            sequence: str):
        self.name: str = name
        self.mutations: List[PMutation] = mutations
        self.sequence: str = sequence


class FastaRecords:
    @staticmethod
    def read_fasta_file(fasta_file) -> List[FastaRecord]:
        fasta_records: List[FastaRecord] = []
        with open(fasta_file, 'r') as fs:
            name = ''
            mutations: List[PMutation] = []
            sequence = ''
            for line in fs:
                l = line.rstrip()
                if l.startswith('>'):
                    if name != '':
                        fasta_records.append(
                            FastaRecord(
                                name,
                                mutations,
                                sequence
                            )
                        )
                    spl = l.split()
                    name = spl[0][1:]
                    sequence = ""
                    mutations: List[PMutation] = []
                    if len(spl) > 1 and spl[1] != "":
                        for mutation in spl[1].split(';'):
                            mutations.append(
                                PMutation_convert(mutation)
                            )
                else:
                    sequence += l
            if name != '':
                fasta_records.append(
                    FastaRecord(
                        name,
                        mutations,
                        sequence
                    )
                )
        return fasta_records


def read_variant(filename: str) -> List[GMutation]:
    results: List[GMutation] = []
    with open(filename) as fs:
        header = fs.readline().rstrip().split('\t')
        id_index = header.index("Id")
        chr_index = header.index("Chr")
        ref_index = header.index("Ref")
        pos_index = header.index("Pos")
        alt_index = header.index("Alt")
        for line in fs:
            spl = line.rstrip().split('\t')
            results.append(
                GMutation(
                    spl[id_index],
                    spl[chr_index],
                    spl[ref_index],
                    int(spl[pos_index]),  # TODO should we -1
                    spl[alt_index]
                )
            )
    return results


class Annotation:
    def __init__(
            self,
            gene_name: str,
            transcript_name: str,
            protein_name: str,
            uniprot_name: str,
            transmembrane: bool):
        self.gene_name: str = gene_name
        self.transcript_name: str = transcript_name
        self.protein_name: str = protein_name
        self.uniprot_name: str = uniprot_name
        self.transmembrane: bool = transmembrane


def read_annotation(filename: str) -> List[Annotation]:
    results: List[Annotation] = []
    with open(filename) as fs:
        header = fs.readline().rstrip().split('\t')
        gene_index = header.index("Gene stable ID")
        transcript_index = header.index("Transcript stable ID")
        protein_index = header.index("Protein stable ID")
        uniprot_index = header.index("UniProtKB/Swiss-Prot ID")
        transmembrane_index = header.index("Transmembrane helices")
        for line in fs:
            spl = line.rstrip().split('\t')
            results.append(
                Annotation(
                    spl[gene_index],
                    spl[transcript_index],
                    spl[protein_index],
                    spl[uniprot_index],
                    spl[transmembrane_index] == "TMhelix"
                )
            )
    return results


class IdentifiedProteinMutation:
    def __init__(
            self,
            protein_name: str,
            mutations: List[PMutation]):
        self.protein_name = protein_name
        self.mutations = mutations


def read_peptide(filename: str) -> List[IdentifiedProteinMutation]:
    results: Dict[str, Dict[str, PMutation]] = {}
    with open(filename) as fs:
        spl = fs.readline().rstrip().split('\t')
        mutated_index = spl.index('Mutated')
        mutation_names_index = spl.index('Mutation names')
        proteins_index = spl.index('Proteins')
        reverse_index = spl.index('Reverse')
        potential_contaminant_index = spl.index('Potential contaminant')
        for line in fs:
            spl = line.rstrip().split('\t')
            if spl[reverse_index] == '+' or \
                    spl[potential_contaminant_index] == '+' or \
                    spl[mutated_index] != 'Yes' or \
                    spl[mutation_names_index] == "":
                continue
            mutations = set(spl[mutation_names_index].split(';'))
            if "ENSP" not in spl[proteins_index]:
                continue
            for x in spl[proteins_index].split('ENSP')[1:]:
                x0 = "ENSP" + x.rstrip(';')
                magic_number = len("ENSP00000345156")
                p = x0[:magic_number]
                ms = x0[magic_number:]
                if ms == "":
                    continue
                for m in ms.split(';'):
                    pmutation = PMutation_convert(m)
                    if pmutation.identity in mutations:
                        if p not in results:
                            results[p] = {}
                        if pmutation.identity not in results[p]:
                            results[p][pmutation.identity] = pmutation
    return [IdentifiedProteinMutation(protein_name, [mutation for _, mutation in id2mutation.items()])
            for protein_name, id2mutation in results.items()]


class ResultRecord:
    def __init__(
            self,
            gmutation: GMutation,
            detected: List[bool],
            transmembrane: bool,
            ensp2pmutation: Dict[str, PMutation]):
        self.gmutation = gmutation
        self.detected = detected
        self.transmembrane = transmembrane
        self.ensp2pmutation = ensp2pmutation


class ResultRecords:
    def __init__(
            self,
            records: Dict[str, ResultRecord],
            ncell_lines: int):
        self.records = records
        self.ncell_lines = ncell_lines

    def write_vcf(self, filename: str):
        with open(filename, "w") as fs:
            fs.write("\t".join(["#CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]) + "\n")

            for name, record in sorted(self.records.items(),
                                       key=lambda item: (item[1].gmutation.chr, item[1].gmutation.nt_pos)):
                if record.transmembrane:
                    continue
                fs.write("\t".join([record.gmutation.chr, f"{record.gmutation.nt_pos}", ".", record.gmutation.nt_before,
                                    record.gmutation.nt_after, ".", ".", f"{sum(record.detected)}"]) + "\n")

    @staticmethod
    def generate_translator(
            origin_proteins_records: List[FastaRecord],
            merge_proteins_records: List[FastaRecord]) -> Dict[str, str]:
        result: Dict[str, str] = {}
        merge_translator: Dict[str, str] = {}
        for record in merge_proteins_records:
            for mutation in record.mutations:
                key = f"{record.name}:{str(mutation)}"
                merge_translator[key] = mutation.identity
        for record in origin_proteins_records:
            for mutation in record.mutations:
                key = f"{record.name}:{str(mutation)}"
                if key in merge_translator:
                    result[mutation.identity] = merge_translator[key]
        return result

    def update(
            self,
            ind,
            proteins_records: List[FastaRecord],
            variants_records: List[GMutation],
            identified_variants: List[IdentifiedProteinMutation],
            proteins_merge_records: List[FastaRecord],
            annotations: List[Annotation]):
        ensp2annotation: Dict[str, Annotation] = {annotation.protein_name: annotation for annotation in annotations}
        gmutations: Dict[str, GMutation] = {variants_record.identity: variants_record for variants_record in
                                            variants_records}
        translator: Dict[str, str] = ResultRecords.generate_translator(proteins_records,
                                                                       proteins_merge_records)
        mutation2identified_proteins: Dict[str, Set[str]] = {}
        for identified_variant in identified_variants:
            for mutation in identified_variant.mutations:
                if mutation.identity not in mutation2identified_proteins:
                    mutation2identified_proteins[mutation.identity] = set()
                mutation2identified_proteins[mutation.identity].add(identified_variant.protein_name)
        for proteins_record in proteins_records:
            for pmutation in proteins_record.mutations:
                gmutation = gmutations[pmutation.identity]
                key = str(gmutation)
                if key not in self.records:
                    self.records[key] = ResultRecord(
                        gmutation,
                        [False for _ in range(self.ncell_lines)],
                        False,
                        {}
                    )
                if translator[pmutation.identity] in mutation2identified_proteins:
                    self.records[key].detected[ind] = True
                if proteins_record.name in ensp2annotation and ensp2annotation[proteins_record.name].transmembrane:
                    self.records[key].transmembrane = True
                self.records[key].ensp2pmutation[proteins_record.name] = pmutation


def prepare_input_ensemble_vep(parameters: Dict):
    annotation: List[Annotation] = read_annotation(parameters["input"]["transmembrane_annotation"])
    proteins_merge_records: List[FastaRecord] = FastaRecords.read_fasta_file(parameters["input"]["proteins"]["merge"])
    cell_lines = list(parameters["input"]["proteins"]["cell_lines"])
    resultsRecord: ResultRecords = ResultRecords({}, len(cell_lines))
    for i in range(len(cell_lines)):
        cell_line = cell_lines[i]
        print(cell_line)
        proteins_filename = parameters["input"]["proteins"]["cell_lines"][cell_line]
        variants_filename = parameters["input"]["variants"]["cell_lines"][cell_line]
        peptides_filename = parameters["input"]["peptides"]["cell_lines"][cell_line]
        proteins_records: List[FastaRecord] = FastaRecords.read_fasta_file(proteins_filename)
        variants_records: List[GMutation] = read_variant(variants_filename)
        identified_variants: List[IdentifiedProteinMutation] = read_peptide(peptides_filename)
        resultsRecord.update(i, proteins_records, variants_records, identified_variants, proteins_merge_records,
                             annotation)
    resultsRecord.write_vcf(parameters["output"]["input.vcf"])


class InputInfo:
    def __init__(
            self,
            mutation: GMutation,
            isDetected: bool):
        self.mutation: GMutation = mutation
        self.isDetected = isDetected


class PredictorResult:
    def __init__(
            self,
            name: str,
            score: float):
        self.name = name
        self.score = score


class OutputInfo:
    def __init__(
            self,
            input_info: InputInfo,
            sift: PredictorResult,
            polyPhen: PredictorResult):
        self.input_info: InputInfo = input_info
        self.sift = sift
        self.polyPhen = polyPhen


def parse_output_ensemble_vep(parameters: Dict):
    input_info: Dict[str, InputInfo] = {}
    with open(parameters["output"]["input.vcf"]) as fs:
        line = fs.readline()
        for line in fs:
            spl = line.rstrip().split('\t')
            mutation = GMutation(
                "",
                spl[0],
                spl[3],
                int(spl[1]),
                spl[4]
            )
            info = InputInfo(
                mutation,
                int(spl[7]) > 0
            )
            input_info[str(mutation)] = info
    output_info: Dict[str, OutputInfo] = {}
    with open(parameters["output"]["output.txt"]) as fs:
        line = fs.readline()
        while line.startswith("##"):
            line = fs.readline()
        for line in fs:
            spl = line.rstrip().split('\t')
            if spl[6] != "missense_variant":
                continue
            spl0 = spl[0].split("_")
            mutation = GMutation(
                "",
                spl0[0],
                spl0[2][0],
                int(spl0[1]),
                spl0[2][-1]
            )
            key = str(mutation)
            if key not in output_info:
                output_info[key] = OutputInfo(
                    input_info[str(mutation)],
                    PredictorResult("unknown", 100),
                    PredictorResult("unknown", -100)
                )
            sift = None
            polyPhen = None
            for info in spl[-1].split(';'):
                spl0 = info.split("=")
                if spl0[0] == "SIFT":
                    i = spl0[1].index('(')
                    sift = PredictorResult(spl0[1][:i], float(spl0[1][(i + 1):-1]))
                if spl0[0] == "PolyPhen":
                    i = spl0[1].index('(')
                    polyPhen = PredictorResult(spl0[1][:i], float(spl0[1][(i + 1):-1]))
            if sift is not None and sift.score < output_info[key].sift.score:
                output_info[key].sift = sift
            if polyPhen is not None and polyPhen.score > output_info[key].polyPhen.score:
                output_info[key].polyPhen = polyPhen
    # with open(parameters["output"]["final.txt"], "w") as fs:
    #     for key, info in output_info.items():
    #         fs.write("\t".join(
    #             [key, str(info.input_info.isDetected), info.sift.name, str(info.sift.score), info.polyPhen.name,
    #              str(info.polyPhen.score)]) + "\n")
    data_sift: Dict[str, List[int]] = {
        "deleterious": [0, 0],
        "deleterious_low_confidence": [0, 0],
        "tolerated_low_confidence": [0, 0],
        "tolerated": [0, 0]
    }
    data_sift_scores = [[], []]
    data_polyphen: Dict[str, List[int]] = {
        "probably_damaging": [0, 0],
        "possibly_damaging": [0, 0],
        "benign": [0, 0]
    }
    data_polyphen_scores = [[], []]
    for key, info in output_info.items():
        if info.sift.name in data_sift:
            data_sift[info.sift.name][int(info.input_info.isDetected)] += 1
            data_sift_scores[int(info.input_info.isDetected)].append(info.sift.score)
        if info.polyPhen.name in data_polyphen:
            data_polyphen[info.polyPhen.name][int(info.input_info.isDetected)] += 1
            data_polyphen_scores[int(info.input_info.isDetected)].append(info.polyPhen.score)
    print(len(data_sift_scores[0]))
    print(len(data_sift_scores[1]))
    print(len(data_polyphen_scores[0]))
    print(len(data_polyphen_scores[1]))
    plot(data_sift, mannwhitneyu(data_sift_scores[0], data_sift_scores[1], alternative='two-sided').pvalue)
    plot(data_polyphen, mannwhitneyu(data_polyphen_scores[0], data_polyphen_scores[1], alternative='two-sided').pvalue)


def plot(data: Dict[str, List[int]], pvalue: float):
    print(data.items())
    print(pvalue)
    items = list(data.items())
    data_bottom = [0, 0]
    data_sum = [sum([value[0] for name, value in items]), sum([value[1] for name, value in items])]
    data_norm = [[value[0] / data_sum[0], value[1] / data_sum[1]] for name, value in items]
    fig, ax = plt.subplots()
    for i in range(len(items)):
        ax.bar([0, 1], data_norm[i], bottom=data_bottom, width=0.8)
        data_bottom[0] += data_norm[i][0]
        data_bottom[1] += data_norm[i][1]
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["false", "true"])
    ax.set_title(f"{','.join([name for name, value in items])}:{pvalue}")
    plt.show()


if __name__ == "__main__":
    with open("parameters.json", 'r') as parameters_fs:
        parameters = json.loads(parameters_fs.read())
    # prepare_input_ensemble_vep(parameters)
    parse_output_ensemble_vep(parameters)
