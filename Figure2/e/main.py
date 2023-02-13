import json
import os
import re

from typing import Dict, List, Set, Tuple, Iterator

parameter_filename = "parameters.json"


class FastaRecord:
    def __init__(
            self,
            name: str,
            sequence: str):
        self.name = name
        self.sequence = sequence


class Peptide:
    def __init__(
            self,
            sequence: str):
        self.sequence = sequence


class PeptideMapped(Peptide):
    def __init__(
            self,
            sequence: str,
            targets: List[Tuple[FastaRecord, int]]):
        super().__init__(sequence)
        self.targets = targets


class Protease:
    def __init__(
            self,
            name: str,
            pattern: str,
            max_missed_cleavage: int):
        self.name = name
        self.pattern = pattern
        self.max_missed_cleavage = max_missed_cleavage


def read_fasta_records(
        file_name: str) \
        -> List[FastaRecord]:
    fasta_records: List[FastaRecord] = []
    fasta_record: FastaRecord = None
    with open(file_name) as file_stream:
        for line in file_stream:
            if line[0] == '>':
                if fasta_record is not None:
                    fasta_records.append(fasta_record)
                fasta_record = FastaRecord(line.split('|')[1], "")
            else:
                fasta_record.sequence += line.rstrip()
        if fasta_record is not None:
            fasta_records.append(fasta_record)
    return fasta_records


def digest(
        fasta_record: FastaRecord,
        protease: Protease,
        length_filter: Tuple[int, int]) \
        -> Iterator[Tuple[int, int]]:
    positions: List[int] = [0, len(fasta_record.sequence)]
    for match in re.finditer(protease.pattern, fasta_record.sequence):
        positions.append(match.start() + 1)
    positions.sort()
    for missed_cleavage in range(0, protease.max_missed_cleavage + 1):
        for i in range(missed_cleavage + 1, len(positions)):
            pos0 = positions[i - missed_cleavage - 1]
            pos1 = positions[i]
            if length_filter[0] <= pos1 - pos0 <= length_filter[1]:
                yield pos0, pos1


def get_peptides_mapped(
        fasta_records: List[FastaRecord],
        proteases: List[Protease],
        length_filter: Tuple[int, int]) \
        -> List[PeptideMapped]:
    peptides: Dict[str, PeptideMapped] = {}
    for protease in proteases:
        for fasta_record in fasta_records:
            for pos0, pos1 in digest(fasta_record, protease, length_filter):
                sequence = fasta_record.sequence[pos0: pos1]
                if sequence in peptides:
                    peptides[sequence].targets.append((fasta_record, pos0))
                else:
                    peptides[sequence] = PeptideMapped(
                        sequence,
                        [(fasta_record, pos0)]
                    )
    return list(peptides.values())


def get_statistics(mapped_peptides: List[PeptideMapped], fasta_records: List[FastaRecord]) -> Tuple[int, float]:
    peptide_count: int = len(mapped_peptides)
    coverage: Dict[str, List[int]] = {}
    for fasta_record in fasta_records:
        coverage[fasta_record.name] = [0 for i in range(len(fasta_record.sequence))]
    for mapped_peptide in mapped_peptides:
        sequence = mapped_peptide.sequence
        targets = mapped_peptide.targets
        for target in targets:
            fasta_record: FastaRecord = target[0]
            position: int = target[1]
            for i in range(position, position + len(sequence)):
                coverage[fasta_record.name][i] += 1
    cov_cnt = 0
    cnt = 0
    for name, value in coverage.items():
        cov_cnt += sum([v != 0 for v in value])
        cnt += len(value)
    print(cov_cnt, cnt)
    return peptide_count, (cov_cnt + 0.0) / cnt


def read_peptides(file_name: str, name_column: str, exclude_columns: List[str]) -> List[str]:
    result: Set[str] = set()
    with open(file_name) as file_stream:
        header = file_stream.readline().rstrip().split('\t')
        name_index = header.index(name_column)
        exclude_indexes = [header.index(exclude_column) for exclude_column in exclude_columns]
        for line in file_stream:
            spl = line.rstrip().split('\t')
            if sum([spl[i] == '+' for i in exclude_indexes]) != 0 or spl[name_index] == "":
                continue
            for name in spl[name_index].split(';'):
                result.add(name)
    return list(result)


def run(
        params: Dict) \
        -> None:
    fasta_records: List[FastaRecord] = read_fasta_records(params["input"]["fasta"]["filename"])
    proteases = {
        protease["name"]: Protease(
            protease["name"],
            protease["pattern"],
            protease["max_missed_cleavage"]
        )
        for protease in params["input"]["proteases"]
    }
    length_filter: Tuple[int, int] = \
        (params["input"]["filter"]["length"]["min"], params["input"]["filter"]["length"]["max"])
    print(get_statistics(get_peptides_mapped(fasta_records, [proteases["Trypsin"]], length_filter), fasta_records))
    print(get_statistics(get_peptides_mapped(fasta_records, list(proteases.values()), length_filter), fasta_records))
    identified_proteins: List[str] = read_peptides(
        params["input"]["peptides"]["filename"],
        params["input"]["peptides"]["columns"]["name"],
        params["input"]["peptides"]["columns"]["exclude"]
    )
    print(len(set(identified_proteins).intersection(set([fasta_record.name for fasta_record in fasta_records]))))
    print(len(fasta_records))
    print(len(proteases))
    pass


if __name__ == "__main__":
    with open(parameter_filename, 'r') as filestream:
        run(json.loads(filestream.read()))
