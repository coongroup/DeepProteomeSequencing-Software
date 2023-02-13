import math
import re

import matplotlib.pyplot as plt

proteaseNames = ["Trypsin", "LysC", "LysN", "GluC", "AspN", "Chymotrypsin"]
proteaseColors = ["#e35d4b", "#f5a42a", "#22968c", "#b7c634", "#81cde9", "#938fc5"]
allColor = "#c8c6c5"


class EdgeData:
    def __init__(self, records, expressionNames, expressionSum):
        self.records = records
        self.expressionNames = expressionNames
        self.expressionSum = expressionSum

    @staticmethod
    def parse(fileName, transcript2exons):
        records = {}
        expressionNames = None
        expressionSum = None
        position2transcripts = {}
        for transcript, exons in transcript2exons.items():
            for exon in exons:
                # if spl[6] == '+':
                #     result[transcriptName].append((int(spl[3]) - 1, int(spl[4]), spl[0]))
                # else:
                #     result[transcriptName].append((-int(spl[4]), -int(spl[3]) + 1, spl[0]))
                if exon[0] > 0:
                    for pos in [(exon[0], exon[2]), (exon[1], exon[2])]:
                        if pos not in position2transcripts:
                            position2transcripts[pos] = set()
                        position2transcripts[pos].add(transcript)
                else:
                    for pos in [(-exon[0], exon[2]), (-exon[1], exon[2])]:
                        if pos not in position2transcripts:
                            position2transcripts[pos] = set()
                        position2transcripts[pos].add(transcript)
        with open(fileName, "r") as fs:
            header = fs.readline().rstrip().split("\t")
            coordinateIndex = header.index("Coordinate")
            isConstitutiveIndex = header.index("IsConstitutive")
            skippedExonStatusIndex = header.index("SkippedExonStatus")
            expressionNames = [h for h in header if h.startswith("Caltech_")]
            expressionSum = [0 for h in header if h.startswith("Caltech_")]
            expressionIndexes = [header.index(expressionName) for expressionName in expressionNames]
            for line in fs:
                spl = line.rstrip().split("\t")
                coordinate = spl[coordinateIndex]
                expression = [int(spl[expressionIndex]) for expressionIndex in expressionIndexes]
                if coordinate not in records:
                    for i in range(len(expression)):
                        expressionSum[i] += expression[i]
                    chr = coordinate.split(':')[0]
                    start = int(coordinate.split(':')[1].split('|')[0])
                    end = int(coordinate.split(':')[1].split('|')[1])
                    if (start, chr) in position2transcripts and (end, chr) in position2transcripts:
                        transcripts = list(
                            position2transcripts[(start, chr)].intersection(position2transcripts[(end, chr)]))
                        records[coordinate] = Edge(
                            chr,
                            start,
                            end,
                            spl[isConstitutiveIndex] == '+',
                            spl[skippedExonStatusIndex],
                            expression,
                            transcripts
                        )
        return EdgeData(list(records.values()), expressionNames, expressionSum)

    def filter(self, robData):
        filterSet = set()
        for chromosome, posFrom, posTo in robData:
            filterSet.add(f"{chromosome}:{posFrom}")
            filterSet.add(f"{chromosome}:{posTo}")
        n = 0
        m = 0
        for record in self.records:
            if record.isConstitutive:
                trigger = False
                m += 1
                for d in [-1, 0, 1]:
                    posFrom = f"{record.chromosome}:{record.posFrom + d}"
                    posTo = f"{record.chromosome}:{record.posTo + d}"
                    if posFrom in filterSet or posTo in filterSet:
                        trigger = True
                if trigger:
                    record.isConstitutive = False
                    n += 1
        print(f"Filter: isConstitutive {n}/{m}")
        # n = 0
        # m = 0
        # for record in self.records:
        #     if record.skippedExonStatus == "Inclusion":
        #         trigger = False
        #         m += 1
        #         for d in [-1, 0, 1]:
        #             posFrom = f"{record.chromosome}:{record.posFrom + d}"
        #             posTo = f"{record.chromosome}:{record.posTo + d}"
        #             if posFrom in filterSet or posTo in filterSet:
        #                 trigger = True
        #         if not trigger:
        #             record.skippedExonStatus = None
        #             n += 1
        # print(f"Filter: skippedExonStatus {n}/{m}")
        pass


class Edge:
    def __init__(self, chromosome, posFrom, posTo, isConstitutive, skippedExonStatus, expression, transcripts):
        self.chromosome = chromosome
        self.posFrom = posFrom
        self.posTo = posTo
        self.isConstitutive = isConstitutive
        self.skippedExonStatus = skippedExonStatus
        self.expression = expression
        self.transcripts = transcripts

    def __str__(self):
        return f"{self.chromosome}:{self.posFrom}|{self.posTo}"


class PeptideData:
    def __init__(self, records, expressionNames):
        self.records = records
        self.expressionNames = expressionNames

    def getSplicedPeptides(self, selectedProteases):
        result = {}  # junction -> [Peptide]
        if selectedProteases == None:
            for record in self.records:
                if len(record.coordinates) > 1:
                    for i in range(1, len(record.coordinates)):
                        coordinate = f"{record.chromosome}:{record.coordinates[i - 1][1]}|{record.coordinates[i][0]}"
                        if coordinate not in result:
                            result[coordinate] = []
                        result[coordinate].append(record)
        else:
            indexes = []
            for i in range(len(self.expressionNames)):
                for protease in selectedProteases:
                    if protease in self.expressionNames[i]:
                        indexes.append(i)
                        continue
            for record in self.records:
                if len(record.coordinates) > 1 and sum([record.msmsCount[index] for index in indexes]) > 0:
                    for i in range(1, len(record.coordinates)):
                        coordinate = f"{record.chromosome}:{record.coordinates[i - 1][1]}|{record.coordinates[i][0]}"
                        if coordinate not in result:
                            result[coordinate] = []
                        result[coordinate].append(record)
        return result

    @staticmethod
    def parse(fileName):
        records = []
        expressionNames = None
        with open(fileName, "r") as fs:
            header = fs.readline().rstrip().split("\t")
            coordinateIndex = header.index("Coordinate")
            peptideSequenceIndex = header.index("PeptideSequence")
            intensityNames = [h for h in header if h.endswith("|Intensity")]
            intensityIndexes = [header.index(intensityName) for intensityName in intensityNames]

            msmsCountNames = [h for h in header if h.endswith("|MsMsCount")]
            msmsCountIndexes = [header.index(msmsCountName) for msmsCountName in msmsCountNames]

            expressionNames = [i.split("|")[0] for i in intensityNames]
            for line in fs:
                spl = line.rstrip().split("\t")
                coordinate = spl[coordinateIndex]
                records.append(
                    Peptide(
                        coordinate.split(':')[0],
                        [(int(i.split('|')[0]), int(i.split('|')[1])) for i in coordinate.split(':')[1].split(';')],
                        spl[peptideSequenceIndex],
                        [int(spl[i]) for i in intensityIndexes],
                        [int(spl[i]) for i in msmsCountIndexes]
                    )
                )
        return PeptideData(records, expressionNames)


class Peptide:
    def __init__(self, chromosome, coordinates, sequence, intensity, msmsCount):
        self.chromosome = chromosome
        self.coordinates = coordinates
        self.sequence = sequence
        self.intensity = intensity
        self.msmsCount = msmsCount


def getScoreArgLys(edge, transcript2exons, transcript2proteinSequence):
    score = 0
    if len(edge.transcripts) == 0:
        return -1
    transcript = edge.transcripts[0]
    if transcript not in transcript2exons:
        return -1
    exons = transcript2exons[transcript]
    if exons[0][0] > 0:
        posFrom = edge.posFrom
    else:
        posFrom = -edge.posTo
    p = 0
    for exon in transcript2exons[transcript]:
        p += exon[1] - exon[0] + 1
        if exon[1] == posFrom:
            break
    p = p // 3
    s = transcript2proteinSequence[transcript]
    for i in range(max(0, p - 1), min(p + 2, len(s))):
        if s[i] == "R" or s[i] == "K":
            score += 1
    return score


def plotProteases(edgeData, peptideData, transcript2proteinSequence, transcript2exons):
    fig, ax = plt.subplots()
    s = sum(edgeData.expressionSum)
    isConstitutive = True
    skippedExonStatus = "None"
    # for selectedProteases, colorProtease in [
    #     ([proteaseNames[0]], proteaseColors[0]),
    #     ([proteaseNames[1]], proteaseColors[1]),
    #     ([proteaseNames[2]], proteaseColors[2]),
    #     ([proteaseNames[3]], proteaseColors[3]),
    #     ([proteaseNames[4]], proteaseColors[4]),
    #     ([proteaseNames[5]], proteaseColors[5]),
    #     (proteaseNames, allColor)
    # ]:
    for protease in proteaseNames:
        print(protease)
        for selectedProteases, colorProtease, alpha, state in [
            ([protease], "grey", 1.0, "all"),
            ([protease], "red", 0.75, "zero"),
            ([protease], "blue", 0.5, "one")
        ]:
            splicedPeptides = peptideData.getSplicedPeptides(selectedProteases)
            xys = []
            selectedEdges = [edge for edge in edgeData.records
                             if edge.isConstitutive == isConstitutive and edge.skippedExonStatus == skippedExonStatus]
            for edge in selectedEdges:
                if sum(edge.expression) <= 10:
                    continue
                x = sum(edge.expression) * 10 ** 6 / s
                y = 0
                if str(edge) in splicedPeptides:
                    y = sum([sum(peptide.msmsCount) for peptide in splicedPeptides[str(edge)]])
                if state == "all":
                    xys.append((x, y))
                else:
                    scoreLysArg = getScoreArgLys(edge, transcript2exons, transcript2proteinSequence)
                    if (scoreLysArg == 0 and state == "zero") or (scoreLysArg > 0 and state == "one"):
                        xys.append((x, y))
            xys.sort()
            xs = []
            ys = []
            slidingBin = 500
            for i in range(slidingBin, len(xys)):
                xs.append(math.log2(xys[i - slidingBin // 2][0] + 1))
                ys.append(sum([xys[i - j][1] > 0 for j in range(slidingBin)]) / slidingBin)
            plt.plot(xs, ys, color=colorProtease, alpha=alpha, label=state)
            # plt.legend()
            plt.xlim(-0.3, 8.3)
            plt.ylim(-0.03, 0.63)
        plt.ylabel("proteomics detection rate, percent")
        plt.xlabel("log2[RPM + 1]")
        plt.legend()
        plt.show()


def plot(edgeData, peptideData):
    splicedPeptides = peptideData.getSplicedPeptides(None)
    fig, ax = plt.subplots()
    s = sum(edgeData.expressionSum)
    for isConstitutive, skippedExonStatus in [
        (True, "None"),
        (False, "Exclusion"),
        (False, "Inclusion")
    ]:
        xys = []
        selectedEdges = [edge for edge in edgeData.records
                         if edge.isConstitutive == isConstitutive and edge.skippedExonStatus == skippedExonStatus]
        # print(isConstitutive, skippedExonStatus, len(selectedEdges))
        for edge in selectedEdges:
            if sum(edge.expression) <= 10:
                continue
            x = sum(edge.expression) * 10 ** 6 / s
            y = 0
            if str(edge) in splicedPeptides:
                y = sum([sum(peptide.msmsCount) for peptide in splicedPeptides[str(edge)]])
            xys.append((x, y))
        xys.sort()
        print(isConstitutive, skippedExonStatus, len(selectedEdges), len(xys))
        xs = []
        ys = []
        slidingBin = 500
        for i in range(slidingBin, len(xys)):
            xs.append(math.log2(xys[i - slidingBin // 2][0] + 1))
            ys.append(sum([xys[i - j][1] > 0 for j in range(slidingBin)]) / slidingBin)
        plt.plot(xs, ys, label=skippedExonStatus)
        plt.xlim(-0.3, 8.3)
        plt.ylim(-0.03, 0.63)
    plt.ylabel("proteomics detection rate, percent")
    plt.xlabel("log2[RPM + 1]")
    plt.legend()
    plt.show()


def plot0(edgeData, peptideData):
    splicedPeptides = peptideData.getSplicedPeptides(None)
    selectedEdges = [edge for edge in edgeData.records
                     if (edge.isConstitutive == True and edge.skippedExonStatus == "None") or
                     (edge.isConstitutive == False and (edge.skippedExonStatus == "Exclusion" or edge.skippedExonStatus == "Inclusion"))]
    protease2indexes = [[] for _ in proteaseNames]
    for i in range(len(proteaseNames)):
        protease2indexes[i] = [j for j in range(len(peptideData.expressionNames))
                               if proteaseNames[i] in peptideData.expressionNames[j]]
    values = []
    n = 0
    s = sum(edgeData.expressionSum)
    for edge in selectedEdges:
        if sum(edge.expression) <= 15:
            continue
        x = sum(edge.expression) * 10 ** 6 / s
        if math.log2(x + 1) < 1:
            continue
        value = 0
        if str(edge) in splicedPeptides:
            for i in range(len(proteaseNames)):
                for splicedPeptide in splicedPeptides[str(edge)]:
                    for j in protease2indexes[i]:
                        if splicedPeptide.msmsCount[j] == 0:
                            continue
                        value |= (1 << i)
        values.append(value)
        n += 1
    print(f"total: {n}")
    cnts = []
    grouped_cnts = {i: [] for i in range(len(proteaseNames) + 1)}
    for i in range(1 << len(proteaseNames)):
        cnt = 0
        for value in values:
            if (i & value) != 0:
                cnt += 1
        # print(bin(i), cnt)
        grouped_cnts[bin(i)[2:].count("1")].append((round(cnt / n, 2), 69.0033 * ((cnt / n) / 0.4), bin(i)))
        cnts.append(cnt)
    print(proteaseNames)
    for i, j in grouped_cnts.items():
        print(i, sorted(j, reverse=True))




def readRobData(filename):
    records = []
    with open(filename, "r") as fs:
        for line in fs:
            spl = line.rstrip().split("\t")
            s0 = spl[1].split(":")
            chrName = s0[0][3:]
            s1 = s0[1].split('-')
            pos0 = int(s1[0])
            pos1 = int(s1[1])
            records.append((chrName, pos0, pos1))
    return records


class GtfAnnotation:
    @staticmethod
    def readFile(filename, nameRegularExpression, transcriptSet):
        result = {}
        with open(filename) as fileStream:
            for line in fileStream:
                if len(line) == 0 or (len(line) > 0 and line[0] == "#"):
                    continue
                spl = line.rstrip().split('\t')
                if spl[2] == "CDS":
                    transcriptName = re.findall(nameRegularExpression, line)[0]
                    if transcriptName not in transcriptSet:
                        continue
                    if transcriptName not in result:
                        result[transcriptName] = []
                    # result[transcriptName].append((int(spl[3]) - 1, int(spl[4]), spl[0]))
                    if spl[6] == '+':
                        result[transcriptName].append((int(spl[3]) - 1, int(spl[4]), spl[0]))
                    else:
                        result[transcriptName].append((-int(spl[4]), -int(spl[3]) + 1, spl[0]))
        return result


class ProteinFastaRecord:
    aminoacids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    @staticmethod
    def readFile(filename, nameRegularExpression):
        fastaRecords = {}
        with open(filename) as fileStream:
            name = ""
            for line in fileStream:
                if line[0] == '>':
                    name = re.findall(nameRegularExpression, line)[0]
                    fastaRecords[name] = ""
                else:
                    fastaRecords[name] += line.rstrip()
        aminoacidsSet = set(ProteinFastaRecord.aminoacids)
        return {name: sequence for name, sequence in fastaRecords.items()
                if sum([c in aminoacidsSet for c in sequence]) == len(sequence)}


if __name__ == "__main__":
    edgesFileName = "D:\\data\\coon\\isoform\\benRequest\\edges.mod3.txt"
    peptidesFileName = "D:\\data\\coon\\isoform\\benRequest\\peptides.mod3.txt"
    alternativeExons = "D:\\data\\coon\\isoform\\benRequest\\Alternative_AS.tab"

    transcript2proteinSequence = ProteinFastaRecord.readFile(
        "D:\\data\\coon\\reference\\Homo_sapiens.GRCh38.pep.all.fa",
        "transcript:([^.]*).[0-9]*"
    )
    transcript2exons = GtfAnnotation.readFile(
        "D:\\data\\coon\\reference\\gtf\\Homo_sapiens.GRCh38.94.gtf",
        "transcript_id \"(ENST[0-9]*)\"",
        set([transcript for transcript in transcript2proteinSequence.keys()])
    )

    edgeData = EdgeData.parse(edgesFileName, transcript2exons)
    print(f"Edge data: {len(edgeData.records)}")

    edgeData.filter(readRobData(alternativeExons))

    peptideData = PeptideData.parse(peptidesFileName)
    print(f"Peptide data: {len(peptideData.records)}")

    print(f"All splicing sites: {len(peptideData.getSplicedPeptides(None))}")
    # plotProteases(edgeData, peptideData, transcript2proteinSequence, transcript2exons)
    # plot(edgeData, peptideData)
    plot0(edgeData, peptideData)
