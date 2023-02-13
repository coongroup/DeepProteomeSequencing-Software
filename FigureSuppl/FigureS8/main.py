import json

import matplotlib.pyplot as plt


class Peptide:
    def __init__(self, sequence,
                 peptideCoordinates,
                 nodeCoordinates,
                 protease2expression):
        self.sequence = sequence
        self.peptideCoordinates = peptideCoordinates
        self.nodeCoordinates = [nodeCoordinates[0]]
        for i in range(1, len(nodeCoordinates)):
            if self.nodeCoordinates[-1][1] == nodeCoordinates[i][0]:
                self.nodeCoordinates[-1][1] = nodeCoordinates[i][1]
            else:
                self.nodeCoordinates.append(nodeCoordinates[i])
        assert len(self.nodeCoordinates) == len(self.peptideCoordinates), sequence
        self.protease2expression = protease2expression

    @staticmethod
    def read(filename, proteases):
        result = {}
        with open(filename, "r") as fs:
            header = fs.readline().rstrip().split('\t')
            sequenceIndex = header.index("Sequence")
            peptideCoordinatesIndex = header.index("PeptideCoordinates")
            nodeCoordinatesIndex = header.index("NodeCoordinates")
            protease2index = {protease: header.index(protease) for protease in proteases}
            for line in fs:
                spl = line.rstrip().split('\t')
                result[spl[sequenceIndex]] = Peptide(
                    spl[sequenceIndex],
                    [[int(coordinate.split('|')[0]), int(coordinate.split('|')[1])]
                     for coordinate in spl[peptideCoordinatesIndex].split(';')],
                    [[int(coordinate.split('|')[0]), int(coordinate.split('|')[1])]
                     for coordinate in spl[nodeCoordinatesIndex].split(';')],
                    {protease: int(spl[index])
                     for protease, index in protease2index.items()}
                )
        return result


def plotSummaryJunctionPeptides(peptide2data, proteases, color, outputFile):
    size = 3
    protease2count = {protease: [0 for _ in range(size)] for protease in proteases}
    for peptide, data in peptide2data.items():
        for protease in proteases:
            if data.protease2expression[protease] != 0:
                protease2count[protease][min(size - 1, len(data.nodeCoordinates) - 1)] += 1
    protease2countNormalized = {protease: [protease2count[protease][i] / sum(protease2count[protease])
                                           for i in range(size)]
                                for protease in proteases}

    width = 0.5

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(5, 2*5), sharex=True, squeeze=True)

    inds = [i for i in range(len(proteases))]
    values = [0 for _ in range(len(proteases))]
    for i in range(len(color)):
        values0 = [protease2count[proteases[j]][i] for j in range(len(proteases))]
        axs[0].bar(inds, values0, width, bottom=values, color=color[i])
        values = [x0 + x1 for x0, x1 in zip(values, values0)]
    axs[0].set_ylabel("Peptides")

    values = [0 for _ in range(len(proteases))]
    for i in range(len(color)):
        values0 = [protease2countNormalized[proteases[j]][i] for j in range(len(proteases))]
        axs[1].bar(inds, values0, width, bottom=values, color=color[i])
        if i == 0:
            for j in range(len(proteases)):
                axs[1].text(j, values0[j] - 0.05, f"{values0[j]}"[:4], color=color[1],
                            verticalalignment='top', horizontalalignment='center')
        values = [x0 + x1 for x0, x1 in zip(values, values0)]
    axs[1].set_ylabel("Percent")
    axs[1].set_xticks(inds)
    proteaseLabels = []
    for i in range(len(proteases)):
        if len(proteases[i]) > 7:
            proteaseLabels.append(proteases[i][:5] + ".")
        else:
            proteaseLabels.append(proteases[i])
    axs[1].set_xticklabels(proteaseLabels)

    plt.savefig(outputFile)


def plotSummaryAroundExon(peptide2data, proteases, color, outputFile, size=10):
    nTerminusCount = {protease: [0 for _ in range(size)] for protease in proteases}
    cTerminusCount = {protease: [0 for _ in range(size)] for protease in proteases}
    nTerminusMsMsCount = {protease: [0 for _ in range(size)] for protease in proteases}
    cTerminusMsMsCount = {protease: [0 for _ in range(size)] for protease in proteases}
    for peptide, data in peptide2data.items():
        for protease in proteases:
            if data.protease2expression[protease] != 0:
                if data.peptideCoordinates[0][0] - data.nodeCoordinates[0][0] < size:
                    nTerminusCount[protease][data.peptideCoordinates[0][0] - data.nodeCoordinates[0][0]] += 1
                    nTerminusMsMsCount[protease][data.peptideCoordinates[0][0] - data.nodeCoordinates[0][0]] += \
                        data.protease2expression[protease]
                if data.nodeCoordinates[-1][1] - data.peptideCoordinates[-1][1] < size:
                    cTerminusCount[protease][data.nodeCoordinates[-1][1] - data.peptideCoordinates[-1][1]] += 1
                    cTerminusMsMsCount[protease][data.nodeCoordinates[-1][1] - data.peptideCoordinates[-1][1]] += \
                        data.protease2expression[protease]

    # data = [[nTerminusCount, cTerminusCount], [
    #     {protease: [nTerminusCount[protease][i] / sum(nTerminusCount[protease]) for i in range(size)]
    #      for protease in proteases},
    #     {protease: [cTerminusCount[protease][i] / sum(cTerminusCount[protease]) for i in range(size)]
    #      for protease in proteases}
    # ]]
    data = [nTerminusCount, cTerminusCount]
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(2*5, 5), squeeze=True)

    inds = [i for i in range(size)]

    for j in range(len(data)):
        for k in range(size):
            for p in range(len(proteases)):
                # print(inds[k]-0.25+p*0.1)
                # print(data[i][j][0][proteases[p]])
                # print(data[i][j][0][proteases[p]][k])
                if j == 0:
                    axs[j].bar(inds[k] - 0.3 + p * 0.12, data[j][proteases[p]][k], 0.11, color=color[p])
                else:
                    axs[j].yaxis.tick_right()
                    axs[j].bar(inds[k] - 0.3 + p * 0.12, data[j][proteases[p]][size - k - 1], 0.11,
                                  color=color[p])
                axs[j].set_ylim([0, 15000])
                axs[j].set_xticks(inds)
                axs[j].set_ylabel("Peptides")
                if j == 0:
                    axs[j].set_xticklabels([r for r in range(size)])
                    axs[j].set_xlabel("Distance to 5' of an exon, nt")
                else:
                    axs[j].set_xticklabels([size - r - 1 for r in range(size)])
                    axs[j].set_xlabel("Distance to 3' of an exon, nt")

    plt.savefig(outputFile)


def plotSummaryExonExon(peptide2data, proteases, color, outputFile, size=11):
    nTerminusCount = {protease: [0 for _ in range(size)] for protease in proteases}
    cTerminusCount = {protease: [0 for _ in range(size)] for protease in proteases}
    nTerminusMsMsCount = {protease: [0 for _ in range(size)] for protease in proteases}
    cTerminusMsMsCount = {protease: [0 for _ in range(size)] for protease in proteases}
    for peptide, data in peptide2data.items():
        if len(data.nodeCoordinates) == 1:
            continue
        for protease in proteases:
            if data.protease2expression[protease] != 0:
                if data.nodeCoordinates[0][1] - data.peptideCoordinates[0][0] < size:
                    nTerminusCount[protease][data.nodeCoordinates[0][1] - data.peptideCoordinates[0][0]] += 1
                    nTerminusMsMsCount[protease][data.nodeCoordinates[0][1] - data.peptideCoordinates[0][0]] += \
                        data.protease2expression[protease]
                if data.peptideCoordinates[-1][1] - data.nodeCoordinates[-1][0] < size:
                    cTerminusCount[protease][data.peptideCoordinates[-1][1] - data.nodeCoordinates[-1][0]] += 1
                    cTerminusMsMsCount[protease][data.peptideCoordinates[-1][1] - data.nodeCoordinates[-1][0]] += \
                        data.protease2expression[protease]
    data = [nTerminusCount, cTerminusCount]
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(2*5, 5), squeeze=True)
    inds = [i for i in range(size)]
    for j in range(len(data)):
        for k in range(size):
            for p in range(len(proteases)):
                axs[j].set_ylim([0, 7500])
                if j == 0:
                    axs[0].bar(inds[k] - 0.3 + p * 0.12, data[j][proteases[p]][k], 0.11, color=color[p])
                    axs[0].set_ylabel("Peptides")
                    axs[0].set_xlabel("Number of nucleotides at 5' donor exon, nt")
                else:
                    axs[1].bar(inds[k] - 0.3 + p * 0.12, data[j][proteases[p]][k], 0.11, color=color[p])
                    axs[1].set_xlabel("Number of nucleotides at 3' acceptor exon, nt")
                axs[j].set_xticklabels([i for i in range(size)])
                axs[j].set_xticks([i for i in range(size)])

    plt.savefig(outputFile)


def plot(parameters):
    peptide2data = Peptide.read(parameters["input"]["peptideFile"],
                                parameters["input"]["proteases"])
    plotSummaryJunctionPeptides(peptide2data, parameters["input"]["proteases"],
                                parameters["input"]["summaryJunctionPeptides"]["color"],
                                parameters["output"]["summaryJunctionPeptides"])
    plotSummaryAroundExon(peptide2data, parameters["input"]["proteases"], parameters["input"]["color"],
                          parameters["output"]["summaryAroundExon"])
    plotSummaryExonExon(peptide2data, parameters["input"]["proteases"], parameters["input"]["color"],
                        parameters["output"]["summaryExonExon"])


if __name__ == "__main__":
    with open("parameters.json", 'r') as parameters_fs:
        plot(json.loads(parameters_fs.read()))
