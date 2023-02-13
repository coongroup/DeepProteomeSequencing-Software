import statistics

import matplotlib.pyplot as plt

minValue = 0
maxValue = 10

cellLines = list(reversed(['HUVEC', 'GM12878', 'hESC', 'HeLaS3', 'HepG2', 'K562']))

geneExpressionFile = "D:\\data\\coon\\22_04_isoforms\\gene_expression.rpkm.protein_encoding.txt"
splicingFile = "D:\\data\\coon\\22_04_isoforms\\splicing.txt"


class OmicsData:
    def __init__(self, path1expression, path1coeff, path2expression, path2coeff):
        self.path1expression = path1expression
        self.path1coeff = path1coeff
        self.path2expression = path2expression
        self.path2coeff = path2coeff

    def getPsi(self):
        a = self.path1expression * self.path1coeff
        b = self.path1expression * self.path1coeff + self.path2expression * self.path2coeff
        if b == 0:
            return -1
        return a / b


class SpliceEvent:
    def __init__(self, id, eventType, geneId, geneName, proteomics, genomics, exonDistance):
        self.id = id
        self.eventType = eventType
        self.geneId = geneId
        self.geneName = geneName
        self.proteomics = proteomics
        self.genomics = genomics
        self.exonDistance = exonDistance

    @staticmethod
    def read(filename):
        # protein2coverage = {}
        # with open("D:\\data\\coon\\standard\\proteomics\\proteinGroups.txt") as fs:
        #     spl = fs.readline().rstrip().split('\t')
        #     proteinIndex = spl.index("Protein IDs")
        #     coverageIndex = spl.index("Sequence coverage [%]")
        #     for line in fs:
        #         spl = line.rstrip().split('\t')
        #         for protein in spl[proteinIndex].split(';'):
        #             if protein.startswith("ENSP"):
        #                 p = protein
        #                 if '.' in p:
        #                     p = p.split('.')[0]
        #                 protein2coverage[p] = float(spl[coverageIndex])
        # gene2proteins = {}
        # with open("D:\\data\\coon\\reference\\Homo_sapiens.GRCh38.pep.all.fa") as fs:
        #     for line in fs:
        #         if line.startswith('>'):
        #             spl = line.rstrip().split(' ')
        #             p = spl[0][1:]
        #             if '.' in p:
        #                 p = p.split('.')[0]
        #             g = ""
        #             for s in spl:
        #                 if s.startswith('gene:'):
        #                     g = s[len("gene:"):]
        #                     if '.' in g:
        #                         g = g.split('.')[0]
        #             if g != "":
        #                 if g not in gene2proteins:
        #                     gene2proteins[g] = []
        #                 gene2proteins[g].append(p)
        # gene2coverage = {}
        # for g, ps in gene2proteins.items():
        #     a = []
        #     for p in ps:
        #         if p in protein2coverage:
        #             a.append(protein2coverage[p])
        #     if len(a) != 0:
        #         gene2coverage[g] = statistics.mean(a)

        genesWithTM = set()
        with open("D:\\data\\coon\\22_04_isoforms\\mart_export.txt") as fs:
            line = fs.readline()
            for line in fs:
                spl = line.split(',')
                if spl[4] == "TMhelix":
                    genesWithTM.add(spl[0])
        cnt = 0
        spliceEvents = {}
        table = {
            "SkippedExon": [0, 0, 0, 0],
            "RetainedIntron": [0, 0, 0, 0],
            "AltDonor": [0, 0, 0, 0],
            "AltAcceptor": [0, 0, 0, 0],
            "Retained Intron": [0, 0, 0, 0],
            "MutuallyExclusiveExons": [0, 0, 0, 0]
        }
        with open(filename) as fs:
            header = fs.readline().rstrip().split('\t')
            omics2dataIndex = {}
            for omics in ["Proteomics", "Genomics"]:
                omics2dataIndex[omics] = {}
                for path in ["Path1", "Path2"]:
                    omics2dataIndex[omics][path] = []
                    for cellLine in cellLines:
                        for h in range(len(header)):
                            if cellLine in header[h] and omics in header[h] and path in header[h]:
                                omics2dataIndex[omics][path].append(h)
            geneIdIndex = header.index("Gene id")
            geneNameIndex = header.index("Gene name")
            positionsPath1Index = header.index("Positions Path1")
            positionsPath2Index = header.index("Positions Path2")
            eventTypeIndex = header.index("Event type")
            for line in fs:
                if line.startswith("#"):
                    continue
                spl = line.rstrip().split('\t')
                geneIdValue = spl[geneIdIndex]
                geneNameValue = spl[geneNameIndex]
                positionsPath1Value = spl[positionsPath1Index]
                positionsPath2Value = spl[positionsPath2Index]
                eventTypeValue = spl[eventTypeIndex]
                if eventTypeValue == "SkippedExon":
                    coeff1 = 2
                    coeff2 = 1
                    k = positionsPath2Value.split(';')[1:-1]
                    d = sum([int(k[i]) - int(k[i - 1]) for i in range(1, len(k), 2)])
                else:
                    continue
                # elif eventTypeValue in ["AltExonStart", "AltExonEnd"]:
                #     coeff1 = 1
                #     coeff2 = 1
                #     d = -1
                # elif eventTypeValue == "AltDonor":
                #     coeff1 = 1
                #     coeff2 = 2
                #     k1 = positionsPath1Value.split(';')
                #     k2 = positionsPath2Value.split(';')
                #     d = int(k2[1]) - int(k1[1])
                # elif eventTypeValue == "AltAcceptor":
                #     coeff1 = 2
                #     coeff2 = 1
                #     k1 = positionsPath1Value.split(';')
                #     k2 = positionsPath2Value.split(';')
                #     d = int(k2[1]) - int(k1[1])
                # elif eventTypeValue == "SkippedExon":
                #     coeff1 = 2
                #     coeff2 = 1
                #     k = positionsPath2Value.split(';')[1:-1]
                #     d = sum([int(k[i]) - int(k[i - 1]) for i in range(1, len(k), 2)])
                # else:
                #     continue
                if d % 3 != 0:
                    continue
                # if geneIdValue not in gene2coverage or gene2coverage[geneIdValue] < 0.8:
                #     cnt += 1
                #     continue
                if geneIdValue in genesWithTM:
                    continue
                event = SpliceEvent(
                    f"{positionsPath1Value}|{positionsPath2Value}",
                    eventTypeValue,
                    geneIdValue,
                    geneNameValue,
                    OmicsData(
                        sum([int(spl[i]) for i in omics2dataIndex["Proteomics"]["Path1"]]),
                        coeff1,
                        sum([int(spl[i]) for i in omics2dataIndex["Proteomics"]["Path2"]]),
                        coeff2
                    ),
                    OmicsData(
                        sum([int(spl[i]) for i in omics2dataIndex["Genomics"]["Path1"]]),
                        coeff1,
                        sum([int(spl[i]) for i in omics2dataIndex["Genomics"]["Path2"]]),
                        coeff2
                    ),
                    d
                )
                pPsi = event.proteomics.getPsi()
                tPsi = event.genomics.getPsi()
                # and :
                # if eventTypeValue in table:
                #     if event.genomics.path1expression + event.genomics.path2expression > 20:
                #         if 0.05 <= tPsi <= 0.95:
                #             table[eventTypeValue][0] += 1
                #             table[eventTypeValue][1] += 1
                #             if 0.0 < pPsi < 1.00:
                #                 table[eventTypeValue][2] += 1
                #                 table[eventTypeValue][3] += 1
                #             if pPsi == 0.0 or pPsi == 1.00:
                #                 table[eventTypeValue][2] += 1
                #         else:
                #             table[eventTypeValue][0] += 1
                if 0.05 <= tPsi <= 0.95 and event.genomics.path1expression + event.genomics.path2expression > 20:
                    # elif eventTypeValue == "AltDonor":
                    #     cntsAltDonor[d % 3] += 1
                    # elif eventTypeValue == "AltAcceptor":
                    #     cntsAltAcceptor[d % 3] += 1
                    spliceEvents[event.id] = event
        # print(cntsAltDonor)
        # print(cntsAltAcceptor)
        print(cnt)
        print(table)
        # exit(0)
        return spliceEvents


def readGene2expression(filename):
    gene2expression = {}
    with open(filename) as fs:
        header = fs.readline().rstrip().split('\t')
        cellLine2index = {cellLine: [h for h in range(len(header)) if cellLine in header[h]] for cellLine in cellLines}
        geneIdIndex = header.index("Gene id")
        for line in fs:
            if line.startswith("#"):
                continue
            spl = line.rstrip().split('\t')
            expression = []
            for cellLine in cellLines:
                for i in cellLine2index[cellLine]:
                    if spl[i] != "NaN":
                        expression.append(float(spl[i]))
            if len(expression) != 0:
                gene2expression[spl[geneIdIndex]] = sum(expression) / len(expression)
    return gene2expression


def plotGeneExpressionBins(gene2expression, spliceEvents):
    # psi2psiPlotExpressionQuantile(spliceEvents, gene2expression)
    # genesBins = [set() for i in range(minValue, maxValue)]
    proteomicsSplicedGenesBins = [set() for _ in range(minValue, maxValue)]
    proteomicsHalfSplicedGenesBins = [set() for _ in range(minValue, maxValue)]
    genomicsSplicedGenesBins = [set() for _ in range(minValue, maxValue)]
    # for geneId, expression in gene2expression.items():
    #     if expression <= minValue or maxValue <= expression:
    #         continue
    #     genesBins[int(expression)].add(geneId)
    cntsSkippedExon = [0, 0, 0]
    cntsSkippedExon0 = [0, 0, 0]
    cntsSkippedExon1 = [0, 0, 0]

    for id, spliceEvent in spliceEvents.items():
        if spliceEvent.geneId not in gene2expression:
            continue
        expression = gene2expression[spliceEvent.geneId]
        if expression <= minValue or maxValue <= expression:
            continue
        # if spliceEvent.proteomics.path1expression > 0 and spliceEvent.proteomics.path2expression > 0:
        #     proteomicsSplicedGenesBins[int(expression)].add(spliceEvent.geneId)
        # if spliceEvent.genomics.path1expression > 20 and spliceEvent.genomics.path2expression > 20:
        #     genomicsSplicedGenesBins[int(expression)].add(spliceEvent.geneId)
        tPsi = spliceEvent.genomics.getPsi()
        pPsi = spliceEvent.proteomics.getPsi()
        if tPsi < 0.05 or tPsi > 0.95 or spliceEvent.genomics.path1expression + spliceEvent.genomics.path2expression <= 20 or spliceEvent.eventType != "SkippedExon":
            continue
        if 0.0 < pPsi < 1.00:
            proteomicsSplicedGenesBins[int(expression)].add((spliceEvent.geneId, spliceEvent.exonDistance))
            cntsSkippedExon0[spliceEvent.exonDistance % 3] += 1
        if pPsi == 0.0 or pPsi == 1.0:
            proteomicsHalfSplicedGenesBins[int(expression)].add((spliceEvent.geneId, spliceEvent.exonDistance))
            cntsSkippedExon1[spliceEvent.exonDistance % 3] += 1
        # if 0.05 <= tPsi <= 0.95 and \
        #         spliceEvent.genomics.path1expression + spliceEvent.genomics.path2expression > 20:
        genomicsSplicedGenesBins[int(expression)].add((spliceEvent.geneId, spliceEvent.exonDistance))
        cntsSkippedExon[spliceEvent.exonDistance % 3] += 1
        #
        # if event.eventType == "SkippedExon":
        #     cntsSkippedExon[event.d % 3] += 1
        # if event.eventType == "SkippedExon" and 0.0 < pPsi < 1.0:
        #     cntsSkippedExon0[d % 3] += 1
        # if event.eventType == "SkippedExon" and (pPsi == 0.00 or pPsi == 1.0):
        #     cntsSkippedExon1[d % 3] += 1
    print(cntsSkippedExon)
    print(cntsSkippedExon0)
    print(cntsSkippedExon1)

    rs = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 10]]
    rsd = {i: 0 for i in range(minValue, maxValue)}
    for ri in range(len(rs)):
        for i in range(rs[ri][0], rs[ri][1]):
            rsd[i] += ri

    fig, ax = plt.subplots()
    width = 1
    labels = [str(r) for r in rs]
    positionsShift = [i + 0.5 for i in range(len(rs))]

    pdata = [0 for _ in range(len(rs))]
    p0data = [0 for _ in range(len(rs))]
    gdata = [0 for _ in range(len(rs))]
    for i in range(minValue, maxValue):
        pdata[rsd[i]] += len(proteomicsSplicedGenesBins[i])
        p0data[rsd[i]] += len(proteomicsHalfSplicedGenesBins[i])
        gdata[rsd[i]] += len(genomicsSplicedGenesBins[i])

    print(pdata)
    print(gdata)
    data = [pdata[i] / (gdata[i] + 1) for i in range(len(rs))]
    ax.bar(positionsShift, data, width=width, color="black")
    data0 = [p0data[i] / gdata[i] for i in range(len(rs))]
    ax.bar(positionsShift, data0, width=width, bottom=data, color="red")
    data1 = [data[i] + data0[i] for i in range(len(rs))]
    data0 = [1.0 - data1[i] for i in range(len(rs))]
    ax.bar(positionsShift, data0, width=width, bottom=data1, color="grey")
    print("data", data)
    print("data0", data0)
    print("data1", data1)

    # data0 = [len(genomicsSplicedGenesBins[i]) - data[i] for i in range(minValue, maxValue)]
    # ax.bar(positionsShift, data0, width=width, bottom=data, color="blue", label="Transcriptomics")

    ax.set_xticks(positionsShift)
    ax.set_xticklabels(labels)
    ax.set_ylabel("# events")
    ax.set_xlabel("Gene expression, RPKM [log2]")
    plt.show()

    ####################################################################################################################

    fig, ax = plt.subplots()
    data = [pdata[i] for i in range(len(rs))]
    ax.bar(positionsShift, data, width=width, color="black")
    data0 = [p0data[i] for i in range(len(rs))]
    ax.bar(positionsShift, data0, width=width, bottom=data, color="red")
    data1 = [data[i] + data0[i] for i in range(len(rs))]
    data0 = [gdata[i] - (data[i] + data0[i]) for i in range(len(rs))]
    ax.bar(positionsShift, data0, width=width, bottom=data1, color="grey")
    print("data", data)
    print("data0", data0)
    print("data1", data1)
    #
    # data0 = [len(genomicsSplicedGenesBins[i]) - data[i] for i in range(minValue, maxValue)]
    # ax.bar(positionsShift, data0, width=width, bottom=data, color="blue", label="Transcriptomics")
    ax.set_xticks(positionsShift)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Percent of fully identifcation events")
    ax.set_xlabel("Gene expression, RPKM [log2]")
    plt.show()


# def plotGeneExpressionSmooth(gene2expression, spliceEvents):
#     values = []
#     for id, spliceEvent in spliceEvents.items():
#         if spliceEvent.geneId not in gene2expression:
#             continue
#         expression = gene2expression[spliceEvent.geneId]
#         if expression <= minValue or maxValue <= expression:
#             continue
#         # if spliceEvent.proteomics.path1expression > 0 and spliceEvent.proteomics.path2expression > 0:
#         #     proteomicsSplicedGenesBins[int(expression)].add(spliceEvent.geneId)
#         # if spliceEvent.genomics.path1expression > 20 and spliceEvent.genomics.path2expression > 20:
#         #     genomicsSplicedGenesBins[int(expression)].add(spliceEvent.geneId)
#         if 0.05 <= spliceEvent.proteomics.getPsi() <= 0.95 and spliceEvent.exonDistance % 3 == 0:
#             values.append((expression, False))
#         if 0.05 <= spliceEvent.genomics.getPsi() <= 0.95 and \
#                 spliceEvent.genomics.path1expression + spliceEvent.genomics.path2expression > 20 and \
#                 spliceEvent.exonDistance % 3 == 0:
#             values.append((expression, True))
#     xs = []
#     ys = []
#     values.sort()
#     n = 500
#     m = 1
#     n0 = 0
#     i0 = 0
#     i1 = 0
#     while i1 != len(values):
#         while n0 != n and i1 != len(values):
#             if values[i1][1]:
#                 n0 += 1
#             i1 += 1
#         if i1 == len(values):
#             break
#         xs.append(values[i0][0])
#         ys.append((i1 - i0 + 1 - n) / (i1 - i0 + 1))
#         m0 = 0
#         while m0 != m and i0 != i1:
#             if values[i0][1]:
#                 n0 -= 1
#                 m0 += 1
#             i0 += 1
#
#     fig, ax = plt.subplots()
#     ax.plot(xs, ys, color="black")
#     ax.set_title(f"Width {n}; Step {m}")
#     ax.set_ylabel("Percent of fully identifcation events")
#     ax.set_xlabel("Gene expression, RPKM [log2]")
#     plt.show()


gene2expression = readGene2expression(geneExpressionFile)
spliceEvents = SpliceEvent.read(splicingFile)
plotGeneExpressionBins(gene2expression, spliceEvents)
# plotGeneExpressionSmooth(gene2expression, spliceEvents)

## task
# fig, ax = plt.subplots()
#
# z = [len([k for k, l in i if l % 3 != 0]) for i in proteomicsSplicedGenesBins]
# xs = [0] + [i + 1 for i in range(minValue, maxValue)]
# ys = [0] + [sum(z[:(i+1)]) / sum(z) for i in range(len(z))]
# ax.plot(xs, ys, color="pink", label="Proteomics12")
#
# z = [len([k for k, l in i if l % 3 == 0]) for i in proteomicsSplicedGenesBins]
# xs = [0] + [i + 1 for i in range(minValue, maxValue)]
# ys = [0] + [sum(z[:(i+1)]) / sum(z) for i in range(len(z))]
# ax.plot(xs, ys, color="red", label="Proteomics0")
#
# z = [len([k for k, l in i if l % 3 != 0]) for i in genomicsSplicedGenesBins]
# xs = [0] + [i + 1 for i in range(minValue, maxValue)]
# ys = [0] + [sum(z[:(i+1)]) / sum(z) for i in range(len(z))]
# ax.plot(xs, ys, color="lightblue", label="Transcriptomics12")
#
# z = [len([k for k, l in i if l % 3 == 0]) for i in genomicsSplicedGenesBins]
# xs = [0] + [i + 1 for i in range(minValue, maxValue)]
# ys = [0] + [sum(z[:(i+1)]) / sum(z) for i in range(len(z))]
# ax.plot(xs, ys, color="blue", label="Transcriptomics0")
#
# z = [len(i) for i in genesBins]
# xs = [0] + [i + 1 for i in range(minValue, maxValue)]
# ys = [0] + [sum(z[:(i+1)]) / sum(z) for i in range(len(z))]
# ax.plot(xs, ys, color="green", label="Genes wo AS")
#
# ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', borderaxespad=0., mode="expand", ncol=3)
# ax.set_ylabel("CDF")
# ax.set_xlabel("Gene expression, RPKM [log2]")
#
# plt.show()


# fig, ax = plt.subplots()
# width = 0.7
# labels = [str(i) for i in range(minValue, maxValue)]
# positions = [i + 0.5 for i in range(minValue, maxValue)]
# data = [len(proteomicsSplicedGenesBins[i]) / len(genomicsSplicedGenesBins[i]) for i in range(minValue, maxValue)]
# ax.bar(positions, data, width=width, color="red")
#
# ax.set_xticks(positions)
# ax.set_xticklabels(labels)
#
# ax.set_ylabel("Percent")
# ax.set_xlabel("Gene expression, RPKM [log2]")
# plt.show()

#
# def psi2psiPlotMod3(spliceEvents):
#     fig, ax = plt.subplots()
#     positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
#     xs = []
#     ys = []
#     colors = []
#     color2label = {
#         "red": "Insertion size % 3 != 0",
#         "grey": "Insertion size % 3 == 0"
#     }
#
#     for id, spliceEvent in spliceEvents.items():
#         if spliceEvent.geneId not in gene2expression:
#             continue
#         expression = gene2expression[spliceEvent.geneId]
#         if expression <= minValue:
#             continue
#         pPsi = spliceEvent.proteomics.getPsi()
#         tPsi = spliceEvent.genomics.getPsi()
#         if 0.05 <= pPsi <= 0.95 and 0.05 <= tPsi <= 0.95 and spliceEvent.genomics.path1expression + spliceEvent.genomics.path2expression > 20:
#             xs.append(pPsi)
#             ys.append(tPsi)
#             if spliceEvent.exonDistance % 3 != 0:
#                 colors.append("red")
#             else:
#                 colors.append("grey")
#
#     for color in set(colors):
#         x0 = [x for x, c in zip(xs, colors) if c == color]
#         y0 = [y for y, c in zip(ys, colors) if c == color]
#         ax.scatter(x0, y0, color=color, alpha=0.6, label=f"{color2label[color]} ({len(x0)})")
#
#     ax.axline((0, 0), (1, 1), color='black')
#     ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', borderaxespad=0., mode="expand", ncol=3)
#     ax.set_xticks(positions)
#     ax.set_yticks(positions)
#
#     ax.set_ylabel("PSI transriptomics")
#     ax.set_xlabel("PSI proteomics")
#     plt.show()
#
#
# def psi2psiPlotSkippedExons(spliceEvents):
#     fig, ax = plt.subplots()
#     positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
#     xs = []
#     ys = []
#     colors = []
#     color2label = {
#         "red": "SkippedExon",
#         "blue": "AltAcceptor",
#         "green": "AltDonor",
#         "grey": "rest"
#     }
#
#     for id, spliceEvent in spliceEvents.items():
#         if spliceEvent.geneId not in gene2expression:
#             continue
#         expression = gene2expression[spliceEvent.geneId]
#         if expression <= minValue:
#             continue
#         pPsi = spliceEvent.proteomics.getPsi()
#         tPsi = spliceEvent.genomics.getPsi()
#         if 0.05 <= pPsi <= 0.95 and 0.05 <= tPsi <= 0.95 and spliceEvent.genomics.path1expression + spliceEvent.genomics.path2expression > 20:
#             xs.append(pPsi)
#             ys.append(tPsi)
#             if spliceEvent.eventType == "SkippedExon":
#                 colors.append("red")
#             elif spliceEvent.eventType == "AltAcceptor":
#                 colors.append("blue")
#             elif spliceEvent.eventType == "AltDonor":
#                 colors.append("green")
#             else:
#                 colors.append("grey")
#
#     for color in set(colors):
#         x0 = [x for x, c in zip(xs, colors) if c == color]
#         y0 = [y for y, c in zip(ys, colors) if c == color]
#         ax.scatter(x0, y0, color=color, alpha=0.6, label=f"{color2label[color]} ({len(x0)})")
#
#     ax.axline((0, 0), (1, 1), color='black')
#     ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', borderaxespad=0., mode="expand", ncol=2)
#     ax.set_xticks(positions)
#     ax.set_yticks(positions)
#
#     ax.set_ylabel("PSI transriptomics")
#     ax.set_xlabel("PSI proteomics")
#     plt.show()
#
#
# def psi2psiPlotExpressionQuantile(spliceEvents, gene2expression):
#     fig, ax = plt.subplots()
#     positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
#     xs = []
#     ys = []
#     colors = []
#     color2label = {
#         "grey": "Q1",
#         "blue": "Q2",
#         "green": "Q3",
#         "red": "Q4"
#     }
#     es = []
#     for gene, expression in gene2expression.items():
#         if 0 < expression < 10:
#             es.append(expression)
#     es.sort()
#     for id, spliceEvent in spliceEvents.items():
#         if spliceEvent.geneId not in gene2expression:
#             continue
#         expression = gene2expression[spliceEvent.geneId]
#         if expression <= minValue:
#             continue
#         pPsi = spliceEvent.proteomics.getPsi()
#         tPsi = spliceEvent.genomics.getPsi()
#         if 0.05 <= pPsi <= 0.95 and 0.05 <= tPsi <= 0.95 and spliceEvent.genomics.path1expression + spliceEvent.genomics.path2expression > 20:
#             xs.append(pPsi)
#             ys.append(tPsi)
#             if spliceEvent.geneId not in gene2expression:
#                 continue
#             if es[0] <= gene2expression[spliceEvent.geneId] < es[int(len(es) * 0.25)]:
#                 colors.append("grey")
#             elif es[int(len(es) * 0.25)] <= gene2expression[spliceEvent.geneId] < es[int(len(es) * 0.50)]:
#                 colors.append("blue")
#             elif es[int(len(es) * 0.50)] <= gene2expression[spliceEvent.geneId] < es[int(len(es) * 0.75)]:
#                 colors.append("green")
#             else:
#                 colors.append("red")
#
#     for color in set(colors):
#         x0 = [x for x, c in zip(xs, colors) if c == color]
#         y0 = [y for y, c in zip(ys, colors) if c == color]
#         ax.scatter(x0, y0, color=color, alpha=0.6, label=f"{color2label[color]} ({len(x0)})")
#
#     ax.axline((0, 0), (1, 1), color='black')
#     ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', borderaxespad=0., mode="expand", ncol=2)
#     ax.set_xticks(positions)
#     ax.set_yticks(positions)
#
#     ax.set_ylabel("PSI transriptomics")
#     ax.set_xlabel("PSI proteomics")
#     plt.show()
