import json
import os
import re

import numpy as np
import pandas as pd
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV

from typing import Dict, List, Tuple
from xgboost import XGBClassifier, plot_importance


# cellLines = list(reversed(['HUVEC', 'GM12878', 'hESC', 'HeLaS3', 'HepG2', 'K562']))


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
                    if spl[6] == '+':
                        result[transcriptName].append((int(spl[3]) - 1, int(spl[4])))
                    else:
                        result[transcriptName].append((-int(spl[4]), -int(spl[3]) + 1))
        return result


class Expression:
    @staticmethod
    def readFile(filename, expressionColumns, transcriptColumn, transcriptSet):
        result = {}
        with open(filename) as fileStream:
            header = fileStream.readline().rstrip().split('\t')
            expressionIndexes = [header.index(expressionColumn) for expressionColumn in expressionColumns]
            transcriptIndex = header.index(transcriptColumn)
            for line in fileStream:
                if len(line) == 0 or (len(line) > 0 and line[0] == "#"):
                    continue
                spl = line.rstrip().split('\t')
                a = [float(spl[i]) for i in expressionIndexes if spl[i] != 'NaN']
                if len(a) == 0:
                    continue
                for transcript in spl[transcriptIndex].split(','):
                    if transcript not in transcriptSet:
                        continue
                    result[transcript] = sum(a) / len(a)
        return result


class TransMembraneInfo:
    def __init__(self,
                 gene_id: str,
                 transcript_id: str,
                 ranges: List[Tuple[int, int]]):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.ranges = ranges

    @staticmethod
    def readTransMembraneInfo(filename, columns):
        transcript2transMembraneInfo = {}
        with open(filename, "r") as fs:
            header = fs.readline().rstrip().split(',')
            indexes = {column: header.index(text) for column, text in columns.items()}
            for line in fs:
                spl = line.rstrip().split(',')
                gene_id = spl[indexes["gene_id"]]
                transcript_id = spl[indexes["transcript_id"]]
                range_start = spl[indexes["range_start"]]
                range_end = spl[indexes["range_end"]]
                if transcript_id not in transcript2transMembraneInfo:
                    transcript2transMembraneInfo[transcript_id] = TransMembraneInfo(gene_id, transcript_id, [])
                if range_start == "":
                    continue
                transcript2transMembraneInfo[transcript_id].ranges.append((int(range_start), int(range_end)))
        return transcript2transMembraneInfo


class SplicingPath:
    def __init__(self, expression, transcripts, positions):
        self.expression = expression
        self.transcripts = transcripts
        self.positions = positions


class SplicingEvent:
    def __init__(self, path1, path2, isMsDetected, isInclusion):
        self.path1 = path1
        self.path2 = path2
        self.isMsDetected = isMsDetected
        self.isInclusion = isInclusion

    def getPsi(self):
        path1expression = sum(self.path1.expression)
        path2expression = sum(self.path2.expression)
        if path1expression == 0 and path2expression == 0:
            return -1
        return 2 * path1expression / (2 * path1expression + path2expression)

    @staticmethod
    def getIndexes(regex, spl):
        indexes = []
        for i in range(len(spl)):
            if len(re.findall(regex, spl[i])):
                indexes.append(i)
        return indexes

    @staticmethod
    def readFile(filename, gColumns, pColumns, transcriptColumns, positionsColumns, eventCodeColumn,
                 eventCodeList, transcriptSet):
        eventCodeSet = set(eventCodeList)
        result = []
        with open(filename) as fileStream:
            header = fileStream.readline().rstrip().split('\t')
            gIndexes = (SplicingEvent.getIndexes(gColumns[0], header),
                        SplicingEvent.getIndexes(gColumns[1], header))
            pIndexes = (SplicingEvent.getIndexes(pColumns[0], header),
                        SplicingEvent.getIndexes(pColumns[1], header))
            peptIndexes = (SplicingEvent.getIndexes("Peptides.([^.]*).Path1", header),
                           SplicingEvent.getIndexes("Peptides.([^.]*).Path2", header))
            transcriptIndexes = (SplicingEvent.getIndexes(transcriptColumns[0], header)[0],
                                 SplicingEvent.getIndexes(transcriptColumns[1], header)[0])
            positionsIndexes = (SplicingEvent.getIndexes(positionsColumns[0], header)[0],
                                SplicingEvent.getIndexes(positionsColumns[1], header)[0])
            eventCodeIndex = SplicingEvent.getIndexes(eventCodeColumn, header)[0]
            for line in fileStream:
                if len(line) == 0 or (len(line) > 0 and line[0] == "#"):
                    continue
                spl = line.rstrip('\n').split('\t')
                if spl[eventCodeIndex] not in eventCodeSet:
                    continue
                path1 = SplicingPath(
                    [int(spl[gIndex]) for gIndex in gIndexes[0]],
                    [transcript for transcript in spl[transcriptIndexes[0]].split(';') if transcript in transcriptSet],
                    [int(position) for position in spl[positionsIndexes[0]].split(';')]
                )
                path2 = SplicingPath(
                    [int(spl[gIndex]) for gIndex in gIndexes[1]],
                    [transcript for transcript in spl[transcriptIndexes[1]].split(';') if transcript in transcriptSet],
                    [int(position) for position in spl[positionsIndexes[1]].split(';')]
                )
                if len(path1.transcripts) == 0 or len(path2.transcripts) == 0:
                    continue
                for i in range(2):
                    for j in peptIndexes[i]:
                        if spl[j] == "":
                            continue
                event = SplicingEvent(path1, path2,
                                      sum([int(spl[pIndex]) for pIndex in pIndexes[0]]) > 0 and
                                      sum([int(spl[pIndex]) for pIndex in pIndexes[1]]) > 0,
                                      sum([int(spl[pIndex]) for pIndex in pIndexes[0]]) > 0)
                if 0.05 <= event.getPsi() <= 0.95 and sum(event.path1.expression) + sum(event.path2.expression) > 20:
                    result.append(event)
        return result


class Point:
    def __init__(self, geneExpression, exonLength, exonMod3, psi, cdsLength, proteasePeptides, scoreLysArg,
                 # isTransMembrane,
                 isMsDetected, event):
        # X
        self.geneExpression = geneExpression
        self.exonLength = exonLength
        self.exonMod3 = exonMod3
        self.psi = psi
        self.cdsLength = cdsLength
        self.proteasePeptides = proteasePeptides
        self.scoreLysArg = scoreLysArg
        # self.isTransMembrane = isTransMembrane
        # Y
        self.isMsDetected = isMsDetected
        # meta
        self.event = event

    @staticmethod
    def __getProteasePeptidesImpl(s0, s1, minLength, maxLength, regExp, missedCleavages):
        result = set()
        pattern = re.compile(regExp)
        starts = [i.start() for i in pattern.finditer(s0)]
        ends = [i.start() for i in pattern.finditer(s1)]
        for i in range(len(starts)):
            for j in range(len(ends)):
                if (len(starts) - i + 1 + j) > missedCleavages:
                    break
                if minLength <= (len(s0) - starts[i] + ends[j]) <= maxLength:
                    result.add(s0[(starts[i] + 1):] + s1[:(ends[j] + 1)])

        return result

    @staticmethod
    def __getProteasePeptides(event, transcript2exons, transcript2proteinSequence, proteases, minLength, maxLength):
        result0 = {protease: [set(), set()] for protease in proteases}
        for ind, path, pathPos in [(0, event.path1, 0), (1, event.path2, 0), (1, event.path2, 2)]:
            for transcript in path.transcripts:
                if transcript not in transcript2exons or transcript not in transcript2proteinSequence:
                    continue
                p = 0
                for exon in transcript2exons[transcript]:
                    p += exon[1] - exon[0] + 1
                    if exon[1] == path.positions[pathPos]:
                        break
                p = p // 3
                s = transcript2proteinSequence[transcript]
                for protease in proteases:
                    result0[protease][ind] |= Point.__getProteasePeptidesImpl(
                        s[max(p - maxLength, 0):p],
                        s[p:min(p + maxLength, len(s))],
                        minLength,
                        maxLength,
                        proteases[protease]["RegExp"],
                        proteases[protease]["MissedCleavages"]
                    )
        result = {}
        for protease in proteases:
            a = result0[protease][0]
            b = result0[protease][1]
            result[protease] = (a - b, b - a)
        return result

    @staticmethod
    def __getScoreArgLys(event, transcript2exons, transcript2proteinSequence):
        score = [0, 0]
        for ind, path, pathPos in [(0, event.path1, 0), (1, event.path2, 0), (1, event.path2, 2)]:
            transcript = path.transcripts[0]
            if transcript not in transcript2exons or transcript not in transcript2proteinSequence:
                continue
            p = 0
            for exon in transcript2exons[transcript]:
                p += exon[1] - exon[0] + 1
                if exon[1] == path.positions[pathPos]:
                    break
            p = p // 3
            s = transcript2proteinSequence[transcript]
            for i in range(max(0, p - 1), min(p + 2, len(s))):
                if s[i] == "R" or s[i] == "K":
                    score[ind] += 1
                    break
        if score[0] == 0 and score[1] == 0:
            return 0
        if score[0] != 0 and score[1] != 0:
            return 2
        return 1

    @staticmethod
    def get(events, transcript2proteinSequence, transcript2expression,
            transcript2exons,
            transcript2transMembraneInfo,
            proteases, minLength, maxLength):
        points = []
        for event in events:
            cdsLengths = [len(transcript2proteinSequence[t]) for t in event.path1.transcripts] + \
                         [len(transcript2proteinSequence[t]) for t in event.path2.transcripts]
            ti = 0
            while ti < len(event.path1.transcripts) and event.path1.transcripts[ti] not in transcript2expression:
                ti += 1
            if ti == len(event.path1.transcripts):
                continue
            isTransMembrane = False
            for t in event.path1.transcripts + event.path2.transcripts:
                if t in transcript2transMembraneInfo and len(transcript2transMembraneInfo[t].ranges) != 0:
                    isTransMembrane = True
            # if isTransMembrane:
            #     continue
            points.append(
                Point(
                    transcript2expression[event.path1.transcripts[0]],
                    event.path2.positions[2] - event.path2.positions[1],
                    (event.path2.positions[2] - event.path2.positions[1]) % 3 == 0,
                    event.getPsi(),
                    sum(cdsLengths) / len(cdsLengths),
                    Point.__getProteasePeptides(
                        event,
                        transcript2exons,
                        transcript2proteinSequence,
                        proteases,
                        minLength,
                        maxLength
                    ),
                    Point.__getScoreArgLys(
                        event,
                        transcript2exons,
                        transcript2proteinSequence
                    ),
                    event.isMsDetected,
                    # isTransMembrane,
                    event
                )
            )
        return points

    @staticmethod
    def convert2pandas(points, proteases):
        cs = ["isMsDetected", "geneExpression", "exonLength", "exonMod3", "psi", "cdsLength", "scoreLysArg"]
        for protease in proteases:
            # cs.append(f"{protease}.MinPeptides")
            # cs.append(f"{protease}.MaxPeptides")
            cs.append(f"{protease} #peptides")
            # cs.append(f"{protease}.Path2")
        ds = []
        for point in points:
            a = [
                int(point.isMsDetected),
                point.geneExpression,
                point.exonLength,
                int(point.exonMod3),
                point.psi,
                point.cdsLength,
                f"{point.scoreLysArg}.0"
            ]
            # sum([int(spl[pIndex]) for pIndex in pIndexes[0]]) > 0 and
            # sum([int(spl[pIndex]) for pIndex in pIndexes[1]]) > 0
            cnt0 = sum(point.event.path1.expression)
            cnt1 = sum(point.event.path2.expression) // 2
            if cnt0 < cnt1:
                for protease in proteases:
                    a.append(len(point.proteasePeptides[protease][0]))
            else:
                for protease in proteases:
                    a.append(len(point.proteasePeptides[protease][1]))

                # x, y = len(point.proteasePeptides[protease][0]), len(point.proteasePeptides[protease][1])
                # a.append(min(x, y))

                # a.append(max(x, y))
                # a.append(len(point.proteasePeptides[protease][0]))
                # a.append(len(point.proteasePeptides[protease][1]))
            ds.append(a)
        return pd.DataFrame(ds, columns=cs)


# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.model_selection import StratifiedKFold
from sklearn.inspection import permutation_importance

from numpy import interp
import seaborn as sb


class BinaryClassification:
    param_grid = {
        "max_depth": [3],
        'min_child_weight': [2],
        "learning_rate": [0.05],
        "gamma": [0.0, 1.5, 2.0, 2.5, 3.0],
        'reg_alpha': [1.15],
        "reg_lambda": [4.0],
        "scale_pos_weight": [4.44],
        "subsample": [0.3],
        "colsample_bytree": [0.65],
        # "max_depth": [2, 3],
        # 'min_child_weight': [3],
        # "learning_rate": [0.1],
        # "gamma": [2.5],
        # 'reg_alpha': [1.0],
        # "reg_lambda": [1.5],
        # "scale_pos_weight": [4.44],
        # "subsample": [0.5],
        # "colsample_bytree": [0.5]
    }

    @staticmethod
    def run(df):
        X = df[df.columns.difference(["isMsDetected"])]
        Y = df["isMsDetected"]

        # sb.heatmap(X.corr(), cmap="YlGnBu", annot=True)
        # plt.show()

        xgb_cl = XGBClassifier(
            use_label_encoder=False,
            verbosity=0, silent=True,
            random_state=42
        )

        # Init Grid Search
        # grid_cv = GridSearchCV(
        #     estimator=xgb_cl,
        #     param_grid=BinaryClassification.param_grid,
        #     n_jobs=5,
        #     cv=5,
        #     scoring='roc_auc'
        # )
        grid_cv = GridSearchCV(
            estimator=xgb_cl,
            param_grid=BinaryClassification.param_grid,
            # n_iter=100,
            n_jobs=8,
            cv=7,
            scoring='roc_auc',
            # random_state=42
        )
        _ = grid_cv.fit(X, Y)
        print(grid_cv.best_params_)
        print(grid_cv.best_score_)
        clf = grid_cv.best_estimator_

        perm_importance = permutation_importance(clf, X, Y, n_repeats=100, random_state=42)
        sorted_idx = perm_importance.importances_mean.argsort()
        plt.barh(X.columns[sorted_idx], perm_importance.importances_mean[sorted_idx],
                 xerr=perm_importance.importances_std[sorted_idx])
        plt.xlabel("Permutation Importance")
        plt.show()

        # explainer = shap.TreeExplainer(clf)
        # shap_values = explainer.shap_values(X)
        # shap.summary_plot(shap_values, X, plot_type="bar")
        #
        # shap.summary_plot(shap_values, X)

        BinaryClassification.plot_importance(clf)

        # show_weights(clf)
        BinaryClassification.plotRocCurve(X, Y, clf)

        # seed = 7
        # test_size = 0.33
        # X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)
        # # fit model no training data
        # model = XGBClassifier()
        # model.fit(X_train, y_train)
        # # plot feature importance
        #
        # print(model)
        #
        # y_pred = model.predict(X_test)
        # predictions = [round(value) for value in y_pred]
        # accuracy = accuracy_score(y_test, predictions)
        # print("Accuracy: %.2f%%" % (accuracy * 100.0))
        #
        # plot_importance(model)
        # pyplot.show()

    @staticmethod
    def plotRocCurve(X, Y, clf):
        fig1 = plt.figure(figsize=[12, 12])
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)
        i = 1
        cv = StratifiedKFold(n_splits=7, shuffle=True, random_state=42)
        for train, test in cv.split(X, Y):
            prediction = clf.fit(X.iloc[train], Y.iloc[train]).predict_proba(X.iloc[test])
            fpr, tpr, t = roc_curve(Y[test], prediction[:, 1])
            tprs.append(interp(mean_fpr, fpr, tpr))
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            plt.plot(fpr, tpr, lw=2, alpha=0.3, color='grey', label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
            i = i + 1

        plt.plot([0, 1], [0, 1], linestyle='--', lw=1, color='black')
        mean_tpr = np.mean(tprs, axis=0)
        mean_auc = auc(mean_fpr, mean_tpr)
        plt.plot(mean_fpr, mean_tpr, color='blue',
                 label=r'Mean ROC (AUC = %0.2f )' % mean_auc, lw=2, alpha=1)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC')
        plt.legend(loc="lower right")
        plt.show()

    @staticmethod
    def plot_importance(clr):
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
        # plot_importance(clr, ax=ax1)
        # ax1.set_title("xgboost defualt importance")

        plot_importance(clr, importance_type="cover", ax=ax1)
        ax1.set_title("xgboost cover")

        plot_importance(clr, importance_type="gain", ax=ax3)
        ax3.set_title("xgboost gain")

        plt.show()

        # plot_tree(clr)


def plotPsiDistribution(df):
    rsx = [[0.00, 0.15], [0.15, 0.25], [0.25, 0.35], [0.35, 0.45], [0.45, 0.55], [0.55, 0.65],
           [0.65, 0.75], [0.75, 0.85], [0.85, 0.95]]
    rsy = [[0, 0] for i in range(len(rsx))]
    for index, row in df.iterrows():
        if row["exonMod3"] == 0 or row["geneExpression"] < 1:
            continue
        val = row["psi"]
        i = 0
        while not (rsx[i][0] < val <= rsx[i][1]):
            i += 1
        if row["isMsDetected"] == 1:
            rsy[i][0] += 1
        else:
            rsy[i][1] += 1
    fig, ax = plt.subplots()
    width = 0.1
    x = [j - 0.05 for i, j in rsx]
    y1 = [i / (i + j) for i, j in rsy]
    ax.bar(x, y1, width=width, color="black")
    print(y1)

    y2 = [j / (i + j) for i, j in rsy]
    ax.bar(x, y2, width=width, bottom=y1, color="grey")

    ax.set_xticks(x)
    # ax.set_xticklabels(x)
    # ax.set_xlim(0.05, 0.95)
    # ax.set_ylim(0, 1)
    plt.show()


if __name__ == "__main__":
    with open("parameters.json", 'r') as parameters_fs:
        params = json.loads(parameters_fs.read())
    if not os.path.isfile(params["Output"]["TmpFile"]):
        transcript2proteinSequence = ProteinFastaRecord.readFile(
            params["Input"]["ProteinFastaFile"]["Filename"],
            params["Input"]["ProteinFastaFile"]["NameRegularExpression"]
        )
        transcript2exons = GtfAnnotation.readFile(
            params["Input"]["GtfAnnotationFile"]["Filename"],
            params["Input"]["GtfAnnotationFile"]["NameRegularExpression"],
            set([transcript for transcript in transcript2proteinSequence.keys()])
        )
        transcript2expression = Expression.readFile(
            params["Input"]["ExpressionFile"]["Filename"],
            params["Input"]["ExpressionFile"]["Columns"]["Expression"],
            params["Input"]["ExpressionFile"]["Columns"]["Transcript"],
            set([transcript for transcript in transcript2proteinSequence.keys()])
        )
        events = SplicingEvent.readFile(
            params["Input"]["SplicingSiteFile"]["Filename"],
            params["Input"]["SplicingSiteFile"]["Columns"]["Genomics"],
            params["Input"]["SplicingSiteFile"]["Columns"]["Proteomics"],
            params["Input"]["SplicingSiteFile"]["Columns"]["Transcript"],
            params["Input"]["SplicingSiteFile"]["Columns"]["Positions"],
            params["Input"]["SplicingSiteFile"]["Columns"]["EventCode"],
            [params["Input"]["SplicingSiteFile"]["EventCode"]],
            set([transcript for transcript in transcript2proteinSequence.keys()])
        )
        transcript2transMembraneInfo = TransMembraneInfo.readTransMembraneInfo(
            params["Input"]["TransMembraneFile"]["Filename"],
            params["Input"]["TransMembraneFile"]["Columns"]
        )
        points = Point.get(
            events,
            transcript2proteinSequence,
            transcript2expression,
            transcript2exons,
            transcript2transMembraneInfo,
            params["Input"]["Digestion"]["Protease"],
            params["Input"]["Digestion"]["Length"]["Minimum"],
            params["Input"]["Digestion"]["Length"]["Maximum"]
        )
        df = Point.convert2pandas(points, list(params["Input"]["Digestion"]["Protease"].keys()))
        df.to_csv(params["Output"]["TmpFile"], index=False)
    df = pd.read_csv(params["Output"]["TmpFile"])
    # plotPsiDistribution(df)
    BinaryClassification.run(df)

    # print(len([point for point in points if point.isMsDetected]))
    # print(len(points))
    # print(len([event.isMsDetected for event in events if event.isMsDetected]))
    # print(len(events))
    # print(transcript2proteinSequence["ENST00000379116"])
    # print(transcript2exons["ENST00000379116"])
    # gColumns, pColumns, transcriptColumns, positionsColumns, eventCodeColumn,
    # eventCodeSet, transcriptSet
    # print(transcript2proteinSequence["ENST00000379116"])
    # print(transcript2exons["ENST00000379116"])
