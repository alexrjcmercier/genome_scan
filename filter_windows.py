#coding:utf-8
import sys

scriptname, inputfile, annotfile, outputdir = sys.argv
inputfile_name = inputfile.split('/')[-1][0:-4] #[0:-4] to remove the '.txt'

annotdict = dict()
with open(annotfile, 'r') as annotf :
    readannot = annotf.readline().rstrip().split('\t')
    while readannot != [''] :
        featuretype = readannot[2]
        if featuretype == 'gene' :
            chrid = readannot[0]
            geneid = readannot[-1].split(';')[0].split('=')[-1]
            featuretype = readannot[2]
            startpos, stoppos = int(readannot[3]), int(readannot[4])
            if chrid not in annotdict.keys() :
                annotdict[chrid] = dict()
            annotdict[chrid][geneid] = [startpos, stoppos]
            readannot = annotf.readline().rstrip().split('\t')
        else :
            readannot = annotf.readline().rstrip().split('\t')

with open(inputfile, 'r') as inputf :
    chrid = inputfile_name.split('.')[-1].split('_')[0]
    chrannots = annotdict[chrid]
    windowsize = int(inputfile_name.rstrip('kb').split('_')[-1])*1000
    with open(''.join(x for x in [outputdir, inputfile_name, '_wf.txt']), 'w') as outputf :
        lineread = inputf.readline().rstrip().split('\t')
        while lineread[0] != 'Position' :
            outputf.write('\t'.join(x for x in lineread + ['\n']))
            lineread = inputf.readline().rstrip().split('\t')
        outputf.write('\t'.join(x for x in lineread + ['Genes/window\n']))
        lineread = inputf.readline().rstrip().split('\t')
        while lineread != [''] :
            startpos, stoppos = int(lineread[0]), int(lineread[0])+(windowsize - 1)
            geneslist = list()
            for k, v in annotdict[chrid].items() :
                genestart, genestop = v[0], v[-1]
                if genestart <= startpos and genestop >= startpos : #gene overlapping the window starting point
                    geneslist.append(k)
                elif genestart <= stoppos and genestop >= stoppos : #gene overlapping the window ending point
                    geneslist.append(k)
                elif genestart >= startpos and genestop <= stoppos : #gene inside the window with no overlapping of the window borders
                    geneslist.append(k)
            genesnb = len(set(geneslist))
            outputf.write('\t'.join(str(x) for x in lineread + [genesnb, '\n']))
            lineread = inputf.readline().rstrip().split('\t')
