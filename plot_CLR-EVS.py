#coding:utf-8
#!/bin/env python2.7
import sys
import os
import matplotlib
matplotlib.use('Agg')
import pylab
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
from matplotlib.backends.backend_pdf import PdfPages

args = sys.argv
inputdirs = args[1:-1] #directories where SweeD_Report.BCINxx..._nnkb.txt and Seqfile_BCINxx..._evolstats_nnkb.txt are stored.
outputdir = args[-1]

#Options
mbconvert = True
scale = 0.2 #in Mbp
windows = ['10kb', '50kb', '100kb'] #Frames in which statistics are calculated. From smaller frame to bigger.
chromosomes = ['BCIN01', 'BCIN02', 'BCIN03', 'BCIN04', 'BCIN05', 'BCIN06', 'BCIN07', 'BCIN08', 'BCIN09', 'BCIN10', 'BCIN11', 'BCIN12', 'BCIN13', 'BCIN14', 'BCIN15']
filestypes = {'CLR':['Position', 'Likelihood', 'Alpha'], 'EVS':['Window', 'Start', 'Stop', '\'Hns\'', '\'S\'', '\'Pi\'', '\'thetaW\'', '\'Hsd\'', '\'D\'']}

#Functions
def Prepline(inputfile) :
    tupline = tuple(inputfile.readline().rstrip().split('\t'))
    return tupline

def GenTreeDict(chromosome, filesdict, filestype, fileslist, windows) :
    filesdict[filestype] = dict()
    for inputfile in fileslist :
        if chromosome in inputfile :
            for window in windows :
                if window in inputfile :
                    filesdict[filestype][window] = inputfile

#Main Script
inputfiles, CLRfiles, EVSfiles = (list() for i in range(3))
for inputdir in inputdirs :
    inputfiles = os.listdir(inputdir)
    for inputfile in inputfiles :
        if 'SweeD_Report' in inputfile :
            CLRfiles.append(''.join([inputdir, inputfile])) #list of valid sweed CLR files
        elif 'evolstats' in inputfile :
            EVSfiles.append(''.join([inputdir, inputfile])) #list of valid evolstats files

filesdict = dict()
for bcin in chromosomes : #creates the arborescence of files for each chromosome studied
    filesdict[bcin] = dict()
    GenTreeDict(bcin, filesdict[bcin], 'CLR', CLRfiles, windows)
    GenTreeDict(bcin, filesdict[bcin], 'EVS', EVSfiles, windows)

datadict = dict()
for chrnb in chromosomes :#['BCIN01'] :
    datadict[chrnb] = dict()

    filetype = 'CLR'
    datadict[chrnb][filetype] = dict()
    for key, value in filesdict[chrnb][filetype].items() :
        with open(value, 'r') as infile :
            datadict[chrnb][filetype][key] = dict()
            intuple = ('0')
            while intuple[0] != filestypes[filetype][0] :
                intuple = Prepline(infile)
            for tupelem in intuple :
                datadict[chrnb][filetype][key][tupelem] = list()
            intuple = Prepline(infile)
            nbcols = len(intuple)
            while intuple != ('',) :
                for i in range(0, nbcols) :
                    if mbconvert == True and i == 0 :
                        datadict[chrnb][filetype][key][filestypes[filetype][i]].append(float(intuple[i])/1000000)
                    elif i in range(1,3) :
                        datadict[chrnb][filetype][key][filestypes[filetype][i]].append(float(intuple[i]))
                intuple = Prepline(infile)
    
    filetype = 'EVS'
    datadict[chrnb][filetype] = dict()
    for key, value in filesdict[chrnb][filetype].items() :
        with open(value, 'r') as infile :
            datadict[chrnb][filetype][key] = dict()
            intuple = ('0')
            while intuple[0] != filestypes[filetype][0] :
                intuple = Prepline(infile)
            for tupelem in intuple :
                datadict[chrnb][filetype][key][tupelem] = list()
            intuple = Prepline(infile)
            while intuple != ('',) :
                for i in range(0, len(intuple)) :
                    if mbconvert == True and i == 1 :
                        datadict[chrnb][filetype][key][filestypes[filetype][i]].append(float(intuple[i])/1000000)
                    elif i in (3, 7, 8) :
                        if intuple[i] == 'None' :
                            datadict[chrnb][filetype][key][filestypes[filetype][i]].append(0.0)
                        else :
                            datadict[chrnb][filetype][key][filestypes[filetype][i]].append(float(intuple[i]))
                intuple = Prepline(infile)

#Graphs generation
with PdfPages('Summary_genome_statistics.pdf') as pdf :
    for bcin in chromosomes : #['BCIN01'] :
        fig = plt.figure(figsize=(8, 6), dpi=1200)
        largestwindow = windows[-1]
        maxpos = datadict[bcin]['CLR'][largestwindow]['Position'][-1]
        linestyles = {'10kb':'-', '50kb':'--', '100kb':'--'}
        linecolors = {'10kb':'darkorange', '50kb':'steelblue', '100kb':'turquoise'}
        plt.suptitle(bcin, fontsize=10)

        #Graph with CLR values along genome
        #(genome position corresponds to the starting point of the sliding window)
        filetype = 'CLR'
        clr1 = plt.subplot(4,1,1)
        clr1.spines['top'].set_color('none')
        clr1.spines['right'].set_color('none')
        clr1.spines['bottom'].set_linewidth(0.5)
        clr1.spines['left'].set_linewidth(0.5)
        clr1.set_xlim(0.0, maxpos)
        clr1.xaxis.set_major_locator(tkr.MultipleLocator(0.2)) #Set main ticks of x axis.
        clr1.xaxis.set_minor_locator(tkr.MultipleLocator(0.05)) #Set secondary tick (without labels) of x axis.
        clr1.tick_params(axis='x', which='major', width=.5, labelsize=4, direction='out')
        clr1.tick_params(axis='x', which='minor', width=.25, direction='out')

        for window in windows :
            clr1.plot(datadict[bcin][filetype][window]['Position'], 
            datadict[bcin][filetype][window]['Likelihood'], 
            color=linecolors[window], linewidth=0.5, linestyle=linestyles[window], zorder=0.1)

        clr1.set_ylabel('CLR', fontsize=6)
        lowlim, uplim = clr1.get_ylim()
        clr1.set_ylim(0.0, uplim)
        clr1.yaxis.set_major_locator(tkr.MultipleLocator(20))
        clr1.yaxis.set_minor_locator(tkr.MultipleLocator(5))
        clr1.tick_params(axis='y', which='major', width=.5, labelsize=4, direction='out')
        clr1.tick_params(axis='y', which='minor', width=.25, direction='out')

        #Graph with Tajima's D values along the genome
        #(genome position corresponds to the starting point of the sliding window)
        filetype = 'EVS'
        evs1 = plt.subplot(4,1,2)
        evs1.spines['top'].set_color('none')
        evs1.spines['right'].set_color('none')
        evs1.spines['bottom'].set_position(('data', 0))
        evs1.spines['bottom'].set_linewidth(0.5)
        evs1.spines['left'].set_linewidth(0.5)
        #evs1.axhline(0, linestyle='-', linewidth=0.5, color='lightgray') #Add a horizontal line at y=0.
        evs1.set_xlim(0.0, maxpos)
        evs1.xaxis.set_major_locator(tkr.IndexLocator(base=0.2, offset=0.2)) #Set main ticks of x axis.
        evs1.xaxis.set_minor_locator(tkr.MultipleLocator(0.05)) #Set secondary tick (without labels) of x axis.
        evs1.tick_params(axis='x', which='major', width=.5, labelsize=4, direction='inout')
        evs1.tick_params(axis='x', which='minor', width=.25, direction='inout')

        for window in windows :
            evs1.plot(datadict[bcin][filetype][window]['Start'], 
            datadict[bcin][filetype][window]['\'D\''], 
            color=linecolors[window], linewidth=0.5, linestyle=linestyles[window], zorder=0.1)

        evs1.set_ylabel('Tajima\'s D', fontsize=6)
        lowlim, uplim = evs1.get_ylim()
        evs1.set_ylim(lowlim, uplim)
        evs1.yaxis.set_major_locator(tkr.MultipleLocator(1))
        evs1.yaxis.set_minor_locator(tkr.MultipleLocator(0.25))
        evs1.tick_params(axis='y', which='major', width=.5, labelsize=4, direction='out')
        evs1.tick_params(axis='y', which='minor', width=.25, direction='out')

        #Graph with Fay & Wu's H (standardized by Zeng) values along the genome
        #(genome position corresponds to the starting point of the sliding window)
        evs2 = plt.subplot(4,1,3)
        evs2.spines['top'].set_color('none')
        evs2.spines['right'].set_color('none')
        evs2.spines['bottom'].set_position(('data', 0))
        evs2.spines['bottom'].set_linewidth(0.5)
        evs2.spines['left'].set_linewidth(0.5)
        evs2.set_xlim(0.0, maxpos)
        evs2.xaxis.set_major_locator(tkr.IndexLocator(base=0.2, offset=0.2)) #Set main ticks of x axis.
        evs2.xaxis.set_minor_locator(tkr.MultipleLocator(0.05)) #Set secondary tick (without labels) of x axis.
        evs2.tick_params(axis='x', which='major', width=.5, labelsize=4, direction='inout')
        evs2.tick_params(axis='x', which='minor', width=.25, direction='inout')

        for window in windows :
            plt.plot(datadict[bcin][filetype][window]['Start'], 
            datadict[bcin][filetype][window]['\'Hsd\''], 
            color=linecolors[window], linewidth=0.5, linestyle=linestyles[window], zorder=0.1)
        
        lowlim, uplim = evs2.get_ylim()
        evs2.set_ylim(lowlim, uplim)
        evs2.yaxis.set_major_locator(tkr.MultipleLocator(1))
        evs2.yaxis.set_minor_locator(tkr.MultipleLocator(0.25))
        evs2.tick_params(axis='y', which='major', width=.5, labelsize=4, direction='out')
        evs2.tick_params(axis='y', which='minor', width=.25, direction='out')
        
        evs2.set_ylabel('Fay and Wu\'s H (Std)', fontsize=6)


        #Graph with Fay & Wu's H (raw) values along the genome
        #(genome position corresponds to the starting point of the sliding window)
        evs3 = plt.subplot(4,1,4)
        evs3.spines['top'].set_color('none')
        evs3.spines['right'].set_color('none')
        evs3.spines['bottom'].set_position(('data', 0))
        evs3.spines['bottom'].set_linewidth(0.5)
        evs3.spines['left'].set_linewidth(0.5)
        evs3.set_xlim(0.0, maxpos)
        evs3.xaxis.set_major_locator(tkr.IndexLocator(base=0.2, offset=0.2)) #Set main ticks of x axis.
        evs3.xaxis.set_minor_locator(tkr.MultipleLocator(0.05)) #Set secondary tick (without labels) of x axis.
        evs3.tick_params(axis='x', which='major', width=.5, labelsize=4, direction='inout')
        evs3.tick_params(axis='x', which='minor', width=.25, direction='inout')

        for window in windows :
            plt.plot(datadict[bcin][filetype][window]['Start'], 
            datadict[bcin][filetype][window]['\'Hns\''], 
            color=linecolors[window], linewidth=0.5, linestyle=linestyles[window], label=window, zorder=0.1)
        
        lowlim, uplim = evs3.get_ylim()
        evs3.set_ylim(lowlim, uplim)
        evs3.yaxis.set_major_locator(tkr.MultipleLocator(100))
        evs3.yaxis.set_minor_locator(tkr.MultipleLocator(20))
        evs3.tick_params(axis='y', which='major', width=.5, labelsize=4, direction='out')
        evs3.tick_params(axis='y', which='minor', width=.25, direction='out')
        
        evs3.set_ylabel('Fay and Wu\'s H (Raw)', fontsize=6)
        #evs3.set_xlabel('Position along chromosome', labelpad=40, fontsize=6)

        #Save figure
        plt.tight_layout()
        fig.text(x=0.5, y=0.015, s='Position along chromosome sequence (Mb)', fontsize=8, horizontalalignment='center')
        fig.legend(loc='lower right', ncol=1, fontsize=6, facecolor='white', framealpha= 0.75, frameon=False)
        pylab.savefig(''.join([bcin, '.png']))
        pdf.savefig()
        plt.close()
