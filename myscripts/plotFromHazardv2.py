#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 16:35:33 2019

@author: rsosaric


USO: python plotFromHazard.py -i pathTo/locations.csv pathTo/fichero.nc
Ejemplo si se ejecuta desde tcrm/myscripts: python plotFromHazard.py -i locations.csv ../output/corrida_min_de_obs/hazard/hazard.nc
"""

from netCDF4 import Dataset
import os
import optparse
import matplotlib.pyplot as plt
import matplotlib
import Utilities.nctools as nctools
import math as math

##Variables globales
output_formats=['.png','.pdf']
legendPosition='lower right'
legendTextSize=11
legendTextSize_Hist=13

## A;os que se plotean en el histograma resumen
years_forHist=[2,10,100,1000]
down_margin=0.38

##Colores para el histograma
colors = ['blue','orange','green','purple','brown','pink','gray','olive','cyan','red']  

#### Funciones aqui ######
#funcion para evaluar G(z), me imagino que z son los a;os (chequea que sea asi la funcion)
def evaluarFuncionG(z,shp,loc,scale):
    return math.exp((1+shp*(z-loc)/scale)-1/shp)

def getEvalGListFromYears(years,shp,loc,scale):
    valList=[]
    for year in years:
        valList.append(evaluarFuncionG(year,shp,loc,scale))
    return valList

def getDictCoord(incoord,dictcoord):
    coord=-100000.00
    id_dict=-1
    for i in range(0,len(dictcoord)):
        ##find coord range
        if dictcoord[i]>incoord:
            if i>0:
                if abs(dictcoord[i]-incoord) < abs(dictcoord[i-1]-incoord):
                    coord=dictcoord[i]
                    id_dict=i
                else:
                    coord=dictcoord[i-1]
                    id_dict=(i-1)
            break
    if coord==-100000.00:
        raise Exception('input coordenate no in dictionarie range. The value was: {}'.format(incoord))
    return id_dict

def getEvalDictFromYears(years,lat,lon,varDict):
    valList=[]
    for yearid in range(0,len(years)):
        #print varDict[yearid,lat,lon]
        valList.append(varDict[yearid,lat,lon])
    return valList
    
def readLocationsFromFile(filePath):
    infile=open(filePath)
    locations={}
    for line in infile.readlines():
        line=line.split(",")
        locations[line[0]]=[float(line[1]),float(line[2])]
        
    infile.close()
    return locations

def plotCurvesFromLocations(locations,windvals,years,logx=True,logy=False):
    names=list(locations)
    #print ("Plotting: ", names)
    plot=plt.figure()
    ax = plot.add_subplot(111)
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    ax.set_xticks((10, 100, 1000))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    x_vals=years
    for name in names:
        lat=locations[name][0]
        lon=locations[name][1]
        
        y_vals=getEvalDictFromYears(years,lat,lon,windvals)
        #print x_vals,y_vals
        plt.plot(x_vals,y_vals,linewidth=0.5,linestyle="-",label=name)
    plt.legend(loc=legendPosition,prop={'size': legendTextSize})
    plt.xlabel(r"""$A\~nos$""", fontsize=20)
    plt.ylabel("Velocidad del viento (m/s)", fontsize=20)
    
    return plot

def summaryHistogram(locations,windvals,windvals_up,windvals_down,years):
    plot=plt.figure()
    ax = plot.add_subplot(111)
   
    x_vals=[]
    x_vals_labels=[]
    delta=1.0
    
    for x in range(1,len(locations)+1):
        x_vals.append(x*delta)
        x_vals_labels.append(x*delta+delta/2)
        
    years_inv=[]
    for i in range(0,len(years)):
        years_inv.append(int(years[len(years)-i-1]))
        
    year_id=len(years)-1
    color_id=0
    for year in years_inv:
        if year in years_forHist:
            y_vals=[]
            for name in locations:
                lat=locations[name][0]
                lon=locations[name][1]
                y_vals.append(windvals[year_id,lat,lon])
            plt.bar(x_vals, y_vals, color=colors[color_id],alpha=0.5,label=str(year))
            color_id+=1
        year_id-=1
    plt.ylabel("Velocidad del viento (m/s)", fontsize=20)
    plt.xticks(x_vals_labels, locations, rotation='vertical')
    plt.legend(title=r"""$A\~nos$""",loc='upper center',prop={'size': legendTextSize_Hist},bbox_to_anchor=(1.02, 1.05),fancybox=True, shadow=True)
    plt.subplots_adjust(bottom=down_margin)
   
    return plot
def savePlot(plot,name,outdir):
    for ext in output_formats:
        path=outdir+"/"+name+ext
        plot.savefig(path)
        print ("Saved plot: "+path)

## Esto es para agregar opciones si hacen falta
p = optparse.OptionParser()
p.add_option("-o", "--outdir", type="string",help="Path to output Dir", dest="outDir", default="output_plotsFromHazard")
p.add_option("-i", "--infile", type="string",help="Path to locations file", dest="inLoc", default="")

## Lavariable options contiene las opciones, y la varaiable args los argumentos que le entras que no sean opciones, por ejemplo: python programa.py -o opcion1 argumento1 argumento2 etc... 
(options, args) = p.parse_args()
##Comprobar que se entre el fichero, si no termina el programa
if len(args)!=1:
    print ("Invalid number of input parameters!!!!!!!")
    exit()
elif options.inLoc=="":
    print ("No input locations file!!!!")
    exit()

try:
    os.makedirs(options.outDir)
except:
    pass


################################################## Parte principal del codigo ##################################################
## Leer data desde .nc
inputNC = nctools.ncLoadFile(args[0])

## importar datos de los lugares desde tu .csv
locations=readLocationsFromFile(options.inLoc)

## Informacion sobre la estructura del .nc y las variables que contiene
nctools.ncFileInfo(args[0])

##leer variables
lons = nctools.ncGetVar(inputNC,"lon")
lats = nctools.ncGetVar(inputNC,"lat")
locs = nctools.ncGetVar(inputNC,"loc")
scales = nctools.ncGetVar(inputNC,"scale")
shps = nctools.ncGetVar(inputNC,"shp")
#valores del viento dict[years,lat,lon]=valor
windvals=nctools.ncGetVar(inputNC,"wspd")
windvals_up=nctools.ncGetVar(inputNC,"wspdupper")
windvals_down=nctools.ncGetVar(inputNC,"wspdlower")
years=nctools.ncGetVar(inputNC,"years")

lat=28.35
lon=268.4
lon=lon-360

print lat,lon,windvals[1,lat,lon]

wspd = inputNC.variables['wspd'][:, j, i]
print wspd
