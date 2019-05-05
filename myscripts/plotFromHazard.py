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

## Formatos en los que saldran los plots
output_formats=['.png','.pdf']

## Ajustes para las curvas
## posicion de la leyenda. Lugares posibles: 'best','upper right','upper left','lower left','lower right','right','center left','center right','lower center','upper center','center'
legendPosition='lower right'
# tama;o de letra de la leyenda (al disminuir el taman; de la letra disminuye automaticamente el espacio que ocupa la leyenda, usa esto para ajustar que la leyenda este arriba de las curvas)
legendTextSize=10


## Ajustes para el histograma
# a;os a plotear
years_forHist=[2,10,100,1000]
# margen inferior (aumentar si los nombres no caben)
down_margin=0.38
# tama;o de letra de la leyenda
legendTextSize_Hist=13
#Colores para el histograma
colors = ['blue','orange','purple','green','brown','pink','gray','olive','cyan','red','navy','firebrick','darkblue']
#Es la transparencia en las barras. Mientras menor el numero mas suave el color
alphaHist=0.3
#include hist with errs
include_hist_errs=True



#### Funciones aqui ######
#funcion para evaluar G(z), me imagino que z son los a;os (chequea que sea asi la funcion)
def evaluarFuncionG(z,shp,loc,scale):
    return math.exp((1+shp*(z-loc)/scale)-1/shp)

def getEvalGListFromYears(years,shp,loc,scale):
    valList=[]
    for year in years:
        valList.append(evaluarFuncionG(year,shp,loc,scale))
    return valList

#Fucion para obtener la posicion de un valor en un arreglo: incoord es el valor y dictcoord es donde se busca su posicion
def getDictCoord(incoord,dictcoord):
    coord=-100000.00
    id_dict=-1
    delta=dictcoord[1]-dictcoord[0]
    
    for i in range(0,len(dictcoord)):
        ##find coord range
        if delta>0:      
            if dictcoord[i]>incoord:
                if i>0:
                    if abs(dictcoord[i]-incoord) < abs(dictcoord[i-1]-incoord):
                        coord=dictcoord[i]
                        id_dict=i
                    else:
                        coord=dictcoord[i-1]
                        id_dict=(i-1)
                break
        elif delta<0:
            if dictcoord[i]<incoord:
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

def getEvalDictFromYears(years,idlat,idlon,varDict):
    valList=[]
    for yearid in range(0,len(years)):
        #print varDict[yearid,lon,lat]
        valList.append(varDict[yearid,idlat,idlon])
    return valList
    
def readLocationsFromFile(filePath):
    infile=open(filePath)
    locations={}
    for line in infile.readlines():
        line=line.split(",")
        locations[line[0]]=[float(line[1]),float(line[2])]
        
    infile.close()
    return locations

def getErrs(varvals,vals):
    errs=[]
    for i in range(0,len(vals)):
        errs.append(abs(vals[i]-varvals[i]))
    return errs

def plotCurvesFromLocations(locations,windvals,windvals_up,windvals_down,years,lats,lons,logx=True,logy=False,plotErrBand=False):
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
    color_id=0
    for name in names:
        lat=locations[name][0]
        lon=locations[name][1]
        idlat=getDictCoord(lat,lats)
        idlon=getDictCoord(lon,lons)
        y_vals=getEvalDictFromYears(years,idlat,idlon,windvals)
        y_vals_up=getEvalDictFromYears(years,idlat,idlon,windvals_up)
        y_vals_down=getEvalDictFromYears(years,idlat,idlon,windvals_down)
        #print x_vals,y_vals
        plt.plot(x_vals,y_vals,linewidth=0.7,linestyle="-",color=colors[color_id],label=name)
        if plotErrBand:
            plt.fill_between(x_vals, y_vals_down, y_vals_up,alpha=0.3, edgecolor=colors[color_id], facecolor=colors[color_id])
            #plt.errorbar(x_vals, y_vals,[getErrs(y_vals_up,y_vals),getErrs(y_vals_down,y_vals)],marker='s',fmt='o',color='black')
        color_id+=1
    plt.legend(loc=legendPosition,prop={'size': legendTextSize},fancybox=True, shadow=True)
    plt.xlabel(r"""$A\~nos$""", fontsize=20)
    plt.ylabel("Velocidad del viento (m/s)", fontsize=20)
    
    return plot

def summaryHistogram(locations,windvals,windvals_up,windvals_down,years,lats,locs,plotHisErrs=False):
    plot=plt.figure()
    ax = plot.add_subplot(111)
    noErrProblem=True
    warnings=[]
    x_vals=[]
    x_vals_labels=[]
    x_vals_middle=[]
    delta=1.0
    
    for x in range(1,len(locations)+1):
        x_vals.append(x*delta)
        x_vals_labels.append(x*delta+delta/2.5)
        x_vals_middle.append(x*delta+delta/2.5)
        
    years_inv=[]
    for i in range(0,len(years)):
        years_inv.append(int(years[len(years)-i-1]))
        
    year_id=len(years)-1
    color_id=0
    for year in years_inv:
        if year in years_forHist:
            y_vals=[]
            err_y_up=[]
            err_y_down=[]
            for name in locations:
                lat=locations[name][0]
                lon=locations[name][1]
                idlat=getDictCoord(lat,lats)
                idlon=getDictCoord(lon,lons)
                wind=windvals[year_id,idlat,idlon]
                y_vals.append(wind)
                err_y_temp=abs(windvals_up[year_id,idlat,idlon]-wind)
                if (err_y_temp/wind>0.8):
                    noErrProblem=False
                    messg="problem detected in errors!!. Found relative err > 80% !!!!. This can happen when there is no error in the .nc file!!! Error bars won't be plotted!"
                    if (messg not in warnings) and plotHisErrs==True:
                        warnings.append(messg)
                        
                err_y_up.append(abs(windvals_up[year_id,idlat,idlon]-wind))
                err_y_down.append(abs(windvals_down[year_id,idlat,idlon]-wind))
                
            plt.bar(x_vals, y_vals, color=colors[color_id],alpha=alphaHist,label=str(year))
            if plotHisErrs and noErrProblem:
                plt.errorbar(x_vals_middle, y_vals,[err_y_up,err_y_down],marker='s',fmt='o',markersize=2.0,color=colors[color_id],markeredgecolor=colors[color_id],capsize=6,elinewidth=1.0,capthick=1.0)
            color_id+=1
        year_id-=1
    plt.ylabel("Velocidad del viento (m/s)", fontsize=20)
    plt.xticks(x_vals_labels, locations, rotation='vertical')
    plt.legend(title=r"""$A\~nos$""",loc='upper center',prop={'size': legendTextSize_Hist},bbox_to_anchor=(1.02, 1.05),fancybox=True, shadow=True)
    plt.subplots_adjust(bottom=down_margin)
    print warnings
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

#####################################################################################################################################################################
################################################## ################## Parte principal del codigo ################# ##################################################
#####################################################################################################################################################################
## Leer data desde .nc
inputFile = Dataset(args[0], mode='r')

## importar datos de los lugares desde tu .csv
locations=readLocationsFromFile(options.inLoc)

## Informacion sobre la estructura del .nc y las variables que contiene
nctools.ncFileInfo(args[0])

##Variable de lons y lats, descomentar si las necesita para algo pero no creo que las necesites porque esos van a ser tus datos de entrada de los lugares
lons = inputFile.variables['lon']
lats = inputFile.variables['lat']

##years
years_nc=inputFile.variables['years']
years=[]
for year in years_nc:
    years.append(int(year))

## Diccionarios de la variables: dict[lat,lon]=valor
locs = inputFile.variables['loc']
scales = inputFile.variables['scale']
shps = inputFile.variables['shp']
#valores del viento dict[years,lon,lat]=valor
windvals=inputFile.variables['wspd']
windvals_up=inputFile.variables['wspdupper']
windvals_down=inputFile.variables['wspdlower']
    
##Obetener unidades de medida 
#dict["var"]=unidad
units={}
for var in ['lon','lat','loc','scale','shp','wspd','wspdupper','wspdlower']:
    units[var]=inputFile.variables[var].units
    print ("Found variable :"+var+" in ["+units[var]+"]")

#########********************* Ficheros ********************##############
##Obetener los datos para cada lugar y plotear las curvas
#Todas las curvas en un solo plot
plot_logx_all=plotCurvesFromLocations(locations,windvals,windvals_up,windvals_down,years,lats,lons)
savePlot(plot_logx_all,"curves_logx_all",options.outDir)
#Todas las curvas en un solo plot (con err band)
plot_logx_err_all=plotCurvesFromLocations(locations,windvals,windvals_up,windvals_down,years,lats,lons,plotErrBand=True)
savePlot(plot_logx_err_all,"curves_logx_errband_all",options.outDir)
## crear un plot de curva para cada location:
#Extraer de locations solo los nombres
names=list(locations)

for name in names:
    temp_locationDict={}
    temp_locationDict[name]=locations[name]
    #temp_curv_plot=plotCurvesFromLocations(temp_locationDict,windvals,windvals_up,windvals_down,years,lats,lons)
    #savePlot(temp_curv_plot,"curves_logx_"+name,options.outDir)
    temp_curv_plot_err=plotCurvesFromLocations(temp_locationDict,windvals,windvals_up,windvals_down,years,lats,lons,plotErrBand=True)
    savePlot(temp_curv_plot_err,"curves_logx_errband"+name,options.outDir)
    
# crear histograma sin errores
plot_summaryHist=summaryHistogram(locations,windvals,windvals_up,windvals_down,years,lats,locs)
savePlot(plot_summaryHist,"summaryHist",options.outDir)

# crear histograma con errores
plot_summaryHist_err=summaryHistogram(locations,windvals,windvals_up,windvals_down,years,lats,locs,plotHisErrs=True)
savePlot(plot_summaryHist_err,"summaryHist_withErrs",options.outDir)




#########********************* Ficheros ********************##############
## sacar informacion del .nc para un fichero
infoFileName="locationsInfo.csv"
infoFile=open(options.outDir+"/"+infoFileName,'w')
loc_unit="["+units["loc"]+"]"
scale_unit="["+units["scale"]+"]"
shp_unit="["+units["shp"]+"]"
infoFile.write("nombre,lat,lon,loc,scale,shp\n")

## Hallar y guardar los valores para cada lugar
for name in names:
    ##Obtener coordenadas del fichero del lugar
    lat=locations[name][0]
    lon=locations[name][1]
    ##Obtener coordenadas del fichero del lugar
    idlat=getDictCoord(lat,lats)
    idlon=getDictCoord(lon,lons)
    ## Escribir en el fichero      latitud y longitud en el .nc                                      parametros de la funcion
    infoFile.write(name+", "+str(lats[idlat])+", "+str(lons[idlon])+", "+str(locs[idlat,idlon])+", "+str(scales[idlat,idlon])+", "+str(shps[idlat,idlon])+"\n")
infoFile.close()

### id para los a;os #### Por ejemplo si los a;os son [2,5,...] El id de 2yr->0;5yr->1;... El id es lo que se pone para llamar al valor de wspd. 
#Por ejemplo si quieres los valores de viento para 2 a;os seria: windvals[0,valor_lat,valor_lon]. Por eso el contador year_id aumenta al final de cada ciclo.
#O sea la variable year_id contiene la posicion en el arreglo de year en years ----->>>>> [2(year_id=0),5(year_id=1),...]
year_id=0
for year in years:
    ##fichero para la informacion del viento
    infoFileName="locationsInfoWind_"+str(int(year))+"yrs.csv"
    infoFile=open(options.outDir+"/"+infoFileName,'w')
    print ("Creating file: "+ str(infoFile).split(" ")[2])
    wind_unit="["+units["wspd"]+"]"
    infoFile.write("nombre,wind"+wind_unit+",wind_up"+wind_unit+",wind_down"+wind_unit+"\n")
    for name in names:
        ##Obtener coordenadas del fichero del lugar
        lat=locations[name][0]
        lon=locations[name][1]
       ##Obtener coordenadas del fichero del lugar
        idlat=getDictCoord(lat,lats)
        idlon=getDictCoord(lon,lons)
        ## Escribir en el fichero               viento                              viento up                                    viento down
        infoFile.write(name+", "+str(windvals[year_id,idlat,idlon])+", "+str(windvals_up[year_id,idlat,idlon])+", "+str(windvals_down[year_id,idlat,idlon])+"\n")
    infoFile.close()
    year_id+=1





