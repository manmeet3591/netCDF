import sys, string
from matplotlib import rc
import numpy
import pylab as pl
import netCDF4
import time as t
import datetime
from dateutil.parser import parse
from pylab import load, meshgrid, title, arange, show
from netcdftime import utime
import scipy.io
import matplotlib as mpl
import argparse
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import datetime as dt
import calendar
import Tkinter as tk

list_evaluate = []
global list_evaluate

def sende1():
    fileLocation = e.get()
    status["text"] = ""
    nc_data = netCDF4.Dataset(fileLocation)
    list_variables = []
    for i in nc_data.variables:
        list_variables.append(str(i))
    for i in range(len(list_variables)):
        #status["text"] = status.cget("text") + str(list_variables[i]) + "\n"
        listbox.insert(tk.END, list_variables[i])
    root.after(1)
    status["text"] = "Welche Variable sollen ausgewertet werden?"
    listbox.grid()

def onselect(evt):
    # Note here that Tkinter passes an event object to onselect()
    w = evt.widget
    index = int(w.curselection()[0])
    value = w.get(index)
    print 'You selected item %d: "%s"' % (index, value)
    list_evaluate.append(value)

def start():
    preData = open('../scripts/test.py','w')
    preData.write('import sys, string \n'+'from matplotlib import rc\n'+'import numpy\n'+'import pylab as pl\n'+'import netCDF4\n'+'import time as t\n'+'import datetime\n'+'from pylab import load, meshgrid, title, arange, show\n'+'from netcdftime import utime\n'+'import scipy.io\n'+'import matplotlib as mpl\n')
    preData.write('nc_data = netCDF4.Dataset('+'/silos/boergel/BREC/files/mean/flux3d_means.nc'+')\n')
    for i in range(len(list_evaluate)):
                preData.write(str(list_evaluate[i])+ " = netCDF4.Dataset['"+ str(list_evaluate[i])+"']\n")
    print list_evaluate
    preData.close()


#### Apperance ####
root = tk.Tk()
status = tk.Label(root, text='Entries')
status.grid()
b = tk.Button(root, text="Get all variables", command=sende1)
c = tk.Button(root, text='Start', command=start)
e = tk.Entry(root)
b.grid()
c.grid()
e.grid()
listbox = tk.Listbox(root)
listbox.bind('<<ListboxSelect>>', onselect)
############################
#### selection ####
status["text"] = listbox

root.mainloop()
