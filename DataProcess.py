#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 07:19:56 2024

@author: tze
"""
from numpy import array,shape, arange
from scipy.interpolate import interp1d

prefixes={'':1,'d':1e-1, 'c':1e-2, 'm':1e-3,'μ':1e-6,'n':1e-9,'p':1e-12,'f':1e-15, 'k':1e3, 'M':1e6}
times={'s':1,'min':60,'h':60*60,'day':24*60*60,'year':365*24*60*60}

def convert_time(timeunit,newtimeunit):
    return times[timeunit]/times[newtimeunit]

def convert_prefix(prefix='',newprefix='n',power=1):
	return prefixes[prefix]**power/prefixes[newprefix]**power

def convert_value(data,prefix,newprefix):
    return data*convert_prefix(prefix,newprefix)

#quantity per unit area
def convert_flux(topprefix='',newtopprefix='m',bottomprefix='',newbottomprefix='m'):
    return convert_prefix(topprefix,newtopprefix)/convert_prefix(bottomprefix,newbottomprefix,2)
def convert_ratio(topprefix='',newtopprefix='m',bottomprefix='',newbottomprefix='m'):
    return convert_prefix(topprefix,newtopprefix)/convert_prefix(bottomprefix,newbottomprefix)

def convert_speed(speed,oldspeed='m/s',newspeed='km/h'):
    [oldlength,oldtime]=oldspeed.split('/')
    [newlength,newtime]=newspeed.split('/')
    return speed*convert_prefix(oldlength[0:-1],newlength[0:-1])/convert_time(oldtime,newtime)
    
#mutable approach
#think about keyword prefix or unit_prefix
def convert_unit(data,newunit):
    data['#data_table'][:,0]=convert_value(data['#data_table'][:,0], data['#data_summary']['x1_unit'][0:-1], newunit[0:-1])
    data['#data_summary']['x1_unit']=newunit


#number of points that you are averaging is 2*n+1 if len(array)-n>i>n
#n is number of surrounding points (left and right) that you are using for averaging
def mov_average(inarray,n):
    #check otherways of smoothing maybe Fourier transform high frequency filter?
    avg_array=[];
    for i in range(0,len(inarray)):
        imin=max(0,i-n)
        imax=min(len(inarray),i+n+1)
        avg_array.append(sum(array[imin:imax]/(imax-imin)))
    return array(avg_array)

#it is expected you get [[xi,yi]] and that wavelegth is with same unit
def numply_at_same_x(A,B):
    xmin=max(min(A[:,0]),min(B[:,0]))
    xmax=min(max(A[:,0]),max(B[:,0]))
    xstep=min(A[:,0][1]-A[:,0][0],B[:,0][1]-B[:,0][0])
    afit = interp1d(A[:,0], A[:,1])
    bfit = interp1d(B[:,0], B[:,1])
    x=arange(xmin,xmax+xstep,xstep)
    y=afit(x)*bfit(x)
    out=array((len(x),2))
    out[:,0]=x
    out[:,1]=y
    return out
#here you work with your dictionary object
def multiply_at_same_x(A,B):
    out={}
    convert_wavelength(A,'nm')
    convert_wavelength(B,'nm')
    out['#data_table']=numply_at_same_x(A['#data_table'],B['#data_table'])
    #you need deep copy
    out['#data_summary']={}
    for keyword,item in A['#data_summary'].items():
        out['#data_summary'][keyword]=item
    out['#data_summary']['tot_col'],out['#data_summary']['tot_row']=shape(out['#data_table'])
    return out

def absolute_reflectance(Data, Rel_reference, Abs_reference):
    return multiply_at_same_x(multiply_at_same_x(Data, Rel_reference),Abs_reference)

"""    

def absolute_reflectance(Abs, Rel, Data):#Abs, Rel, Data [#data_table] col1 wavelenght col2 values
	wmin=max(min(Abs['#data_table'][:,0]),min(Rel['#data_table'][:,0]),min(Data['#data_table'][:,0]))
	wmax=min(max(Abs['#data_table'][:,0]),max(Rel['#data_table'][:,0]),max(Data['#data_table'][:,0]))
    wstep=min(Abs.wlength[1]-Abs.wlength[0],Rel.wlength[1]-Rel.wlength[0])
    afit = ipl.interp1d(A.wlength, A.data)
    bfit = ipl.interp1d(B.wlength, B.data)
    out.wlength=np.arange(wmin,wmax+wstep,wstep)
    out.data=afit(out.wlength)*bfit(out.wlength)
    return out

def average_curves(self,data):#data.list in which A.wlength A.data
        wstep=0
        wmin=0
        wmax=np.inf
        out=container()
        for item in data:
            wmin=max(min(item.wlength),wmin)
            wmax=min(max(item.wlength),wmax)
            wstep=max((item.wlength[1]-item.wlength[0]),wstep)
        out.wlength=np.arange(wmin,wmax+wstep,wstep)
        sum=np.zeros(len(out.wlength))
        for item in data:
            fititem=ipl.interp1d(item.wlength, item.data)
            sum=sum+fititem(out.wlength)
        out.data=sum/len(data)
        out.wlength_units=data[0].wlength_units
        return out
        #interpolate reference to measured data find Si Sio2 data from 190 to 1100
        #measured.wlegnthmin:measured.wlength_step:measured.wlength_max #see in TMM
        #after that just multiply the data
"""