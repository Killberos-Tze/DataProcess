#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 07:19:56 2024

@author: tze
"""
from numpy import empty,shape, arange, array
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
def convert_unit_IHTM(data,newprefix,column):
    if data['#data_summary'][f'{column}_col']=='#data_table':
        data['#data_table']=convert_value(data['#data_table'], data['#data_summary'][f'{column}_prefix'], newprefix)
    else:
        data['#data_table'][:,data['#data_summary'][f"{column}_col"]]=convert_value(data['#data_table'][:,data['#data_summary'][f"{column}_col"]], data['#data_summary'][f'{column}_prefix'], newprefix)
    data['#data_summary'][f'{column}_prefix']=newprefix

#number of points that you are averaging is 2*n+1 if len(array)-n>i>n
#n is number of surrounding points (left and right) that you are using for averaging
def mov_average(inarray,n):
    #check otherways of smoothing maybe Fourier transform high frequency filter?
    avg_array=[];
    for i in range(0,len(inarray)):
        imin=max(0,i-n)
        imax=min(len(inarray),i+n+1)
        avg_array.append(sum(inarray[imin:imax])/(imax-imin))
    return array(avg_array)

def average_IHTM(A,n,column):#n points to average, column to average
    out=copy_IHTM(A)
    out['#data_table'][:,out["#data_summary"][f"{column}_col"]]=mov_average(out['#data_table'][:,out["#data_summary"][f"{column}_col"]],n)
    out['#data_summary'][f"{column}_label"]=out['#data_summary'][f"{column}_label"]+'_avg'
    out['#data_summary'][f"{column}_smooth_points"]=n
    return out

#it is expected you get [[xi,yi]] and that wavelegth is with same unit
def numply_at_same_x(A,B):
    xmin=max(min(A[:,0]),min(B[:,0]))
    xmax=min(max(A[:,0]),max(B[:,0]))
    xstep=min(A[:,0][1]-A[:,0][0],B[:,0][1]-B[:,0][0])
    afit = interp1d(A[:,0], A[:,1])
    bfit = interp1d(B[:,0], B[:,1])
    x=arange(xmin,xmax+xstep,xstep)
    y=afit(x)*bfit(x)
    out=empty((len(x),2))
    out[:,0]=x
    out[:,1]=y
    return out

def absorb_at_same_x(A,B):
    xmin=max(min(A[:,0]),min(B[:,0]))
    xmax=min(max(A[:,0]),max(B[:,0]))
    xstep=min(A[:,0][1]-A[:,0][0],B[:,0][1]-B[:,0][0])
    afit = interp1d(A[:,0], A[:,1])
    bfit = interp1d(B[:,0], B[:,1])
    x=arange(xmin,xmax+xstep,xstep)
    y=1-afit(x)-bfit(x)
    out=empty((len(x),2))
    out[:,0]=x
    out[:,1]=y
    return out

def numdiv_at_same_x(A,B):
    xmin=max(min(A[:,0]),min(B[:,0]))
    xmax=min(max(A[:,0]),max(B[:,0]))
    xstep=min(A[:,0][1]-A[:,0][0],B[:,0][1]-B[:,0][0])
    afit = interp1d(A[:,0], A[:,1])
    bfit = interp1d(B[:,0], B[:,1])
    x=arange(xmin,xmax+xstep,xstep)
    y=afit(x)/bfit(x)
    out=empty((len(x),2))
    out[:,0]=x
    out[:,1]=y
    return out
#here you work with your dictionary object
def multiply_2col_IHTM(A,B):#ihtm A and B have to have same x units and same y units
    out=copy_IHTM(A)
    out['error']=''
    if A['#data_summary']['tot_col']!=B['#data_summary']['tot_col'] and B['#data_summary']['tot_col']!=2:
        out['error']='wrong number of columns'
    if A['#data_summary']['x1_prefix']!=B['#data_summary']['x1_prefix'] or A['#data_summary']['y1_prefix']!=B['#data_summary']['y1_prefix']:
        out['error']='prefixes do not match'
    if A['#data_summary']['x1_prefix']==B['#data_summary']['x1_prefix'] and A['#data_summary']['y1_prefix']==B['#data_summary']['y1_prefix']: 
        out['#data_table']=numply_at_same_x(A['#data_table'],B['#data_table'])
        out['#data_summary']['tot_col'],out['#data_summary']['tot_row']=shape(out['#data_table'])
    return out

def divide_2col_IHTM(A,B):
    out=copy_IHTM(A)
    out['error']=''
    if A['#data_summary']['tot_col']!=B['#data_summary']['tot_col'] and B['#data_summary']['tot_col']!=2:
        out['error']='wrong number of columns'
    if A['#data_summary']['x1_prefix']!=B['#data_summary']['x1_prefix'] or A['#data_summary']['y1_prefix']!=B['#data_summary']['y1_prefix']:
        out['error']='prefixes do not match'
    if A['#data_summary']['x1_prefix']==B['#data_summary']['x1_prefix'] and A['#data_summary']['y1_prefix']==B['#data_summary']['y1_prefix']: 
        out['#data_table']=numdiv_at_same_x(A['#data_table'],B['#data_table'])
        out['#data_summary']['tot_col'],out['#data_summary']['tot_row']=shape(out['#data_table'])
    return out

def copy_IHTM(A):
    out={}
    out['#data_table']=array(A['#data_table'])
    out['#data_summary']={}
    for keyword,item in A['#data_summary'].items():
        out['#data_summary'][keyword]=item
    out['#data_summary']['tot_col'],out['#data_summary']['tot_row']=shape(out['#data_table'])
    return out

def absolute_reflectance_IHTM(Data, Rel_reference, Abs_reference):
    return multiply_2col_IHTM(divide_2col_IHTM(Data, Rel_reference),Abs_reference)

def absorbance_IHTM(A,B):
    out=copy_IHTM(A)
    out['error']=''
    if A['#data_summary']['tot_col']!=B['#data_summary']['tot_col'] and B['#data_summary']['tot_col']!=2:
        out['error']='wrong number of columns'
    if A['#data_summary']['x1_prefix']!=B['#data_summary']['x1_prefix'] or A['#data_summary']['y1_prefix']!=B['#data_summary']['y1_prefix']:
        out['error']='prefixes do not match'
    if A['#data_summary']['x1_prefix']==B['#data_summary']['x1_prefix'] and A['#data_summary']['y1_prefix']==B['#data_summary']['y1_prefix']:
        out['#data_table']=absorb_at_same_x(A['#data_table'],B['#data_table'])
        out['#data_summary']['tot_col'],out['#data_summary']['tot_row']=shape(out['#data_table'])
        out['#data_summary']['y1_name']='Absorbance'
        if 'R01' in out['#data_summary']['y1_label']:
            out['#data_summary']['y1_label']=out['#data_summary']['y1_label'].replace('R01','A')
        elif 'T01' in out['#data_summary']['y1_label']:
            out['#data_summary']['y1_label']=out['#data_summary']['y1_label'].replace('T01','A')
        else:
            out['#data_summary']['y1_label']=out['#data_summary']['y1_label']+'_A'
    return out