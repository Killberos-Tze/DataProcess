#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 07:19:56 2024

@author: tze
"""

class Conversion:
    def __init__(self):
        self.prefixes={'':1,'d':1e-1, 'c':1e-2, 'm':1e-3,'μ':1e-6,'n':1e-9,'p':1e-12,'f':1e-15, 'k':1e3, 'M':1e6}
    
    def convert_prefix(self,prefix='',newprefix='n',power='1'):
        return self.prefixes[prefix]**power/self.prefixes[newprefix]**power
    
    def convert_flux(self,topprefix='',newtopprefix='m',bottomprefix='',newbottomprefix='m'):
        return self.convert_prefix(topprefix,newtopprefix)/self.convert_prefix(bottomprefix,newbottomprefix,2)
    