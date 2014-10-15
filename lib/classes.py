import sys
import numpy as np
from lib.unitconv import *
from lib.constants import mp, kb

class cl_star():
        def __init__(self,passvar=[]):
                self.Ms = c_Msun2kg(1.44)    #1.4-0.6 # from Msun to kg
#               Rs = c_Rsun2m(0.0045)   #0.0045-0.02 # from Rsun to m
                self.Rs = 10600.     # m
                self.Temp = 1.e+5    # K
                self.densStars = 5.e-4       # NumNS/pc3
                self.rhoc = self.Ms*1.e+3/ (4./3.*m.pi*m.pow(self.Rs*1.e+2,2.)) # g/cm3
                self.vs = 22000.    # m/s
                self.An = 1
                self.mn = self.An * mp    # kg
                self.nb = self.rhoc*1e-3 / self.mn # numb.bary / cm3
#                self.pF = 426.e+6 * 5.344286e-28 # from eV to kg m/s, Fermi momentum 
                self.Eth= (3./2.)*kb*self.Temp
                self.dist = 200	# pc 

        def prints(self):
                print "--- Stellar characteristics"
                print "    Ms = %.2f Msun,  Rs = %.5f Rsun,  rhoc = %.2e g/cm3,  Temp = %.2e K,   An = %i,  vs = %.2e m/s "  % (c_kg2Msun(self.Ms),c_m2Rsun(self.Rs),self.rhoc,self.Temp,self.An,self.vs)

class cl_dm():
	def __init__(self, sigx, anncs):
	# inputs:
        	self.rhox  	= 0.4     					# GeV/cm3
#       	self.mx_v 	= [ 10.**i for i in range(4,8+1)]  		# GeV
        	self.mx_v	= [ 10.**(2*i-2) for i in range(4,4+1)]# (1,3+1)]  	# GeV
                self.mxkg_v     = [ c_GeV2kg(mx) for mx in self.mx_v ]    				# kg
                self.mxg_v  	= [ c_GeV2g(mx) for mx in self.mx_v ]      				# g
        	self.anncs 	= anncs # 1.e-50  						# cm3/s
        	self.sigx  	= sigx    ; self.sigx_m  = self.sigx/1.e+4  	# m2
        	self.sigxx 	= 1.e-30  ; self.sigxx_m = self.sigxx/1.e+4    	# m2
        	self.vx    	= 27000.  					# m/s
        	self.SELFCAPT 	= False
        	self.BOSON    	= False
		self.DEG 	= False		# ?????????????????????
        	if self.BOSON==False:
          		self.BEC = False
	# initializing other variables:
		self.Nsg_v 	= np.zeros(len(self.mx_v))
		self.tNsg_v	= np.zeros(len(self.mx_v))
        	self.notes_v	= ['' for k in range(0,len(self.mx_v))]
#		self.Flux_v	= np.zeros(len(self.mx_v))
		self.FluxEmit_v = np.zeros(len(self.mx_v))
                self.FluxEarth_v= np.zeros(len(self.mx_v))
                self.Capt_v	= np.zeros(len(self.mx_v))
        # more:
		self.Nxinit = 0.	# np.array([0.])                   # initial values N(t=0)
	
	def prints(self):
		print "--- Dark matter:"
        	print "    rhox = %.2e GeV/cm3,  anncs = %.2e cm3/s,  sigx = %.2e cm2,  vx = %.2e m/s"   % (self.rhox,self.anncs,self.sigx,self.vx)
        	print "    selfinteractions:",self.SELFCAPT," (sigxx = %.2e cm2)" % self.sigxx
        	print "    bosonic DM:",self.BOSON

        def printOutput(self):
            print ""
            print "# DM collapse and destruction of NS or self-annihilation burst:"
            print "# mx(GeV)  period(yr)  Nx(collapse)  Flux@Earth(GeV/m2/yr)  notes"
            for it1,it2,it3,it4,it5 in zip(self.mx_v,self.tNsg_v,self.Nsg_v,self.FluxEarth_v,self.notes_v):
                print " %.2e  %.2e    %.2e      %.2e              %s" % (it1,c_s2yr(it2),it3,it4,it5)
            print """\nNotes:
            1:      Nsg not reached within the age of the Universe  
            2:      Degenerate Dark Star will be formed (Nsg>Nde and >Nfe) only if 4 holds
            4:      Degenerate Dark Star definitely formed (Nsg>Nde and >Nfe)
            5:      Ngeoxx was reached, and Nx_C_sC_A_geoxx is not correct
            6:      Nbec   was reached, check if BEC was correctly taken into account
            7:      BEC occurs before geoxx limit, check if implementation is correct in this case
             """
            f = open('NeutrinoBursts.dat', 'a')
            for it1,it2,it3,it4,it5,it6 in zip(self.mx_v, self.Nsg_v, self.tNsg_v, self.FluxEmit_v, self.FluxEarth_v, self.notes_v):
#                f.write("%.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %s\n" % (self.sigx, self.anncs, 
#                                                         it1, it2, c_s2yr(it3), it4, it5, it6) )
                f.write("%.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e\n" % (self.sigx, self.anncs, it1, it2, c_s2yr(it3), it4, it5) )
            f.close()


class cl_timevars():
	def __init__(self):
		self.time     	= [ (10.**(k/3.))-1. for k in range(0,54)]  # s
#       	time  = np.linspace(0.,5.e+12,20)     # numpy.linspace(start, stop, num=50, endpoint=True, retstep=False)
        	self.rx		= np.zeros(len(self.time))
        	self.rxtag    	= ['yes' for k in range(0,len(self.time))]
        	self.ann	= np.zeros(len(self.time))
        	self.annirate 	= np.zeros(len(self.time))
		self.Nx		= np.zeros(len(self.time))
		self.Nsg	= np.zeros(len(self.time))
		self.Nx22     	= np.zeros(len(self.time))
        	self.Nx3      	= np.zeros(len(self.time))
        	self.Nx31     	= np.zeros(len(self.time))
        	self.Nx32     	= np.zeros(len(self.time))
