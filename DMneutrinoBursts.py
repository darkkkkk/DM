#!/usr/bin/python
import numpy as np
from scipy.integrate import odeint
import math as m
#from lib.DMcapture import F_captGould,F_captBert,F_captBert2,F_selfCapt
from lib.DMcapture import F_capture
from lib.functions import *
from lib.unitconv import *
from lib.classes import *
from lib.constants import *

def calcALL(sigx, anncs):

    print " ********************************** "
    print "               INPUT " 
    print " ********************************** "
    # Input and initialize dark matter:
    dm = cl_dm(sigx, anncs)
    dm.prints()

    # Input star:
    st = cl_star()
    st.prints()

    print " ********************************** "
    print "            CALCULATIONS "
    print " ********************************** "
    for i in xrange(len(dm.mx_v)):
	  print " ********************************** "
	  print "          mx = %.2e " % dm.mx_v[i]
	  print " ********************************** "
	  mr = st.mn*dm.mxkg_v[i]/(st.mn+dm.mxkg_v[i])	# reduced mass, kg

          # Initializes lists that evolve with time:
          t = cl_timevars()

	  # check bullet cluster constraints:
          F_checksigxxconstraints(i, dm)

          # Characteristic time scales
          t1, tth = F_calctimescales(i, st, dm)

          # Characteristic radius of dark matter distribution, sets t.rx:
          rxco, rxth = F_calcDMradius(i, t, st, dm, t1, tth)

	  # Characteristic DM Numbers:
          Nsg,Nsgvect,Nfe,Nde,Ngeoxx,Ngeoxxvect,Nbec,Nchadeg = F_calcDMnumbers(i, t, st, dm, rxth)

          # Capture rates calculation:
          capt, selfcapt = F_capture(i, st, dm)

	  # Calc Nx(t) numerically (solves dN/dt = capt + selfCapt*N(t) - ann*(N(t)**2)):
          t.Nx, Nx_C_A, Nx_C_sC_A = F_calcNxnum(t, dm, capt, selfcapt)

          # Check if numerical solution agrees with analytical:
#          F_checkNxanalytical(t, dm, capt, rxth, Nx_C_A)

	  # Check if self-collapse is produced, if geometrical limit in selfcapture is reached, if degenerate core and Bose-Einstein condensate are formed:
          SELFCOLL, kscoll, kgeoxx, kbec  = F_checkSelfColl_Deg_Bec(i, t, dm, t.Nx, Nsgvect, Nde, Nfe, Ngeoxxvect, Nbec)

	  # Time when self-collapse starts and time of BEC formation:
          tNsg, tNbec = F_calctimeSelfColl_Bec(t, dm, t.Nx, Nsg, Nbec, kbec)

	  # In the special case when BEC is formed, calculates rxbec and annbec, and updates t.rx :
	  if dm.BEC: rxbec, annbec = F_calcBECcharact(i, t, dm, kbec)

	  # Recalcs Nx_C_sC_A(t) if selfcapture geoxx is reached (NOT CORRECT)
#          if dm.SELFCAPT:
#              Nx_C_sC_A = F_recalcNx_C_sC_A_with_geoxx(t, dm, capt, selfcapt, Nx_C_sC_A, kgeoxx )
#              t.Nx = Nx_C_sC_A

          # Updates kbec and t.rx with rxbec (may have changed with new Nx_C_sC_A)
          if dm.BEC and dm.SELFCAPT: kbec = F_recalckbec(i, t, dm, Nx_C_sC_A, Nbec, rxbec)

	  # Calcs Nx_C_A_BEC(t) (capture and annihilation, with BEC after t>tNbec):
          if dm.BEC and not dm.SELFCAPT:
              Nx_C_A_BEC = F_calcNx_C_A_BEC(t, dm, capt, kbec, annbec)
              t.Nx = Nx_C_A_BEC

          # Calcs Nx_C_sC_A_BEC(t) (capture, selfcapture and annihilation, with BEC after for t>tNbec):
          if dm.BEC and dm.SELFCAPT:
              Nx_C_sC_A_BEC = F_calcNx_C_sC_A_BEC(t, dm, Nx_C_sC_A, capt, selfcapt, kgeoxx, kbec, annbec)
              t.Nx = Nx_C_sC_A_BEC

	  # Calcs DM annihilation rate:
	  for k in xrange(len(t.time)):
              try:
                  t.annirate[k] = F_ann(dm.anncs,t.rx[k]) * m.pow(t.Nx[k],2.)
              except OverflowError:
                  t.annirate[k] = 9.e+99

          # Recheck if self-collapse is produced with new Nx(t):
          if not (not dm.SELFCAPT and not dm.BEC): SELFCOLL, kscoll = F_recheckSelfColl(i, t, dm, t.Nx, Nsgvect, Nde, Nfe)

	  # Computes the collapse of the DM cloud:
          if SELFCOLL: timecoll, rxcoll ,vxcoll, Nxcoll, annrate, tff, kend = F_computeCollapse(i, t, dm, t.Nx, capt, kscoll)

          # Prints time evolution of parameters
          print" time(s)   time(yr)     Nx        Nsg      rx(type)          ann      flux(GeV/m2/yr)"
          for k in xrange(0,kscoll+1):
              print "%.2e   %.2e   %.2e   %.2e   %.2e (%s)   %.2e   " % (t.time[k],c_s2yr(t.time[k]),t.Nx[k],Nsgvect[k],t.rx[k],t.rxtag[k],t.annirate[k])
	  if SELFCOLL:
	      print" now DM self collapses in a timescale < %.2e s" % tff
	      print " t         rx(cm)       vx/c       Nx      annrate(part/s)"
              for k in xrange(0,kend+1):
                  print "%.2e   %.2e   %.2e   %.2e   %.2e" % (timecoll[k],rxcoll[k],vxcoll[k]/clight,Nxcoll[k],annrate[k])
   
          if SELFCOLL: 
              # Flux of energy emitted:
              dm.FluxEmit_v[i]  = dm.Nsg_v[i] * dm.mx_v[i] / c_s2yr(dm.tNsg_v[i])		# GeV/yr
              # Flux of neutrinos at Earth:
              dm.FluxEarth_v[i] = dm.FluxEmit_v[i] /4./m.pi/m.pow(c_pc2m(st.dist),2.)    	# GeV/m2/yr
              dm.Capt_v[i] = capt


    # OUTPUT:
    dm.printOutput()

def main():
    for sigx in [1.e-30, 1.e-32, 1.e-34, 1.e-36, 1.e-38, 1.e-40, 1.e-42, 1.e-44, 1.e-48, 1.e-50 ]:
        for anncs in [1.e-26, 1.e-30, 1.e-35, 1.e-40, 1.e-45, 1.e-50, 1.e-55, 1.e-60]:
            calcALL(sigx, anncs)
  
if __name__ == "__main__":
    main()
