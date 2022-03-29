# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 19:05:09 2014

 Solve pressure screen given Inlet pressure, Dilution Pressure, Valve Positions,
 Inlet consistency, and Rejects Thickening Factor

@author: joshua
"""
import numpy as np
from scipy.optimize import minimize

#
# Objective function for screen numerical solution
#
def objScreen(F,Fin,CV,P):
    """ F = Flow 0, 1 , and 2 from mesh network full screen model,
    Fin = Inlet flow to screen, EU = GPM...    
    CV =  N/A, Pipe, Reject Internal, Reject Valve,...
    ...Screen, Accepts Valve, Dil Valve, EU = GPM/sqrt(PSI)
    P = Inlet, Dilution Pump, Rejects Tank Head,...
    ...Accepts Tank Head, EU = PSIA """
    for i in range(np.shape(CV)[0]):
        CV[i] = max(1.,CV[i])
    
    # LOOP 0
    X0 = (Fin)**2*CV[1]**-2 + (CV[4]**-2 + CV[5]**-2)*F[0]**2 -\
    P[0] + P[3]
    
    # LOOP 1
    X1 = (Fin)**2*CV[1]**-2 + F[1]**2*CV[2]**-2 + \
    (F[1] + F[2])**2*CV[3]**-2 - P[0] +P[2]
    
    # LOOP 2
    X2 = (F[1] + F[2])**2*CV[3]**-2 + F[2]**2*CV[6]**-2 - P[1] + P[2]
    
    # Inlet flow is known
    X3 = F[0] + F[1] - Fin
    
    return X0**2 + X1**2 + X2**2 + X3**2

#
# Mass balance equality constraint for screen numerical solution
#
def consMB(F):
    """ Densities from screen model lb/gal
    F = flows from mesh network """
    # Constants 
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    DEN = np.array([7.42, 8.03, 7.62, 7.49])
    
    # Mass balance
    eq = DEN[INLET]*(F[0] + F[1]) + DEN[DIL]*F[2] - DEN[ACCPT]*F[0] -\
    DEN[REJ]*(F[1] + F[2])
    
    return eq
    
#
# Pressure Screen Model Solution Given CV's, Pressures and Inlet Flow
#
def solveScreen(CV, P, Fin, F0):
    """ 
    F0 = Flow 0, 1 , and 2 from mesh network screen model initial guesses
    Fin = Inlet flow GPM
    CV =  N/A, Pipe, Reject Internal, Reject Valve,...
    ...Screen, Accepts Valve, Dil Valve, EU = GPM/sqrt(PSI)
    P = Inlet, Dilution Pump, Rejects Tank Head,...
    ...Accepts Tank Head, EU = PSIA
    """
    for i in range(np.shape(CV)[0]):
        CV[i] = max(1.,CV[i])
    
    # Build tuple of systems values
    SYS = (Fin,CV,P)
    
    # Create boundries 
    lim = (1.,10000.) # 1 to 10,000 gpm
    bnds = (lim,lim,lim) # set all the same
    
    # setup equality contraint
    #con1 = {'type': 'eq', 'fun':consMB} #constraints=con1,

    # solve nonlinear set of equations    
    sol = minimize(objScreen, F0, args=SYS, method='SLSQP',\
    bounds=bnds, options={'ftol':1e-6})
    print 'Screen Solve: ', sol.message, '  ' , sol.fun
    
    # Debug
    # print sol
    
    Inlet = (sol.x[0] + sol.x[1])
    Accepts = sol.x[0]
    Rejects = sol.x[1] + sol.x[2]
    Dilution = sol.x[2]
    
    # return flows in GPM
    return np.array([Inlet, Dilution, Rejects, Accepts])

#
# Consistency Model
#
def solveCons(F,Cf,RTF):
    """F = Inlet flow, Dilution Flow, Rejects flow, Accepts
    Cf = Feed Consistency in %, RTF = Rejects thickening factor unitless,
    Densities from screen model lb/gal
    """
    
    # Constants 
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    DEN = np.array([7.42, 8.03, 7.62, 7.49])
    
    # calculate mass flows
    Mf = max(1.,F[INLET]*DEN[INLET])
    Md = max(1.,F[DIL]*DEN[DIL])
    Mr = max(1.,F[REJ]*DEN[REJ])
    
    # Calc Consistencies
    Cf = Cf/100.   #convert from percent to decimal   
    # calculate rejects consistency base on rejects thickening factor
    Cr = RTF*Cf    
    Cr = min(Cf*Mf/Mr,Cr)
    #solve sequentially
    CaMa = Cf*Mf - Cr*Mr
    Ma = (1.-Cf)*Mf+Md-(1.-Cr)*Mr+CaMa
    Ca = CaMa/Ma
    
    # return Consistency for Accepts and Rejects in %
    return np.array([Ca*100.,Cr*100.])

#
# Calculate Rejects and Accepts Pressures given inlet..
# ..dilution pump pressures and flows
#
def solvePress(F,CV,P):
    """ F = Inlet, Dilution, Rejects, Accepts flow
    CV =  N/A, Pipe, Reject Internal, Reject Valve,...
    ...Screen, Accepts Valve, Dil Valve, EU = GPM/sqrt(PSI)
    P = Inlet, Dilution Pump, Rejects Tank Head,...
    ...Accepts Tank Head, EU = PSIA """
    for i in range(np.shape(CV)[0]):
        CV[i] = max(1.,CV[i])
    # Constants 
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    
    # Accepts Head
    P[ACCPT] = P[INLET] - (F[INLET])**2*CV[1]**-2 -\
    (CV[4]**-2 + CV[5]**-2)*F[REJ]**2
    
    # LOOP 2
    P[REJ] = P[DIL] - (F[REJ])**2*CV[3]**-2 - F[DIL]**2*CV[6]**-2

    return P

#
# Objective Function for Thevenin Equvilant Calcs
#
def objThev(F,CV,P):
    """
    F = flow 0 and 1 from Thevenin Model with sources energized
    CV =  N/A, Pipe, Reject Internal, Reject Valve,...
    ...Screen, Accepts Valve, Dil Valve, EU = GPM/sqrt(PSI)
    P = Inlet, Dilution Pump, Rejects Tank Head,...
    ...Accepts Tank Head, EU = PSIA"""
    for i in range(np.shape(CV)[0]):
        CV[i] = max(1.,CV[i])
    # Constants 
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    
    # LOOP 0
    X0 = P[REJ] - P[DIL] + F[0]**2/CV[6]**2 + (F[0] - F[1])**2/CV[3]
        
    # LOOP 1
    X1 = P[ACCPT] - P[REJ] + (F[1] - F[0])**2/CV[3]**2 + F[1]**2*(CV[2]**-2 + \
          CV[4]**-2 + CV[5]**-2)
    
    
    return X0**2 + X1**2

#
# Calculate Thevenin Equvilant Cv
#
def theveninCV(CV):
    """
    CV =  N/A, Pipe, Reject Internal, Reject Valve,...
    ...Screen, Accepts Valve, Dil Valve, EU = GPM/sqrt(PSI)
    """
    for i in range(np.shape(CV)[0]):
        CV[i] = max(1.,CV[i])    
    # CV3 and CV6 in parallel
    CV36 = CV[3]+CV[6]
    # CV36 in series with CV2
    CV236 = CV36*CV[2]/(CV36**2+CV[2]**2)**0.5
    # CV 4 and 5 in series
    CV45 = CV[4]*CV[5]/(CV[4]**2+CV[5]**2)**0.5
    # CV45 and CV236 in parallel
    CV23456 = CV45 + CV236
    # CV1 in series with CV23456
    Gth = CV[1]*CV23456/(CV[1]**2+CV23456**2)**0.5
    
    return Gth

#
# solve Thevenin Equvalant Circuit
#
def solveThevenin(CV,P,F0):
    """
    F0 = Flow 0, & 1 from mesh network for thevenin screen model
    ..with P0 removed from the circuit initial guesses
    CV =  N/A, Pipe, Reject Internal, Reject Valve,...
    ...Screen, Accepts Valve, Dil Valve, EU = GPM/sqrt(PSI)
    P = Inlet, Dilution Pump, Rejects Tank Head,...
    ...Accepts Tank Head, EU = PSIA  
    """
    for i in range(np.shape(CV)[0]):
        CV[i] = max(1.,CV[i])
        
    # Constants 
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    
    # Build tuple of systems values
    SYS = (CV,P)
    
    # Create boundries 
    lim = (1.,10000.) # 1 to 10,000 gpm
    bnds = (lim,lim) # set all the same

    # solve nonlinear set of equations    
    sol = minimize(objThev, F0, args=SYS, method='SLSQP',\
    bounds=bnds, options={'ftol':1e-6})
    print 'Thevenin Solve: ', sol.message, '  ' , sol.fun
    
    # Return Overall Thevenin Cv and Pressure   
    return theveninCV(CV), sol.x[1]**2*(CV[4]**-2+CV[5]**-2) + P[ACCPT]

#
#  Classes
#
import LiquidValve as lv
class PressureScreen:
    """Pressure Screen Model"""
    #
    # Class Constants
    #
    CVREJ, CVACPT, CVDIL = 3, 5, 6
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    
    #
    # Class Varaiables
    #
    # Inlet, Dilution Pump, Rejects Tank Head, Accepts Tank Head, PSIA
    P = np.zeros(4)
    # Flows Inlet, Dilution, Rejects, Accepts flow in GPM
    F = np.zeros(4)
    # Consisency Inlet, Dilution, Rejects and Accepts in %
    C = np.array([6.,0.,12.,2.])
    # Thevenin Circuit Equivalant CV [GPM/sqrt(PSI)], Press [PSIA] 
    Th = np.array([100.,10.])
    # Delta Press accross Screen [PSI]
    DP = 0.0

    #
    # class methods
    #
    def __init__(self, MaxCVs, P, ConsFeed):
        """Initialize Screen Model"""

        # N/A, Pipe, Reject Internal, Reject Valve, Screen, Accepts Valve, Dil Valve
        self._MaxCV = MaxCVs
        # Inlet, Dilution Pump, Rejects Tank Head, Accepts Tank Head, PSIA
        self.P = P
        self.C[self.INLET] = ConsFeed
        self._DilVLV = lv.LiquidValve(MaxCVs[self.CVDIL],curve = 'Linear',\
                                      rangeability = 10000.)
        self._AcceptsVLV = lv.LiquidValve(MaxCVs[self.CVACPT],curve = 'Linear',\
                                          rangeability = 10000.)
        self._RejectsVLV = lv.LiquidValve(MaxCVs[self.CVREJ],curve = 'Linear',\
                                          rangeability = 10000.)
        self._CV = MaxCVs
        self.C[self.INLET] = ConsFeed # Feed Consistency %


    #
    # Build CV array from Valve Positions
    #
    def UpdateCV(self, dVLV, rVLV, aVLV):
        """Build new CV array from Valve Positions in %"""

        self._CV[self.CVDIL] = self._DilVLV.GetCV(dVLV)
        self._CV[self.CVACPT] = self._AcceptsVLV.GetCV(aVLV)
        self._CV[self.CVREJ] = self._RejectsVLV.GetCV(rVLV)
        return self._CV
    
    #
    # Get Inlet, Dilution, Accepts and Rejects flows
    # given valve percentages open and Pressures
    #    
    def UpdateFlows(self, dVLV, rVLV, aVLV, P,\
                 F0 = np.array([2000.,2000.,2000.])):
        """Valve Pos in % and P = Inlet, Dilution Pump,
        Rejects Tank Head, Accepts Tank Head, PSIA"""

        CV = self.UpdateCV(dVLV, rVLV, aVLV)
        Gth, Pth = self.UpdateThevenin(dVLV, rVLV, aVLV,P)
        Fin = Gth*np.sqrt(abs(P[self.INLET]-Pth))
        self.F = solveScreen(CV,P,Fin,F0)
        return self.F
    #
    # Calculate Rejects and Accepts Pressures given inlet..
    # ..dilution pump pressures and flows
    #
    def UpdatePress(self, dVLV, rVLV, aVLV, P, F):
        """Flows 0-2 array, Valve Pos in % and 
           Inlet, Dilution Pump Pressures in PSIA"""
        CV = self.UpdateCV(dVLV, rVLV, aVLV)
        for i in range(np.shape(CV)[0]):
            CV[i] = max(1.,CV[i])
        self.P = self.P #solvePress(F,CV,P)
        self.DP = F[self.ACCPT]**2/CV[4]**2
        return self.P       

    #
    # Solve Consistencies
    #
    def UpdateCons(self, F, ConsFeed, RTF = 2.0,):
        """Solve for Accepts and Rejects Consistency %
        given inlet consistency, reject thickening factor and 
        flows = Feed, Dilution, Accepts, Rejects"""

        self.C[self.INLET] = ConsFeed
        self.C[self.ACCPT], self.C[self.REJ] = solveCons(F,ConsFeed,RTF)
        return self.C
    
    #
    # Solve for Thevenin CV and Pressure
    #
    def UpdateThevenin(self, dVLV, rVLV, aVLV, P,\
                    Fth = np.array([2000.,2000.])):
        """Valve Pos in % and P = Inlet, Dilution Pump,
        Rejects Tank Head, Accepts Tank Head, PSIA"""

        CV = self.UpdateCV(dVLV, rVLV, aVLV)
        self.Th = solveThevenin(CV,P,Fth)
        return self.Th     
    
    #
    # Iterate the screen one solution
    #
    def Iterate(self,Pin, Pdil, dVLV, rVLV, aVLV, ConsF, RTF):
        """ Iterate the screen given Pin = P Inlet, Pdil = P Dilution
            Dilution, Rejects and Accepts Valve position in %
            Feed Consistency in % and Rejects Thickening Factor unitless """
            
        # Build initial guess
        F0 = np.array([self.F[self.ACCPT],self.F[self.INLET] -\
                       self.F[self.ACCPT],self.F[self.DIL]])
        self.P[self.INLET],self.P[self.DIL] = Pin, Pdil
        F = self.UpdateFlows(dVLV,rVLV,aVLV,self.P,F0)
        C = self.UpdateCons(F,ConsF,RTF)
        P = self.P# self.UpdatePress(dVLV,rVLV,aVLV,self.P,F)
        return F,P,self.DP,C
#
# if run as a script create test object
#
import matplotlib.pyplot as plt
import time as tm
if __name__ == "__main__":
    # Constants
    INLET, DIL, REJ, ACCPT = 0, 1, 2, 3
    
    # Inlet, Dilution Pump, Rejects Tank Head, Accepts Tank Head, PSIA
    P = np.array([61., 70., 20., 18.])
    # N/A, Pipe, Reject Internal, Reject Valve, Screen, Accepts Valve, Dil Valve
    MaxCV = np.array([0., 1000., 200., 600., 1200., 500., 40.])
    # Initial flow guesses
    F0 = np.array([2000.,2000.,2000.])
    # Feed Consistency
    ConsF = 6. #%
    # Rejects Thickening Factor
    RTF = 1.9
    
    
    Screen1B = PressureScreen(MaxCV,P,ConsF)
    
    # null, dilution, rej and accpt valve positions
    VLV = np.array([0.,100.,100.,100.])

    iterations = 10
    Fs = np.zeros((iterations,4))
    Ps = np.zeros_like(Fs)
    Cs = np.zeros_like(Fs)
    DPs = np.zeros(iterations)
    
    for i in range(iterations-1):
        t1 = tm.clock()
        Fs[i+1], Ps[i+1], DPs[i+1], Cs[i+1] = Screen1B.Iterate(P[INLET],\
         P[DIL],VLV[DIL], VLV[REJ],VLV[ACCPT],ConsF,RTF)
        t2 = tm.clock()
        dt = t2 - t1
        print 'Iteration Time ' + str(dt)
    
    t = range(iterations)
    plt.figure()
    plt.subplot(911)
    plt.plot(t,Ps[:,INLET],'r-',label='Inlet Press')
    plt.legend(loc='best')
    plt.subplot(912)
    plt.plot(t,Ps[:,DIL],'b-',label='Dil Press')
    plt.legend(loc='best')
    plt.subplot(913)
    plt.plot(t,Ps[:,REJ],'g-',label='Rej Press')
    plt.legend(loc='best')
    plt.subplot(914)
    plt.plot(t,Ps[:,ACCPT],'k-',label='Accept Press')
    plt.legend(loc='best')
    plt.subplot(915)
    plt.plot(t,DPs,'m-',label='Screen DP')
    plt.legend(loc='best')
    #======================================
    plt.subplot(916)
    plt.plot(t,Fs[:,INLET],'r.',label='Inlet Flow')
    plt.legend(loc='best')
    plt.subplot(917)
    plt.plot(t,Fs[:,DIL],'b.',label='Dil Flow')
    plt.legend(loc='best')
    plt.subplot(918)
    plt.plot(t,Fs[:,REJ],'g.',label='Rej FLow')
    plt.legend(loc='best')
    plt.subplot(919)
    plt.plot(t,Fs[:,ACCPT],'k.',label='Accept Flow')
    plt.legend(loc='best')
    plt.show()
        
        
        


    
#    Flows, P, DP, Cons = Screen1B.UpdateScreen(P[INLET],P[DIL],VLV[DIL], VLV[REJ],VLV[ACCPT],ConsF,RTF)
    
#    print 'Inlet Flow      = % 8.2f' % (Flows[INLET]), "GPM"
#    print 'Dilution Flow   = % 8.2f' % (Flows[DIL]), "GPM"
#    print 'Accepts Flow    = % 8.2f' % (Flows[ACCPT]), "GPM"
#    print 'Rejects Flow    = % 8.2f' % (Flows[REJ]), "GPM"
#    print 'Screen DP       = % 8.2f' % (DP), "PSI"
#    print " "
#    print 'Feed Cons       = % 8.2f' % (Cons[INLET]), "%"
#    print 'Accepts Cons    = % 8.2f' % (Cons[ACCPT]), "%"
#    print 'Rejects Cons    = % 8.2f' % (Cons[REJ]), "%"
#    print ' ' 
#    print 'Thevenin CV     = % 8.2f' % (Screen1B.Th[0]), "GPM/PSI^0.5"
#    print 'Thevenin Press  = % 8.2f' % (Screen1B.Th[1]), "PSIA"
#    print ' '
#    print 'Inlet Press     = % 8.2f' % (P[INLET]), 'PSIA'
#    print 'Dilution Press  = % 8.2f' % (P[DIL]), 'PSIA'
#    print 'Rejects Press   = % 8.2f' % (P[REJ]), 'PSIA'
#    print 'Accepts Press   = % 8.2f' % (P[ACCPT]), 'PSIA'
