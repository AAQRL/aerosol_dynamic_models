# Written by Girish Sharma under the instruction of Dr. Pratim Biswas in the Washington University in St. Louis 
# in March, 2016. 
from tkinter import *
import math as ma # maths operations
import matplotlib.pyplot as plt # plotting
import numpy as np # for arrays
from scipy.integrate import ode # for integration

COAG = 0
REAC = 0
COND = 0
NUCL = 0

############################################################ GUI SECTION ####################################################
class Application(Frame):
    
    def __init__(self,master):
        """Initialise the Frame"""
        Frame.__init__(self,master)
        self.grid()
        self.create_widgets()
    
    def create_widgets(self):
        """Create button, text and entry widgets"""
        # all the labels for the GUI are written here
        Label(self, text = "SELECT THE SYSTEM PROPERTIES").grid(row = 0)
        Label(self, text = "Temperature ").grid(row = 1, column =0, columnspan = 2, sticky = W)
        Label(self, text = " Kelvin ").grid(row = 1, column =3, sticky = W)
        Label(self, text = "Pressure ").grid(row = 2, column =0, columnspan = 2, sticky = W)
        Label(self, text = " atm ").grid(row = 2, column =3, sticky = W)
        
        Label(self, text = "SELECT THE SPECIES PROPERTIES").grid(row =3)
        Label(self, text = "Saturation Pressure ").grid(row = 4, column =0, columnspan = 2, sticky = W)
        Label(self, text = " Pa ").grid(row = 4, column =3, sticky = W)
        Label(self, text = "Molecular Weight ").grid(row = 5, column =0, columnspan = 2, sticky = W)
        Label(self, text = " g/mol ").grid(row = 5, column =3, sticky = W)
        Label(self, text = "Density ").grid(row = 6, column =0, columnspan = 2, sticky = W)
        Label(self, text = " Kg/m3 ").grid(row = 6, column =3, sticky = W)
        Label(self, text = "Surface Tension ").grid(row = 7, column =0, columnspan = 2, sticky = W)
        Label(self, text = " dyne/cm ").grid(row = 7, column =3, sticky = W)
        
        Label(self, text = "TICK FOR THE PHENOMENON (SELECT ALL THAT APPLY)").grid(row = 8)
        Label(self, text = "Dimensionless Reaction Rate ").grid(row = 12, column =1)
        Label(self, text = " 0.1(low) - 1(high) ").grid(row = 12, column =3)
        
        Label(self, text = "INITIAL CONDITIONS").grid(row = 14)
        Label(self, text = "Particle Number Concentration").grid(row = 15, column =0, columnspan = 2, sticky = W)
        Label(self, text = " #/cm3 ").grid(row = 15, column =3, sticky = W)
        Label(self, text = "Particle Geometric Average Diameter").grid(row = 16, column =0, columnspan = 2, sticky = W)
        Label(self, text = " nm ").grid(row = 16, column =3, sticky = W)
        Label(self, text = "Geometric Standard Deviation (>=1)").grid(row = 17, column =0, columnspan = 2, sticky = W)
        Label(self, text = " dimensionless ").grid(row = 17, column =3, sticky = W)
        Label(self, text = "Saturation Ratio ").grid(row = 18, column =0, columnspan = 2, sticky = W)
        Label(self, text = " dimensionless ").grid(row = 18, column =3, sticky = W)
        Label(self, text = "Simulation time ").grid(row = 19, column =0, columnspan = 2, sticky = W)
        Label(self, text = " sec ").grid(row = 19, column =3, sticky = W)
        
        def_T = StringVar(self,value = '298')
        self.T = Entry(self, textvariable = def_T)
        self.T.grid(row=1, column =2, sticky =W)
        def_P = StringVar(self,value = '1')
        self.P = Entry(self, textvariable = def_P)
        self.P.grid(row=2, column =2, sticky =W)
        def_Ps = StringVar(self,value = '0.000001')
        self.Ps = Entry(self, textvariable = def_Ps)
        self.Ps.grid(row=4, column =2, sticky =W)
        def_M = StringVar(self,value = '100')
        self.M = Entry(self, textvariable = def_M)
        self.M.grid(row=5, column =2, sticky =W)
        def_rho = StringVar(self,value = '1000')
        self.rho = Entry(self, textvariable = def_rho)
        self.rho.grid(row=6, column =2, sticky =W)
        def_gam = StringVar(self,value = '25')
        self.gam = Entry(self, textvariable = def_gam)
        self.gam.grid(row=7, column =2, sticky =W)
        
        self.COAG = BooleanVar()
        Checkbutton(self, text="Coagulation", variable = self.COAG, command = self.coag).grid(row = 9, column = 0, sticky =W)
        self.COND = BooleanVar()
        Checkbutton(self, text="Condensation", variable = self.COND, command = self.cond).grid(row = 10, column = 0, sticky =W)
        self.NUCL = BooleanVar()
        Checkbutton(self, text="Nucleation", variable = self.NUCL, command = self.nucl).grid(row = 11, column = 0, sticky =W)
        self.REAC = BooleanVar()
        Checkbutton(self, text="Reaction", variable = self.REAC, command = self.reac).grid(row = 12, column = 0, sticky =W)
        def_r = StringVar(self,value = '0.1')
        self.r = Entry(self, textvariable = def_r)
        self.r.grid(row=12, column =2, sticky =W)
        
        def_num_zero = StringVar(self,value = '0e0')
        self.num_zero = Entry(self, textvariable = def_num_zero)
        self.num_zero.grid(row=15, column =2, sticky =W)
        def_dpg_zero = StringVar(self,value = '0e0')
        self.dpg_zero = Entry(self, textvariable = def_dpg_zero)
        self.dpg_zero.grid(row=16, column =2, sticky =W)
        def_sigma_g_zero = StringVar(self,value ='1')
        self.sigma_g_zero = Entry(self, textvariable = def_sigma_g_zero)
        self.sigma_g_zero.grid(row=17, column =2, sticky =W)
        def_S_zero = StringVar(self,value = '0e0')
        self.S_zero = Entry(self, textvariable = def_S_zero)
        self.S_zero.grid(row=18, column =2, sticky =W)
        def_t_final_input = StringVar(self,value = '9000')
        self.t_final_input = Entry(self, textvariable = def_t_final_input)
        self.t_final_input.grid(row=19, column =2, sticky =W)        
        
        self.submit_button = Button(self, text ="Submit", command = self.input_param)
        self.submit_button.grid(row = 20,column = 1)
        # this button will call input_param function and this is how we more to main program 
        
    def coag(self):
        global COAG
        if self.COAG.get():
            COAG = 1
        else:
            COAG = 0
        
    def cond(self):
        global COND
        if self.COND.get():
            COND = 1
        else:
            COND = 0
            
    def nucl(self):
        global NUCL
        if self.NUCL.get():
            NUCL = 1
        else:
            NUCL = 0
                
    def reac(self):
        global REAC
        if self.REAC.get():
            REAC = 1
        else:
            REAC = 0
    
    def input_param(self):
        """Display messages based on the input passed in"""
        # brings all the input values together and are made global after converting string to float
        global T,P,Ps,M,rho,gam,r_rate # all SI except r_rate which is dimensionless 
        T = float(self.T.get())
        P = 101325.0*float(self.P.get()) # atm to Pa
        Ps= float(self.Ps.get())       
        M = 1.0e-3*float(self.M.get()) # g/mol to Kg/mol
        rho = float(self.rho.get())   
        gam = 1.0e-3*float(self.gam.get()) # dyne/cm to N/m
        r_rate = float(self.r.get())
        
        # initial conditions (convert all to SI)
        num_zero = 1.0e6*float(self.num_zero.get()) # #/cm3 to #/m3
        dpg_zero = 1.0e-9*float(self.dpg_zero.get()) # nm to m
        sigma_g_zero = float(self.sigma_g_zero.get()) 
        S_zero = float(self.S_zero.get())   
        vg_zero = (22.0/7)/6*(dpg_zero)**(3.0)  
        t_final_input = float(self.t_final_input.get())
        initial = [num_zero,vg_zero,sigma_g_zero,S_zero, t_final_input]

        main(initial)
        





############################################ FUNCTIONS FOR THE MAIN FUNCTION ###########################################################

# MOMENT FUNCTION
# # Inputs are :
#             Nv  : Number conentration (#/m3)
#            vgv  : geometric mean volume (m3)
#            sdv  : standard deviation (dimensionless)
#            kv   : kth moment
def moment(Nv,vgv,sdv,kv):
    return Nv*(vgv)**(kv)*ma.exp(9.0/2*(kv*kv)*ma.log(sdv)*ma.log(sdv))

# COAGULATION FUNCTION
# Inputs are :
#            sd  : standard deviation (dimensionless)
#            rg  : radius corresponding to vg (m)
#            flag: takes 0 or 1 value (if 0: no coagulation, if 1: coagulation takes place)
# Outputs is :
#            com : com has com0 (for the coagulation coefficient of zeroth moment m3/s) and 
#                          com2 (for the coagulation coefficient of second moment m3/s)

def coag(sd,rg,flag):
    if flag == 0 or sd <= 1e-20 or rg <= 1e-40:
        com = [0.0 , 0.0]
        return com
    
    else:
        b0 = 0.633 + 0.092*sd*sd - 0.022*sd*sd*sd
        b2 = 0.39 + 0.5*sd - 0.214*sd*sd + 0.029*sd*sd*sd
        c1 = (6.0*kb*T*rg/rho)**(1.0/2)
        c2 = 2.0*kb*T/(3.0*mu)
        
        lnsig = ma.log(sd)*ma.log(sd)
        b5 = 1.257
        
        com0fm = c1 * b0 * ( ma.exp(25.0/8* lnsig) + 2.0*ma.exp(5.0/8*lnsig) + ma.exp(1.0/8*lnsig) ) 
        com0cn = c2 * ( 1.0 + ma.exp(lnsig) + b5 * (lam/rg) * ma.exp(1.0/2 * lnsig) * (1.0 + ma.exp(2.0*lnsig)) )
        
        com2fm = 2.0 * c1 * b2 * ma.exp(3.0/2*lnsig) * ( ma.exp(25.0/8* lnsig) + 2.0*ma.exp(5.0/8*lnsig) + ma.exp(1.0/8*lnsig) )
        com2cn = 2.0 * c2 * ( 1.0 + ma.exp(lnsig) + b5 * (lam/rg) * ma.exp(-1.0/2 * lnsig) * (1.0 + ma.exp(-2.0*lnsig)))
        
        com0 = (com0fm * com0cn) / (com0fm + com0cn) # com0fm #
        com2 = (com2fm * com2cn) / (com2fm + com2cn) # com2fm #
        com  = [com0 , com2]
        return com


# CONDENSATION FUNCTION
# Inputs are :
#            sd  : standard deviation (dimensionless)
#            vg  : geometric mean volume (m^3)
#            flag: takes 0 or 1 value (if 0: no condensation, if 1: condensation takes place)
# Outputs is :
#            con: con has coefm1  (for the condensation coefficient of first moment m3/s) and 
#                         coefsx  (for the condensation coefficient of monomer balance m3/s)
#                         coefm2  (for the condensation coefficient of second moment m3/s)


def cond(sd,vg,flag):
    lnsig = ma.log(sd)*ma.log(sd)
    A  = 1.0/3
    F1 = (36.0*PI)**A * v1 * ns * ma.sqrt(kb*T/2/PI/m1)
    F2 = ma.sqrt(8.0*kb*T/(PI*m1))
    C1 = (48.0*PI*PI)**A*ns*v1*lam/3
    if flag == 0 or vg <= 0:
        con = [ 0.0 , 0.0 , 0.0 ]
        return con
    else:
        etac1 = v1**A /tau * vg**(2.0*A) * ma.exp(2.0*lnsig)
        delc1 = etac1
        psic1 = 2.0 * v1**A /tau *vg**(2.0*A) * ma.exp(8.0 * lnsig)
        etac2 = 8.0/3/tau/(2.0*r1) *lam * v1**(2.0*A) * vg**A * ma.exp(lnsig/2)
        delc2 = etac2
        psic2 = 16.0/3/tau/(2.0*r1) * lam * v1 **(2.0*A) * vg**(A) * ma.exp(3.5 * lnsig)
        coefm1 = etac1 * etac2 / (etac1 + etac2) # etac1 #
        coefsx = delc1 * delc2 / (delc1 + delc2) #  delc1 #
        coefm2 = psic1 * psic2 / (psic1 + psic2) # psic1 #
        con = [ coefm1 , coefsx , coefm2 ]
        return con

    
    
# NUCLEATION FUNCTION
# Inputs are :
#            S   : saturation ratio (dimensionless)
#            flag: takes 0 or 1 value (if 0: no nucleation, if 1: nucleation takes place
# Outputs are:
#            nuc : nuc has kstar (dimensionless critical particle size) and 
#                          x11    nucleation rate ( #/m3/s)

def nucl(S,flag):
    if S < 1.001 or flag == 0:
        nuc = [ 0.0 , 0.0 ]
        return nuc
    else:
        sigmad = gam * v1 ** (2.0/3) / (kb * T)
        xks = PI/6 * (4.0 * sigmad)**(3.0)
        kstar = xks / (ma.log(S))**(3.0)
        coef1 = (2.0 / 9.0 / PI) ** (1.0/3)
        pep = kstar * (ma.log (S)/2)
        if pep < 60000.0:
            coef2 = ma.exp(-pep)
            x11 = ns / tau * S**(2.0) * coef1 * sigmad **(0.5) * coef2
            nuc = [ kstar , x11 ]
            return nuc
        else:
            coef2 = 0.0
            x11 = 0.0
            nuc = [ kstar , x11 ]
            return nuc
        

# REACTION FUNCTION
# Inputs are : 
#             (currently no inputs)
# Outputs are:
#           reac : rate of the reaction (#/m3/s) 

def reac():
    if REAC == 0:
        R = 0.0
    else:
        R  = r_rate*(ns/tau)                                                  
    return R


def vectorfield(t , y):
    """
    Defines the differential equations for the 4 parameters to fins the lognormal size distribution and vapor phase conc.

    Arguments:
        w :  vector of the state variables:
                  w = [N, V, V2, S]
        t :  time (dimensionless)
        p :  parameter
                  p = K
    """
    M0 = y[0]*n_d
    M1 = y[1]*n_d*v1
    M2 = y[2]*n_d*v1*v1
    S  = y[3]
    
    #pdb.set_trace()
    
    if y[0] > 1e-99 and y[1] > 1e-99 and y[2] > 1e-99:
        vg = M1*M1 / ( M0 ** (1.5) * M2 ** (0.5))
        if vg > 0:
            rg = (3.0*vg/(4.0*PI))**(1.0/3)
            linsig = 1.0/9 * ma.log(M0*M2/(M1*M1))
            if linsig > 0:
                sd = ma.exp(linsig**0.5)
            else:
                sd = 1.0
        else:
            rg = 0.0
            linsig = 0.0
            sd = 1.0
    else:
        vg = 0.0
        rg = 0.0
        linsig = 0.0
        sd = 1.0
    
    
    # Call all the functions that are required here and store the values in appropriate variables
    # unpack
    com0 , com2                   =   coag(sd,rg,COAG)
    coefm1 , coefsx , coefm2      =   cond(sd,vg,COND)
    kstar , x11                   =   nucl(S,NUCL)
    R                             =   reac()

    
    ######################################## NUCLEATION ############################################ 
    dxNUC_S =   x11 / ns * kstar * tau
    dxNUC_N =   x11 / n_d * tau
    dxNUC_V =   x11 / n_d * kstar * tau
    dxNUC_W =   x11 / n_d * kstar * kstar * tau
    
    
    ######################################## REACTION ############################################ 
    dxREA_S =   R / (ns / tau)
    dxREA_N =   0.0
    dxREA_V =   0.0
    dxREA_W =   0.0
    
    
    ######################################## CONDENSATION ############################################ 
    if y[0] > 1e-50 and y[1] > 1e-50:
        dxCON_S =    coefsx * y[0] * ( S - 1.0 ) * tau /v1       
        dxCON_V =    coefm1 * y[0] * ( S - 1.0 ) * tau /v1 * ns/n_d
        dxCON_W =    coefm2 * y[1] * ( S - 1.0 ) * tau /v1 * ns/n_d
    else:
        dxCON_S =    0.0       
        dxCON_V =    0.0
        dxCON_W =    0.0
    
    
    ######################################## COAGULATION ################################################## 
    if y[0] > 1e-50 and y[1] > 1e-50 and y[2] > 1e-50:
        dxCOA_N =    com0 * y[0] ** 2.0 * n_d * tau
        dxCOA_W =    com2 * y[1] ** 2.0 * n_d * tau
    else:
        dxCOA_N =    0.0
        dxCOA_W =    0.0
     
    
    

    ############ Gathering all the phenomenon together into f(t,y) #########################################
    f = [dxNUC_N - dxCOA_N,
         dxNUC_V + dxCON_V,
         dxNUC_W + dxCON_W + dxCOA_W,
         dxREA_S - dxNUC_S - dxCON_S]
    return f


# Use ODEINT to solve the differential equations defined by the vector field


# The ``driver`` that will integrate the ODE(s):
def main(initial):
    global kb,Na,PI
    global ns,lam,mu
    global m1,v1,r1,s1,sig,tau,K,Kn1
    global n_d
    
    kb  = 1.38064852e-23                                                     # Boltzmann constant (J/K)
    Na  = 6.022140857e23                                                     # Avagadro's Number (mol-1)
    PI  = 3.14159                                                            # pi 
                                                       
    
    # DERIVED CONSTANTS
    Mg = 28.97e-3                                                            # Molecular Weight of air (Kg/m3)
    ns = Ps / (kb*T)        
    mu=1.716e-5 * (T/273)**(2.0/3)                                           # viscosity of the medium
    lam = mu/P * ma.sqrt(PI*kb*T/(2.0*Mg/Na))                                # mean free path of air (m)
    
    m1 = M/Na                                                                # mass of a monomer (Kg)
    v1 = m1/rho                                                              # volume of a monomer (m3)
    r1 = (1.0/2) * (6.0*v1/PI)**(1.0/3)                                      # radius of a monomer (m)
    s1 = 4.0 * PI * r1 * r1                                                  # area of a monomer (m2)
    sig= gam * v1**(2.0/3) / (kb*T)                                          # surface tension group (dimensionless)
    tau= (ns*s1*(kb*T/(2.0*PI*m1))**(1.0/2))**(-1.0)                         # Characterstic time for particle growth (s)
    K  = (2.0*kb*T/(3.0*mu))*ns*tau                                          # Coagulation coefficient (dimensionless)



    ###################################################### PLACE WHERE I NEED TO MAKE CHANGES ################################

 
    num_zero = initial[0]                                                     # Particle Number Concentration (#/m3)  at t=0
    v_g_zero = initial[1]                                                     # Geometric mean volume ( in m3) at t=0
    sigma_g_zero = initial[2]                                                 # Geometric standard deviation at t=0
    n1_zero = initial[3] * ns                                                 # Vapor Concentration (#/m3)  at t=0 
    t_final_input = initial[4]
    
    if num_zero > 1e9:
        n_d = num_zero
    else:
        n_d = ns

    # Set the time range
    t_start = 0.0
    t_final = t_final_input/tau
    delta_t = 0.1

    # Start by specifying the integrator:
    # use ``vode`` with "backward differentiation formula"
    r = ode(vectorfield).set_integrator('vode', method='bdf', nsteps = '100000') # atol= '1e-9', rtol= '1e-9' 
 

    # Number of time steps: 1 extra for initial condition
    num_steps = np.floor((t_final - t_start)/delta_t) + 1
    

    # Set initial condition(s): for integrating variable and time!
    N_t_zero = num_zero/n_d
    V_t_zero = moment(num_zero,v_g_zero,sigma_g_zero,1.0)/(v1*n_d)
    V2_t_zero  = moment(num_zero,v_g_zero,sigma_g_zero,2.0)/(v1*v1*n_d)
    S_t_zero = n1_zero/ns
    r.set_initial_value([N_t_zero, V_t_zero, V2_t_zero, S_t_zero], t_start)
    

    
 
    # create vectors to store variables
    t = np.zeros((num_steps, 1))    # time dimensionless   
    N_t = np.zeros((num_steps, 1))  # Number Concentration Dimensionless
    V_t = np.zeros((num_steps, 1))  # Moment 1 Dimensionless
    V2_t = np.zeros((num_steps, 1)) # Moment 2 Dimensionless
    S_t = np.zeros((num_steps, 1))  # Saturation ratio
    W_t = np.zeros((num_steps, 1))  # polydispersity index
    dp_t = np.zeros((num_steps, 1))
    
    # initialise the vectors for integration
    t[0] = t_start
    N_t[0] = N_t_zero
    V_t[0] = V_t_zero
    V2_t[0] = V2_t_zero
    S_t[0] = S_t_zero
    W_t[0] = 0.0
    
 
    # Integrate the ODE(s) across each delta_t timestep
    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + delta_t)
 
        # Store the results to plot later
        t[k] = r.t
        N_t[k] = r.y[0]
        V_t[k] = r.y[1]
        V2_t[k] = r.y[2]
        S_t[k] = r.y[3]
        k += 1
    
    

    

    # print the solution in dimensional form 
    with open('Moment_ode', 'w') as f:
        print('[t(in s)]', '[M0(#/m3)]', '[M1(m3/m3)]', '[M2(m6/m3)]', '[S(unitless)]', file=f)
        for i in range(0,k):
            M0 = N_t[i]*n_d
            M1 = V_t[i]*n_d*v1
            M2 = V2_t[i]*n_d*v1*v1
            S  = S_t[i]
            print(tau*t[i], M0, M1, M2, S, file=f)
            
    # create vectors to save the size distribution data
    Nc = np.zeros((num_steps, 1))
    dpmg = np.zeros((num_steps, 1))
    sigmag = np.zeros((num_steps,1))
    
    # print the size distribution
    with open('Size_distribution', 'w') as f2:
        # Print & save the solution.
        print('[time(min)]',' [Number Conc ( /cm3)]', '[Geomteric mean diameter(nm)]', '[standard deviation]', file=f2)
        for i in range(0,k):
            Nc[i] = n_d * N_t[i]* 1e-6
            # if else to avoid division by zero
            if N_t[i] < 1e-80 or V2_t[i] < 1e-80 or V_t[i] < 1e-80:
                vg = 0.0
                temp = 0.0
                dp_t[i] = 0.0
            else:
                vg = v1 * V_t[i]*V_t[i]/(N_t[i]**(3.0/2)*V2_t[i]**(1.0/2))
                pop = v1 * V_t[0]*V_t[0]/(N_t[0]**(3.0/2)*V2_t[0]**(1.0/2))
                temp = N_t[i] * V2_t[i] / ( V_t[i] * V_t[i] )
                dp_t[i] = (vg/v1)**(1.0/3)
            
            # check if dp_t is less than the diameter for 1 molecule
            if dp_t[i] < 1:
                dp_t[i] = 1.0
            # check if the dp_t value is nan
            elif ma.isnan(dp_t[i]):
                dp_t[i] = 0.0
            else:
                pass
            
            # convert dp_t to dimensional form
            dpmg[i] = 2.0*r1*dp_t[i]*1e9
            
            # check the value of temp
            if temp >  1.0:
                sigmag[i] = ma.exp(ma.sqrt(((1.0/9)*ma.log(temp))))
                W_t[i] = ma.sqrt(ma.exp(1.0/9 * ma.log(temp))-1)
            else:
                sigmag[i] = 1.0
                W_t[i] = 0.0
                
            print(t[i]*tau/60, Nc[i], dpmg[i], sigmag[i], file=f2)
            # if number concentration is less than 1e2 #/m3, put it equal to zero in plots
            if N_t[i]*n_d <= 1e2: 
                N_t[i] = 1e2/n_d
                
    with open('check', 'w') as f:
        print(sigma_g_zero, v_g_zero, num_zero, N_t_zero, n_d, pop, file=f)
        print(moment(num_zero,v_g_zero,sigma_g_zero,1) , moment(num_zero,v_g_zero,sigma_g_zero,2), file=f)
        print(V_t[0], V_t[0] * (v1*n_d), (V_t[0]/N_t[0])**(1.0/3) * 2.0 * r1, file=f)
    
    # import data from Prastinis paper for comparison
#     N = 393
#     file = open('dp_C_R1.txt', 'r')
#     i = 0
#     t_paper = range(N)
#     dp_t_paper = range(N)
#     for line in file:
#         line = line.split(", ")
#         t_paper[i] = float(line[0])
#         dp_t_paper[i] = float(line[1])
#         i = i+1
#     N = 395
#     file = open('dp_nC_R1.txt', 'r')
#     i = 0
#     t_paper = range(N)
#     dp_t_paper = range(N)
#     for line in file:
#         line = line.split(", ")
#         t_paper[i] = float(line[0])
#         dp_t_paper[i] = float(line[1])
#         i = i+1
#     N = 379
#     file = open('dp_nC_R01.txt', 'r')
#     i = 0
#     t_paper = range(N)
#     dp_t_paper = range(N)
#     for line in file:
#         line = line.split(", ")
#         t_paper[i] = float(line[0])
#         dp_t_paper[i] = float(line[1])
#         i = i+1
#     N = 360
#     file = open('dp_nC_R01_FM_paper.txt', 'r')
#     i = 0
#     t_paper = range(N)
#     dp_t_paper = range(N)
#     for line in file:
#         line = line.split(", ")
#         t_paper[i] = float(line[0])
#         dp_t_paper[i] = float(line[1])
#         i = i+1    

    print(ns*tau*r_rate)



#     plot the results
    plt.figure(figsize=(15,15))   
    plt.subplot(311)
    plt.grid(True)
    plt.ylabel('Number Concentration N (#/cm3)')
    plt.semilogy(t*tau/60, N_t*n_d*1e-6,color="blue", linewidth = 1.0, linestyle ="-")
    plt.title('Lognormal Size Distribution Parameters')
    plt.subplot(312)
    plt.grid(True)
    plt.ylabel('Average diameter dp (nm)')
    plt.plot(t*tau/60, dp_t*2.0*r1*1e9,color="blue", linewidth = 1.0, linestyle ="-")
    plt.subplot(313)   
    plt.grid(True)
    plt.xlabel('time (in minutes)')
    plt.ylabel('sigma_g')
    plt.plot(t*tau/60, sigmag,color="blue", linewidth = 1.0, linestyle ="-")
    plt.show()
#     plt.figure(2)
#     plt.grid(True)
#     plt.xlabel('time (in minutes)')
#     plt.ylabel('Number Concentration N (#/cm3)')
#     plt.semilogy(t*tau/60, N_t*n_d*1e-6,color="blue", linewidth = 1.0, linestyle ="-")
#     plt.show()
        
if __name__ == "__main__":
    root = Tk()
    root.title("Simulator for Moment Method")
    root.geometry("900x500")

    app= Application(root)
    root.mainloop()
