import numpy as np
from polychrom.starting_conformations import grow_cubic

class cell:
    """
    Container class that stores information about a bacterial system, including the confinement (nucleoid) size and loop-extrusion parameters.


    Length units are in simulation monomer size,
    (approximate) time units are in minutes
    """

    def __init__(self, name, N, monomerSize, radius, L_0, growthRate, LE_lifetime, LE_stepRate, LE_backstepRate, parS_degrees, parS_strength, terRegion, terStrength, tReplication, timeTraverse=7/60):
        
        assert len(parS_degrees)==len(parS_strength)
        self.name=name
        self.N=N
        self.terStrength=terStrength
        self.terLength=len(terRegion)
        #lengths in lattice units
        self.monomer_size=monomerSize
        self.radius=radius/monomerSize
        self.L_0=L_0/monomerSize

        #rates and times using minutes
        self.rateGrowth=growthRate
        self.rateReplication=N/tReplication/2 #step rate for forks, monomers per minute

        #Calculate arrays for LE simulations
        self.loopSize=LE_lifetime*2*LE_stepRate
        self.offloadingRates=np.ones(N)*1/LE_lifetime
        self.offloadingRates[terRegion]*=terStrength
    
        parSsites=[int((parS%360)*N/360) for parS in parS_degrees]
        self.loadingRates=np.ones(N)
        self.loadingRates[parSsites]=parS_strength
        if terStrength>1:
            #if we have enhanced off-loading at ter, also don't load SMCs there
            self.loadingRates[terRegion]=0

        self.stepRate=LE_stepRate
        self.backstepRate=LE_backstepRate
        self.timeTraverse=timeTraverse #seconds

    def multiply_length(self, length_factor):
        self.L_0*=length_factor

    def multiply_radius(self, radius_factor):
        self.radius*=radius_factor

    #for steady state, only circular chromosome and then linear replicated segment
    def start_point_unreplicating(self,topMonomer):
        height=self.R_to_height(0)
        start_data = grow_cubic(self.N, (int(height)-1), method="extended")  # creates a compact conformation that fills the height
    
        for i in range(self.N):
            start_data[i][:]=(start_data[i][:]-int(height-1)/2) #center it

        #roll the array so that top monomer gets the position of monomer 1
        if topMonomer>0:
            start_data=np.roll(start_data, self.N-topMonomer)
        
        return start_data

    def inferred_time(self, R):
        inferred_t=R/self.rateReplication #in minutes
        return inferred_t

    def t_to_height(self, t):
        height=self.L_0*np.exp(self.rateGrowth*t) #in lattice units
        return height

    def R_to_height(self, R):
        t=self.inferred_time(R)
        return self.t_to_height(t)

#E coli
def e_coli(N, monomer_size, loop_size, ter_size_kb=800, terStrength=100, no_bypass=False):
    t_replication=50 #minutes
    t_doubling=60
    growth_rate=np.log(2)/t_doubling #exponential growth rate for cells, per minute
    loop_extruder_speed = 46 #kb/min
    radius=1000/2 #nm
    
    bp_per_monomer=N//4600
    cut_from_ter=800-ter_size_kb
    total_kb=4600-cut_from_ter
    N_cut=N-cut_from_ter//bp_per_monomer

    half_ter_region=ter_size_kb//2//bp_per_monomer
    mid_ter_region=N_cut//2
    
    ter_sites = np.arange(mid_ter_region-half_ter_region,mid_ter_region+half_ter_region) #ter region centered around ter
    parS=[] #no parS sites
    parS_strength=[]#no parS sites
    L_0=1600*np.exp(growth_rate*(t_doubling-t_replication))#cell size at start of replication, nm. cells grow 10 minutes before replication starts
    lifetime=loop_size/(2*loop_extruder_speed)
    stepRate=loop_extruder_speed/(4600/N)
    backStepRate=0
    if no_bypass:
        return cell("e_coli_no_bypass",N_cut, monomer_size, radius,L_0,growth_rate,lifetime,stepRate,backStepRate,parS,parS_strength,ter_sites,terStrength,t_replication, 0)
    else:
        return cell("e_coli",N_cut, monomer_size, radius,L_0,growth_rate,lifetime,stepRate,backStepRate,parS,parS_strength,ter_sites,terStrength,t_replication)
