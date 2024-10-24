import sys
import numpy as np
from os.path import exists
from numba import jit


# ## Continuous potentials

# In[3]:


# Parameter delta controls the assimetry:
# delta < 0 -> powerstroke
# delta > 0 -> brownian ratchet
# if particle moves to positive cycles

def saw_delta(n = 30, delta = -.6, ymin = -1):
    x = np.linspace(-10, 10, n)
    y = np.zeros(n)
    
    i = int((0.5 - delta/2)*n)
    
    y[:i] = np.linspace(0, ymin, i+1)[:i]
    y[i:] = np.linspace(ymin, 0, len(x) - i)
    return [x, y]

@jit(nopython = True)
def saw_langevin(x, L = 1, delta = -.6, ymin = -1):
    valley = (0.5 - delta/2)*L
    
    # y = ax + b -> y = (ymin/valley) x
    if x <= valley:
        return ymin/valley * x
    else:
        return ymin/(valley - L) * (x - L)

@jit(nopython = True)
def dsawdx(x, L = 1, delta = -.6, ymin = -1):
    valley = (0.5 - delta/2)*L

    # y = ax + b -> y = (ymin/valley) x
    if x <= valley:
        return ymin/valley
    else:
        return ymin/(valley - L)


# ## Discrete potentials

# In[4]:


def flat(n = 30):
    x = np.linspace(0, 1, n+1)[:n]
    x = x + x[1]/2
    y = np.zeros(x.shape)
    
    return [x, y]

def saw_gillespie(n = 30, delta = -.6, ymin = -1):
    x = np.linspace(0, 1, n+1)[:n]
    x = x + x[1]/2
    y = np.zeros(n)
    
    i = int((0.5 - delta/2)*n)
    
    y[:i] = np.linspace(0, ymin, i+1)[:i]
    y[i:] = np.linspace(ymin, 0, len(x) - i)
    return [x, y]


# ## Coupling terms

# In[5]:


@jit(nopython = True)
def coupling(x, xp1, xm1, L, spring_k, spring_l):
    d = xp1 - x - spring_l
    E = spring_k/2*(1 - np.cos(2*np.pi*d/L))
    
    d = x - xm1 - spring_l
    E += spring_k/2*(1 - np.cos(2*np.pi*d/L))
    
    return E

@jit(nopython = True)
def dcouplingdx(x, xp1, xm1, L, spring_k, spring_l):
    d = xp1 - x - spring_l
    dEdx = -spring_k*np.sin(2*np.pi*d/L)*np.pi/L
    
    d = x - xm1 - spring_l
    dEdx += spring_k*np.sin(2*np.pi*d/L)*np.pi/L
    
    return dEdx

@jit(nopython = True)
def asym_coupling(state, M, N, spring_k, spring_l):
    E_k = np.zeros(M)
    for m in range(M):
        dist = state[(m+1)%M] - state[m] - spring_l
        
        E_k[m] = spring_k/2*(1 - np.cos(dist*np.pi/(N/2)))
        
    E = np.sum(E_k)/M
    return E

@jit(nopython = True)
def asym_coupling_split(state, M, N, spring_k, spring_l):
    Ef = np.zeros(M)
    Eb = np.zeros(M)
    for m in range(M):
        dist = state[(m+1)%M] - state[m] - spring_l
        Ef[m] = spring_k/2*(1 - np.cos(dist*np.pi/(N/2)))
        
        dist = state[m] - state[(m-1)%M] - spring_l
        Eb[m] = spring_k/2*(1 - np.cos(dist*np.pi/(N/2)))
        
    return Ef/(2*M), Eb/(2*M)


# ## Non-equilibrium energy

# In[6]:


@jit(nopython = True)
def perturb_langevin(x, energy, pos, spread, pulse = False):
    if pulse:
        if x < pos + spread / 2 and x > pos - spread / 2:
            return energy/spread
        return 0
    
    return np.exp(-(x - pos)**2/(2*spread**2))*energy/(spread*np.sqrt(2*np.pi))

@jit(nopython = True)
def perturb_gillespie(coord, N, M, offset, perturb_index, perturb_energy):    
    ks = np.ones(2*M)
    
    for i in range(M):
        offset_index = (coord[i] + offset*i) % N
        if perturb_index[offset_index, 0]:
            ks[2*i] = np.exp(perturb_energy[offset_index, 0])
        if perturb_index[offset_index, 1]:
            ks[2*i + 1] = np.exp(-perturb_energy[offset_index, 1])
            
    return ks
        


# ## Auxiliary functions

# In[7]:


@jit(nopython = True)
def energy_state_langevin(state, M, L, delta, spring_k, spring_l, flat = False, V_min = -1):
    m = np.arange(M)
    mp1 = (m + 1)%M
    mm1 = (m - 1)%M
    
    E = sum([coupling(state[i], state[mp1[i]], state[mm1[i]], L, spring_k, spring_l) for i in m])/(2*M)
    V = 0
    if not flat:
        V = sum([saw(state[i], L, delta, V_min) for i in m])
    
    return E + V

@jit(nopython = True)
def energy_state_split(state, M, L, delta, spring_k, spring_l, flat = False, V_min = -1):
    m = np.arange(M)
    mp1 = (m + 1)%M
    mm1 = (m - 1)%M
    
    E = sum([coupling(state[i], state[mp1[i]], state[mm1[i]], L, spring_k, spring_l) for i in m])/(2*M)
    V = 0
    if not flat:
        V = sum([saw(state[i], L, delta, V_min) for i in m])
    
    return E, V

@jit(nopython = True)
def energy_particle(x, xp1, xm1, M, L, delta, spring_k, spring_l, flat = False, V_min = -1):
    E = coupling(x, xp1, xm1, L, spring_k, spring_l) / M
    V = 0
    if not flat:
        V = saw(x, L, delta, V_min)
    
    return E + V

@jit(nopython = True)
def denergy_particledx(x, xp1, xm1, M, L, delta, spring_k, spring_l, flat = False, V_min = -1):
    dEdx = dcouplingdx(x, xp1, xm1, L, spring_k, spring_l) / M
    dVdx = 0
    if not flat:
        dVdx = dsawdx(x, L, delta, V_min)
    
    return dEdx + dVdx

def occupancy_langevin(traj, nbins, L = 1):
    M = traj.shape[0]
    freq = np.zeros((nbins, M))
    
    for m in range(M):
        hist = np.histogram(traj[m], bins = nbins, range = (0, L))[0]
        hist = hist / traj.shape[1]
        freq[:, m] = hist
    
    return freq

@jit(nopython = True)
def prob_matrix_langevin(nbins, M, L, spring_k, spring_l = 0, delta = -.6, flat = False, kT = 1, V_min = -1):
    E_tot = np.zeros(nbins**M, dtype = np.float64)
    x = np.linspace(0, L, nbins+1)[:nbins]
    x = x + x[1]/2
    
    idx = np.zeros(M, dtype = np.int_)

    for k in range(nbins**M):
        K = k
        for m in range(M):
            idx[m] = K % nbins
            K = K // nbins
            
        E_tot[k] = energy_state(x[idx], M, L, delta, spring_k, spring_l, flat, V_min)

    p = np.exp(-E_tot/kT)
    return p/np.sum(p)

def prob_marginal(p, nbins, M, m = 0):
    p = p.reshape(np.array([nbins for i in range(M)]))
    
    if M != 1:
        ms = np.arange(M)
        ms = tuple(ms[ms != m])
        marginal_p = np.sum(p, axis = ms)
    else:
        marginal_p = p
    
    return marginal_p/np.sum(marginal_p)

@jit(nopython = True)
def energy_state_gillespie(state, energy, offset, spring_k, spring_l):
    M = len(state)
    N = len(energy)
    
    E = asym_coupling(state, M, N, spring_k, spring_l)
    V = sum([energy[(state[m] + offset*m)%N] for m in range(M)])
    
    return E+V


@jit(nopython = True)
def prob_matrix_gillespie(x, energy, M, kT = 1, offset = 0, spring_k = 10, spring_l = 0):
    N = len(x)
    E_tot = np.zeros(N**M, dtype = np.float64)
    
    idx = np.zeros(M, dtype = np.int_)

    for k in range(N**M):
        K = k
        for m in range(M):
            idx[m] = K % N
            K = K // N
        
        E_tot[k] = energy_state_gillespie(idx, energy, offset, spring_k, spring_l)

    p = np.exp(-E_tot/kT)
    return p


@jit(nopython = True)
def W_matrix(prob_matrix, N, M):
    W = np.zeros((N**M, N**M))
    
    idx = np.zeros(M, dtype = np.int_)
    factor = np.array([N**m for m in range(M)])
    for k in range(N**M):
        K = k
        for m in range(M):
            idx[m] = K % N
            K = K // N
        
        for m in range(M):
            i = idx[m]
            for i_ in [(i-1)%N, (i+1)%N]:
                idx_ = np.copy(idx)
                idx_[m] = i_
                k_ = int(np.sum(idx_*factor))
                
                delta = (prob_matrix[k]/prob_matrix[k_])**0.5
                W[k, k_] = delta
                
            
    for k in range(N**M):
        W[k,k] = - np.sum(W[:,k])

    return W


@jit(nopython = True)
def occupancy_gillespie(trajectory, time, N, M):
    freq = np.zeros(N**M)
    factor = np.array([N**m for m in range(M)])
    
    t = time[1:] - time[:-1]
    for i in range(len(time)-1):
        k = int(np.sum(trajectory[:, i]*factor))
        freq[k] += t[i]
    freq /= np.sum(freq)
    return freq


# ## Langevin trajectories

# In[8]:


@jit(nopython = True)
def langevin(x0, L, t_max, dt = 1e-2, 
             flat = False, delta = -.6, V_min = -1, 
             spring_k = 0, spring_l = 0, 
             kT = 1, seed = 14, 
             perturb_E = 0, perturb_x = 0, 
             perturb_s = 1, perturb_pulse = False):
    
    M = len(x0)
    time = np.arange(0, t_max, dt)
    x = np.zeros((len(time), M))
    cycle = np.zeros((len(time), M))
    
    np.random.seed(seed)
    x[0] = x0
    
    mp1 = (np.arange(M) + 1) % M
    mm1 = (np.arange(M) - 1) % M
    
    
    for t in range(1, len(time)):
        for m in range(M):
            v = -denergy_particledx(x[t-1, m], x[t-1, mp1[m]], x[t-1, mm1[m]], 
                                    M, L, delta, spring_k, spring_l, flat, V_min)/kT
            
            v += perturb_langevin(x[t-1, m], perturb_E, perturb_x, perturb_s, perturb_pulse)
            
            d = np.random.normal()*2**0.5
            dx = v*dt + d*dt**0.5
            
            if abs(dx) > L:
                print("Bad choice of dt, a step was larger than L")
                return x[:t].T, time[:t], cycle[:t].T
                
            xnew = x[t-1, m] + dx
            
            if xnew > L:
                xnew -= L
                cycle[t, m] = cycle[t-1, m] + 1
            elif xnew < 0:
                xnew += L
                cycle[t, m] = cycle[t-1, m] - 1
            else:
                cycle[t, m] = cycle[t-1, m]
            
            x[t, m] = xnew
                
    return x.T, time, cycle.T


# ## Gillespie trajectories

# In[9]:


@jit(nopython = True)
def next_step_rates(old_coord, r, kT, N, M, offset,
                    energy, spring_k, spring_l, tau,
                    perturb_index, perturb_energy):  
    
    new_coord = np.zeros((2*M, M), dtype = np.int_)
    effective_tau = np.zeros((2*M), dtype = np.float_)
    
    for m in range(M):
        i = old_coord[m]
        
        # Assign same position to the possible next states
        new_coord[2*m] = old_coord
        new_coord[2*m + 1] = old_coord
        
        # Change only the position of mth particle
        new_coord[2*m, m] = (i+1)%N
        new_coord[2*m + 1, m] = (i-1)%N
        
        # Calculate the effective tau
        effective_tau[2*m] = (tau[i] + tau[(i+1)%N]) / 2
        effective_tau[2*m + 1] = (tau[i] + tau[(i-1)%N]) / 2
            
    # Calculate rates from conservative energies
    ks = np.zeros(2*M)
    old_energy = energy_state_gillespie(old_coord, energy, offset, spring_k, spring_l)
    
    for m2 in range(2*M):
        new_energy = energy_state_gillespie(new_coord[m2], energy, offset, spring_k, spring_l)
        delta = np.exp((old_energy - new_energy)/(2*kT))
        
        ks[m2] = delta / (delta + 1/delta)
        ks[m2] = ks[m2] / effective_tau[m2]
        
    # Energetic contribution from  perturbation
    if np.sum(perturb_index**2) != 0:
        non_eq_rate = perturb_gillespie(old_coord, N, M, offset, 
                              perturb_index, perturb_energy)
        ks *= non_eq_rate
    
    
    t = -np.log(r[0])/(np.sum(ks))
    
    c = 0
    cycle = np.zeros(M)
    for m2 in range(2*M):
        c += ks[m2]/np.sum(ks)
        if r[1] < c:
            if new_coord[m2, m2//2] == 0 and old_coord[m2//2] == N-1:
                cycle[m2//2] = 1
            elif new_coord[m2, m2//2] == N-1 and old_coord[m2//2] == 0:
                cycle[m2//2] = -1
            
            return new_coord[m2], t, cycle
        

def trajectory(x, energy, rep, i0, precalculate_W = True,
               offset = 0, spring_k = 10, spring_l = 0,
               kT = 1, seed = 14, tau = None,
               perturb_index = None, perturb_energy = None):
    
    M = len(i0)
    N = len(energy)
    traj = np.zeros((rep, M), dtype = np.int_)
    time = np.zeros(rep)
    cycle = np.zeros((rep, M))
    
    if type(perturb_index) == type(None):
        perturb_index = np.zeros((len(energy), 2))
        perturb_energy = np.zeros((len(energy), 2))
    
    if type(tau) == type(None):
        tau = np.ones(len(energy))
    
    np.random.seed(seed)
    random = np.random.uniform(size = (rep, 2))
    traj[0] = i0
    
    @jit(nopython = True)
    def helper(traj, time, cycle):
        for r in range(1, rep):
            i, t, c = next_step_rates(traj[r-1, :], random[r], kT, N, M, offset,
                                      energy, spring_k, spring_l, tau,
                                      perturb_index, perturb_energy)

            time[r] = time[r-1] + t
            for m in range(M):
                cycle[r, m] = cycle[r-1, m] + c[m]
            traj[r, :] = i
            
        return traj.T, time.T, cycle.T
    
    return helper(traj, time, cycle)


# argument list:

#01. N
N = int(sys.argv[1])
#02. M
M = int(sys.argv[2])
#03. niter
niter = int(sys.argv[3])
#04. seed
seed = int(sys.argv[4])*10000
#05. init => 1: zeros, 2: ordered, 3: shuffled
init_opt = int(sys.argv[5])
#06. kT
kT = float(sys.argv[6])
#07. spring_k
k = float(sys.argv[7])
#08. spring_l (-1: sequential spring_l)
l = int(sys.argv[8])
#09. offset
off = int(sys.argv[9])
#10. Vmax (-1: original Vmax)
Vmax = float(sys.argv[10])
#11. Drive (kT)
mu_ATP = float(sys.argv[11])
#12. Tau free (s) DNA free
tauF = float(sys.argv[12])
#13. Tau bound (s) DNA bound
tauB = float(sys.argv[13])

#14. file
output = sys.argv[14]
#15. complete output ('y' = yes, 'n' = no)
complete = sys.argv[15]
#16. output trajectory ('y' = yes, 'n' = no)
output_traj = sys.argv[16]
#17. output time ('y' = yes, 'n' = no)
output_time = sys.argv[17]
#18. output cycle ('y' = yes, 'n' = no)
output_cycl = sys.argv[18]
#19. number of repetitions (only for complete = 'n')
nrep = int(sys.argv[19])

# Processing
if l == -1:
    l = N/M

if init_opt == 1:
    init = np.zeros(M)
elif init_opt == 2:
    init = np.linspace(0, N, M+1)[:M]
    init = np.int_(init)
elif init_opt == 3:
    np.random.seed(seed)
    init = np.random.randint(0, N+1, M)

x = np.linspace(0, 1, N+1)[:N]
y = np.load("../energy_landscape.npy")
if Vmax != -1:
    y = y/np.max(y)*Vmax

# Drive
perturb_index = np.zeros((N, 2))

# Pi release
#perturb_index[27, 0] = 1
#perturb_index[28, 1] = 1
# ADP removal
perturb_index[10, 0] = 1
perturb_index[11, 1] = 1
# ATP binding
perturb_index[13, 0] = 1
perturb_index[14, 1] = 1

perturb_energy = perturb_index * mu_ATP / np.sum(perturb_index)

# Tau
tau = np.zeros(N)
tau[0:10] = tauF
tau[10:N] = tauB

# Distance (DNA translocation, in nm)
dist = 2 * 0.34

# Calc
if complete != "n":
    traj, time, cycle = trajectory(x, y, niter, init, False,
                                   off, k, l, kT, seed, tau,
                                   perturb_index, perturb_energy)
    speed = np.sum(cycle[:,-1]) * dist / time[-1]

elif complete == "n":
    traj = np.zeros((nrep, M))
    time = np.zeros(nrep)
    cycle = np.zeros((nrep, M))
    speed = np.zeros(nrep)

    for rep in range(nrep):
        if init_opt == 3:
            np.random.seed(seed + rep)
            init = np.random.randint(0, N+1, M)

        traj_, time_, cycle_ = trajectory(x, y, niter, init, False,
                                          off, k, l, kT, seed + rep, tau,
                                          perturb_index, perturb_energy)
        traj[rep] = traj_[:,-1]
        time[rep] = time_[-1]
        cycle[rep] = cycle_[:,-1]
        speed[rep] = np.sum(cycle[rep]) * dist / time[rep]


# Output
if output_traj != "n":
    with open(output + "_traj.npy", "wb") as f:
        np.save(f, traj)

if output_time != "n":
    with open(output + "_time.npy", "wb") as f:
        np.save(f, time)
    with open(output + "_speed.npy", "wb") as f:
        np.save(f, speed)

if output_cycl != "n":
    with open(output + "_cycle.npy", "wb") as f:
        np.save(f, cycle)
