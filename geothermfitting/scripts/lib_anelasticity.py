import numpy as np
from scipy.optimize import minimize_scalar, brent
from os.path import join
from os import getenv, sys

class anelasticity_model_cammarano_type_1(object):
   # This class is used to build anelasticity model by Cammarano et al 2003
   # but only for Q1 to Q4 models
   from scipy.spatial import cKDTree
   u_l_bord=660e3; # [upper and lower mantle boundary]
   omega = 1;  # seismic frequency
   def __init__(self, B, g, a):
      if len(B)!=2 or len(g)!= 2 or len(a)!=2:
         raise('All input parameters should be of length 2, for upper and lower mantle!')
      self.B = B;
      self.g = g;
      self.a = a;
   def create_Q_radial(self, depths, adiabat_temp):
      from scipy.interpolate import interp1d   

      if len(adiabat_temp) !=len(depths):
         raise('input adiabat_temp is not consistent with size')

      self.dpths = depths; 
      B_arr = np.zeros(len(depths));
      a_arr = np.zeros(len(depths));
      g_arr = np.zeros(len(depths));

      B_arr[:]                        = self.B[1];
      B_arr[self.dpths<self.u_l_bord] = self.B[0];
      g_arr[:]                        = self.g[1];
      g_arr[self.dpths<self.u_l_bord] = self.g[0];
      a_arr[:]                    = self.a[1];
      a_arr[self.dpths<self.u_l_bord] = self.a[0];

      fi_rad, fi_solid = np.loadtxt('/Users/fredrichards/space1/LLSVP_Density/DATA/SOLIDUS/SiaG_solidus.csv',\
              comments='#', delimiter=',', usecols=(0,2), unpack=True);
      my_intepr1d = interp1d(fi_rad, fi_solid,fill_value='extrapolate');
      solid_temp = my_intepr1d(self.dpths)
      self.my_solid = solid_temp; 
      Q_s =B_arr*(self.omega)**a_arr*np.exp((a_arr*g_arr*solid_temp)/(adiabat_temp));
      return Q_s

   def create_Q_2d(self, depths, temps):
       
      B_arr = np.zeros((len(depths), len(temps)));
      a_arr = np.zeros((len(depths), len(temps)));
      g_arr = np.zeros((len(depths), len(temps)));

      B_arr[:,:]                      = self.B[1];
      B_arr[depths<self.u_l_bord,:]   = self.B[0];
      g_arr[:,:]                      = self.g[1];
      g_arr[depths<self.u_l_bord,:]   = self.g[0];
      a_arr[:,:]                      = self.a[1];
      a_arr[depths<self.u_l_bord,:]   = self.a[0];

      # make a copy of a matrix to be used in libseis_min lib
      self.a_matrix = a_arr;

      # read in the solidus
      fi_rad, fi_solid = np.loadtxt('/Users/fredrichards/space1/LLSVP_Density/DATA/SOLIDUS/SiaG_solidus.csv',\
               comments='#', delimiter=',', usecols=(0,2), unpack=True);

      my_tree = self.cKDTree(fi_rad[:,np.newaxis], leafsize=10);
      w_idx = my_tree.query(depths[:,np.newaxis], k=1)
      #solid_temp = np.sum((fi_solid[w_idx[1]]*(1./w_idx[0])),axis=1)/np.sum(1./w_idx[0],axis=1)
      solid_temp = fi_solid[w_idx[1]];

      temps_x, solid_temp_x = np.meshgrid(temps, solid_temp);
      
      Q_s =B_arr*(self.omega)**a_arr*np.exp((a_arr*g_arr*solid_temp_x)/(temps_x));

      return Q_s

class anelasticity_model_cammarano_type_2(object):
   # This class is used to build anelasticity model by Cammarano et al 2003
   # but only for Q1 to Q4 models
   from scipy.spatial import cKDTree
   u_l_bord=660e3; # [upper and lower mantle boundary]
   omega = 1;# seismic frequency
   R=8.3145
   def __init__(self, A, d, E, V, a):
      if len(A)!=2 or len(d)!= 2 or len(E)!=2 or len(V)!=2 or len(a)!=2:
         raise('All input parameters should be of length 2, for upper and lower mantle!')
      self.A = A;
      self.d = d;
      self.E = E;
      self.V = V;
      self.a = a;
   def create_Q_radial(self, depths, adiabat_temp):
      if len(adiabat_temp) != len(depths):
         raise('input adiabat_temp is not consistent with size')

      self.dpths = depths; 
      _, _, self.pres = dens_grav_pres_prof(6371.e3-self.dpths);

      B_arr = np.zeros(len(depths));
      a_arr = np.zeros(len(depths));
      E_arr = np.zeros(len(depths));
      V_arr = np.zeros(len(depths));

      B_arr[:]                         = self.A[1]*(self.d[1]**self.a[1]);
      B_arr[self.dpths<self.u_l_bord]  = self.A[0]*(self.d[0]**self.a[0]);

      a_arr[:]                         = self.a[1];
      a_arr[self.dpths<self.u_l_bord]  = self.a[0];

      E_arr[:]                         = self.E[1];
      E_arr[self.dpths<self.u_l_bord]  = self.E[0];

      V_arr[:]                         = self.V[1];
      V_arr[self.dpths<self.u_l_bord]  = self.V[0];

      Q_s =B_arr*(self.omega)**a_arr*np.exp((a_arr*(E_arr+self.pres*V_arr))/(self.R*adiabat_temp));

      return Q_s

def load_q_model(q_target=None, dpths=None,  q_mu_lb=20.0, q_mu_ub=4900.0,\
                                             q_p_lb=20.0, q_p_ub=57827.0):
   """
   This routine loads a requested Q model on the given depths

   These values will be used to limit the variations of 
       reference model Q factors:
         q_mu_lb  = 20.0; q_mu_ub = 4900.0; q_p_lb = 20.0; q_p_ub = 57827.0;

   Dpths for which a Q model should be given out:
         dpths = np.linspace(0,2890., 257);

   Which model should be loaded
         q_target = 'PREM';
   """
   import h5py
   import numpy as np
   from os.path import join
   from os import getenv
   from scipy.interpolate import interp1d
   
   # Q_MODELS is the hdf5 collection of all the available Q models
   q_hdf5_name = '/Users/fredrichards/space1/LLSVP_Density/DATA/Q_collection.hdf5';
   # open the hdf5 file
   hdf5_fi = h5py.File(q_hdf5_name, mode='r');
   # see which models are available
   avail_Qs = list(hdf5_fi.keys());
   # see if the model is available
   if q_target in avail_Qs:
         q_model = hdf5_fi.get(q_target);
         mod_depth   = np.array(q_model.get('depth'));
         mod_Q_mu    = np.array(q_model.get('Q_mu'));
         mod_Q_p     = np.array(q_model.get('Q_kappa'));
         hdf5_fi.close()
   else:
         hdf5_fi.close()
         raise ValueError(\
               'The requested Q model is not available. '\
               +'Available Q models:\n'+', '.join(avail_Qs))
   
   # limit the variations of Q factors
   mod_Q_mu[mod_Q_mu <= q_mu_lb] = q_mu_lb; mod_Q_mu[mod_Q_mu >= q_mu_ub] = q_mu_ub;
   mod_Q_p [mod_Q_p  <= q_p_lb ] = q_p_lb;  mod_Q_p [mod_Q_p  >= q_p_ub ] = q_p_ub;
   
   # interpolate the model q factor on the requested dpths 
   interp_func = interp1d(mod_depth, mod_Q_mu, kind='slinear');
   q_mu = interp_func(dpths);
   interp_func = interp1d(mod_depth, mod_Q_p , kind='slinear');
   q_p  = interp_func(dpths);

   return q_mu, q_p


class anelasticity_model(object):
   """
      This class creates a model of anelasticity 
      as done by Richards, Hoggard, and White (2018), 
      followed by a formalism introduced in Takei 2017
   """
   from scipy.spatial import cKDTree 
   def __init__(self, eta_0, mu0, dmudT, dmudp, act_eng, act_vol,\
                  solgrad, sol50, alpha_B, A_B, delphi, temps, dens_flg): 
   
      # backgraound viscosity
      self.eta_0 = eta_0;
      # unrelaxed compliance parameters
      self.mu0      = mu0;
      self.dmudT    = dmudT;
      self.dmudp    = dmudp;
      # activation energy
      self.act_eng  = act_eng;
      # activation volume
      self.act_vol  = act_vol;
      # solidus temperature gradient
      self.solgrad  = solgrad;
      # solidus temperature gradient
      self.sol50  = sol50;
      # Anelasticity model parameters
      self.A_B      = A_B;
      #  frequency dependence of anelasticity
      self.alpha_B  = alpha_B;
      #  ratio of melt
      self.delphi   = delphi

      # seismic period
      self.period = 100;
      self.tau_p  = 6.e-5;
      self.T_eta  = 0.94;

      # what is temperature resoution
      # an array of temperature
      self.temps = temps; 
      self.dens_flg = dens_flg; 

   def temp_2_dens(self, temp_fld):
      """
         Convert a temperature field previously computed by vs_2_temp
         to density using the already known depths in the system.

         We reshape all the input temp_fld to a 1D array.
         make a Tree out of the self.temps array, so we can easily find 
         all the neighbouring points in the field. and convert 
         all of them to density! 
      """
      # make sure temp_fld size is consistent with lengths of depths and temps
      if len(temp_fld) != len(self.dpths):
         raise('Ooops! temp_fld should be of length (%i)' %(len(self.dpths)));

      # define the output dens_fld based on temp_fld
      dens_fld = [None]*len(temp_fld);

      # average the two neighbouring points and reshape the field at each depth 
      for i in range(len(temp_fld)):
         orig_shape = np.shape(temp_fld[i]);
         ts=np.reshape(temp_fld[i],np.size(temp_fld[i]))
         print('Converting layer at %4.0f km to density.' %(self.dpths[i]/1.e3));
         # from here, we want to find the neighbouring points for temperatures
         my_temp_tree = self.cKDTree(self.temps[:,np.newaxis], leafsize=10);
         # here we query the neighbouring point for all the temperatures
         temp_tree_res = my_temp_tree.query(\
                  temp_fld[i].reshape((np.prod(np.shape(temp_fld[i])),1)),k=2);
         # make sure that no distance is zero, so we do not end up with 1./0
         temp_tree_res[0][:,0 ][(temp_tree_res[0]<=1e-4)[:,0]] =  1.0;
         temp_tree_res[0][:,1:][(temp_tree_res[0]<=1e-4)[:,0]] = 1e10;
         dens_fld[i] = \
            np.sum((1./temp_tree_res[0])*self.rho[i,temp_tree_res[1]],axis=1)/\
                  np.sum(1./temp_tree_res[0],axis=1);
         dens_fld[i] = dens_fld[i].reshape(orig_shape);
         

      return dens_fld

   def vs_2_temp(self, depths, seis_model):
      """
         This routine converts a seismic model (list of length len(depths) to temperature)
      """

      if True in depths > 400.0e3:
         print("Warning! You are trying to use this model at depths larger than 400 km");

      if len(seis_model) != len(depths):
         raise('seis_model should be an array of length %i!\n' %len(depths));
      # set the depths for the whole class
      self.dpths = depths;
      # output temperature fld has the same dimension as the input seis_model
      temp_fld = [None]*len(seis_model);
      
      # compute seismic velocities
      self.vs, _, _, _ = self.compute_vs_vsu_anel_qual\
                                 (self.dpths, self.temps);

      # go through every layer of the seismic model and convert it 
      # to temperature
      for i in range(len(seis_model)):
         print('Converting layer at %4.0f km to temperature.' %(self.dpths[i]/1.e3));
         orig_shape = np.shape(seis_model[i]);
         # make a tree of seismic values at that depth
         seismic_layer_tree = self.cKDTree(self.vs[i,:].reshape((len(self.vs[i,:]),1)), leafsize=10);
         # get the two neighbouring points for each seismic value
         seis_tree_res = seismic_layer_tree.query(\
                  seis_model[i].reshape((np.prod(np.shape(seis_model[i])),1)), k=2);
         # make sure that no distance is zero, so we do not end up with 1./0
         seis_tree_res[0][:,0 ][(seis_tree_res[0]<=1e-4)[:,0]] = 1.0;
         seis_tree_res[0][:,1:][(seis_tree_res[0]<=1e-4)[:,0]] = 1e10;
         # make a weighted arithmetic average for the temperature of each point
         temp_fld[i] = np.sum((1./seis_tree_res[0])*self.temps[seis_tree_res[1]],axis=1)/\
                  np.sum(1./seis_tree_res[0],axis=1);
         # reshape each field to have the same shape as the input seis_model
         temp_fld[i] = temp_fld[i].reshape(orig_shape);
      return temp_fld

   def compute_vs_vsu_anel_qual(self, depths, temps):
      from scipy.special import erf
      """
         input: 
            depths [km]
            temps  [K]
            anel_prms typew
      """
      # Thermodynamic gas constant
      R = 8.3145   
      # Earth's radius
      R_earth = 6371.e3; 
   
      # model radii
      rshl = R_earth - depths;
      
      # a profile of density, gravity, and pressure for the given depths(radii)
      # dens, grav, pres = dens_grav_pres_prof(rshl);
      pres = depths*(1./30.)*1e6

      # load a model of mantle material solidus
      T_solidus = load_solidus(dpths=depths*1e-3, mode='const_grad', solgrad=self.solgrad, T50=self.sol50)
      
      # load a 2d field of viscosity for different temperatures and pressures(dpths) 
      viscosity = self.function_viscosity(temps, T_solidus, pres,\
                        T_eta=self.T_eta, back_visc=self.eta_0,\
                        act_eng=self.act_eng, act_vol=self.act_vol)
      # load a 2d field of unrelaxed compliances for different temperatures and pressures(dpths)
      unrlx_compl = self.unrelaxed_compl(temps=temps, pres=pres, mu0=self.mu0,\
                                    dmudT=self.dmudT, dmudp=self.dmudp, \
                                    T_ref = 273, pres_ref = 0)
      ## Just to check values of A_eta 
      #A_eta = self.function_A_eta(temps, T_solidus, T_eta);
      
      # load a 2d field of A_p for different temperatures and solidus temperatures(dpths)
      A_p = self.function_A_p(temps=temps, T_solid=T_solidus, phi=self.delphi);
      
      # load a 2d field of sigma_p for different temperatures and solidus temperatures(dpths)
      sigma_p = self.function_sigma_p(temps, T_solid=T_solidus);
      
      # compute 2d maxwellian time, given viscosity and unrelaxed compliances
      tau_maxwell = viscosity*unrlx_compl
      
      # include depth dependence of period
      period=np.zeros(np.shape(tau_maxwell))
      for i in np.arange(0,np.size(depths),1):
          period[i,:]=(3.*(depths[i]/1.e3))/4.2

      # non-dimensional period: see Equation (21) of Takei 2017
      # period_non = self.period/(2*np.pi*tau_maxwell)
      period_non = period/(2*np.pi*tau_maxwell)
      
      # real compliance, see Equation (26) of Takei 2017
      J_1 = unrlx_compl*(1+self.A_B*(period_non)**self.alpha_B/self.alpha_B\
               + np.sqrt(2*np.pi)/2*A_p*sigma_p*\
               (1.-erf(np.log(self.tau_p/period_non)/(np.sqrt(2)*sigma_p))))
      # imaginary compliance, see Equation(26) of Takei 2017
      J_2 = unrlx_compl*np.pi/2*(self.A_B*(period_non**self.alpha_B)\
               + A_p*np.exp(-(np.log(self.tau_p/period_non)**2)/(2*(sigma_p**2))))\
               + unrlx_compl*period_non
   
      # compute attenuation
      qual_fac = J_1/J_2;

      # here we compute densities based on temperature/pressure variations
      # Bulk modulus Gpa, and thermal expansion
      
      if self.dens_flg==0:
          # calculate the density variations
          bmod     = 115.2e9
          alpha_T  = 3.59e-5
          
          # calculate the density variations
          self.rho = density_fld(temps, pres, alpha=alpha_T, bmod=bmod,\
                   dens_ref= 3291., T_ref=(600.+273.), pres_ref=0)
                   
      elif self.dens_flg==1:
          a0=2.832e-5
          a1=0.758e-8
          K0=130.e9
          KT=4.8
          grun=6.
          AY=1.0
          CY=2.0
          rhom=3330.

          # calculate the density variations
          self.rho = density_fld_ga(temps, pres, grun=grun, a0=a0, a1=a1,\
                 dens_ref=rhom, K0=K0, KT=KT, A=AY, C=CY);
      else:
          print ('Warning! Density flag %i not recognised.' %(self.dens_flg));

      # compute shear wave velocity, frequency dependent
      Vs = 1./(np.sqrt(self.rho*J_1)*1.e3); # km/s
   
      # compute shear wave velocity, unharmonic
      Vsu = 1./(np.sqrt(self.rho*unrlx_compl)*1.e3); # km/s
   
      # compute the anelastic correction
      anel_corr = Vsu - Vs
   
      return Vs, Vsu, anel_corr, qual_fac

   # def unrelaxed_compl
   def unrelaxed_compl(self, temps, pres, mu0, dmudT, dmudp, T_ref = 273, pres_ref = 0):
      """
         Compute the unrelaxed compliance base on Equation 27 of Takei 2017
         temps:   range of temperatures
         pres:    range of pressures with changing depth (depths)
         mu0, dmudT, dmudp are coefficients in Table 2 of Takei 2017
         T_ref, and pres_ref: are the reference temperatures at which these values are obtained
      """
   
      temps_x, pres_x = np.meshgrid(temps, pres);
      J_u = 1. / (mu0 + dmudp*(pres_x-pres_ref) + dmudT*(temps_x-T_ref));
      return J_u

   # def function_A_eta
   def function_A_eta(self, temps, T_solid, T_eta, gamma=5.0, phi=0.0, lambda_fac=0.0):
      """
      function_A_eta(): Equation (29) of Takei 2017,
      inputs:
         temps:   an array of the possible tempereature
         T_solid: solidus temperature at each depth 
         t_eta:   normalized temperature above which the activation energy is H+\Delta H
         gamma:   the factor of extra reduction 
         phi:     melt fraction
      """
      A_eta = np.zeros((len(T_solid), len(temps)))
   
      temps_x, t_solid_x = np.meshgrid(temps, T_solid);
      T_holo = temps_x/t_solid_x;
      A_eta[T_holo<T_eta] = 1.0;
      A_eta[T_holo>=T_eta] =\
                  np.exp((-1.*((T_holo[T_holo>=T_eta]-T_eta)/\
                  (T_holo[T_holo>=T_eta]-(T_holo[T_holo>=T_eta]*T_eta))))*np.log(gamma));
      A_eta[T_holo>=1.0] =\
                  (1./gamma)*np.exp(-lambda_fac*phi);
      return A_eta 
   
   # def func_beta
   def func_beta(self, phi=0.0):
      """
         This function is supposed to model the impact of melt 
         at this moment is set to zero
      """
      beta = 0.0;
      return beta
   
   # def function_A_p
   def function_A_p(self, temps, T_solid, phi):
      """
      function_A_p(): Equation (24) of Takei 2017,
      inputs:
         temps:   an array of the possible tempereature
         T_solid: solidus temperature at each depth 
         phi:     melt fraction
      """
      A_p   = np.zeros((len(T_solid),len(temps)))
   
      temps_x, t_solid_x = np.meshgrid(temps, T_solid);
      T_holo = temps_x/t_solid_x;
   
      A_p[T_holo< 0.91] = 1e-2;
      A_p[T_holo>=0.91] =\
               1e-2+0.4*(T_holo[T_holo>=0.91]-0.91);
      A_p[T_holo>=0.96] =\
               3e-2;
      A_p[T_holo>=1.0] =\
               3e-2+self.func_beta(phi);
   
      return A_p
   
   # def function_sigma_p
   def function_sigma_p(self, temps, T_solid):
      """
      function_sigma_p(): Equation (25) of Takei 2017,
      inputs:
         temps:   an array of the possible tempereature
         T_solid: solidus temperature at each depth 
      """
      sigma_p   = np.zeros((len(T_solid),len(temps)))
   
      temps_x, t_solid_x = np.meshgrid(temps, T_solid);
      T_holo = temps_x/t_solid_x;
   
      sigma_p[T_holo< 0.92] = 4;
      sigma_p[T_holo>=0.92] = 4 + 37.5*(T_holo[T_holo>=0.92]-0.92);
      sigma_p[T_holo>=1.0 ] = 7;
   
      return sigma_p
   def function_viscosity(self, temps, T_solid, pres, T_eta, back_visc, act_eng, act_vol):
      """
      function_viscosity(): works out the viscosity given Equation (28) of Takei 2017,
      inputs:
         temps:      an array of the possible tempereature
         T_solid:    solidus temperature at each depth 
         pres:       Pressure, which in this case can be translated into depth
         T_eta:      Threshold as given by Takei 2017, for computing A_eta (see function_A_eta)
         back_visc:  a radial backgroud viscosity, should have the same dimensions as T_solid
         act_eng,    act_vol: activation energy and volume
      """
   
      # Gas constant
      R = 8.3145
      # Reference pressure
      P_ref = 1.5e9
      # Reference temperature
      T_ref = 1473.
      
      rad_visc = np.zeros(len(pres));
      if type(back_visc) == float:
         rad_visc[:] = back_visc;
      elif type(back_visc) == np.ndarray and \
                  type(back_visc)==list and len(back_visc)==len(pres):
         rad_visc = back_visc;
      else:
        raise ValueError(\
              'Length of the inpute back_visc should'+\
              ' either be %i or a constant viscosity for all depths' %len(pres));
   
      # temps are the same at different depths
      # pres and rad_visc are changing with depth, but have a single value at each depth
      temps_x, pres_x = np.meshgrid(temps, pres);
      _, rad_visc_x = np.meshgrid(temps, rad_visc);
   
      # compute A_eta (eq 29 of Takei 2017)
      A_eta = self.function_A_eta(temps, T_solid, T_eta);
   
      # compute viscosity (eq 28 of Takei 2017)
      visc = (rad_visc_x*np.exp((act_eng/R)*(1./temps_x-1./T_ref))*\
               np.exp((act_vol/R)*(pres_x/temps-P_ref/T_ref)))*A_eta

      return visc
# END of class   
   
# def mckenzie_bickle_solidi
def mckenzie_bickle_solidi(dpths):
   """   
      Computes the solidus temperature profile
         as modeled by McKenzie and Bickle, can be used for a single depth
      dpths: Depths at which the solidus temperature is computed
   """
   if np.prod(np.shape(dpths)) > 1: 
        raise ValueError(\
              'MCKENZIE routine works only for a single value of depth.')
      
   # Solidus parameters:
   T0=1373.1;
   Bs=1.2e-2;
   As=4.968e-4;
   Cs=136.05;
   #
   P=dpths/30.0;
   T1=(Cs*P)+T0
   Pd=(T1-T0)/Cs+As*np.exp(Bs*(T1-T0))-P
   dPdT=1./Cs+As*Bs*np.exp(Bs*(T1-T0))
   dT=-Pd/dPdT
   T1=T1+dT
   err=np.abs(dT/T1)
   for N in range(1,10):
      if np.abs(dT) <= 0.1:
         break;
      Pd=(T1-T0)/Cs+As*np.exp(Bs*(T1-T0))-P
      dPdT=1./Cs+As*Bs*np.exp(Bs*(T1-T0))
      dT=-Pd/dPdT
      T1=T1+dT
   Temp_solid=T1
   dPdT=1./Cs+As*Bs*np.exp(Bs*(Temp_solid-T0))
   dTs=1./dPdT
   return Temp_solid

# def load_solidus
def load_solidus(dpths=None, mode=None, Tp_katz=None, XH2O_katz=None, solgrad=None, T50=None):
   """   This routine loads in the solidus temperature array, 
         based on the input mode   
      dpths:   depths at which solidus will be loaded
      mode:    which model/dataset to be used
         """
   from scipy.interpolate import interp1d
   # case: const_grad
   #     Solidus with a constant gradient
   if mode == 'const_grad':
      # Check if we have the necessary parameters
      if T50==None or solgrad==None:
         raise ValueError(\
           'load_solidus: for constant gradient solidus T50 and solgrad should be given!');
      # Computation of the solidus
      T_sol = T50 + (solgrad*(dpths-50.0))
   # case: MB
   #     Solidus based on McKenzie and Bickle 2003
   elif mode == 'MB':
      T_sol = np.zeros(len(dpths));
      for i in range(len(dpths)):
        T_sol[i] = mckenzie_bickle_solidi(dpths[i]);
   # case: katz
   #     Solidus of Katz et al. 2013
   #     possible inputs for Tp_katz & XH2O_katz
   #                             1306., 500
   #                             1322., 112
   #                             1334., 0
   elif mode == 'katz':
      if Tp_katz == None or XH2O_katz == None:
         raise ValueError(\
           'load_solidus: Either of Tp_katz and XH2O_katz are empty!');
      # load a profile based on Tp_katz and XH2O_katz
      katz_dpth, katz_T = np.loadtxt(join('/Users/fredrichards/space1/LLSVP_Density/DATA/katz_solidi/',\
                     'zTF_'+str('%i' %Tp_katz)+'_'+str('%i' %XH2O_katz)+'.txt'),\
                     usecols=(0,1), unpack=True);
      # interpolate the profile onto the given depths
      T_sol = interp1d(\
            katz_dpth[katz_dpth.argsort()], katz_T[katz_dpth.argsort()],\
                        fill_value='extrapolate')(dpths);
   # case: raise an error if input is not clear 
   else:
     raise ValueError(\
           'Available Q models: const_grad, MB, STH, given input:'+str(mode))

   T_sol = T_sol + 273.# conversion to kelvin
   return T_sol


def density_fld(temps, pres, alpha, bmod, dens_ref= 3291., T_ref=(600+273), pres_ref=0):
    
   """
      This routine generates a 3d field of density 
      based on pressure and temperature variations
      over a reference density profile(let's say PREM)
   """
   # 
   temps_x, pres_x  = np.meshgrid(temps, pres);
   density = dens_ref*((1.-alpha*(temps_x - T_ref))+((pres_x - pres_ref)/bmod));
   return density
 
def dV(x,K0,KT,pres_x):
    return np.abs((K0*(3./2.)*(x**(7./3.)-x**(5./3.))*(1.+(((3./4.)*(KT-4.))*(x**(2./3.)-1.))))-pres_x)
      
def density_fld_ga(temps, pres, grun, a0, a1, dens_ref, K0, KT, A, C):
   """
      This routine generates a 3d field of density 
      based on pressure and temperature variations
      over a reference density profile(let's say PREM)
   """
   temps_x, pres_x  = np.meshgrid(temps, pres);
   dV_x=np.zeros(np.shape(pres_x))
   for i in np.arange(0,np.size(pres_x[:,0]),1):
       for j in np.arange(0,np.size(pres_x[0,:]),1):
           res = minimize_scalar(dV,bounds=(A, C), args=(K0,KT,pres_x[i,j]), method='brent')
           dV_x[i,j]=res.x
   alphaP0=dV_x*np.exp((grun+1.)*((dV_x**(-1.))-1.))
   rhoP0=dens_ref*dV_x
   intalphaT=(a0*(temps_x-273.))+((a1/2.)*((temps_x**2.)-(273.**2.)))
   density=rhoP0*(1.-(alphaP0*intalphaT))
   return density

def dens_grav_pres_prof(rad):
   """ 
      This subroutine generates a profile of
      the Earth's density, gravity, and pressure
      at the given radii
         inputs:
            rad: input radii
   """
   from scipy.interpolate import interp1d   
   # Universal gravitational constant
   grav_const = 6.67408e-11;
   # load the Earth's density from PREM
   rshl_prem, dens_prem = np.loadtxt('/Users/fredrichards/space1/LLSVP_Density/DATA/PREM_1s.csv',delimiter=',',\
           usecols=(0,2), unpack=True)
   # Converting to SI units
   rshl_prem = rshl_prem*1.e3; dens_prem = dens_prem*1.e3;
   # Generate arrays for mass, gravity, and pressure
   mass_prem = np.zeros(np.shape(rshl_prem));
   grav_prem = np.zeros(np.shape(rshl_prem));
   pres_prem = np.zeros(np.shape(rshl_prem));
   # Compute the inner Mass of the Earth at each radius
   for i in range(len(rshl_prem)-2,0,-1):
      mass_prem[i] = mass_prem[i+1] +\
               (4./3.)*np.pi*dens_prem[i]*(((rshl_prem[i]+rshl_prem[i-1])/2.)**3.\
               -((rshl_prem[i]+rshl_prem[i+1])/2.)**3.)
   mass_prem[0] = mass_prem[1] +\
               (4./3.)*np.pi*dens_prem[i]*((rshl_prem[i])**3.-\
               ((rshl_prem[i]+rshl_prem[i+1])/2.)**3.);
   # Compute the Earth's gravity at each radius
   grav_prem[rshl_prem!=0] = mass_prem[rshl_prem!=0]*grav_const/rshl_prem[rshl_prem!=0]**2.;
   # Compute the Earths pressure profile at each radius
   for i in range(1,len(rshl_prem)-1):
      if rshl_prem[i]==0.0:
         del_pres = 0.0
      else:
         del_pres = 1./3.*grav_prem[i]*dens_prem[i]*\
               ((rshl_prem[i-1]**3.-rshl_prem[i]**3.)/(rshl_prem[i]**2.));
      pres_prem[i] = pres_prem[i-1] + del_pres; 
   # interpolate density gravity and pressure onto the input radii
   dens = interp1d(rshl_prem, dens_prem)(rad);
   grav = interp1d(rshl_prem, grav_prem)(rad);
   pres = interp1d(rshl_prem, pres_prem)(rad);

   return dens, grav, pres