import h5py
import numpy as np
import sys
from numpy import loadtxt
from os.path import join
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

def load_grd_file(fname):
   """
      load_grd_file loads the content of a grd file,
         knowing that we have lon lat and z
   """
   file_id = h5py.File(fname, 'r')
   data = np.array(file_id.get('z'))
   lat  = np.array(file_id.get('lat'))
   lon  = np.array(file_id.get('lon'))
   file_id.close()

   return lat, lon, data

def plot_grd_file(name_file, plot_title=None, clb_title=None):
   """
      plot a given grid file
         name_file: name of the file
         plot_title: if given, plot_title will appear as the title
         clb_title: if given, appears as the title of the colorbar
   """
   import matplotlib.pyplot as plt
   from mpl_toolkits.basemap import Basemap

   # open the grd file and read lat, lon and data:z
   file_set = h5py.File(name_file,'r')
   data     = file_set.get('z')
   np_vel   = np.array(data)
   data     = file_set.get('lat')
   np_lat   = np.array(data)
   data     = file_set.get('lon')
   np_lon   = np.array(data)
   file_set.close()


   h = plt.figure(figsize=(16,9))
   grd_map = Basemap(projection='moll',lon_0=0.0,lat_0=0.0);
   np_lon,np_lat = np.meshgrid(np_lon,np_lat);
   x,y = grd_map(np_lon,np_lat)
   palette = plt.cm.seismic_r
   palette.set_bad('k')
   #abs_val = np.max(abs(np_vel))
   cntr_map = grd_map.contourf(x,y,np_vel,64,cmap=palette,extend="both")
   grd_map.drawcoastlines()
   cb = grd_map.colorbar(cntr_map, "bottom")
   if clb_title==None:
      cb.set_label('[%]')
   else:
      cb.set_label(clb_title)
   if plot_title !=None:
      plt.title(plot_title)
   return h;

class mineralogy_table:
   """
      Class to convert seismic data with specified depths to
      temperatures and/or density

   """
   models_avail = ['SLB_16'];
   compositions_avail = ['pyrolite', 'basalt', 'solar', 'fepyrolite', 'harzburgite', 'bridgmanite'];
   data_avail = ['v_p', 'v_s', 'shear_mod', 'bulk_mod', 'rho']
   def __init__(self, model=None, composition=None, depths=None,\
                     q_model=None, k_rad=2,\
                     E=None, alpha=None, Q_cammarano=None, Q_lu=None, Q_matas=None,\
                     linearisation=False, geo_temp_var=None, geotherm=None):
      """
      This class initilizes the following:
         mineralogy_table.intrp_miner_data: interpolated minerlogy tables

      Input parameters:
         model          which model to use? e.g., SLB_16
         composition    which composition to look up? pyrolite or basalt
         depths         the depths at which the mineralogy tables are computed
         k_rad          how many neighboring layers should be interpolated to one?
                        default k_rad=2
      """
      # h5py to read in the mineralogy tables
      import h5py
      # for 1D interpolation of radial profiles
      from scipy.interpolate import interp1d
      # Load the anelasticity_model for Cammarano in case we need to use it
      from lib_anelasticity import anelasticity_model_cammarano_type_1,anelasticity_model_lu_type

      # Gas constant
      self.R = 8.3144598;

      # Seismic frequency
      self.omega=1.;

      # Here we define the parameters for the existing tables
      # modify this section if you have any new data
      self.path = '/Users/fredrichards/space1/LLSVP_Density/DATA/'
      # depths should be in kms
      if len(str("%1.0f" %(depths[0]))) > 4:
         print("Warning: The input depths are in meters! Let's convert it to kilometers");
         depths = depths/1e3;
      if True in depths>=2890.:
         raise('depth values are not suitable for mantle!')

      self.dpths = depths;
      self.k_rad = k_rad;

      # check if the input parameters are consistent with our available data
      if not model in self.models_avail:
         raise  ValueError('Models available:'+\
               str(self.models_avail)+'\n Please give in a valid model name')
      self.model = model;
      # if anelastic correction is required, ask for E and \alpha
      if q_model=='Cammarano' and (Q_cammarano==None or (type(Q_cammarano) != anelasticity_model_cammarano_type_1)):
         raise ValueError('For a Cammarano Anelastic correction Q_cammarano should be given by function anelasticity_model_cammarano_type_1')
      elif (q_model!=None and q_model!='Lu' and q_model!='Matas' and q_model!='Cammarano' and (E==None or alpha==None)):
         raise ValueError('For an anelastic correction based on a seismic profile E and alpha are necessary.'+\
                        'One of them is empty.');
      elif (q_model!=None and q_model!='Cammarano' and (E!=None and alpha!=None)):
          self.alpha=alpha
          self.E=E
      # We store q_model and Q_cammarano
      #  so that we know which parameters we have chosen for this
      self.q_model = q_model;
      self.Q_cammarano=Q_cammarano;
      self.Q_lu=Q_lu;
      self.Q_matas=Q_matas;

      # See if it is a mixed model or not
      # If there we have a mixed model, then we make self.composition and self.vol_fracs array
      #  if this is a pure model, these lists are having a length of one
      if '!' in composition:
         if np.mod(len(composition.split('!')),2) != 0:
            raise ValueError('Error with "Composition" format. Example: Pyrolite!20!Basalt!80!')
         self.composition = [None]*int(len(composition.split('!'))/2);
         self.vol_fracs   = np.zeros(int(len(composition.split('!'))/2));
         print('Note: You have chosen a mixed model with:')
         for i in range(int(len(composition.split('!'))/2)):
            print(str('%15s |  %s %%')\
               %(composition.split('!')[2*i], composition.split('!')[2*i+1]))
            if not composition.split('!')[2*i] in self.compositions_avail:
               raise  ValueError(\
                  'Composition '+composition.split('!')[2*i]+' is not available!')
            self.composition[i] = composition.split('!')[2*i];
            self.vol_fracs[i] = composition.split('!')[2*i+1];
      else:
         # Next in case of having one composition
         #  here we still make a list of compositions with 1 member
         if not composition in self.compositions_avail:
            raise  ValueError(\
            'Please give in a valid compositions.\nCompositions available:'\
            +str(self.compositions_avail))
         self.composition = [composition];
         self.vol_fracs   =np.array([100]);

      # Check if the sum of vol_fracs is correct
      if np.sum(self.vol_fracs)!= 100:
         raise ValueError(\
         'Some of the volumetric fractions should be 100!')

      # Do we have depths?
      if len(depths)==0:
         raise  ValueError('"depths" is empty!')

      # If we have linearisation, we should also have the geotherm
      if linearisation and (len(geotherm)==0 or geo_temp_var==None):
         raise ValueError('Either of the "geotherm" or "geo_temp_var" is not given')

      if linearisation:
         self.geotherm = geotherm;
         self.geo_temp_var = geo_temp_var;
         self.rad_grad = [None]*len(self.dpths);
         self.bvals = [None]*len(self.dpths);

      # Up until here it should be only reading in parameters
      # open the hdf5 file and read in the required dataset
      self.load_data(linearisation_flag=linearisation);

   #--------------------------------------------------------------------------------
   #  Apply the anelasticity correction
   #--------------------------------------------------------------------------------
   def apply_anelasticity(self, vel_data_set, my_data_type):
      """
         apply_anelasticity:
               loads a specific model of anelasticity for mantle material
               and corrects the given velocity values
         my_data_type: could be either v_s or v_p;
                  Right now we only restrict it to v_s

      """
      # At the moment we only correct shear wave velocities
      if my_data_type!='v_s':
         return None
      # load a Cammarano type of annelasticity model, so we can determine if the input
      #  function for Cammarano anelasticity is correct or not
      from lib_anelasticity import anelasticity_model_cammarano_type_1,anelasticity_model_lu_type
      from scipy.interpolate import interp1d

      # Contruct the Q and Alpha matrices
      if self.q_model =='Cammarano' and my_data_type=='v_s':
         if type(self.Q_cammarano) != anelasticity_model_cammarano_type_1:
            raise('Q_cammarano should be of type: lib_anelasticity.anelasticity_model_type_1!');
         # Create a Cammarano Q matrix
         # Q cammarano is already initiated outside this routine
         #     the ready-to-be-used can be called here for alpha and q matrices
         mat_q_3d,solid_temp = self.Q_cammarano.create_Q_2d(self.dpths*1e3, self.temp);
         alpha = self.Q_cammarano.a_matrix;
      elif self.q_model =='Lu' and my_data_type=='v_s':
         # Create a Lu Q matrix
         # Q lu is already initiated outside this routine
         #     the ready-to-be-used can be called here for alpha and q matrices
         mat_q_3d = self.Q_lu.create_Q_2d(self.dpths*1e3, self.temp);
         alpha = self.Q_lu.a_matrix;

      elif self.q_model =='Matas' and my_data_type=='v_s':
         # Create a Lu Q matrix
         # Q lu is already initiated outside this routine
         #     the ready-to-be-used can be called here for alpha and q matrices
         mat_q_3d = self.Q_matas.create_Q_2d(self.dpths*1e3, self.temp);
         alpha = self.Q_matas.a_matrix;

      elif my_data_type == 'v_s':
         # load either q_mu or q_kappa from datasets
         #     right now: I am not sure if q_kappa values can be used
         #     for p wave anelastic correction
         if my_data_type=='v_s':
            q_rad, _ = self.load_q_model(q_target=self.q_model, dpths=self.dpths);
         elif my_data_type=='v_p':
             _, q_rad = self.load_q_model(q_target=self.q_model, dpths=self.dpths);
         # q_rad[self.dpths>670.]=q_rad[self.dpths>670.]-212.
         # load the radial adiabat
         #  Note: Ideally, one would use a terra geotherm for these conversions
         # SLB_dpths, SLB_temp_ref = \
         #          np.loadtxt('/Users/fredrichards/space1/LLSVP_Density/DATA/Q_MODELS/1600K_adi_SLB2011.txt', skiprows=1, unpack=True);
         SLB_dpths, SLB_temp_ref = \
                  np.loadtxt('/Users/fredrichards/space1/LLSVP_Density/DATA/TerraMT512vs.dat', usecols=(0,1), unpack=True);
         # interpolate the SLB 2011 adiabat on the models depths
         temp_adiabat = interp1d(SLB_dpths, SLB_temp_ref,\
                  kind='slinear', fill_value='extrapolate')(self.dpths);

         # make a matrix out of temperatures, adiabat and q model, to perform the operations
         #  via matrices
         R=self.R
         mat_temps, mat_temp_adiabat = np.meshgrid(self.temp, temp_adiabat)
         _, mat_q_rad = np.meshgrid(self.temp, q_rad)
         mat_q_3d = mat_q_rad*\
               (np.exp((self.alpha*self.E/R)*(1./mat_temps - 1./mat_temp_adiabat)))
         alpha=self.alpha


      # # Store unchanged velocities
      vel_data_set_old=vel_data_set

      # The following is the physical equations we use here
      # v_{s,p}_corrected = v_{s,p}*(1-(F/pi*alpha)*(1/Q_cur)
      vel_data_set = \
            vel_data_set*(1.- \
            self.calculate_F(alpha)/\
            (np.pi*alpha)*(1./mat_q_3d));

      # Check with plot
      # X,Y=np.meshgrid(self.temp,self.dpths)
      #
      # plt.pcolormesh(X,np.flipud(Y),np.flipud(((vel_data_set-vel_data_set_old)/vel_data_set_old)*100.),cmap=plt.cm.Reds, vmin=-3,vmax=0)
      # plt.colorbar()
      # plt.plot(self.geotherm[:],self.dpths[:],'-k')
      # plt.xlabel('Temperature (K)')
      # plt.ylabel('Depth (km)')
      # plt.gca().invert_yaxis()
      # plt.savefig(repr(self.q_model)+'.png')

      return vel_data_set
   #--------------------------------------------------------------------------------
   #
   #--------------------------------------------------------------------------------
   # loads and computes averages based on given self.composition, and self.vol_fracs
   #  out_data_type, is what we are in the end interested in
   #  depending on out_data_type, different dataset would be read in
   #     and the way to compute averages may be different
   #--------------------------------------------------------------------------------
   def load_data(self, linearisation_flag):
      # for cKDTree routines
      import scipy.spatial as scisp
      for i in range(len(self.composition)):
         data_fi = h5py.File(join(self.path,self.model+'_'+self.composition[i]+'.hdf5'),'r')
         if i==0:
            # read in temperature and depths from the first file
            # depths and temperatures are supposed to be identical
            self.temp  = np.array(data_fi.get('Temperatures'));
            miner_dpth = np.array(data_fi.get('Depths'));
            # We are going to make principle arrays of bulk
            bulk_voigt_ave  = (self.vol_fracs[i]/100.)*np.array(data_fi.get('bulk_mod'));
            bulk_reuss_ave  = np.divide(self.vol_fracs[i]/100.,np.array(data_fi.get( 'bulk_mod')));
            shear_voigt_ave = (self.vol_fracs[i]/100.)*np.array(data_fi.get('shear_mod'));
            shear_reuss_ave = np.divide(self.vol_fracs[i]/100.,np.array(data_fi.get('shear_mod')));
            dens_ave        = (self.vol_fracs[i]/100.)*np.array(data_fi.get('rho'));
         else:
            # From here on we had the contributions from other components to the
            #     already defined arrays
            bulk_voigt_ave  = bulk_voigt_ave  + (self.vol_fracs[i]/100.)*np.array(data_fi.get('bulk_mod'));
            shear_voigt_ave = shear_voigt_ave + (self.vol_fracs[i]/100.)*np.array(data_fi.get('shear_mod'));
            bulk_reuss_ave  = bulk_reuss_ave  + np.divide(self.vol_fracs[i]/100.,np.array(data_fi.get( 'bulk_mod')));
            shear_reuss_ave = shear_reuss_ave + np.divide(self.vol_fracs[i]/100.,np.array(data_fi.get('shear_mod')));
            dens_ave        = dens_ave        + (self.vol_fracs[i]/100.)*np.array(data_fi.get('rho'));
         data_fi.close();

      bulk_mod    = bulk_voigt_ave*0.5    + np.divide(1.,bulk_reuss_ave)*0.5;
      shear_mod   = shear_voigt_ave*0.5   + np.divide(1.,shear_reuss_ave)*0.5;
      # to interpolate the shear, bulk moduli and densities onto the input dpths
      #     we make a tree of the available dpths in the model, and find the
      #     two neighbor points in miner_dpths for each depth of the input dpths
      tree_miner_dpth = scisp.cKDTree(\
                  miner_dpth.reshape(len(miner_dpth),1), leafsize=100)

      # find out the neighbouring layers for each seismic depths
      tree_res_miner_dpth = tree_miner_dpth.query(\
                  self.dpths.reshape(len(self.dpths),1), k=self.k_rad)

      # make sure that no distance is zero, so we do not end up with 1/0
      tree_res_miner_dpth[0][:,0 ][(tree_res_miner_dpth[0]<=1e-4)[:,0]]=1.0
      tree_res_miner_dpth[0][:,1:][(tree_res_miner_dpth[0]<=1e-4)[:,0]]=1e10

      print ("Removing NaNs")
      for z in range(np.shape(bulk_mod)[0]):
          for x in range(np.shape(bulk_mod)[1]):
              if np.isnan(bulk_mod[z,0])==True:
                  bulk_mod[z,0]=bulk_mod[z,1]
                  shear_mod[z,0]=shear_mod[z,1]
                  dens_ave[z,0]=dens_ave[z,1]
              elif np.isnan(bulk_mod[z,x])==True:
                  if np.isnan(bulk_mod[z,x-1])==False and np.isnan(bulk_mod[z,x+1])==False:
                      bulk_mod[z,x]=(bulk_mod[z,x-1]+bulk_mod[z,x+1])/2
                      shear_mod[z,x]=(shear_mod[z,x-1]+shear_mod[z,x+1])/2
                      dens_ave[z,x]=(dens_ave[z,x-1]+dens_ave[z,x+1])/2
                  else:
                      bulk_mod[z,x]=(bulk_mod[z,x-2]+bulk_mod[z,x+2])/2
                      shear_mod[z,x]=(shear_mod[z,x-2]+shear_mod[z,x+2])/2
                      dens_ave[z,x]=(dens_ave[z,x-2]+dens_ave[z,x+2])/2
              elif np.isnan(bulk_mod[z,-1])==True:
                  bulk_mod[z,-1]=bulk_mod[z,-2]
                  shear_mod[z,-1]=shear_mod[z,-2]
                  dens_ave[z,-1]=dens_ave[z,-2]

      # a weighted average to map the mineralogy model to the input depths
      self.bulk_mod = \
            np.dot(np.diag(1./np.sum(1./tree_res_miner_dpth[0],axis=1)),\
                  np.einsum('km,kmj->kj',1./tree_res_miner_dpth[0], \
                    bulk_mod[tree_res_miner_dpth[1]]))
      self.shear_mod = \
            np.dot(np.diag(1./np.sum(1./tree_res_miner_dpth[0],axis=1)),\
                  np.einsum('km,kmj->kj',1./tree_res_miner_dpth[0], \
                    shear_mod[tree_res_miner_dpth[1]]));
      self.dens = \
            np.dot(np.diag(1./np.sum(1./tree_res_miner_dpth[0],axis=1)),\
                  np.einsum('km,kmj->kj',1./tree_res_miner_dpth[0], \
                    dens_ave[tree_res_miner_dpth[1]]));

      # if we are to linearise the table
      if linearisation_flag:
         print("Linearising the input models");
         self.bulk_mod  = self.linearise_table(self.bulk_mod );
         self.shear_mod = self.linearise_table(self.shear_mod);
         self.dens      = self.linearise_table(self.dens     );

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   #
   #--------------------------------------------------------------------------------
   def shear_wave(self):
      """
         Calculates shear wave velocity based on loaded bulk, shear, and dens
      """
      # This routine generates shear wave-velocities out of the loaded densy and shear modulus
      results =  np.sqrt(np.divide(self.shear_mod,self.dens));
      # In case we need to apply an anelasticity correction
      if self.q_model !=None:
         results = self.apply_anelasticity(results, my_data_type='v_s');
      return results

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   #
   #--------------------------------------------------------------------------------
   def primary_wave(self):
      """
         This routine generates primary wave-velocities out of the loaded
                  densy and shead+bulk moduli
      """
      results = np.sqrt(np.divide(self.bulk_mod+(4./3.)*self.shear_mod,self.dens));
      return results

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   #
   #--------------------------------------------------------------------------------
   # calculate the F factor in anelastic correction
   def calculate_F(self, alpha):
      """

      """
      import numpy as np
      F = ((np.pi*alpha)/2.)*(1./np.tan(np.pi*alpha/2.))
      return F

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   # linearises input array "non_lin_table", previously loaded,
   #     around a thermal profile given as the self.geotherm, by \pm self.geo_temp_var
   #--------------------------------------------------------------------------------
   def linearise_table(self, non_lin_table):
      bvals    = np.zeros(np.shape(non_lin_table)[0]);
      rad_grad = np.zeros(np.shape(non_lin_table)[0]);
      lin_table = np.zeros(np.shape(non_lin_table));
      for ir in range(np.shape(non_lin_table)[0]):
         ir_mean = abs(self.temp-(self.geotherm[ir])).argmin();
         ir_min  = abs(self.temp-(self.geotherm[ir]-self.geo_temp_var)).argmin();
         ir_max  = abs(self.temp-(self.geotherm[ir]+self.geo_temp_var)).argmin();

         if ir_min==0 and (ir_mean == ir_min):
            ir_min=0; ir_mean=1;
         rad_grad[ir] = np.nanmedian(np.gradient(\
            non_lin_table[ir][ir_min:ir_max], self.temp[ir_min:ir_max]\
                                                   ));
         bvals[ir] =  non_lin_table[ir][ir_mean]-rad_grad[ir]*self.temp[ir_mean];
         lin_table[ir] = rad_grad[ir]*self.temp + bvals[ir]

      return lin_table

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   # load a specific q model
   #--------------------------------------------------------------------------------
   def load_q_model(self, q_target=None, dpths=None,  q_mu_lb=20.0, q_mu_ub=4900.0,\
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
      q_hdf5_name = '/Users/fredrichards/space1/LLSVP_Density/DATA/Q_collection.hdf5'
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
      interp_func = interp1d(mod_depth, mod_Q_mu, kind='slinear', fill_value='extrapolate');
      q_mu = interp_func(dpths);
      interp_func = interp1d(mod_depth, mod_Q_p , kind='slinear', fill_value='extrapolate');
      q_p  = interp_func(dpths);

      return q_mu, q_p

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   #
   #--------------------------------------------------------------------------------
   def convert2temperature(self, data_type, seis_model, k_lateral=4):
      """
      This class initilizes the following:
         data_set interpolated minerlogy tables

      Input parameters:
         model          which model to use? e.g., SLB_16
         composition    which composition to look up? pyrolite or basalt
         data_type      which dataset are we converting? v_s, v_p, shear_mod?
         depths         the depths at which the mineralogy tables are computed
         k_rad          how many neighboring layers should be interpolated to one?
                        default k_rad=1
      """
      # for cKDTree routines
      import scipy.spatial as scisp
      temp_fld = [None]*len(self.dpths);

      # Check if the dataset is given and is available
      if not data_type in self.data_avail:
         raise  ValueError('data_type should be one of',str(self.data_avail))

      # based on the input data we decide which
      #     data set to use as intrp_miner_data
      if    data_type == 'v_p':
         mydataset = self.primary_wave();
      elif  data_type == 'v_s':
         mydataset = self.shear_wave();
      elif  data_type == 'shear_mod':
         mydataset = self.shear_mod
      elif  data_type == 'bulk_mod':
         mydataset = self.bulk_mod
      else:
         mydataset = self.dens


      # go through a do loop an convert seismic values to temperature
      for i in range(len(self.dpths)):

         # Get mean vs
         meanobsvs=np.mean(seis_model[i])

         # Get interpolated Vs for geotherm
         tree_temperature = scisp.cKDTree(self.temp.reshape(len(self.temp),1),leafsize=100);

         tree_temperature_res = tree_temperature.query([self.geotherm[i]], k=2)
         res_conv = \
            np.sum((1./tree_temperature_res[0])*mydataset[i][tree_temperature_res[1]])/\
                        np.sum(1./tree_temperature_res[0])

         # Vs correction to get mean Vs equal to geotherm
         seis_model[i]=seis_model[i]-(meanobsvs-res_conv)

         # print('%s to temperature layer at %5.0i km' %(data_type,self.dpths[i]));
         # build the coordinates, by creating longitudes and latitudes
         lon = np.linspace(-180,180,int(np.shape(seis_model[i])[1]))
         lat = np.linspace(-90,90,int(np.shape(seis_model[i])[0]))
         # mesh out of longitude and latitude so we can have the
         grd_lon, grd_lat = np.meshgrid(lon, lat);
         # here we make a column stack of longitudes and latitudes as coordinates
         # coords = np.column_stack(( grd_lon.reshape(np.prod(np.shape(grd_lon))),\
         #                            grd_lat.reshape(np.prod(np.shape(grd_lat)))));
         # make a tree out of the mineral table at that specific depth
         # mydataset[i] = mineral table values at that depth
         tree_miner_at_depths = scisp.cKDTree(\
                  mydataset[i].reshape(len(mydataset[i]),1),leafsize=100);
         # find the two neighboring points in the mineral table for each entry in the seismic model
         tree_res = tree_miner_at_depths.query(\
             seis_model[i].reshape(np.prod(np.shape(seis_model[i])),1), k=k_lateral)
         # We find the closes point, and then label the next index in the table as the
         #  second closest point in the table
         neighbor1 = tree_res[1][:,0];
         neighbor2 = tree_res[1][:,0]+1;

         # we calculate the distance of all the point
         dist1 = seis_model[i].reshape(np.prod(np.shape(seis_model[i])))\
                  - mydataset[i][neighbor1]

         # based on the distances from neighbor 1, we my change the second closest neighbors
         neighbor2[dist1<0.0] += -2;
         # adjust the members equaling the minimum or maximum indices
         neighbor2[neighbor2>=len(mydataset[i])] = \
               len(mydataset[i])-1;
         neighbor2[neighbor2<0] = 0
         # what is the distance from the second closes member
         dist2 = seis_model[i].reshape(np.prod(np.shape(seis_model[i])))\
                  - mydataset[i][neighbor2]
                  #- np.take(mydataset[i],neighbor2,mode='wrap')
         temp_fld[i] = (1./dist1*self.temp[neighbor1]+\
                       1./dist2*self.temp[neighbor2])/\
                        ((1./dist1)+(1./dist2))

         ## Here is the original method that was a bit more complicated
         #tree_res[0][:,0 ][(tree_res[0]<=1e-4)[:,0]] = 1.0  ;
         #tree_res[0][:,1:][(tree_res[0]<=1e-4)[:,0]] = 1e10 ;
         ## interpolate between the two neighboring points
         #temp_fld[i] = \
         #   np.sum((1./tree_res[0][:,0:2])*self.temp[tree_res[1][:,0:2]],axis=1)/\
         #            np.sum(1./tree_res[0][:,0:2],axis=1)
         ## if the two neighboring points are not sequential data from the table, it means we have multiple
         ##  candidates, so better to find the value from the actual neighboring points in the grid
         #temp_fld[i][np.var(tree_res[1],axis=1)>np.var(range(k_lateral))] = float('NaN')
         ## make a tree out of the coordinates of the nodes with a value
         #tree_coords = scisp.cKDTree(coords[temp_fld[i]==temp_fld[i]], leafsize=100);
         ## find the neighboring points of the nodes with "nan"
         #tree_coords_res = tree_coords.query(coords[temp_fld[i]!=temp_fld[i]],k=k_lateral)
         #tree_coords_res[0][:,0 ][(tree_coords_res[0]<=1e-8)[:,0]] = 1.0  ;
         #tree_coords_res[0][:,1:][(tree_coords_res[0]<=1e-8)[:,0]] = 1e10 ;
         ## set the values by averaging the four nearest points
         #temp_fld[i][temp_fld[i]!=temp_fld[i]]=\
         #      np.sum((1./tree_coords_res[0])*temp_fld[i][temp_fld[i]==temp_fld[i]][tree_coords_res[1]],axis=1)/\
         #         np.sum(1./tree_coords_res[0],axis=1)
         # reshape data based on the input data
         temp_fld[i] = temp_fld[i]+(self.geotherm[i]-np.mean(temp_fld[i]))
         # print np.mean(temp_fld[i]), self.geotherm[i]
         temp_fld[i] = temp_fld[i].reshape(np.shape(seis_model[i]))
      return temp_fld

   #--------------------------------------------------------------------------------
   #--------------------------------------------------------------------------------
   #
   #--------------------------------------------------------------------------------
   def temperature2convert(self, temp_fld, data_type=None):
      """
      This class initilizes the following:
         mineralogy_table.intrp_miner_data: interpolated minerlogy tables

      Input parameters:
         model          which model to use? e.g., SLB_16
         composition    which composition to look up? pyrolite or basalt
         data_type      which dataset are we converting? v_s, v_p, shear_mod?
         depths         the depths at which the mineralogy tables are computed
         k_rad          how many neighboring layers should be interpolated to one?
                        default k_rad=1
      """
      # for cKDTree routines
      import scipy.spatial as scisp

      # check if temperature is having the right dimensions
      if np.shape(temp_fld)[0] != len(self.dpths):
         raise ValueError(str('input temperature should have %i rows,'\
                  +' as in the defined dpths array' %len(self.dpths)));

      # check if the data_type is valid
      if not data_type in self.data_avail:
         raise  ValueError('data_type should be one of',str(self.data_avail));

      # based on the input data we decide which data set to use as intrp_miner_data
      if    data_type == 'v_p':
         mydataset = self.primary_wave();
      elif  data_type == 'v_s':
         mydataset = self.shear_wave();
      elif  data_type == 'shear_mod':
         mydataset = self.shear_mod
      elif  data_type == 'bulk_mod':
         mydataset = self.bulk_mod
      else:
         mydataset = self.dens

      # initiating the array for the results
      res_conv = [None]*np.shape(temp_fld)[0];

      # make a tree out of the available temperatures in the mineralogy model
      tree_temperature = scisp.cKDTree(self.temp.reshape(len(self.temp),1),leafsize=100);

      # go through every layer and convert it to the corresponding data_set
      for i in range(np.shape(temp_fld)[0]):
         # print('T to %s layer at %5.0i km' %(data_type,self.dpths[i]));
         tree_temperature_res = tree_temperature.query(\
                     temp_fld[i].reshape(np.prod(np.shape(temp_fld[i])),1), k=2)
         # make sure that no distance is zero, so we do not end up with 1/0
         tree_temperature_res[0][:,0 ][(tree_temperature_res[0]<=1e-4)[:,0]]=1.0
         tree_temperature_res[0][:,1:][(tree_temperature_res[0]<=1e-4)[:,0]]=1e10
         res_conv[i] = \
            np.sum((1./tree_temperature_res[0])*mydataset[i][tree_temperature_res[1]],axis=1)/\
                        np.sum(1./tree_temperature_res[0],axis=1)
         res_conv[i] = res_conv[i].reshape(np.shape(temp_fld[i]));
         if data_type=='rho':
            res_conv[i] = res_conv[i]*1e3;
      return res_conv

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
   rshl_prem, dens_prem = np.loadtxt('/Users/fredrichards/space1/LLSVP_Density/DATA/PREM_1s.csv',\
            delimiter=',',usecols=(0,2), unpack=True)
   # Converting to SI units
   rshl_prem = rshl_prem*1e3; dens_prem = dens_prem*1e3;
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
   grav_prem = mass_prem*grav_const/rshl_prem**2.;
   # Compute the Earths pressure profile at each radius
   for i in range(1,len(rshl_prem)):
      pres_prem[i] = pres_prem[i-1] + grav_prem[i]*dens_prem[i]*(rshl_prem[i-1]-rshl_prem[i]);
   # interpolate density gravity and pressure onto the input radii
   dens = interp1d(rshl_prem, dens_prem)(rad);
   grav = interp1d(rshl_prem, grav_prem)(rad);
   pres = interp1d(rshl_prem, pres_prem)(rad);

   return dens, grav, pres



class lay_seis_model(object):
   def __init__(self, model_path):
      self.path  = model_path;
      self.depth = float(model_path[model_path.rfind('/')+1:model_path.rfind('km.grd')]);
      self.lon, self.lat, self.data = self.import_grd();
      self.temp = None;
      self.rho = None;
      self.sph = None;

   def import_grd(self):
      file_set = h5py.File(self.path,'r')
      data     = np.array(file_set.get('z'))
      lat      = np.array(file_set.get('lat'))
      lon      = np.array(file_set.get('lon'))
      file_set.close()
      return lon, lat, data

   def export_grd(self, output_name):
      h5py_file = h5py.File(name=output_name,mode='w')
      h5py_file.create_dataset(name='z'  , data=self.data)
      h5py_file.create_dataset(name='lat', data=self.lat)
      h5py_file.create_dataset(name='lon', data=self.lon)
      h5py_file.close()

   def export_grd_rho(self, output_name):
      h5py_file = h5py.File(name=output_name,mode='w')
      h5py_file.create_dataset(name='z'  , data=self.rho)
      h5py_file.create_dataset(name='lat', data=self.lat)
      h5py_file.create_dataset(name='lon', data=self.lon)
      h5py_file.close()

   def export_grd_temp(self, output_name):
      h5py_file = h5py.File(name=output_name,mode='w')
      h5py_file.create_dataset(name='z'  , data=self.temp)
      h5py_file.create_dataset(name='lat', data=self.lat)
      h5py_file.create_dataset(name='lon', data=self.lon)
      h5py_file.close()

   def export_grd_temp_ascii(self, output_name):
      ascii_id = open(file=output_name, mode='w')
      lon_x, lat_x = np.meshgrid(self.lon, self.lat);
      lon_x = lon_x.reshape(np.prod(np.shape(lon_x)))
      lat_x = lat_x.reshape(np.prod(np.shape(lat_x)))
      my_temp = self.temp.reshape(np.prod(np.shape(lat_x)))
      for i in range(len(lon_x)):
         ascii_id.write(str('%8.2f %8.2f %12.8e\n' %(lon_x[i], lat_x[i], my_temp[i])))
      ascii_id.close()
