/* -------------------------------------------------
    3PHMLT
   -------------------------------------------------  */

/* 
 SOLVE EQUATIONS PERTAINING TO PRODUCTIVITY AND PRESSURE TEMPERATURE
 FOR MELTING OF A THREE LITHOLOGY MANTLE 

 PAREMTERSIATIONS OF KATZ 2003 FOR PERIDOTITE MELTING AND LAMBART
 FOR KG MELTING WITH PRODUCTIVITY AFTER CPX OUT SET BY KATZ HARZBURGITE MELTING EQUATION
 ALSO INCLUDE A HARZBURGITE LITHOLOGY WHICH DOES NOT MELT
*/
/* COMPILE WITH: g++ -o 3phmlt.o 3phmlt_1.cc */
/* -------------------------------------------------
    INCLUDES 
   -------------------------------------------------  */

// std
#include <cstdio>
#include <string>
#include <iostream>
#include <getopt.h>

// nrh
#include "nr3.h"
// this overloads rk4 between the version I modified and the original
#include "rk4_mod.h"
#include "rk4.h"
#include "roots.h"


/* -------------------------------------------------
    STRUCTURES 
   -------------------------------------------------  */

struct lith {
  // initialise variables
  Doub A[3];	 // solidus parameterisation, degC:degC/GPa:degC/GPa^2
  Doub B[3];     // lhrz liquidus parameterisation, degC:degC/GPa:degC/GPa^2
  Doub C[3];     // hrz liquidus parameterisation, degC:degC/GPa:degC/GPa^2
  Doub r[2];	 // cpx out parameterisation, cpx/melt:cpx/melt/GPa
  Doub M;	 // weight fraction of cpx in the unmelted rock
  Doub bt[2];	 // exponent for F(T) parameterisation
  Doub Cp;	 // heat capacity, J/(KgK)
  Doub alph[2];  // thermal expansivity, solid-liquid, /K
  Doub rho[2];   // density, solid-liquid, Kg/m^3
  Doub DS;

};

/* -------------------------------------------------
    FUNCTION AND FUNCTORS 
   -------------------------------------------------  */

// PERIDOTITE
// parameterisations from Katz 2003
// T_liquidus_lherzolite
Doub Tlhrz (const Doub& P, const lith& lith) {
    return lith.B[0] + lith.B[1]*P + lith.B[2]*pow(P,2);
}
// T_liquidus
Doub Tl_pd (const Doub& P, const lith& lith) {
  return lith.C[0] + lith.C[1]*P + lith.C[2]*pow(P,2);
}
// T_solidus
Doub Ts_pd (const Doub& P, const lith& lith) {
  return lith.A[0] + lith.A[1]*P + lith.A[2]*pow(P,2);
}
// D(T_liquidus_lherzolite)/DP
Doub dtlhrz_dp (const Doub& P, const lith& lith) {
  return lith.B[1]+2*lith.B[2]*P;
}
// D(T_liquidus)/DP
Doub dtl_dp (const Doub& P, const lith& lith) {
  return lith.C[1]+2*lith.C[2]*P;
}
// D(T_solidus)/DP
Doub dts_dp (const Doub& P, const lith& lith)  {
  return lith.A[1]+2*lith.A[2]*P;
}

// HARZ - these equations are not used presently
//        need modification and inclusion when the depleted lithology is real.
// T_liquidus
Doub Tl_hrz (const Doub& P, const lith& lith) {
  return lith.C[0] + lith.C[1]*P + lith.C[2]*pow(P,2);
}
// T_solidus
Doub Ts_hrz (const Doub& P, const lith& lith) {
  return lith.A[0] + lith.A[1]*P + lith.A[2]*pow(P,2);
}
// D(T_liquidus_harzburgite)/DP
Doub dthrz_dp (const Doub& P, const lith& lith) {
  return lith.B[1]+2*lith.B[2]*P;
}
// D(T_liquidus)/DP
Doub dtl_dp_hrz (const Doub& P, const lith& lith) {
  return lith.C[1]+2*lith.C[2]*P;
}
// D(T_solidus)/DP
Doub dts_dp_hrz (const Doub& P, const lith& lith)  {
  return lith.A[1]+2*lith.A[2]*P;
}


// PYROXENITE : KG1
// parameterisations from Lambart
// T_solidus
Doub Ts_px (const Doub& P, const lith& lith) {
  return lith.A[0] + lith.A[1]*P + lith.A[2]*pow(P,2);
}
// T_cpx_out
Doub Tco_px (const Doub& P, const lith& lith) {
  return lith.B[0] + lith.B[1]*P + lith.B[2]*pow(P,2);
}
// T_liquidus
Doub Tl_px (const Doub& P, const lith& lith) {
  return lith.C[0] + lith.C[1]*P + lith.C[2]*pow(P,2);
}
// D(T_solidus)/DP
Doub dts_dp_px (const Doub& P, const lith& lith)  {
  return lith.A[1]+2*lith.A[2]*P;
}
// D(T_co)/DP
Doub dtco_dp_px (const Doub& P, const lith& lith) {
  return lith.B[1]+2*lith.B[2]*P;
}
// D(T_liquidus)/DP
Doub dtl_dp_px (const Doub& P, const lith& lith) {
  return lith.C[1]+2*lith.C[2]*P;
}
// F_cpx_out
Doub Fco_px () {
  const Doub rho_0=0.415405808583069;
  const Doub rho_1=0.318786433253217;

  return rho_0 + rho_1;
}

// PERIDOTITE
// F_cpx-out, pressure dependance accounts for the varying mass fraction of cpx in source
// as a function of the phase's stability
Doub Fcpx_o (const Doub& P, const lith& lith) {
  return lith.M/( lith.r[0]+lith.r[1]*P );
}
// T_cpx_out, temperature at which cpx is exhausted from the source
Doub Tcpx_o (Doub Fcpx_o, const Doub& P, const lith& lith) {
  return pow(Fcpx_o,1/lith.bt[0])*( Tlhrz(P,lith)-Ts_pd(P,lith) ) + Ts_pd(P,lith);
}
// D(Fcpx_o^(1/beta1))/DP
Doub dfcpxobt_dp (const Doub& P, const lith& lith) {
  return -1*lith.r[1]*(lith.M/lith.bt[0])*				\
    ( pow(Fcpx_o(P,lith),(1-lith.bt[0])/lith.bt[0]) )/pow(lith.r[0]+lith.r[1]*P,2);
}
// D(Fcpx_o))/DP
Doub dfcpxo_dp (const Doub& P, const lith& lith) {
  return -lith.r[1]*lith.M/pow(lith.r[0]+lith.r[1]*P,2);
}
// D(T_cpx_out)/DP
Doub dtcpxo_dp (Doub Fcpx_o, const Doub& P, const lith& lith) {
  return pow(Fcpx_o,1/lith.bt[0])*(dtlhrz_dp(P,lith)-dts_dp(P,lith)) +	\
    (Tlhrz(P,lith)-Ts_pd(P,lith))*						\
    dfcpxobt_dp(P,lith) + dts_dp(P,lith);
}


// PERIDOTITE
// dT/dF|p
Doub dt_df_p(const Doub& F, const Doub& P, const lith& lith) {
  if ( F<=Fcpx_o(P,lith) && F>0 ) {
    return pow( F,(1-lith.bt[0])/lith.bt[0] )*\
      ( Tlhrz(P,lith)-Ts_pd(P,lith) )/lith.bt[0];
  } else if (F>0) {
    return 1/lith.bt[1]*						\
      pow( (F-Fcpx_o(P,lith))/(1-Fcpx_o(P,lith)), (1-lith.bt[1])/lith.bt[1] )* \
      ( Tl_pd(P,lith)-Tcpx_o(Fcpx_o(P,lith),P,lith) )/( 1-Fcpx_o(P,lith) );
  } else {				// fix to not break when F=0
    return 1e6;}
}
// dT/dP|f
Doub dt_dp_f(const Doub& F, const Doub& P, const lith& lith) {
  if ( F<=Fcpx_o(P,lith) ) {
    return pow( F,1/lith.bt[0] )*( dtlhrz_dp(P,lith)-dts_dp(P,lith) ) + dts_dp(P,lith);
  } else {
    return dtcpxo_dp(Fcpx_o(P,lith),P,lith) +				\
      (dtl_dp(P,lith)-dtcpxo_dp(Fcpx_o(P,lith),P,lith))*		\
      pow( (F-Fcpx_o(P,lith))/(1-Fcpx_o(P,lith)), 1/lith.bt[1]) +	\
      ( Tl_pd(P,lith)-Tcpx_o(Fcpx_o(P,lith),P,lith) )/lith.bt[1] *	\
      pow( (F-Fcpx_o(P,lith))/(1-Fcpx_o(P,lith)), (1-lith.bt[1])/lith.bt[1] )* \
      ( (F-Fcpx_o(P,lith))*dfcpxo_dp(P,lith) - dfcpxo_dp(P,lith)*(1-Fcpx_o(P,lith)) )/ \
      pow(1-Fcpx_o(P,lith),2);
  }
}

// HARZBURGITE
Doub dt_df_p_hrz(const Doub& F, const Doub& P, const lith& lith) {
  // **************SPOOFED FOR NOW
  return 1;
}

Doub dt_dp_f_hrz(const Doub& F, const Doub& P, const lith& lith) {
  // **************SPOOFED FOR NOW
  return 1;
}

// PYROXENITE
// dT/dF|p
Doub dt_df_p_px(const Doub& F, const Doub& P, const lith& lith) {
  // these constants are from the lambart parameterisation of F(P,T)
  const Doub rho_0=0.415405808583069;
  const Doub rho_1=0.318786433253217;

  if ( F<Fco_px() && F>0 ) {
    // this is the melting regime defined by the Lambart parametrisation

    const Doub T0=Ts_px(P,lith);
    const Doub Tco=Tco_px(P,lith);
    const Doub Dt=(Tco-T0);

    return Dt*pow(rho_0*rho_0+4*rho_1*F,-0.5);
  } else if (F>0) {
    return ( Tl_px(P,lith)-Tco_px(P,lith) )/( (1-Fco_px())*lith.bt[1] )*	\
      pow( (F-Fco_px())/(1-Fco_px()), (1-lith.bt[1])/lith.bt[1] );
  } else {			// fix to not break when F=0
    return 1e6;
  }
}
// dT/dP|f
Doub dt_dp_f_px(const Doub& F, const Doub& P, const lith& lith) {
  // these constants are from the lambart parameterisation of F(P,T)
  const Doub rho_0=0.415405808583069;
  const Doub rho_1=0.318786433253217;

  if ( F<Fco_px() ) {

    const Doub T0=Ts_px(P,lith);
    const Doub Tco=Tco_px(P,lith);
    const Doub Dt=(Tco-T0);

    return dts_dp_px(P,lith) + (dtco_dp_px(P,lith)-dts_dp_px(P,lith))*	\
      (-rho_0+pow(rho_0*rho_0+4*rho_1*F,0.5))/(2*rho_1);
  } else if (F<1) {
    return dtco_dp_px(P,lith)+						\
      pow((F-Fco_px())/(1-Fco_px()),1/lith.bt[1])*			\
      (dtl_dp_px(P,lith)-dtco_dp_px(P,lith));
  } else {
    return 0;
  }
}

// DERIVATIVES FUNCTION FOR RK4 - df
void derivs_pf(const Doub P, VecDoub_I& y, VecDoub_O& dydx) {
  dydx[0]=P;
}

void derivs_av(const Doub P, VecDoub_I& y, VecDoub_O& dydx, const VecDoub pin) {
  // dPbardP=P'F, where P' is the streamline average P
  // px
  if (pin[8]>0)
    dydx[0]=pin[0]*(pin[4]/(1-pin[7]));
  else 
    dydx[0]=0;
  // pd
  if (pin[9]>0)
    dydx[1]=pin[1]*(pin[5]/(1-pin[7]));
  else
    dydx[1]=0;
  // hz
  if (pin[10]>0)
    dydx[2]=pin[2]*(pin[6]/(1-pin[7]));
  else
    dydx[2]=0;
  // bulk
  if (pin[11]>0)
    dydx[3]=pin[3]*(pin[7]/(1-pin[7]));
  else
    dydx[3]=0;
}

// DERIVATIVES FUNCTION FOR RK4 - dp
// peridotite first to melt integrals
void derivs_pd(const Doub P, VecDoub_I& y, VecDoub_O& dydx, const VecDoub phi) {

  // peridotite katz 2003
  lith pd;
  pd.A[0]=1085.7; pd.A[1]=132.9; pd.A[2]=-5.1; // degC:degC/GPa:degC/GPa^2
  pd.B[0]=1475; pd.B[1]=80; pd.B[2]=-3.2;      // degC:degC/GPa:degC/GPa^2
  pd.C[0]=1780; pd.C[1]=45; pd.C[2]=-2;	       // degC:degC/GPa:degC/GPa^2
  pd.r[0]=0.5; pd.r[1]=0.08;		       // cpx/melt:cpx/melt/GPa
  pd.M=0.15;				       // weight fraction cpx
  pd.bt[0]=1.5; pd.bt[1]=1.5;		       // 
  pd.Cp=1187;				       // J/(KgK)
  pd.alph[0]=30e-6; pd.alph[1]=68e-6;	       // /K
  pd.rho[0]=3.3e3; pd.rho[1]=2.9e3;	       // Kg/m^3
  pd.DS=407;				       // entropy of fusion, J/Kg/K

  // pyroxenite (KG1) Lambart
  lith px;
  px.A[0]=1095.36654354154; px.A[1]=124.061186889434; px.A[2]=-4.70657895065642;
  px.B[0]=1179.58449528221; px.B[1]=157.15408596868; px.B[2]=-11.0582998232629;
  // this final parameterisation of the true liquidus taken from Katz 2003
  px.C[0]=1780; px.C[1]=45; px.C[2]=-2;	  // degC:degC/GPa:degC/GPa^2
  px.Cp=1140;				  // J/(KgK)
  px.bt[0]=1.5; px.bt[1]=1.5;		  // exponent, only use for hrz section for kg1/2
  px.alph[0]=30e-6; px.alph[1]=68e-6;	  // /K
  px.rho[0]=3.3e3; px.rho[1]=2.9e3;	  // Kg/m^3
  px.DS=380;				  // entropy of fusion, J/Kg/K

  lith hrz;
  hrz.A[0]=1385.6; hrz.A[1]=132.9; hrz.A[2]=-5.1;
  hrz.B[0]=1780; hrz.B[1]=45; hrz.B[2]=-2;	// degC:degC/GPa:degC/GPa^2
  hrz.Cp=1000;				// J/(KgK)
  hrz.alph[0]=30e-6; hrz.alph[1]=68e-6;	// /K
  hrz.rho[0]=3.25e3; hrz.rho[1]=2.9e3;	// Kg/m^3
  hrz.DS=407;

  const Doub t0=273.15;
  const Doub g=9.81;

  // perform test to check that current T is below peridotite solidus
  // Doub Tpds=pd.A[0]+pd.A[1]*P+pd.A[2]*pow(P,2);
  // check for px solidus
  Doub Tpxs=px.A[0]+px.A[1]*P+px.A[2]*pow(P,2);
  // check for cpx out px
  // Doub Tpx_co=px.B[0] + px.B[1]*P + px.B[2]*pow(P,2);
  Doub Thrz=hrz.A[0]+hrz.A[1]*P+hrz.A[2]*pow(P,2);

  // track phi in the solid as the pyoxenite, peridotite and harzburgite fractions 
  // vary from melt generation
  VecDoub phi_tmp(3);
  phi_tmp[0]=(1-y[0])*phi[0]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);
  phi_tmp[1]=(1-y[1])*phi[1]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);
  phi_tmp[2]=(1-y[2])*phi[2]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);

  // Structure of differntials in y vector
  // y: Fpx Fpd Fhrz T tcpx tcpd tchrz 

  // four melting regimes:
  // 1. pd only
  // 2. px + pd
  // 3. px_hrz + pd
  // 4. px_hrz + pd_hrz + hrz (if the Tp is high)

  // pow(10,9) present in the differential terms with alpha and Cp because
  // these were from ds/dp in Pa and the dp need to be in GPa for these equations

  if (y[3]<Tpxs) {		// PERIDOTITE ONLY MELTING
    
    //// df/dP|s - f(P,T)
    // px
    dydx[0]=0;
    // pd
    dydx[1]=-1*( phi_tmp[0]*(px.Cp/(y[3]+t0)*dt_dp_f(y[1],P,pd)	\
			     - pow(10.0,9)*px.alph[0]/px.rho[0]) +	\
		 (phi_tmp[1])*(pd.Cp/(y[3]+t0)*dt_dp_f(y[1],P,pd) -	\
			       pow(10.0,9)*pd.alph[0]/pd.rho[0]) +	\
		 (phi_tmp[2])*(hrz.Cp/(y[3]+t0)*dt_dp_f(y[1],P,pd) -	\
			       pow(10.0,9)*hrz.alph[0]/hrz.rho[0]))/	\
      ( phi_tmp[1]*(pd.DS + pd.Cp/(y[3]+t0)*dt_df_p(y[1],P,pd)) +	\
	phi_tmp[0]*(px.Cp/(y[3]+t0)*dt_df_p(y[1],P,pd)) +		\
	phi_tmp[2]*(hrz.Cp/(y[3]+t0)*dt_df_p(y[1],P,pd)) );
    // hrz
    dydx[2]=0;
    // bulk df/dp
    dydx[4]=dydx[0]*phi_tmp[0]+dydx[1]*phi_tmp[1]+dydx[2]*phi_tmp[2];

    // df/dp integrals for each lithology
    dydx[5]=0;
    dydx[6]=phi_tmp[1]*dydx[1];
    dydx[7]=0;

    // dT/dP
    // all lithologies - assume thermal equilibrium maintained
    dydx[3]=dydx[1]*dt_df_p(y[1],P,pd) + dt_dp_f(y[1],P,pd);

  } else if ( y[3]<Thrz ) {	// PERIDOTITE + PYROXENITE MELTING
                                // either could be melting as hrz

    // df/dP|s - f(P,T)
    // px
    dydx[0]=-1*(							\
           ( phi_tmp[0]*px.Cp/(y[3]+t0)+phi_tmp[1]*pd.Cp/(y[3]+t0) + \
	     phi_tmp[2]*hrz.Cp/(y[3]+t0) )*dt_dp_f_px(y[0],P,px) -	\
	   pow(10.0,9)*( phi_tmp[0]*px.alph[0]/px.rho[0]+phi_tmp[1]*pd.alph[0]/pd.rho[0]+ \
		       phi_tmp[2]*hrz.alph[0]/hrz.rho[0] ) +		\
	   phi_tmp[1]*pd.DS*						\
	   (dt_dp_f_px(y[0],P,px)-dt_dp_f(y[1],P,pd))			\
	   /dt_df_p(y[1],P,pd))/					\
      ( phi_tmp[0]*px.DS + phi_tmp[1]*pd.DS*dt_df_p_px(y[0],P,px)/dt_df_p(y[1],P,pd) + \
	(phi_tmp[0]*px.Cp/(y[3]+t0)+phi_tmp[1]*pd.Cp/(y[3]+t0)+phi_tmp[2]*hrz.Cp/(y[3]+t0) \
	 )*dt_df_p_px(y[0],P,px));
    // pd
    dydx[1]=dt_df_p_px(y[0],P,px)/dt_df_p(y[1],P,pd) * dydx[0] +	\
      (dt_dp_f_px(y[0],P,px)-dt_dp_f(y[1],P,pd))/dt_df_p(y[1],P,pd);
    // hrz
    dydx[2]=0;
    // bulk df/dp
    dydx[4]=dydx[0]*phi_tmp[0]+dydx[1]*phi_tmp[1]+dydx[2]*phi_tmp[2];

    // df/dp integrals for each lithology
    dydx[5]=phi_tmp[0]*dydx[0];
    dydx[6]=phi_tmp[1]*dydx[1];
    dydx[7]=0;

    // dT/dP
    // all lithologies - assume thermal equilibrium maintained
    dydx[3]=dt_dp_f_px(y[0],P,px) + dt_df_p_px(y[0],P,px) * dydx[0];

  }  else {			// PERIDOTITE + HARZBURGITE MELTING
    // **************SPOOFED FOR NOW, i.e. no harzburgite melting
    dydx[0]=0;
    dydx[1]=0;
    dydx[2]=0;

    dydx[4]=0;

    dydx[5]=0;
    dydx[6]=0;
    dydx[7]=0;

    dydx[3]=0;
  }

  // Pressure at base of crust
  // PYROXENITE
  dydx[8]= -1*( y[5]/(1-y[4]) );
  // PERIDOTITE
  dydx[9]= -1*( y[6]/(1-y[4]) );
  // HARZBURGITE
  dydx[10]=0;
  // ALL
  dydx[11]= -1*( y[4]/(1-y[4]) );

  // - I[F/(1-Fbulk)]dP without phi scaling
  // PYROXENITE
  dydx[12]= -1*( y[0]/(1-y[4]) );
  // PERIDOTITE
  dydx[13]= -1*( y[1]/(1-y[4]) );
  // HARZBURGITE
  dydx[14]=0;

}

// pyroxenite first to melt integrals
void derivs_px(const Doub P, VecDoub_I& y, VecDoub_O& dydx, const VecDoub phi) {

  // peridotite katz 2003
  lith pd;
  pd.A[0]=1085.7; pd.A[1]=132.9; pd.A[2]=-5.1; // degC:degC/GPa:degC/GPa^2
  pd.B[0]=1475; pd.B[1]=80; pd.B[2]=-3.2;      // degC:degC/GPa:degC/GPa^2
  pd.C[0]=1780; pd.C[1]=45; pd.C[2]=-2;	       // degC:degC/GPa:degC/GPa^2
  pd.r[0]=0.5; pd.r[1]=0.08;		       // cpx/melt:cpx/melt/GPa
  pd.M=0.15;				       // weight fraction cpx
  pd.bt[0]=1.5; pd.bt[1]=1.5;		       // 
  pd.Cp=1187;				       // J/(KgK)
  pd.alph[0]=30e-6; pd.alph[1]=68e-6;	       // /K
  pd.rho[0]=3.3e3; pd.rho[1]=2.9e3;	       // Kg/m^3
  pd.DS=407;			               // entropy of fusion, J/Kg/K

  // pyroxenite (KG1) Lambart
  lith px;
  px.A[0]=1095.36654354154; px.A[1]=124.061186889434; px.A[2]=-4.70657895065642;
  px.B[0]=1179.58449528221; px.B[1]=157.15408596868; px.B[2]=-11.0582998232629;
  // this final parameterisation of the true liquidus taken from Katz 2003
  px.C[0]=1780; px.C[1]=45; px.C[2]=-2;	  // degC:degC/GPa:degC/GPa^2
  px.Cp=1140;				  // J/(KgK)
  px.bt[0]=1.5; px.bt[1]=1.5;             // exponent, only use for hrz section for kg1/2
  px.alph[0]=30e-6; px.alph[1]=68e-6;	  // /K
  px.rho[0]=3.3e3; px.rho[1]=2.9e3;	  // Kg/m^3
  px.DS=380;				  // entropy of fusion, J/Kg/K

  lith hrz;
  hrz.A[0]=1385.6; hrz.A[1]=132.9; hrz.A[2]=-5.1;
  hrz.B[0]=1780; hrz.B[1]=45; hrz.B[2]=-2;	// degC:degC/GPa:degC/GPa^2
  hrz.Cp=1000;				// J/(KgK)
  hrz.alph[0]=30e-6; hrz.alph[1]=68e-6;	// /K
  hrz.rho[0]=3.25e3; hrz.rho[1]=2.9e3;	// Kg/m^3
  hrz.DS=407;

  const Doub t0=273.15;
  const Doub g=9.81;

  // perform test to check that current T is below peridotite solidus
  Doub Tpds=pd.A[0]+pd.A[1]*P+pd.A[2]*pow(P,2);
  // check for px solidus
  // Doub Tpxs=px.A[0]+px.A[1]*P+px.A[2]*pow(P,2);
  // check for cpx out px
  // Doub Tpx_co=px.B[0] + px.B[1]*P + px.B[2]*pow(P,2);
  Doub Thrz=hrz.A[0]+hrz.A[1]*P+hrz.A[2]*pow(P,2);

  // track phi in the solid as the pyoxenite, peridotite and harzburgite fractions 
  // vary from melt generation
  VecDoub phi_tmp(3);
  phi_tmp[0]=(1-y[0])*phi[0]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);
  phi_tmp[1]=(1-y[1])*phi[1]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);
  phi_tmp[2]=(1-y[2])*phi[2]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);

  // Structure of differntials in y vector
  // y: Fpx Fpd Fhrz T tcpx tcpd tchrz 

  // four melting regimes:
  // 1. pd only
  // 2. px + pd
  // 3. px_hrz + pd
  // 4. px_hrz + pd_hrz + hrz (if the Tp is high)

  // pow(10,9) present in the differential terms with alpha and Cp because
  // these were from ds/dp in Pa and the dp need to be in GPa for these equations

  if (y[3]<Tpds) {		// PYROXENITE ONLY MELTING
    
    //// df/dP|s - f(P,T)
    // px
    dydx[0]=-1*( phi_tmp[0]*(px.Cp/(y[3]+t0)*dt_dp_f_px(y[0],P,px)	\
			     - pow(10.0,9)*px.alph[0]/px.rho[0]) +	\
		 (phi_tmp[1])*(pd.Cp/(y[3]+t0)*dt_dp_f_px(y[0],P,px) -	\
			       pow(10.0,9)*pd.alph[0]/pd.rho[0]) +	\
		 (phi_tmp[2])*(hrz.Cp/(y[3]+t0)*dt_dp_f_px(y[0],P,px) -	\
			       pow(10.0,9)*hrz.alph[0]/hrz.rho[0]))/	\
      ( phi_tmp[0]*(px.DS + px.Cp/(y[3]+t0)*dt_df_p_px(y[0],P,px)) +	\
	phi_tmp[1]*(pd.Cp/(y[3]+t0)*dt_df_p_px(y[0],P,px)) +		\
	phi_tmp[2]*(hrz.Cp/(y[3]+t0)*dt_df_p_px(y[0],P,px)) );
    // pd
    dydx[1]=0;
    // hrz
    dydx[2]=0;
    // bulk df/dp
    dydx[4]=dydx[0]*phi_tmp[0]+dydx[1]*phi_tmp[1]+dydx[2]*phi_tmp[2];

    // df/dp integrals for each lithology
    dydx[5]=phi_tmp[0]*dydx[0];
    dydx[6]=0;
    dydx[7]=0;

    // dT/dP
    // all lithologies - assume thermal equilibrium maintained
    dydx[3]=dydx[0]*dt_df_p_px(y[0],P,px) + dt_dp_f_px(y[0],P,px);

  } else if ( y[3]<Thrz ) {	// PERIDOTITE + PYROXENITE MELTING
                                // either could be melting as hrz

    // df/dP|s - f(P,T)
    // px
    dydx[0]=-1*(							\
           ( phi_tmp[0]*px.Cp/(y[3]+t0)+phi_tmp[1]*pd.Cp/(y[3]+t0) + \
	     phi_tmp[2]*hrz.Cp/(y[3]+t0) )*dt_dp_f_px(y[0],P,px) -	\
	   pow(10.0,9)*( phi_tmp[0]*px.alph[0]/px.rho[0]+phi_tmp[1]*pd.alph[0]/pd.rho[0]+ \
		       phi_tmp[2]*hrz.alph[0]/hrz.rho[0] ) +		\
	   phi_tmp[1]*pd.DS*						\
	   (dt_dp_f_px(y[0],P,px)-dt_dp_f(y[1],P,pd))			\
	   /dt_df_p(y[1],P,pd))/					\
      ( phi_tmp[0]*px.DS + phi_tmp[1]*pd.DS*dt_df_p_px(y[0],P,px)/dt_df_p(y[1],P,pd) + \
	(phi_tmp[0]*px.Cp/(y[3]+t0)+phi_tmp[1]*pd.Cp/(y[3]+t0)+phi_tmp[2]*hrz.Cp/(y[3]+t0)\
	 )*dt_df_p_px(y[0],P,px));
    // pd
    dydx[1]=dt_df_p_px(y[0],P,px)/dt_df_p(y[1],P,pd) * dydx[0] +	\
      (dt_dp_f_px(y[0],P,px)-dt_dp_f(y[1],P,pd))/dt_df_p(y[1],P,pd);
    // hrz
    dydx[2]=0;
    // bulk df/dp
    dydx[4]=dydx[0]*phi_tmp[0]+dydx[1]*phi_tmp[1]+dydx[2]*phi_tmp[2];

    // df/dp integrals for each lithology - scaled
    dydx[5]=phi_tmp[0]*dydx[0];
    dydx[6]=phi_tmp[1]*dydx[1];
    dydx[7]=0;

    // dT/dP
    // all lithologies - assume thermal equilibrium maintained
    dydx[3]=dt_dp_f_px(y[0],P,px) + dt_df_p_px(y[0],P,px) * dydx[0];

  }  else {			// PERIDOTITE + HARZBURGITE MELTING
    // **************SPOOFED FOR NOW, i.e. no harzburgite melting
    dydx[0]=0;
    dydx[1]=0;
    dydx[2]=0;

    dydx[4]=0;

    dydx[5]=0;
    dydx[6]=0;
    dydx[7]=0;

    dydx[3]=0;
  }

  // Pressure at base of crust
  // PYROXENITE
  dydx[8]= -1*( y[5]/(1-y[4]) );
  // PERIDOTITE
  dydx[9]= -1*( y[6]/(1-y[4]) );
  // HARZBURGITE
  dydx[10]=0;
  // ALL
  dydx[11]= -1*( y[4]/(1-y[4]) );

  // - I[F/(1-Fbulk)]dP without phi scaling
  // PYROXENITE
  dydx[12]= -1*( y[0]/(1-y[4]) );
  // PERIDOTITE
  dydx[13]= -1*( y[1]/(1-y[4]) );
  // HARZBURGITE
  dydx[14]=0;

}


// ROOT FINDING ALGORITHM FOR TEMPERATURE PRESSURE OF SOLIDUS INTERSECTION
// 
struct root_func_px {

  Doub Tp;			// K
  VecDoub phi;
  struct lith px, pd, hrz;
  root_func_px(const Doub tp, const struct lith& Px, const struct lith& Pd, const struct lith& Hrz, const VecDoub& Phi) : px(Px), pd(Pd), hrz(Hrz), Tp(tp), phi(Phi) {}

  Doub operator() (const Doub P) {
    const Doub t0=273.15;
    return (Tp+t0)*exp( ( phi[0]*px.alph[0]/px.rho[0] + phi[1]*pd.alph[0]/pd.rho[0] + \
			  phi[2]*hrz.alph[0]/hrz.rho[0] )*P*pow(10.0,9)/	\
			(phi[0]*px.Cp + phi[1]*pd.Cp + phi[2]*hrz.Cp) ) \
      - px.A[0] - px.A[1]*P - px.A[2]*pow(P,2) - t0;
  }

};

struct root_func_pd {

  Doub Tp;			// K
  VecDoub phi;
  struct lith px, pd, hrz;
  root_func_pd(const Doub tp, const struct lith& Px, const struct lith& Pd, const struct lith& Hrz, const VecDoub& Phi) : px(Px), pd(Pd), hrz(Hrz), Tp(tp), phi(Phi) {}

  Doub operator() (const Doub P) {
    const Doub t0=273.15;
    return (Tp+t0)*exp( ( phi[0]*px.alph[0]/px.rho[0] + phi[1]*pd.alph[0]/pd.rho[0] + \
			  phi[2]*hrz.alph[0]/hrz.rho[0] )*P*pow(10.0,9)/	\
			(phi[0]*px.Cp + phi[1]*pd.Cp + phi[2]*hrz.Cp) ) \
      - pd.A[0] - pd.A[1]*P - pd.A[2]*pow(P,2) - t0;
  }

};


Doub tp_px_i(const Doub& Tp, const Doub& Pi, const Doub& Pf, const struct lith& px, const struct lith& pd, const struct lith& hrz, const VecDoub& phi) {

  // user end function for calculating the point of solidus intersection
  Doub xacc=1e-6;		// x accuracy to aim for

  struct root_func_px Tcalc(Tp,px,pd,hrz,phi); // define function

  // run root bisection algorithm
  return rtbis(Tcalc,Pi,Pf,xacc);

}

Doub tp_pd_i(const Doub& Tp, const Doub& Pi, const Doub& Pf, const struct lith& px, const struct lith& pd, const struct lith& hrz, const VecDoub& phi) {

  // user end function for calculating the point of solidus intersection
  Doub xacc=1e-6;		// x accuracy to aim for

  struct root_func_pd Tcalc(Tp,px,pd,hrz,phi); // define function

  // run root bisection algorithm
  return rtbis(Tcalc,Pi,Pf,xacc);

}


/* -------------------------------------------------
   MAIN LOOP 
   -------------------------------------------------  */

int main (int argc, char *argv[]) {

  // OPTARGS
  int c;
  bool tc_mode=false,phi_mode=false,tp_set=false,dy_pmin=false;
  VecDoub tps(10); // max number of potential temperatures that can be accepted
  char *tp_tmp;
  VecDoub phi(3);		// fraction of lithologies 1) px 2) pd 3) hrz
  Doub Pmin=1.;			// pressure to decompress to, GPa

  while ((c = getopt(argc, argv, "m:Stx:p:h:T:")) != -1) {
    switch(c) {
    case 't':			// crustal thickness mode
      tc_mode = true;
      break;
    case 'm':			// set depth of maximum upwelling
      Pmin = atof(optarg);
      break;
    case 'S':			// set dynamic Pmin
      // should be 1 for true 0 for false, also initialise Pmin
      cerr << "Pmin set to dynamic mode" << endl;
      dy_pmin = true;
      Pmin=0.;
      break;
    case 'x':			// set fracttion of pyroxenite if run in fixed source mode
      phi_mode = true;
      phi[0]=atof(optarg);
      break;
    case 'p':			// set fracttion of peridotite if run in fixed source mode
      phi_mode = true;
      phi[1]=atof(optarg);
      break;
    case 'h':			// set fracttion of harzburgite if run in fixed source mode
      phi_mode = true;
      phi[2]=atof(optarg);
      break;
    case 'T':			// set mantle potential temperature
      tp_tmp = optarg;
      tp_set=true;
      break;
    }
  }

  if (tc_mode == false && phi_mode == false) {
    cout << "Set tc calculation or phi variation mode, -t or -p." << endl;
    return 1;
  }
  
  int tp_count=0;
  if (tp_set==true) {
    char *pick=strtok (tp_tmp, ",");
    while (pick != NULL)
      {
  	tps[tp_count]=atof(pick);
  	pick = strtok (NULL, ",");
  	tp_count++;
      }
  } else {
    cout << "Set potential temperature(s) for run using -T flag." << endl;
    return 1;
  }


  // LITHOLOGICAL INFORMATION
  // peridotite katz 2003
  lith pd;
  pd.A[0]=1085.7; pd.A[1]=132.9; pd.A[2]=-5.1; // degC:degC/GPa:degC/GPa^2
  pd.B[0]=1475; pd.B[1]=80; pd.B[2]=-3.2;      // degC:degC/GPa:degC/GPa^2
  pd.C[0]=1780; pd.C[1]=45; pd.C[2]=-2;	       // degC:degC/GPa:degC/GPa^2
  pd.r[0]=0.5; pd.r[1]=0.08;	               // cpx/melt:cpx/melt/GPa
  pd.M=0.15;				       // weight fraction cpx
  pd.bt[0]=1.5; pd.bt[1]=1.5;		       // 
  pd.Cp=1187;				       // J/(KgK)
  pd.alph[0]=30e-6; pd.alph[1]=68e-6;	       // /K
  pd.rho[0]=3.3e3; pd.rho[1]=2.9e3;	       // Kg/m^3
  pd.DS=407;                                   // entropy of fusion, J/Kg/K

  // pyroxenite (KG1) Lambart
  lith px;
  px.A[0]=1095.36654354154; px.A[1]=124.061186889434; px.A[2]=-4.70657895065642;
  px.B[0]=1179.58449528221; px.B[1]=157.15408596868; px.B[2]=-11.0582998232629;
  // this final parameterisation of the true liquidus taken from Katz 2003
  px.C[0]=1780; px.C[1]=45; px.C[2]=-2;	  // degC:degC/GPa:degC/GPa^2
  px.Cp=1140;				  // J/(KgK)
  px.bt[0]=1.5; px.bt[1]=1.5;		  // exponent, only use for hrz section for kg1/2
  px.alph[0]=30e-6; px.alph[1]=68e-6;	  // /K
  px.rho[0]=3.3e3; px.rho[1]=2.9e3;	  // Kg/m^3
  px.DS=380;				  // entropy of fusion, J/Kg/K

  lith hrz;
  hrz.A[0]=1385.6; hrz.A[1]=132.9; hrz.A[2]=-5.1;
  hrz.B[0]=1780; hrz.B[1]=45; hrz.B[2]=-2;	// degC:degC/GPa:degC/GPa^2
  hrz.Cp=1000;				// J/(KgK)
  hrz.alph[0]=30e-6; hrz.alph[1]=68e-6;	// /K
  hrz.rho[0]=3.25e3; hrz.rho[1]=2.9e3;	// Kg/m^3
  hrz.DS=407;

  Doub Tp;
  const Doub Pmax=15.;
  const Doub t0=273.15;		// K
  const int size_ar=15; // length of the input output arrays for kr4
  
  VecDoub y(size_ar), dydx(size_ar), yout(size_ar);
  VecDoub pf_px(1), pf_pd(1), pf_hz(1), pf_bulk(1), dpfdf_px(1), dpfdf_pd(1), \
    dpfdf_hz(1), dpfdf_bulk(1);
  VecDoub pf_px_out(1), pf_pd_out(1), pf_hz_out(1), pf_bulk_out(1);
  VecDoub pin(12), dpdp(4), pout(4), pbar(4);
  const Doub h=-1e-5;	// step size, GPa
  
  /* -------------------------------------------------
     PHI MODE

     Set a constant phi and loop over mantle potential
     temperatures, outputting a continuous record of P
     T F tc and dF/dP
     -------------------------------------------------  */

  if (phi_mode==true) {
    cerr << "Starting phi mode" << endl;

    // first make sure lithologies sum to one
    if (phi[0]+phi[1]+phi[2]-1 > 0.001) {
      cerr << "ERROR: lithology fractions must sum to 1" << endl;
      return 1;
    }

    // loop over potential temperatures to calculate at
    cerr << "Starting Tp loop" << endl;
    for (int i=0;i<tp_count;i++) {
      cerr << "T: " << "\r";
      cout << ">" << endl;
      Tp=tps[i];			// set Tp for inside loop

      // BOUNDARY CONDITION ASSIGNMENT
      /* calculate P and T of solidus intersection
	 given mantle potential temperature
	 pass potential temperature and pmin pmax for root search */
   
      Doub p1=tp_pd_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
      Doub p2=tp_px_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
      Doub p=0;
      Doub p_init=0;		// initial P of intersection
      if (p1>p2) {
	p=p1;
	y[3]=pd.A[0] + pd.A[1]*p + pd.A[2]*pow(p,2); // T initial, degC
	cerr << "pd intersects first" << endl;
      } else {
	p=p2;
	y[3]=px.A[0] + px.A[1]*p + px.A[2]*pow(p,2); // T initial, degC
	cerr << "px intersects first" << endl;
      }
      y[0]=0;			 // F initial, px
      y[1]=0;			 // F initial, pd
      y[2]=0;			 // F initial, hz
      y[4]=0;			 // mean F
      y[5]=0;			 // weighted dF/dP intrgral, px
      y[6]=0;			 // weighted dF/dP intrgral, pd
      y[7]=0;			 // weighted dF/dP intrgral, hz
      y[8]=0;			 // P at base of crust, px
      y[9]=0;			 // P at base of crust, pd
      y[10]=0;			 // P at base of crust, hz
      y[11]=0;			 // P at base of crust, total
      y[12]=0;			 // P at base of crust, pd - not modified for phi
      y[13]=0;			 // P at base of crust, hz - not modified for phi
      y[14]=0;			 // P at base of crust, total - not modified for phi
      // initialise P, F integral
      pf_px[0]=0;
      pf_pd[0]=0;
      pf_hz[0]=0;
      // initialise PF, pressure integral
      pbar[0]=0;pbar[1]=0;pbar[2]=0;pbar[3]=0;
      // assign initial P of intersection
      p_init=p;
      // print initial intersection of solidus to stderr
      fprintf(stderr, "%.0f: %.2f %.2f %.2f %.2f %.2f\n",
	      Tp, p, y[3], y[0], y[1], y[2]);

      // initialise dynamic Pmin, if selected
      if (dy_pmin==true)
	Pmin=0;

      // PRINT THE ADIABAT BETWEEN BETWEEN MAX P AND PX SOLIDUS INTERSECTION
      Doub p_tmp=Pmax;
      while ( p_tmp>=p ) {
	Doub t_tmp = (Tp+t0)*						\
	  exp( (phi[0]*px.alph[0]/px.rho[0]+phi[1]*pd.alph[0]/pd.rho[0]	\
		+ phi[2]*hrz.alph[0]/hrz.rho[0])*p_tmp*pow(10.0,9)/		\
	       (phi[0]*px.Cp + phi[1]*pd.Cp + phi[2]*hrz.Cp) ) - t0;
	
	printf("%.9f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f 0 0 0 0 0 0 0 0\n",	\
	       p_tmp, t_tmp,						\
	       yout[0], yout[1], yout[2], yout[4],			\
	       yout[8], yout[9], yout[10], yout[11],			\
	       Pmin, phi[0], phi[1], phi[2]);

	p_tmp=p_tmp+h;
      }

      // INTEGRATION LOOP
      // calculate intial derivatives
      // df is calcualted from the dp integrals for the df integrals
      //  - px pd hz bulk
      Doub df_px, df_pd, df_hz, df_bulk;
      while ( p>=Pmin ) {
	// perform integration
	if (p1>p2) {
	  // peridotite melts first
	  derivs_pd(p,y,dydx,phi);
	  rk4(y,dydx,p,h,yout,phi,derivs_pd);

	  // df integrals
	  df_px=yout[0]-y[0];
	  df_pd=yout[1]-y[1];
	  df_hz=yout[2]-y[2];
	  df_bulk=yout[4]-y[4];

	  // set y to now incremented
	  y=yout;

	  // base derivs takes: x y dydx
	  // base rk4 takes: y dydx x h yout derivs
	  if (yout[0]>0) {
	    derivs_pf(p,pf_px,dpfdf_px);
	    rk4(pf_px, dpfdf_px, p, df_px, pf_px_out, derivs_pf);
	    pf_px[0]=pf_px_out[0];	  }
	  if (yout[1]>0) {
	    derivs_pf(p,pf_pd,dpfdf_pd);
	    rk4(pf_pd, dpfdf_pd, p, df_pd, pf_pd_out, derivs_pf);
	    pf_pd[0]=pf_pd_out[0];	  }
	  if (yout[2]>0) {
	    derivs_pf(p,pf_hz,dpfdf_hz);
	    rk4(pf_hz, dpfdf_hz, p, df_hz, pf_hz_out, derivs_pf);
	    pf_hz[0]=pf_hz_out[0];	  }
	  if (yout[4]>0) {
	    derivs_pf(p,pf_bulk,dpfdf_bulk);
	    rk4(pf_bulk, dpfdf_bulk, p, df_bulk, pf_bulk_out, derivs_pf); 
	    pf_bulk[0]=pf_bulk_out[0];	  }

	  // calculate P'*F
	  // -give the P' info to pin
	  pin[0]=pf_px[0]/yout[0];pin[1]=pf_pd[0]/yout[1];
	  pin[2]=pf_hz[0]/yout[2];pin[3]=pf_bulk[0]/yout[4];
	  // -give the melt fractions to pin
	  pin[4]=yout[0];pin[5]=yout[1];pin[6]=yout[2];pin[7]=yout[4];
	  // -give the crustal pressure integrals to pin
	  // --these are only used to test whether the integration should be performed
	  pin[8]=yout[8];pin[9]=yout[9];pin[10]=yout[10];pin[11]=yout[11];

	  // modified derivs takes: x y dydx extra
	  // modified rk4 takes: y dydx x h yout extra derivs
	  // now go back to pressure integral including P' (streamline average P)
	  derivs_av(p,pbar,dpdp,pin);
	  rk4(pbar,dpdp,p,h,pout,pin,derivs_av);
	} else {
	  // pyroxenite melts first
	  derivs_px(p,y,dydx,phi);
	  rk4(y,dydx,p,h,yout,phi,derivs_px);

	  // df integrals
	  df_px=yout[0]-y[0];
	  df_pd=yout[1]-y[1];
	  df_hz=yout[2]-y[2];
	  df_bulk=yout[4]-y[4];

	  // set y to now incremented
	  y=yout;

	  // base derivs takes: x y dydx
	  // base rk4 takes: y dydx x h yout derivs
	  if (yout[0]>0) {
	    derivs_pf(p,pf_px,dpfdf_px);
	    rk4(pf_px, dpfdf_px, p, df_px, pf_px_out, derivs_pf);
	    pf_px[0]=pf_px_out[0];	  }
	  if (yout[1]>0) {
	    derivs_pf(p,pf_pd,dpfdf_pd);
	    rk4(pf_pd, dpfdf_pd, p, df_pd, pf_pd_out, derivs_pf);
	    pf_pd[0]=pf_pd_out[0];	  }
	  if (yout[2]>0) {
	    derivs_pf(p,pf_hz,dpfdf_hz);
	    rk4(pf_hz, dpfdf_hz, p, df_hz, pf_hz_out, derivs_pf);
	    pf_hz[0]=pf_hz_out[0];	  }
	  if (yout[4]>0) {
	    derivs_pf(p,pf_bulk,dpfdf_bulk);
	    rk4(pf_bulk, dpfdf_bulk, p, df_bulk, pf_bulk_out, derivs_pf); 
	    pf_bulk[0]=pf_bulk_out[0];	  }

	  // calculate P'*F
	  // -give the P' info to pin
	  pin[0]=pf_px[0]/yout[0];pin[1]=pf_pd[0]/yout[1];
	  pin[2]=pf_hz[0]/yout[2];pin[3]=pf_bulk[0]/yout[4];
	  // -give the melt fractions to pin
	  pin[4]=yout[0];pin[5]=yout[1];pin[6]=yout[2];pin[7]=yout[4];
	  // -give the crustal pressure integrals to pin
	  // --these are only used to test whether the integration should be performed
	  pin[8]=yout[8];pin[9]=yout[9];pin[10]=yout[10];pin[11]=yout[11];
	  // modified derivs takes: x y dydx extra
	  // modified rk4 takes: y dydx x h yout extra derivs
	  // now go back to pressure integral including P' (streamline average P)
	  derivs_av(p,pbar,dpdp,pin);
	  rk4(pbar,dpdp,p,h,pout,pin,derivs_av);
	}

	// set pbar to now incremented
	pbar=pout;

	// PRINT imcremental P, T, F_bar information
	// Note, y: Fpx Fpd Fhrz T tcpx tcpd tchrz 
	// 1:P_cur     2:T_cur
	// 3:F_px      4:F_pd       5:F_hrz       6: bulk F
	// 7:P_px      8:P_pd       9:P_hrz      10:P_bulk
	// 11: Pmin
	// 12/13/14: px/pd/hz remaining solid fraction
	// 15/16/17/18: mean melt fraction for each lithology and bulk, px/pd/hz/bulk
	// 19/20/21/22: mean P FRACTIONAL melting, each lithology and bulk, px/pd/hz/bulk
	// note that both mean calculations apply modification to output from integrals
	// - so that the result is appropriate for a trapezoid melting region
	// - without this modication it would be like calculating the integral
	// - for increasingly tall triangular melt regions.
	// NOTES:
	// need -1 factor for crustal thickness as integral returns negative tc
	// --- a consequence of the reversed limits
	VecDoub phi_tmp(3);
	phi_tmp[0]=(1-y[0])*phi[0]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);
	phi_tmp[1]=(1-y[1])*phi[1]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);
	phi_tmp[2]=(1-y[2])*phi[2]/((1-y[0])*phi[0] + (1-y[1])*phi[1] + (1-y[2])*phi[2]);

	printf("%.9f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",	\
	       p, yout[3],						\
	       yout[0], yout[1], yout[2], yout[4],			\
	       yout[8], yout[9], yout[10], yout[11],			\
	       Pmin,							\
	       phi_tmp[0], phi_tmp[1], phi_tmp[2],			\
	       (-1*yout[12]+yout[0]*(Pmin-p))/(Pmin-p_init),		\
	       (-1*yout[13]+yout[1]*(Pmin-p))/(Pmin-p_init),		\
	       (-1*yout[14]+yout[2]*(Pmin-p))/(Pmin-p_init),		\
	       (-1*yout[11]+yout[4]*(Pmin-p))/(Pmin-p_init),		\
	       -1*(pout[0]+pin[0]*yout[0]*(Pmin-p))/(yout[12]-yout[0]*(Pmin-p)), \
	       -1*(pout[1]+pin[1]*yout[1]*(Pmin-p))/(yout[13]-yout[1]*(Pmin-p)), \
	       -1*(pout[2]+pin[2]*yout[2]*(Pmin-p))/(yout[14]-yout[2]*(Pmin-p)), \
	       -1*(pout[3]+pin[3]*yout[4]*(Pmin-p))/(yout[11]-yout[4]*(Pmin-p)));

	// increment pressure
	p=p+h;

	if ( dy_pmin==true ) {
	  Pmin=yout[11];
	}

      }
    }
  }
  
  /* ---------------------------------------------------
     TC MODE MODE

     scan over range of lithology mass fractions (phi), 
     predicting crustal thickness for each peridotite-
     pyroxenite-harzburgite source mixture, Tp is fixed for
     each iteration across phi
     ---------------------------------------------------  */

  if (tc_mode==true) {
    // amount to increment lithologies
    const Doub alph_inc=0.025;

    // POTENTIAL TEMPERATURE LOOP
    for (int i=0;i<tp_count;i++) { // loop over potential temperatures to calculate at
      cout << ">" << endl;
      Tp=tps[i];		// set Tp for inside loop

      // LITHOLOGY LOOP
      // use alpha to represent phi[2], mass fraction harzburgite
      // use beta to represent relative masses of phi[0] and phi[1] 
      // (equivalent to two lithology original phi)
      // 
      //    hz    ^
      //   /  \   | alpha, 0-1
      //  pd--px
      //    -> beta, 0-1
      //    
      Doub alph=0.;
      while (alph<= 1) {
	// scale down bet inc as alph increases
	Doub bet=0.;
	Doub bet_inc=1/(1+(1-alph)/alph_inc);
	while (bet<1) {
	  // track progress
	  cerr << alph << " " << bet << "\r";
	  
	  // set phi
	  phi[0]=bet*(1-alph);
	  phi[1]=(1-bet)*(1-alph);
	  phi[2]=alph;

	  // set dynamic pmin
	  if ( dy_pmin==true )
	    Pmin=0;
	  
	  // Boundary condition assignment
	  /* calculate P and T of solidus intersection
	     given mantle potential temperature
	     pass potential temperature and pmin pmax for root search */
	  Doub p1=tp_pd_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
	  Doub p2=tp_px_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
	  Doub p=0;
	  if (p1>p2) {
	    p=p1;
	    y[3]=pd.A[0] + pd.A[1]*p + pd.A[2]*pow(p,2); // T initial, degC
	  } else {
	    p=p2;
	    y[3]=px.A[0] + px.A[1]*p + px.A[2]*pow(p,2); // T initial, degC
	  }

	  y[0]=0;  // F initial, px
	  y[1]=0;  // F initial, pd
	  y[2]=0;  // F initial, hz
	  y[4]=0;  // mean F
	  y[5]=0;  // weighted dF/dP intrgral, px
	  y[6]=0;  // weighted dF/dP intrgral, pd
	  y[7]=0;  // weighted dF/dP intrgral, hz
	  y[8]=0;  // P at base of crust, px
	  y[9]=0;  // P at base of crust, pd
	  y[10]=0; // P at base of crust, hz
	  y[11]=0; // P at base of crust, total
	  y[12]=0; // P at base of crust, pd - not modified for phi
	  y[13]=0; // P at base of crust, hz - not modified for phi
	  y[14]=0; // P at base of crust, total - not modified for phi

	  // INTEGRATION LOOP
	  Doub Ppx5, Plz5, Phz5;
	  while ( p>=Pmin ) {	// calculate intial derivatives
	    // perform integration
	    if (p1>p2) {
	      // peridotite melts first
	      derivs_pd(p,y,dydx,phi);
	      rk4(y,dydx,p,h,yout,phi,derivs_pd);
	    } else {
	      // pyroxenite melts first
	      derivs_px(p,y,dydx,phi);
	      rk4(y,dydx,p,h,yout,phi,derivs_px);
	    }

	    // set y to now incremented
	    y=yout;
	    
	    if (yout[1]<=0.05) {
	      Ppx5=yout[8];
	      Plz5=yout[9];
	      Phz5=yout[10];
	    }

	    // increment pressure
	    p=p+h;
	    if (dy_pmin==true) {
	      Pmin=y[11];
	    }
	  }
	  
	  // print: incremental phi_px, phi_pd, phi_hrz,
	  //                    P_px, P_pd, P_hz,
	  //                    P_tot, P_tot, Pmin
	  printf("%.9f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",	\
                 phi[0], phi[1], phi[2],				\
		 yout[8], yout[9], yout[10],				\
		 yout[8]+yout[9]+yout[10], yout[11], Pmin,		\
		 Ppx5, Plz5, Phz5, yout[1]);
	  
	  // increment counters
	  bet=bet+bet_inc;
	}

	// FINAL RUN OF BETA WITH BET=1 ------------------------------------
	bet=1;
	// track progress
	cerr << alph << " " << bet << "\r";
	
	// set phi
	phi[0]=bet*(1-alph);
	phi[1]=(1-bet)*(1-alph);
	phi[2]=alph;
	
	// Boundary condition assignment
	/* calculate P and T of solidus intersection
	   given mantle potential temperature
	   pass potential temperature and pmin pmax for root search */
	Doub p1=tp_pd_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
	Doub p2=tp_px_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
	Doub p=0;
	if (p1>p2) {
	  p=p1;
	  y[3]=pd.A[0] + pd.A[1]*p + pd.A[2]*pow(p,2); // T initial, degC
	} else {
	  p=p2;
	  y[3]=px.A[0] + px.A[1]*p + px.A[2]*pow(p,2); // T initial, degC
	}
	
	y[0]=0;	y[1]=0;	y[2]=0;
	y[4]=0;	y[5]=0;	y[6]=0; y[7]=0; y[8]=0; y[9]=0;
	y[10]=0; y[11]=0; y[12]=0; y[13]=0; y[14]=0;

	// re-initialise Pmin for dynamic Pmin mode
	if (dy_pmin==true)
	  Pmin=0;
	
	// INTEGRATION LOOP
	Doub Ppx5, Plz5, Phz5;
	while ( p>=Pmin ) {	// calculate intial derivatives
	  // perform integration
	  if (p1>p2) {
	    // peridotite melts first
	    derivs_pd(p,y,dydx,phi);
	    rk4(y,dydx,p,h,yout,phi,derivs_pd);
	  } else {
	    // pyroxenite melts first
	    derivs_px(p,y,dydx,phi);
	    rk4(y,dydx,p,h,yout,phi,derivs_px);
	  }
	  
	  // set y to now incremented
	  y=yout;

	  if (yout[1]<=0.05) {
	    Ppx5=yout[8];
	    Plz5=yout[9];
	    Phz5=yout[10];
	  }
	  
	  // increment pressure
	  p=p+h;
	  
	  if (dy_pmin==true) {
	    Pmin=y[11];
	  }
	}
	
	// print: incremental phi_px, phi_pd, phi_hrz,
	//                    P_px, P_pd, P_hz,
	//                    P_tot, P_tot, Pmin
	printf("%.9f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",	\
	       phi[0], phi[1], phi[2],					\
	       yout[8], yout[9], yout[10],				\
	       yout[8]+yout[9]+yout[10], yout[11], Pmin,		\
	       Ppx5, Plz5, Phz5, yout[1]);
	
	alph=alph+alph_inc;
      }
      
      // FINAL RUN OF BETA WITH ALPH=1 ---------------------------------------
      alph=1;Doub bet=0;
      // track progress
      cerr << alph << " " << bet << "\r";
      
      // set phi
      phi[0]=bet*(1-alph);
      phi[1]=(1-bet)*(1-alph);
      phi[2]=alph;
      
      // Boundary condition assignment
      /* calculate P and T of solidus intersection
	 given mantle potential temperature
	 pass potential temperature and pmin pmax for root search */
      Doub p1=tp_pd_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
      Doub p2=tp_px_i(Tp, 0., Pmax, px, pd, hrz, phi); // P initial, GPa
      Doub p=0;
      if (p1>p2) {
	p=p1;
	y[3]=pd.A[0] + pd.A[1]*p + pd.A[2]*pow(p,2); // T initial, degC
      } else {
	p=p2;
	y[3]=px.A[0] + px.A[1]*p + px.A[2]*pow(p,2); // T initial, degC
      }

      y[0]=0;  y[1]=0;  y[2]=0;
      y[4]=0;  y[5]=0;  y[6]=0; y[7]=0; y[8]=0; y[9]=0;
      y[10]=0; y[11]=0; y[12]=0; y[13]=0; y[14]=0;
      
      // re-initialise Pmin for dynamic Pmin mode
      if (dy_pmin==true)
	Pmin=0;
      
      // INTEGRATION LOOP
      Doub Ppx5, Plz5, Phz5;
      while ( p>=Pmin ) {	// calculate intial derivatives
	// perform integration
	if (p1>p2) {
	  // peridotite melts first
	  derivs_pd(p,y,dydx,phi);
	  rk4(y,dydx,p,h,yout,phi,derivs_pd);
	} else {
	  // pyroxenite melts first
	  derivs_px(p,y,dydx,phi);
	  rk4(y,dydx,p,h,yout,phi,derivs_px);
	}
	
	// set y to now incremented
	y=yout;

	if (yout[1]<=0.05) {
	  Ppx5=yout[8];
	  Plz5=yout[9];
	  Phz5=yout[10];
	}
	
	// increment pressure
	p=p+h;
	
	if (dy_pmin==true) {
	  Pmin=y[11];
	}
      }
      
      // print: incremental phi_px, phi_pd, phi_hrz,
      //                    P_px, P_pd, P_hz,
      //                    P_tot, P_tot, Pmin
      printf("%.9f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", \
	     phi[0], phi[1], phi[2],					\
	     yout[8], yout[9], yout[10],				\
	     yout[8]+yout[9]+yout[10], yout[11], Pmin,			\
	     Ppx5, Plz5, Phz5, yout[1]);
      
    }
  }    
  return 0;
}
  
