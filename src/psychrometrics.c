/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 

#include <stdio.h>
#include <math.h>
#include "psychrometrics.h"
#include "svp.h"



static double 
c_pa( void )
{
   return 1.006;
}

/*
static double
v_sp( double T, double P, double W)      // specific volume 
{
    return  0.287055*T*( 1.0 + 1.6078*W ) / P   ;
}
*/

static double                 /* E4a */
h_W_db(double W, double t_db) /* enthalpy given W and dry bulb */
{
  return c_pa()*t_db + W*(2501 + 1.805*t_db);
}

static double                  /* E4c */
W_h_t(double h, double t_db)   /* humidity ratio given enthalpy and dry bulb */
{
   return  (h - c_pa() * t_db) / ( 2501 + 1.805 * t_db) ;
}

static double                  /* E5 */ 
W_p_pw( double p, double p_w)  /* humidity ratio given pressure and vp */
{
    return 0.62198*( p_w / ( p - p_w ));
}



/* Find the wet bulb temperature given dry bulb, dew point and pressure.
   This procedure is needed several times and so it is given its own function
*/
   /*  The T_wb_zero is used by T_wet_bulb, it is the equation we are trying 
    *  to find the root of. This function seems to have more than one root but
    *  we want to use only the first one. 
   */

static int wb_NO_ROOT = 0; /* error checking hack */

static double
T_wb_zero(double t_wb, double t_dp, double t_db, double P)
{
   double Wsat_wb, W, retval;

   W = W_p_pw( P, p_sat(t_dp + 273.15 ));
   Wsat_wb = W_p_pw( P, p_sat(t_wb + 273.15));
   
   retval = (2501*Wsat_wb - (2501 + 1.805*t_db)*W - c_pa()*t_db)/
                      (2.381*Wsat_wb - 4.186*W - c_pa())  -  t_wb;
 
   return retval;
}
static double
T_wet_bulb( double t_dp, double t_db, double P)
{
   int i ;
   double t1, t2, tm;
   double  f1, f2, fm;
  
   /* try to bracket root  */
   t1 = t_dp;         /* dew point must be lower than wet bulb */     
                      /* so this is a good lower value */
   f1 = T_wb_zero(t1, t_dp, t_db, P);
   /* go in steps until sign changes */
   t2 = t1;
   do
   {
     t2 = t2 + 1;
     f2 = T_wb_zero(t2, t_dp, t_db, P);
     if( t2 > t_db ) wb_NO_ROOT=1; 
   }while( f1*f2 > 0.0 ); 

   /* home in on root using bisection */
   i = 0;
   do{  
      tm = (t1 + t2) / 2;
      fm = T_wb_zero(tm, t_dp, t_db, P);
      if( f1*fm < 0.0) 
      {
          t2 = tm;
          f2 = fm;
      }else
      {
          t1 = tm;
          f1 = fm;
      }
      i = i + 1;
    }while( fabs(t2-t1) > 0.0001 && i < 100);
   return tm;
}
/* end T_wet_bulb */


/*   range checking values */

float min_db = -10;
float max_db = 110;

float min_RH = 0;
float max_RH = 100;

float min_W = 0;
float max_W = 0.1; /* ?? */

/*
    Start the definitions of the psycrometric state equations 
 */

int 
P_db_dp( PsyState *psy)                    /* ref:p2 */
{
   double p_ws, p_w, W, t_wb;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->T_dew >= psy->T_db ) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   p_w = p_sat(psy->T_dew + 273.15);   
   W = W_p_pw( psy->P, p_w);
   t_wb = T_wet_bulb(psy->T_dew, psy->T_db, psy->P);

   psy->h = h_W_db(W, psy->T_db);  
   psy->RH = 100 * p_w / p_ws ;
   psy->T_wb = t_wb;
   psy->W = W;

   return NO_ERROR;
}

int 
P_db_W (PsyState *psy)                      /* ref:p3 */
{
   double p_w, p_ws, t_dp, t_wb;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
                 // invalid combinations possible

   p_ws = p_sat(psy->T_db + 273.15);
   p_w = psy->P*psy->W/(psy->W + 0.62198);    
   t_dp = T_sat(p_w) - 273.15;
   t_wb = T_wet_bulb(t_dp, psy->T_db, psy->P);
 
   psy->h = h_W_db(psy->W, psy->T_db);
   psy->RH = 100 * p_w / p_ws ;
   psy->T_dew = t_dp;
   psy->T_wb = t_wb;

   return NO_ERROR;
}

int
P_db_wb(PsyState *psy)                     /* ref:p6 */               
{
   double p_ws, p_sat_wb, Wsat_wb, W, p_w;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->T_wb >= psy->T_db ) return VAR_OUT_OF_RANGE;
                                            
   p_ws = p_sat(psy->T_db + 273.15);         
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   Wsat_wb = W_p_pw( psy->P, p_sat_wb );
   W = (Wsat_wb*(2501-2.381*psy->T_wb) - c_pa()*(psy->T_db - psy->T_wb))/
                               (2501 + 1.805*psy->T_db - 4.186*psy->T_wb);
   p_w = psy->P*W/(W + 0.62198);
 
   psy->h = h_W_db(W, psy->T_db);
   psy->RH = 100*p_w/p_ws;
   psy->T_dew = T_sat(p_w) - 273.15;
   psy->W = W;

   return NO_ERROR;
}

int 
P_db_rh(PsyState *psy)                     /* ref:p8 */
{
   double p_ws, p_w, W, t_dp, t_wb;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   p_w = p_ws * psy->RH/100 ;
   W = W_p_pw( psy->P, p_w ); 
   t_dp = T_sat(p_w) - 273.15;
   t_wb = T_wet_bulb( t_dp, psy->T_db, psy->P);
   
   psy->h = h_W_db(W, psy->T_db);
   psy->W = W;
   psy->T_dew = t_dp;
   psy->T_wb = t_wb;

   return NO_ERROR;
}

int 
P_db_h (PsyState *psy)              /* ref:p12 */
{
   double p_ws, W, p_w, t_dp, t_wb;
   
   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   // h can be too high 

   p_ws = p_sat(psy->T_db + 273.15);       
   W = W_h_t( psy->h, psy->T_db);   
   p_w = psy->P*W / (W + 0.62198);
   t_dp = T_sat(p_w) - 273.15;
   t_wb = T_wet_bulb( t_dp, psy->T_db, psy->P);

   psy->RH = 100 * p_w / p_ws;
   psy->W = W;
   psy->T_wb = t_wb;
   psy->T_dew = t_dp;
 
   return NO_ERROR;
}

int 
P_wb_dp(PsyState *psy)                     /* ref:p16 */
{
   double p_w, p_sat_wb, Wsat_wb, W, t_db, p_ws;

   if( psy->T_wb < psy->T_dew ) return VAR_OUT_OF_RANGE; 

   p_w = p_sat(psy->T_dew + 273.15);
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   Wsat_wb =  W_p_pw( psy->P, p_sat_wb);     
   W     =    W_p_pw( psy->P, p_w);
   t_db = ((2501 - 2.381*psy->T_wb)*Wsat_wb 
               - (2501 - 4.186*psy->T_wb)*W + c_pa()*psy->T_wb)
               / ( c_pa() + 1.805*W );
   p_ws = p_sat(t_db + 273.15);

   psy->h = h_W_db( W, t_db);
   psy->RH = 100*p_w/p_ws;
   psy->T_db = t_db;
   psy->W = W;

   return NO_ERROR;
}

int 
P_wb_W (PsyState *psy)              /* ref:p17 */
{
   double p_w, t_dp, p_sat_wb, Wsat_wb, t_db, p_ws;
                                                      // wb can be too low
   p_w = psy->P*psy->W / (psy->W + 0.62198);
   t_dp = T_sat(p_w) - 273.15;
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   Wsat_wb = W_p_pw( psy->P, p_sat_wb);
   t_db = ((2501 - 2.381*psy->T_wb)*Wsat_wb 
           - (2501 - 4.186*psy->T_wb)*psy->W + c_pa()*psy->T_wb )
           / (c_pa() + 1.805*psy->W);
   p_ws = p_sat(t_db + 273.15);

   psy->h = h_W_db( psy->W, t_db);
   psy->RH = 100*p_w/p_ws;
   psy->T_db = t_db;
   psy->T_dew = t_dp ;

   return NO_ERROR;
}


/* P_wb_rh is complicated; we need to solve for the dry bulb temerature
   by finding a root. This function is called zero1, params contains 
   various numbers needed by zero1 */
struct params{
    double Wsat_wb;
    double t_wb;
    double rh;
    double P;  };

static double 
zero1( double t_db, struct params pms)
{
   double p_ws, p_w, W, z;

   p_ws = p_sat( t_db +273.15 );
   p_w = pms.rh * p_ws;
   W = 0.62198*p_w / (pms.P - p_w);
   
   z = ((2501 - 2.381*pms.t_wb)*pms.Wsat_wb - c_pa()*(t_db - pms.t_wb))
        / (2501 + 1.805*t_db - 4.186*pms.t_wb) -  W;
   return z;
}
int 
P_wb_rh(PsyState *psy)                       /* ref:p21 */
{
   double p_sat_wb, Wsat_wb, t_db, p_ws, p_w, W;
   struct params pms;

   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;
   
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   Wsat_wb = W_p_pw( psy->P, p_sat_wb) ;
   
   pms.Wsat_wb = Wsat_wb;
   pms.t_wb = psy->T_wb;
   pms.rh  = psy->RH / 100;
   pms.P  = psy->P;
   
   double t_l, t_r, t_m, f_l, f_r, f_m;
   int i = 0;
  
   t_l = psy->T_wb;   /* wb less than db, good guess for lower limit */
   t_r = t_l + 50;      /* 50 arbitrary number, UGH */  
   f_l = zero1( t_l, pms);      
   f_r = zero1( t_r, pms);
   if( f_l*f_r > 0.0 ) return NO_ROOT;
   do{
      t_m =  (t_r + t_l) / 2;
      f_m = zero1( t_m, pms);
      if( f_m*f_l < 0 ){
         t_r = t_m;
         f_r = f_m;
      }
      else{
         t_l = t_m;
         f_l = f_m;  
      }
      i = i + 1;
    }while( fabs( (t_l-t_r)/2 ) > 0.001 && i < 100);
    t_db = t_l;
   
   p_ws = p_sat(t_db + 273.15);
   p_w = psy->RH*p_ws/100;
   W = 0.62198*p_w/(psy->P - p_w);

   psy->h = h_W_db(W, t_db);
   psy->T_db = t_db;
   psy->W = W;
   psy->T_dew =  T_sat(p_w) - 273.15;

   return NO_ERROR;
}
/* end of P_wb_rh */

int 
P_dp_rh(PsyState *psy)                   /* ref:p29 */
{
   double p_w, p_ws, t_db, W, t_wb;

   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   p_w = p_sat(psy->T_dew + 273.15);
   p_ws = 100*p_w/psy->RH;
   t_db = T_sat(p_ws) - 273.15;
   W = W_p_pw( psy->P, p_w);
   t_wb = T_wet_bulb( psy->T_dew, t_db, psy->P);

   psy->T_db = t_db;
   psy->T_wb = t_wb;
   psy->h = h_W_db(W, t_db);
   psy->W = W;

   return NO_ERROR;
}

int 
P_dp_h (PsyState *psy)              /* ref:p31 */
{
   double p_w, W, t_db, p_ws, t_wb;
                                       // invalid combinations possible
				       // h can be too low
   p_w = p_sat(psy->T_dew + 273.15);
   W = W_p_pw( psy->P, p_w);  
   t_db = (psy->h - 2501*W)/(c_pa() + 1.805*W);
   p_ws = p_sat(t_db + 273.15);
   t_wb = T_wet_bulb( psy->T_dew, t_db, psy->P);

   psy->T_db = t_db;
   psy->T_wb = t_wb;
   psy->RH = 100*p_w/p_ws;
   psy->W = W;

   return NO_ERROR;
}

int 
P_rh_W (PsyState *psy)                 /* ref:p32 */
{
   double p_w, t_db, p_ws, t_dp, t_wb;

   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   p_w = psy->P*psy->W/(psy->W + 0.62198);
   p_ws = 100*p_w/psy->RH;  
   t_dp = T_sat(p_w) - 273.15;
   t_db = T_sat(p_ws) - 273.15;
   t_wb = T_wet_bulb( t_dp, t_db, psy->P);
   
   psy->h = h_W_db(psy->W, t_db);
   psy->T_db = t_db;
   psy->T_wb = t_wb;
   psy->T_dew = t_dp;

   return NO_ERROR;
}

int 
P_W_h  (PsyState *psy)                 /* ref:p34 */
{
   double t_db, p_ws, p_w, t_dp;
                                               // invalid combinations possible
					       // h can be too low
   t_db = (psy->h - 2501*psy->W)/(c_pa() + 1.805*psy->W);
   p_ws = p_sat(t_db + 273.15);
   p_w = psy->P*psy->W/(psy->W + 0.62198);
   t_dp = T_sat(p_w) - 273.15;

   psy->T_wb = T_wet_bulb( t_dp, t_db, psy->P);
   psy->T_db = t_db;
   psy->T_dew = t_dp;
   psy->RH = 100*p_w/p_ws;

   return NO_ERROR;
}

/* P_rh_h needs to find dry bulb temperature by the root 
   of a function called zero2 */
static double
zero2(double t_db, double h, double RH, double P)
{  
   double p_ws, p_w, W, z;

   p_ws = p_sat(t_db + 273.15);
   p_w = p_ws * RH/100;
   W = W_p_pw( P, p_w);
   
   z = c_pa()*t_db + W*(2501 + 1.805*t_db) - h;
   
   return z;
}

int P_rh_h (PsyState *psy)            /* ref:p33 */
{
   double t_db, p_w, W, t_dp;
   
   double t_l, t_r, t_m, f_l, f_r, f_m;
   
   if( psy->RH < min_RH || psy->RH > max_RH ) return VAR_OUT_OF_RANGE;
   
   int i = 0;
  
   /* NOTE: I have to guess limits with no clues. How to improve?? */
   t_l = -40;
   t_r = 100;
   f_l = zero2( t_l, psy->h, psy->RH, psy->P);      
   f_r = zero2( t_r, psy->h, psy->RH, psy->P);
   if( f_l*f_r > 0.0 ) printf("ERROR root in zero2") ;
   do{
      t_m =  (t_r + t_l) / 2;
      f_m = zero2( t_m, psy->h, psy->RH, psy->P);
      if( f_m*f_l < 0 ){
         t_r = t_m;
         f_r = f_m;
      }
      else{
         t_l = t_m;
         f_l = f_m;  
      }
      i = i + 1;
    }while( fabs( (t_l-t_r)/2 ) > 0.001 && i < 100);
   t_db = t_l;
      
   p_w = p_sat(t_db + 273.15)*psy->RH/100;
   W = 0.62198*p_w/(psy->P - p_w);
   t_dp = T_sat(p_w) - 273.15;
   
   psy->T_db = t_db;
   psy->T_wb = T_wet_bulb( t_dp, t_db, psy->P);
   psy->T_dew = t_dp;
   psy->W = W;

   return NO_ERROR;
}

int 
W_rh_h  (PsyState *psy)                 /* ref:p35 */
{
   double t_db, p_ws, p_w;

   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   t_db = (psy->h - 2501*psy->W) / (c_pa()  + 1.805*psy->W);
   p_ws = p_sat( t_db + 273.15);
   p_w = p_ws * psy->RH/100;

   psy->T_db = t_db;
   psy->T_dew = T_sat(p_w) - 273.15;
   psy->P = 0.62198*p_w / psy->W + p_w;
   psy->T_wb = T_wet_bulb( psy->T_dew, t_db, psy->P);

   return NO_ERROR;
}

int 
db_h_dp(PsyState *psy)                 /* ref:p14 */
{
   double p_ws, W, p_w, p_tot;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->T_dew > psy->T_db ) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   W = (psy->h - c_pa()*psy->T_db) / (2501 + 1.805*psy->T_db);
   p_w = p_sat(psy->T_dew + 273.15);
   p_tot = 0.62198*p_w / W + p_w;
   
   psy->T_wb = T_wet_bulb( psy->T_dew, psy->T_db, p_tot);
   psy->RH = 100*p_w / p_ws;
   psy->P = p_tot;
   psy->W = W;
 
   return NO_ERROR;
}

int
wb_dp_W(PsyState *psy)                       /* ref:p20 */
{
   double p_w, p_tot, p_sat_wb, Ws_wb, t_db, p_ws;

   p_w = p_sat(psy->T_dew + 273.15);
   p_tot = 0.62198*p_w / psy->W + p_w;
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   Ws_wb = W_p_pw( p_tot, p_sat_wb);
   t_db = ((2501 - 2.381*psy->T_wb)*Ws_wb - (2501 - 4.186*psy->T_wb)*psy->W 
       + c_pa()*psy->T_wb) / (c_pa() + 1.805*psy->W);
   p_ws = p_sat(t_db + 273.15);

   psy->h = h_W_db( psy->W, t_db);
   psy->RH = 100*p_w/p_ws;
   psy->T_db = t_db;
   psy->P = p_tot;

   return NO_ERROR;
}

int
dp_rh_h(PsyState *psy)                /* ref:p28 */
{
   double p_w, p_ws, t_db, W;

   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   p_w = p_sat(psy->T_dew + 273.15);
   p_ws = 100*p_w / psy->RH;
   t_db = T_sat(p_ws) - 273.15;
   W = (psy->h - c_pa()*t_db) / (2501 + 1.805*t_db);

   psy->T_db = t_db;
   psy->W = W;
   psy->P = 0.62198*p_w / W + p_w;
   psy->T_wb = T_wet_bulb( psy->T_dew, t_db, psy->P);

   return NO_ERROR;
}

int
dp_h_W(PsyState *psy)              /* ref:p30 */
{
   double p_w, p_tot, t_db, p_ws;
   
   p_w = p_sat(psy->T_dew + 273.15);
   p_tot = 0.62198*p_w / psy->W + p_w;
   t_db = (psy->h - 2510*psy->W) / (c_pa() + 1.805*psy->W);
   p_ws = p_sat(t_db + 273.15);


   psy->P = p_tot;
   psy->T_db = t_db;
   psy->T_wb = T_wet_bulb( psy->T_dew, t_db, p_tot);
   psy->RH = 100*p_w/p_ws;

   return NO_ERROR;
}

int 
db_W_rh (PsyState *psy)                      /* ref:p4 */
{
   double p_ws, p_w, p_tot;
  
   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   p_w = psy->RH * p_ws / 100;
   p_tot = p_w + 0.62198 * p_w / psy->W ; 
    
   psy->T_dew = T_sat(p_w) - 273.15;
   psy->T_wb = T_wet_bulb( psy->T_dew, psy->T_db, p_tot);
   psy->P = p_tot;
   psy->h = h_W_db(psy->W, psy->T_db);

   return NO_ERROR;
}
   
int
dp_W_rh(PsyState *psy)                /* ref:p27 */
{
   double p_w, p_ws, t_db;

   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;


   p_w = p_sat(psy->T_dew + 273.15);
   p_ws = 100*p_w / psy->RH;
   t_db = T_sat(p_ws) - 273.15;

   psy->T_db = t_db;
   psy->P = 0.62198*p_w / psy->W + p_w;
   psy->h = c_pa()*t_db + psy->W*(2501 + 1.805*t_db);
   psy->T_wb = T_wet_bulb( psy->T_dew, t_db, psy->P);

   return NO_ERROR;
}

int 
db_rh_h (PsyState *psy)                      /* ref:p5 */
{
   double p_ws, p_w, W, p_tot;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->RH < min_RH || psy->RH > max_RH) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   p_w = psy->RH * p_ws / 100;
   W = (psy->h - c_pa()*psy->T_db) / (2510 + 1.805*psy->T_db);
   p_tot = 0.62198 * p_w / W + p_w; 
   
   psy->T_dew = T_sat(p_w) - 273.15;
   psy->T_wb = T_wet_bulb( psy->T_dew, psy->T_db, p_tot);
   psy->P = p_tot;
   psy->W = W;

   return NO_ERROR;
}

int
wb_W_h(PsyState *psy)                    /* ref:p22 */
{  
   double t_db, p_ws, Ws_wb, p_sat_wb, p_tot, p_w;

   t_db = (psy->h - 2501*psy->W) / (c_pa() + 1.805*psy->W);
   p_ws = p_sat(t_db + 273.15);
   Ws_wb = ((2501 + 1.805*t_db - 4.186*psy->T_wb)*psy->W 
                   + c_pa()*(t_db - psy->T_wb)) 
                   / (2501 - 2.381*psy->T_wb);
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   p_tot = 0.62198*p_sat_wb / Ws_wb + p_sat_wb;
   p_w = p_tot*psy->W / (psy->W + 0.62198);

   psy->T_db = t_db;
   psy->T_dew = T_sat(p_w) - 273.15;
   psy->RH = 100*p_w / p_ws;
   psy->P = p_tot;   

   return NO_ERROR;
}

int
db_dp_W(PsyState *psy)             /* ref:p7 */
{
   double p_ws, p_w, p_tot;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if( psy->T_dew > psy->T_db ) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   p_w = p_sat(psy->T_dew + 273.15);
   p_tot = 0.62198*p_w / psy->W + p_w;

   psy->T_wb = T_wet_bulb( psy->T_dew, psy->T_db, p_tot);
   psy->P = p_tot;
   psy->RH = 100 * p_w / p_ws;
   psy->h = h_W_db(psy->W, psy->T_db);

   return NO_ERROR;
}

int 
db_wb_h(PsyState *psy)                 /* ref:p15 */
{
   double p_ws, W, p_sat_wb, Ws_wb, p_tot, p_w;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if(psy->T_wb > psy->T_db ) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   W = (psy->h - c_pa()*psy->T_db) / (2501 + 1.805*psy->T_db);
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   Ws_wb = ((2501 + 1.805*psy->T_db - 4.186*psy->T_wb)*W 
             + c_pa()*(psy->T_db - psy->T_wb))
             / (2501 - 2.381*psy->T_wb);
   p_tot = 0.62198*p_sat_wb/Ws_wb + p_sat_wb;
   p_w = p_tot*W / (W + 0.62198);
   
   psy->T_dew = T_sat(p_w) - 273.15 ;
   psy->RH = 100*p_w / p_ws;
   psy->P = p_tot;
   psy->W = W;
 
   return NO_ERROR;
}

int
db_wb_W(PsyState *psy)             /* ref:p10 */
{
   double p_ws, p_sat_wb, W_s_wb, p_tot, p_w;

   if( psy->T_db < min_db  || psy->T_db > max_db ) return VAR_OUT_OF_RANGE;
   if(psy->T_wb > psy->T_db ) return VAR_OUT_OF_RANGE;

   p_ws = p_sat(psy->T_db + 273.15);
   p_sat_wb = p_sat(psy->T_wb + 273.15);
   W_s_wb = ((2501 + 1.805*psy->T_db - 4.186*psy->T_wb)*psy->W 
             + c_pa()*(psy->T_db - psy->T_wb))
             / (2501 - 2.381*psy->T_wb);
   p_tot = 0.62198*p_sat_wb / W_s_wb + p_sat_wb;
   p_w = p_tot*psy->W / (psy->W + 0.62198);

   psy->T_dew = T_sat(p_w) - 273.15;
   psy->P = p_tot;
   psy->RH = 100*p_w / p_ws;
   psy->h = h_W_db(psy->W, psy->T_db);

   return NO_ERROR;
}
