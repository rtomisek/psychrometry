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
#include "svp.h"

/* The saturated vapor partial pressure of water is a fundamental part
   of all psychrometric calcualtions. Several functions are available in
   the literature and some are given here.   
   
   All these functions have been adjusted so that temperature T is 
   in degrees Kelvin and vapor pressure p_sat is in Pascals
*/

double p_sat( double T )   /* selects the equation to be used in following calculations */
{                          /* TO DO: this selection should be in software NOT compiled in */
    return p_ASH( T );     
}

/* This function was given by ASHRAE. 
   
   Range given as 213.15 to 473.15  K ;
   agrees with NIST data within 0.02% over this range */
double p_ASH(double T)
{
   double a, A, B, C, D;

   if( T < 213.15 )
   {
      return -1;
    } else if( T < 273.15 )
   {
      A = -0.7297593707E-5;
      B =  0.5397420727E-2;
      C =  0.2069880620E2;
      D = -0.6042275128E4;
   } else if( T < 322.15 )
   {
      A =  0.1255001965E-4;
      B = -0.1923595289E-1;
      C =  0.2705101899E2;
      D = -0.6344011577E4;
   } else if( T < 373.15 )
   {
      A =  0.1246732157E-4;    
      B = -0.1915465806E-1;
      C =  0.2702388315E2;
      D = -0.6340941639E4; 
      
   } else if( T < 423.15 )
   {
      A =  0.1204507646E-4;
      B = -0.1866650553E-1;
      C =  0.2683629403E2;
      D = -0.6316972063E4;
   }  else if( T <= 473.15 )
   {
      A =  0.1069730183E-4;
      B = -0.1698965754E-1;
      C =  0.2614073298E2;
      D = -0.6220781230E4;
   } else
   {
      return -1;
   }
    
   a = ( (A*T) + B )*T + C + D/T;
   
   return 1000 * exp( a );
}


/* Equation from IAPWS
   Range given 273.15 to 646.15 K
 */   
   
double p_IAPWS(double T)
{	    
        double C[] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502 };
	double v, A, P_ws;

	v = 1 - T/647.096;
        A = (647.096/T)*(C[0]*v + C[1]*pow(v, 1.5) + C[2]*pow(v, 3) + C[3]*pow(v, 3.5) + C[4]*pow(v, 4) + C[5]*pow(v, 7.5))  ;
        P_ws = 220640 *exp(A);

	return P_ws*100;
}


/*  equation from Buck     
    ??? significant errors ???
    range unknown   */
double p_BUCK(double T)
{
    double kPa;
    double T_C;
    
    T_C = T - 273.15;  /*given equation is for degrees C */
    
    if( T_C > 0.0)
    {
        kPa = 0.61121 * exp((18.678 - T_C/234.5) * (T_C/(257.14 + T_C)));
    } else if( T_C <= 0.0)
    {
        kPa = 0.61115 * exp((23.036 - T_C/333.7) * (T_C/(279.82 + T_C)));
    }
    return (kPa * 1000) ;   /*return as Pa*/
}

/*   equation from Sonntag  
 
     range unknown, but compared to ASHRAE
     0.1 % difference up to 424
     0.5 % difference up to 464K    */
double p_SONN(double T)
{
    double res;
   
    if( T >= 273.15 )
    {
        res = exp( -6.0969385E3/T + 2.1240964E1 - 2.711193E-2*T + 1.673952E-5*pow(T, 2) + 2.433502*log(T) );
    }else if( T < 273.15 )
    {
        res = exp( -6.0245282E3/T + 2.932707E1 + 1.0613868E-2*T - 1.3198825E-5*pow(T, 2) - 4.9382577e-1*log(T) );
    }
 
    return res;   
}


/* equation from Teten
 
   range 252 to 364 K for less than 0.5 difference from ASHRAE
*/
double p_TET(double T)
{
    double kPa;
    double T_C;
    
    T_C = T - 273.15;
    
    if( T_C > 0.0 )
    {
        kPa = 0.61078 * exp( 17.27*T_C/(T_C + 237.3) );
    }
    else if( T_C < 0.0 )
    {
        kPa = 0.61078 * exp( 21.875*T_C/(T_C + 265.5) );
    }
    
    return kPa*1000;
}
    
    
    
/* --------------------------------------------------------------- */

/* calculate the water vapor saturation temperature for a given partial 
   pressure , this is the reverse of the above equation and is also
   from ASHRAE. */
/* !!!!!!!!!does not work !!!!!!*/
double 
T_ASH(double p_wsat)
{
    double b;
    double E, F, G, H, K;
    
    if( p_wsat < 1 )
    {
        return 0;
    } else if( p_wsat < 611 )
    {
        E = 0.1004926534E-2;
        F = 0.1392917633E-2;
        G = 0.2815151574;
        H = 0.7311621119E1;
        K = 0.125893734E3;
    } else if( p_wsat < 12350 )
    {
        E = 0.5031062503E-2;
        F = -0.8826779380E-1;
        G = 0.1243688446E1;
        H = 0.3388534296E1;
        K = 0.2150077993E3;
    } else if( p_wsat < 101420 )
    {
        E = 0.1209512517E-4;
        F = -0.3545542105;
        G = 0.2125893734E3;
        H = -0.2050301050E2;
        K = 0.2718585432E3;
    } else if( p_wsat < 476207 )
    {
         E = 0.2467291016E-1;
         F = -0.9367112883;
         G = 0.1514142334E2;
         H = -0.9882417501E2;
         K = 0.4995092948E3;
    } else if( p_wsat < 1555099)
    {
          E = 0.2748402484E-4;
          F = -0.1068661307E1;
          G = 0.1742964962E2;
          H = -0.1161208532E3;
          K = 0.5472618120E3;
    } else
    {
          return 0;
    }
        
    b = log( p_wsat );
    
    return (   (((E*b + F)*b + G)*b + H)*b + K );
    
}


/* or we can just invert p_sat numerially,
 * using binary search */
/* pressure in Pascals 
 * T_sat in Kelvins */
double 
T_sat(double p_wsat)
{
   unsigned int i=0;
   double tol=0.001;
   double t_m, t_1, t_2, f_m, f_1, f_2;

   t_1 = 213.15;
   f_1 = p_sat(t_1) - p_wsat;
   t_2 = 473.15;
   // f_2 = p_sat(t_2) - p_wsat;

   do{
      t_m =  (t_2 + t_1) / 2;
      f_m = p_sat(t_m) - p_wsat;
      if( f_m*f_1 < 0 ){
         t_2 = t_m;
         f_2 = f_m;
      }
      else{
         t_1 = t_m;
         f_1 = f_m;
      }
      i = i + 1;
    }while( fabs( (t_1-t_2)/2 ) > tol && i < 100);

   return t_m;
}

