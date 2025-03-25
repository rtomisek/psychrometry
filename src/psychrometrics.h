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



#ifndef _PSYCHROMETRICS_H
#define _PSYCHROMETRICS_H

enum errors
{  NO_ERROR,
   ERROR,
   NO_ROOT,
   VAR_OUT_OF_RANGE,
   UNRESOLVABLE
};

/*   The structure PsyState contains seven psychrometric variables.
 *   They are:
 *     P      total pressure          Pascals
 *     T_db   dry bulb temperature    centigrade
 *     T_wb   wet bulb temperature    centigrade
 *     T_dew  dew point temperature   centigrade
 *     RH     relative humidity       %
 *     W      humidity ratio          KG/KG           
 *     h      enthalpy                kJoules/KG
 *
 *   All units should be in SI.
 */

/* NOTE: keep all temperatures below 100C, range check not yet implemented */
 
typedef struct {
   double T_db;
   double T_wb;
   double T_dew;
   double P;
   double RH;
   double W;
   double h;
} PsyState ;


/* cases where P is specified */
int P_db_dp( PsyState *);
int P_db_W ( PsyState *);
int P_db_wb( PsyState *);
int P_db_rh( PsyState *);
int P_db_h ( PsyState *);
int P_wb_dp( PsyState *);
int P_wb_W ( PsyState *);
int P_wb_rh( PsyState *);
int P_dp_rh( PsyState *);
int P_dp_h ( PsyState *);
int P_rh_W ( PsyState *);
int P_W_h  ( PsyState *);
int P_rh_h ( PsyState *);
/* cases where P is unknown */
int db_wb_h( PsyState *);
int W_rh_h ( PsyState *);
int db_h_dp( PsyState *);
int wb_dp_W( PsyState *);
int dp_rh_h( PsyState *);
int dp_h_W ( PsyState *); 
int db_W_rh( PsyState *);
int dp_W_rh( PsyState *);
int db_rh_h( PsyState *);
int wb_W_h ( PsyState *);
int db_dp_W( PsyState *);
int db_wb_W( PsyState *);
/*   as yet not implemented
int db_wb_rh(PsyState *);  
int wb_dp_rh(PsyState *);  
int db_wb_dp(PsyState *);  
int dp_w_rh(PsyState *);
int wb_w_rh(PsyState *); 
int wb_rh_h(PsyState *);
int wb_dp_h (PsyState *); */
    
/*
   Cannot be Implemented:
 
   int P_wb_h ( PsyState *);   
     this one has problems with accuracy because the wet bulb and 
     enthapy lines are nearly parallel. It is effectivley unresolvable.
     
   int P_dp_W( PsyState *);
      this is unresolvable because dew point and W basically redundant;
      they both specify the same horizontal line on the pchart.
      
   int db_W_h (PsyState *);   
       unresolvable because of redundacy.  Any two of the three
       variables determine the third one uniquely.
     
   int db_dp_rh (PsyState *);
        same problem as db_w_h above.
        
*/

#endif /* _PSYCHROMETRICS_H */
