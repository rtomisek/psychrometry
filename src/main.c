
#include <stdio.h>

#include "psychrometrics.h"
#include "svp.h"


void print_state( PsyState ps)
{
    printf("Pressure   %f\n", ps.P);
    printf("Dry bulb   %f\n", ps.T_db);
    printf("Wet bulb   %f\n", ps.T_wb);
    printf("Dew point  %f\n", ps.T_dew);
    printf("RH         %f\n", ps.RH);
    printf("Ratio      %f\n", ps.W);
    printf("Enthalpy   %f\n", ps.h);
}

int
main(void)
{
    PsyState ps;
  
    printf("by IAPWS - %f\n", p_IAPWS(350));
    printf("by ASHRE - %f\n", p_ASH(350));    
   
    return 0;
}
