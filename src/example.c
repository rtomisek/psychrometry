
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
  
    /* Example, figure dew point given dry bulb and RH
     * We need to specify three variables so pressure is also given */

    ps.P = 101325;  /* 1 atm in pascals */
    ps.T_db = 27;   /* dry bulb should be below 100C */
    ps.RH = 80;

    P_db_rh( &ps );  /* note: the fuction argument is a pointer to PsyState */

    printf("The dew point is %f\n\n", ps.T_dew);
    printf("or print the whole thing\n");
    print_state( ps );    

    return 0;
}
