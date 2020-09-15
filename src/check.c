
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
    float dbt;

    ps.P = 101325;  /* 1 atm in pascals */
    ps.RH = 80;

    for(dbt = 20 ; dbt < 120 ; dbt = dbt+20){
       ps.T_db = dbt;
       P_db_rh( &ps );

       print_state( ps );    }

    return 0;
}
