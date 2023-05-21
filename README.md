# building_kinetics_model
Automating the construction of kinetics models with symbolic python. 
Currently working for the following reaction scheme under continuous stirring reactor:

HCOOH(g) + * --> HCOOH* 
HCOOH* + * ---> HCOObi* + H*
HCOObi* --> HCOOmono* 
HCOOmono* + * --> CO2* + H*
2H* --> H2* + *
H2* --> H2(g) + *
CO2* --> CO2(g) + *

For continuous stirring reactor, the three reactions are:

HCOOH_in --> HCOOH(g)
CO2_in --> CO2(g)
H2_in --> H2(g)
