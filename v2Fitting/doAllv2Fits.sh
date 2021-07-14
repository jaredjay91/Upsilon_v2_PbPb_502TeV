#!/bin/bash
<<COMMENT
root -l -b -q "GetAllYieldsvsPhi.C(1,0)"
#root -l -b -q "GetAllYieldsvsPhi.C(1,1)"
#root -l -b -q "GetAllYieldsvsPhi.C(1,2)"
root -l -b -q "GetAllYieldsvsPhi.C(1,3)"
root -l -b -q "GetAllYieldsvsPhi.C(1,4)"
#root -l -b -q "GetAllYieldsvsPhi.C(1,5)"
COMMENT
root -l -b -q "Get_v2_vs_var.C(1,0)"
#root -l -b -q "Get_v2_vs_var.C(1,1)"
#root -l -b -q "Get_v2_vs_var.C(1,2)"
root -l -b -q "Get_v2_vs_var.C(1,3)"
root -l -b -q "Get_v2_vs_var.C(1,4)"
#root -l -b -q "Get_v2_vs_var.C(1,5)"
<<COMMENT
root -l -b -q "GetAllYieldsvsPhi.C(2,0)"
#root -l -b -q "GetAllYieldsvsPhi.C(2,1)"
#root -l -b -q "GetAllYieldsvsPhi.C(2,2)"
root -l -b -q "GetAllYieldsvsPhi.C(2,3)"
root -l -b -q "GetAllYieldsvsPhi.C(2,4)"
#root -l -b -q "GetAllYieldsvsPhi.C(2,5)"
COMMENT
root -l -b -q "Get_v2_vs_var.C(2,0)"
#root -l -b -q "Get_v2_vs_var.C(2,1)"
#root -l -b -q "Get_v2_vs_var.C(2,2)"
root -l -b -q "Get_v2_vs_var.C(2,3)"
root -l -b -q "Get_v2_vs_var.C(2,4)"
#root -l -b -q "Get_v2_vs_var.C(2,5)"

