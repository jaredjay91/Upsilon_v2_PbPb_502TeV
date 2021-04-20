#!/bin/bash

root -b -l -q "draw_v2.C(1)"
root -b -l -q "draw_v2.C(2)"
root -b -l -q draw_v2_ratios.C
root -b -l -q draw_v2_theory_comp.C
