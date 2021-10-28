#!/usr/bin/python3
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


def naca_profile(x,t):
    # for now only for sharp TE, change when need a blunt TE
    y = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 +0.2843*x**3 -0.1036*x**4)
    return y

def rescale_by_chord(x,xb,y_top,y_bot,chord):
    x = x * chord
    xb = xb * chord
    y_top = y_top * chord
    y_bot = y_bot * chord
    return x,xb,y_top,y_bot

chord = 1.0
# chord = 0.2
# N = 100
N = 500
point_type = 3
# split_lines = False
split_lines = True
ref = 12
create_mesh = True
file_name = "naca2.geo"

if point_type == 1:
    #equally space
    t = np.linspace(0,1.0,N)
    x = t
    xb = x
    y_top = naca_profile(t,ref/100.0)
    y_bot = -1.0 * y_top
    x,xb,y_top,y_bot = rescale_by_chord(x,xb,y_top,y_bot,chord)

elif point_type ==2:
    #refined le
    t = np.linspace(0,1.0,N)
    x = ((1 - t) ** 2)
    xb = t ** 2
    y_top = naca_profile((1-t)**2,ref/100.0)
    y_bot = -1.0 * naca_profile(t**2,ref/100.0)
    x,xb,y_top,y_bot = rescale_by_chord(x,xb,y_top,y_bot,chord)

else:
    # refine in each edge
    t = np.linspace(0,np.pi,N)
    x = 0.5 * (1 - np.cos(t))
    y_gen = naca_profile(x,ref/100.0)
    xb = x
    y_ini = np.zeros_like(y_gen)
    y_top = y_ini + y_gen
    y_bot = y_ini - y_gen
    x,xb,y_top,y_bot = rescale_by_chord(x,xb,y_top,y_bot,chord)

def export_geo(x,xb,y_top,y_bot):
    # with open("naca.geo","w") as file:
    with open(file_name,"w") as file:
        file.write("c1 = 1.0;\n")
        file.write("c2 = 10.0;\n")
        file.write("c3 = 12.0;\n")
        i = 0
        # write points
        file.write("\n// airfoil top\n")
        for xx,yy in zip(x,y_top):
            i+=1
            file.write(f"Point({i}) = {{{xx}, {yy}, 0.0, c1}};\n")
        file.write("\n// airfoil bottom\n")
        # for now only for sharp TE, change when need a blunt TE
        for xx,yy in zip(xb[1:-1],y_bot[1:-1]):
            i+=1
            file.write(f"Point({i}) = {{{xx}, {yy}, 0.0, c1}};\n")
        final_point = i

        # write splines
        file.write("\n// airfoil lines\n")
        if not split_lines:
            points_top = range(1,N+1)
            points_top = list(map(str,points_top))
            points_bot = [1]+list(range(N+1,2*N-1))+[N]
            points_bot = list(map(str,points_bot))
            file.write(f"Spline(1) = {{{','.join(points_top)}}};\n")
            file.write(f"Spline(2) = {{{','.join(points_bot)}}};\n")

            final_line = 2
        else:
            first_point_top = np.where(x>0.02)[0][0] + 1 # add 1 since array starts at 0 and Points numbered from 1
            first_point_bot = first_point_top + N -1
            second_point_top = np.where(x>0.65)[0][0] + 1 # add 1 since array starts at 0 and Points numbered from 1
            second_point_bot = second_point_top + N -1

            points_top = range(1,first_point_top+1)
            points_top = list(map(str,points_top))
            points_bot = [1]+list(range(N+1,first_point_bot+1))
            points_bot = list(map(str,points_bot))
            file.write(f"Spline(1) = {{{','.join(points_top)}}};\n")
            file.write(f"Spline(2) = {{{','.join(points_bot)}}};\n")

            points_top = range(first_point_top,second_point_top+1)
            points_top = list(map(str,points_top))
            points_bot = range(first_point_bot,second_point_bot+1)
            points_bot = list(map(str,points_bot))
            file.write(f"Spline(3) = {{{','.join(points_top)}}};\n")
            file.write(f"Spline(4) = {{{','.join(points_bot)}}};\n")

            points_top = range(second_point_top,N+1)
            points_top = list(map(str,points_top))
            points_bot = list(range(second_point_bot,2*N-1)) + [N]
            points_bot = list(map(str,points_bot))
            file.write(f"Spline(5) = {{{','.join(points_top)}}};\n")
            file.write(f"Spline(6) = {{{','.join(points_bot)}}};\n")

            final_line = 6

    if create_mesh:
        export_mesh(final_point, final_line, first_point_top, second_point_top, first_point_bot, second_point_bot)

def export_mesh(last_point, last_line, ft, st, fb, sb):
    with open(file_name,"a") as file:
        file.write("\n// points for arc of cmesh\n")
        j = last_point + 1
        file.write(f"Point({j}) = {{c2*Cos(Pi*100/180), c2*Sin(Pi*100/180), 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{c2*Cos(Pi*100/180), -c2*Sin(Pi*100/180), 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{-c2, 0.0, 0.0, c1 }};\n")
        j += 1

        file.write("// points for rest of cmesh\n")
        distance_a = 2.5
        distance_b = 1.0
        distance_c = 0.4
        factor_a = 1.015
        factor_b = 1.01
        file.write(f"Point({j}) = {{{distance_a}, c2, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{c3, c2, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{{distance_a}, -c2, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{c3, -c2, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{c3, 0.0, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{-c2*Cos(Pi/4), c2*Sin(Pi/4), 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{{distance_b}, c2*{factor_a}, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{-{distance_c}, c2*{factor_b}, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{-c2*Cos(Pi/4), -c2*Sin(Pi/4), 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{{distance_b}, -c2*{factor_a}, 0.0, c1 }};\n")
        j += 1
        file.write(f"Point({j}) = {{-{distance_c}, -c2*{factor_b}, 0.0, c1 }};\n")
        j += 1
        file.write("//+\n")

        k = last_line + 1
        file.write(f"Circle({k}) = {{{last_point+9}, 1, {last_point+1}}};\n")
        k += 1
        file.write(f"Circle({k}) = {{{last_point+12}, 1, {last_point+2}}};\n")
        k += 1
        file.write(f"Circle({k}) = {{{last_point+3}, 1, {last_point+9}}};\n")
        k += 1
        file.write(f"Circle({k}) = {{{last_point+3}, 1, {last_point+12}}};\n")
        k += 1
        file.write("//+\n")
        file.write(f"Line({k}) = {{{last_point+10}, {last_point+4}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+4}, {last_point+5}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+5}, {last_point+8}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+13}, {last_point+6}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+6}, {last_point+7}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+7}, {last_point+8}}};\n")
        k += 1
        file.write(f"Line({k}) = {{1, {last_point+3}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{N}, {last_point+8}}};\n")
        # file.write(f"Line({k}) = {{{last_point}, {last_point+8}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+4}, {N}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+6}, {N}}};\n")
        # file.write(f"Line({k}) = {{{last_point+6}, 1, {last_point}}};\n")
        k += 1
        file.write("//+\n")
        file.write(f"Line({k}) = {{{last_point+9}, {ft}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+10}, {st}}};\n")
        k += 1
        file.write(f"Spline({k}) = {{{last_point+1}, {last_point+11}, {last_point+10}}};\n")
        k += 1
        file.write("//+\n")
        file.write(f"Line({k}) = {{{last_point+12}, {fb}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+13}, {sb}}};\n")
        k += 1
        file.write(f"Spline({k}) = {{{last_point+2}, {last_point+14}, {last_point+13}}};\n")
        k += 1
        file.write("\n// Surfaces\n")
        file.write(f"Curve Loop(1) = {{{last_line+11}, {last_line+3}, {last_line+15}, -1}};\n")
        file.write(f"Plane Surface(1) = {{1}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(2) = {{{last_line+15}, 3, -{last_line+16}, -{last_line+17}, -{last_line+1}}};\n")
        file.write(f"Plane Surface(2) = {{2}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(3) = {{{last_line+16}, 5, -{last_line+13}, -{last_line+5}}};\n")
        file.write(f"Plane Surface(3) = {{3}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(4) = {{{last_line+13}, {last_line+12}, -{last_line+7}, -{last_line+6}}};\n")
        file.write(f"Plane Surface(4) = {{4}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(5) = {{{last_line+11}, {last_line+4}, {last_line+18}, -2}};\n")
        file.write(f"Plane Surface(5) = {{5}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(6) = {{{last_line+18}, 4, -{last_line+19}, -{last_line+20}, -{last_line+2}}};\n")
        file.write(f"Plane Surface(6) = {{6}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(7) = {{{last_line+19}, 6, -{last_line+14}, -{last_line+8}}};\n")
        file.write(f"Plane Surface(7) = {{7}};\n")
        file.write("//+\n")
        file.write(f"Curve Loop(8) = {{{last_line+14}, {last_line+12}, -{last_line+10}, -{last_line+9}}};\n")
        file.write(f"Plane Surface(8) = {{8}};\n")
        file.write("//+\n")
        file.write(f"Transfinite Surface {{1}};\n")
        file.write(f"Transfinite Surface {{2}} = {{{ft},{st},{last_point+10},{last_point+9}}};\n")
        file.write(f"Transfinite Surface {{3}};\n")
        file.write(f"Transfinite Surface {{4}};\n")
        file.write(f"Transfinite Surface {{5}};\n")
        file.write(f"Transfinite Surface {{6}} = {{{fb},{sb},{last_point+13},{last_point+12}}};\n")
        file.write(f"Transfinite Surface {{7}};\n")
        file.write(f"Transfinite Surface {{8}};\n")

        file.write("\n// Curves divisions, change to custom\n")
        file.write(f"Transfinite Curve {{{last_line+13}, -{last_line+11}, {last_line+7}, {last_line+14}, {last_line+10}, {last_line+15}, {last_line+16}, {last_line+18}, {last_line+19}}} = 50 Using Progression 0.865;\n")
        file.write("//+\n")
        file.write(f"Transfinite Curve {{1, 2, {last_line+3}, {last_line+4}}} = 11 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{{last_line+1},{last_line+2}}} = 95 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{{last_line+17}, {last_line+20}}} = 46 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{-3,-4}} = 140 Using Progression 0.998;\n")
        file.write(f"Transfinite Curve {{{last_line+8}, {last_line+5}, 5, 6}} = 90 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{-{last_line+12}, -{last_line+6}, -{last_line+9}}} = 60 Using Progression 0.91;\n")

        file.write("//+\n")
        file.write(f"Recombine Surface {{1, 2, 3, 4, 5, 6, 7, 8}};\n")
        file.write("//+\n")

def plot_airfoil(x,xb,y_top,y_bot,chord):
    fig, ax = plt.subplots()
    plt.plot(x,y_top,'kd')
    plt.plot(xb,y_bot,'kd')
    lim = 0.2*chord
    ax.set_ylim([-lim,lim])
    plt.grid()
    plt.show()

if __name__ == "__main__":
    export_geo(x,xb,y_top,y_bot)
    # plot_airfoil(x,xb,y_top,y_bot,chord)
