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
N = 100
point_type = 3
# split_lines = False
split_lines = True
ref = 12

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
    with open("naca2.geo","w") as file:
        file.write("c1 = 1.0;\n")
        file.write("c2 = 0.8;\n")
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
        # write splines
        file.write("\n// airfoil lines\n")
        if not split_lines:
            points_top = range(1,N+1)
            points_top = list(map(str,points_top))
            points_bot = [1]+list(range(N+1,2*N-1))+[N]
            points_bot = list(map(str,points_bot))
            file.write(f"Spline(1) = {{{','.join(points_top)}}};\n")
            file.write(f"Spline(2) = {{{','.join(points_bot)}}};\n")
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


def plot_airfoil(x,xb,y_top,y_bot,chord):
    fig, ax = plt.subplots()
    plt.plot(x,y_top,'kd')
    plt.plot(xb,y_bot,'kd')
    lim = 0.2*chord
    ax.set_ylim([-lim,lim])
    plt.grid()
    plt.show()

if __name__ == "__main__":
    # export_geo(x,xb,y_top,y_bot)
    plot_airfoil(x,xb,y_top,y_bot,chord)
