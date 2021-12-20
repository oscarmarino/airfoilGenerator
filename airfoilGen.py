#!/usr/bin/python3
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import os


def naca_profile_4_digits(x,t,blunt):
    if blunt:
        a4 = -0.1015
    else:
        a4 = -0.1036
    y = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 +0.2843*x**3 + a4*x**4)
    return y

def read_naca456():
    df_temp = pd.read_csv('naca.out', skiprows=7, sep='\s+', skipfooter=1)
    return np.array(df_temp['y'])

def create_naca456_input_file(x,t):
    with open('naca456/myNaca63016.nml',"w") as file:
        file.write("&NACA\n")
        file.write("  name    = 'NACA 63-018',\n")
        file.write("  profile = '63',\n")
        file.write(f"  toc     = {t},  \n")
        file.write("  camber  = '0'\n")
        if len(x)<90:
            file.write("  denCode = 3/\n")
        else:
            file.write("  denCode = 0\n")
            xx = list(map(str,x))
            file.write(f"  xTable  = {','.join(xx)}/\n")

def naca_profile_6_digits(x,t,blunt):
    if blunt:
        # not changes needed here, it'll be cut and rescaled later
        pass
    create_naca456_input_file(x,t)
    os.system('naca456/naca456 > /dev/null')
    y = read_naca456()
    os.system('rm ./naca.gnu ./naca.dbg ./naca.out')
    return y

def rescale_by_chord(x,xb,y_top,y_bot,chord,resc_y=True):
    x = x * chord
    xb = xb * chord
    if resc_y:
        chord_y = chord
    else:
        chord_y = 1
    y_top = y_top * chord_y
    y_bot = y_bot * chord_y
    return x,xb,y_top,y_bot

def naca_profile(x,t,type,blunt):
    if type == 4:
        y = naca_profile_4_digits(x,t,blunt)
    elif type == 6:
        y = naca_profile_6_digits(x,t,blunt)
    else:
        print("using 4 digit naca es default")
        y = naca_profile_4_digits(x,t,blunt)

    return y

def get_straight_te(x,xb,y_top,y_bot,te_w,chord):
    N = len(x)
    star_index = int(N/2)
    for i,(yt,yb) in enumerate(zip(y_top[star_index:],y_bot[star_index:])):
        width = np.abs(yt - yb)
        if width <= te_w:
            break
    i = i + star_index
    print(f"The geometry was cut and rescaled to have a straight TE, it has now {i} points")

    symmetric = y_top[i] == -y_bot[i]
    if not symmetric:
        raise Exception("Straight TE rescale only implemented for simmetric airfoils")

    x_new = x[:i]
    xb_new = xb[:i]
    yt_new = y_top[:i]
    yb_new = y_bot[:i]
    yt_new[-1] = te_w / 2.0
    yb_new[-1] = -te_w / 2.0

    def x_linear_inter(xv,yv,y):
        xx = (y - yv[0]) * (xv[1]-xv[0]) / (yv[1]-yv[0]) + xv[0]
        return xx

    if width == te_w:
        final_length = x[i-1]
    else:
        final_length = x_linear_inter(x[i:i+2],y_top[i:i+2],yt_new[-1])

    x_new[-1] = final_length
    xb_new[-1] = final_length

    partial_chord = chord / final_length
    x,xb,y_top,y_bot = rescale_by_chord(x_new,xb_new,yt_new,yb_new,partial_chord,resc_y=False)

    return x, xb, y_top, y_bot

def get_airfoil_points(**kwargs):

    # get naca 0012 y default with 100 points ref at le and te
    point_type = kwargs.get('point_type',3)
    N = kwargs.get('N',100)
    ref = kwargs.get('ref',12)
    type = kwargs.get('type',4)
    blunt = kwargs.get('blunt',False)
    te_w = kwargs.get('te_w',None)
    chord = kwargs.get('chord',1.0)

    if point_type == 1:
        #equally space
        t = np.linspace(0,1.0,N)
        x = t
        xb = x
        y_top = naca_profile(t,ref/100.0,type,blunt)
        y_bot = -1.0 * y_top
        x,xb,y_top,y_bot = rescale_by_chord(x,xb,y_top,y_bot,chord)

    elif point_type ==2:
        #refined le
        t = np.linspace(0,1.0,N)
        x = ((1 - t) ** 2)
        xb = t ** 2
        y_top = naca_profile((1-t)**2,ref/100.0,type,blunt)
        y_bot = -1.0 * naca_profile(t**2,ref/100.0,type,blunt)
        index_sort = x.argsort()
        x = x[index_sort]
        y_top = y_top[index_sort]
        x,xb,y_top,y_bot = rescale_by_chord(x,xb,y_top,y_bot,chord)

    else:
        # refine in each edge
        t = np.linspace(0,np.pi,N)
        x = 0.5 * (1 - np.cos(t))
        y_gen = naca_profile(x,ref/100.0,type,blunt)
        xb = x
        y_ini = np.zeros_like(y_gen)
        y_top = y_ini + y_gen
        y_bot = y_ini - y_gen
        x,xb,y_top,y_bot = rescale_by_chord(x,xb,y_top,y_bot,chord)

    if te_w is not None:
        x,xb,y_top,y_bot = get_straight_te(x,xb,y_top,y_bot,te_w,chord)

    return x, xb, y_top, y_bot

def export_geo(x,xb,y_top,y_bot,file_name,**kwargs):
    blunt = kwargs.get('blunt',False)
    create_mesh = kwargs.get('create_mesh',False)
    split_lines = kwargs.get('split_lines',False)
    create_zones = kwargs.get('create_zones',False)
    N = len(x)
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
        if blunt:
            bot_final_point = len(xb)
        else:
            bot_final_point = -1
        for xx,yy in zip(xb[1:bot_final_point],y_bot[1:bot_final_point]):
            i+=1
            file.write(f"Point({i}) = {{{xx}, {yy}, 0.0, c1}};\n")
        final_point = i

        # write splines
        file.write("\n// airfoil lines\n")
        if blunt:
            bot_final_point = final_point
        else:
            bot_final_point = N
        if not split_lines:
            points_top = range(1,N+1)
            points_top = list(map(str,points_top))
            points_bot = [1]+list(range(N+1,2*N-1))+[bot_final_point]
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
            points_bot = list(range(second_point_bot,2*N-1)) + [bot_final_point]
            points_bot = list(map(str,points_bot))
            file.write(f"Spline(5) = {{{','.join(points_top)}}};\n")
            file.write(f"Spline(6) = {{{','.join(points_bot)}}};\n")

            final_line = 6

            if blunt:
                last_points = [N,bot_final_point]
                last_points = list(map(str,last_points))
                file.write(f"Line(7) = {{{','.join(last_points)}}};\n")
                final_line += 1

    if create_mesh:
        if not split_lines:
            raise Exception("Mesh can only be created with spline of airfoil in parts")
        else:
            y_position_last_point = y_top[-1]
            export_mesh(N, final_point, final_line, first_point_top, second_point_top, first_point_bot, second_point_bot, create_zones, y_position_last_point, blunt)

def export_mesh(N, last_point, last_line, ft, st, fb, sb, create_zones, yl, blunt):
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
        file.write(f"Point({j}) = {{c3, {yl}, 0.0, c1 }};\n")
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
        if blunt:
            file.write(f"Point({j}) = {{c3, -{yl}, 0.0, c1 }};\n")
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
        if blunt:
            file.write(f"Line({k}) = {{{last_point+7}, {last_point+15}}};\n")
        else:
            file.write(f"Line({k}) = {{{last_point+7}, {last_point+8}}};\n")
        k += 1
        file.write(f"Line({k}) = {{1, {last_point+3}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{N}, {last_point+8}}};\n")
        k += 1
        file.write(f"Line({k}) = {{{last_point+4}, {N}}};\n")
        k += 1
        if blunt:
            file.write(f"Line({k}) = {{{last_point+6}, {last_point}}};\n")
        else:
            file.write(f"Line({k}) = {{{last_point+6}, {N}}};\n")
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
        if blunt:
            file.write("//+\n")
            file.write(f"Line({k}) = {{{last_point+8}, {last_point+15}}};\n")
            k += 1
            file.write(f"Line({k}) = {{{last_point}, {last_point+15}}};\n")
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
        if blunt:
            file.write(f"Curve Loop(8) = {{{last_line+14}, {last_line+22}, -{last_line+10}, -{last_line+9}}};\n")
            file.write(f"Plane Surface(8) = {{8}};\n")
            file.write(f"Curve Loop(9) = {{{last_line}, {last_line+22}, -{last_line+21}, -{last_line+12}}};\n")
            file.write(f"Plane Surface(9) = {{9}};\n")
        else:
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
        if blunt:
            file.write(f"Transfinite Surface {{9}};\n")

        file.write("\n// Curves divisions, change to custom\n")
        file.write(f"Transfinite Curve {{{last_line+13}, -{last_line+11}, {last_line+7}, {last_line+14}, {last_line+10}, {last_line+15}, {last_line+16}, {last_line+18}, {last_line+19}}} = 50 Using Progression 0.865;\n")
        file.write("//+\n")
        file.write(f"Transfinite Curve {{1, 2, {last_line+3}, {last_line+4}}} = 11 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{{last_line+1},{last_line+2}}} = 95 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{{last_line+17}, {last_line+20}}} = 46 Using Progression 1;\n")
        file.write(f"Transfinite Curve {{-3,-4}} = 140 Using Progression 0.998;\n")
        file.write(f"Transfinite Curve {{{last_line+8}, {last_line+5}, 5, 6}} = 90 Using Progression 1;\n")
        if blunt:
            file.write(f"Transfinite Curve {{-{last_line+12}, -{last_line+6}, -{last_line+9}, -{last_line+22}}} = 60 Using Progression 0.91;\n")
            file.write(f"Transfinite Curve {{{last_line}, {last_line+21}}} = 5 Using Progression 1;\n")
        else:
            file.write(f"Transfinite Curve {{-{last_line+12}, -{last_line+6}, -{last_line+9}}} = 60 Using Progression 0.91;\n")

        file.write("//+\n")
        if blunt:
            file.write(f"Recombine Surface {{1, 2, 3, 4, 5, 6, 7, 8, 9}};\n")
        else:
            file.write(f"Recombine Surface {{1, 2, 3, 4, 5, 6, 7, 8}};\n")
        file.write("//+\n")
    # zones may change!!!!
    if create_zones:
        with open(file_name,"a") as file:
            if blunt:
                #extrusion
                file.write("Extrude {0, 0, 0.1} {\n")
                file.write("  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Layers{20}; Recombine;\n")
                file.write("}\n")
                file.write("//+\n")
                #surface zones
                file.write('Physical Surface("airfoil") = {50, 65, 91, 184, 158, 143};\n')
                file.write('Physical Surface("airfoil_te") = {224};\n')
                file.write('Physical Surface("outflow") = {42, 77, 73, 99, 121, 117, 232, 210, 214, 192, 166, 170, 135};\n')
                file.write('Physical Surface("sidea") = {3, 2, 1, 5, 6, 7, 9, 4, 8};\n')
                file.write('Physical Surface("sideb") = {51, 78, 100, 122, 237, 215, 193, 171, 144};\n')
                file.write("//+\n")
                #volume zone
                file.write('Physical Volume("fluid") = {1, 2, 3, 4, 5, 6, 7, 8, 9};\n')
            else:
                #extrusion
                file.write("Extrude {0, 0, 0.1} {\n")
                file.write("  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Layers{20}; Recombine;\n")
                file.write("}\n")
                file.write("//+\n")
                #surface zones
                file.write('Physical Surface("airfoil") = {140, 155, 181, 88, 62, 47};\n')
                file.write('Physical Surface("outflow") = {39, 74, 70, 96, 118, 114, 207, 211, 189, 163, 167, 132};\n')
                file.write('Physical Surface("sidea") = {4, 3, 2, 1, 5, 6, 7, 8};\n')
                file.write('Physical Surface("sideb") = {48, 75, 97, 119, 212, 190, 168, 141};\n')
                file.write("//+\n")
                #volume zone
                file.write('Physical Volume("fluid") = {1, 2, 3, 4, 5, 6, 7, 8};\n')

def export_txt(x,xb,y_top,y_bot,file_name,le_start=True):
    if le_start:
        k = 0
    with open(file_name,"w") as file:
        file.write(f"x          y\n")
        for a,b in zip(x,y_top):
            file.write(f"{a:6.4e} {b:6.4e}\n")
        for a,b in zip(xb,y_bot):
            file.write(f"{a:6.4e} {b:2.6e}\n")


def plot_airfoil(x,xb,y_top,y_bot,**kwargs):
    chord = kwargs.get('chord',1.0)
    prop = kwargs.get('prop','kd')
    fig, ax = plt.subplots()
    plt.plot(x,y_top,prop)
    plt.plot(xb,y_bot,prop)
    lim = 0.2*chord
    ax.set_ylim([-lim,lim])
    plt.grid()

if __name__ == "__main__":
    # chord = 0.2
    # N = 100
    N = 500
    ref = 18
    file_name = "naca2.geo"
    type = 6
    blunt = True
    trailing_edge_width = 0.001 / 0.2 # 1 mm scaled by chord
    # x, xb, y_top, y_bot = get_airfoil_points(N=N)
    # x, xb, y_top, y_bot = get_airfoil_points(N=N,blunt=blunt)
    # x, xb, y_top, y_bot = get_airfoil_points(N=N,te_w=trailing_edge_width)
    x, xb, y_top, y_bot = get_airfoil_points(N=N,blunt=True, te_w=trailing_edge_width)
    plot_airfoil(x,xb,y_top,y_bot,prop='k-')
    export_geo(x,xb,y_top,y_bot,file_name,blunt=blunt,create_mesh=True,split_lines=True,create_zones=True)
    # fn = "naca5str.txt"
    # # fn = "naca2b.txt"
    # export_txt(x,xb,y_top,y_bot,fn)

    trailing_edge_width = 0.31 / 1000 / 0.2 # 0.31 mm scaled by chord
    file_name = "naca63.geo"
    x, xb, y_top, y_bot = get_airfoil_points(N=N, ref=ref, type=type, blunt=blunt,te_w=trailing_edge_width)

    plot_airfoil(x,xb,y_top,y_bot,prop='k-')
    # plot_airfoil(x,xb,y_top,y_bot)
    # export_geo(x,xb,y_top,y_bot,file_name)
    # export_geo(x,xb,y_top,y_bot,file_name,blunt=blunt)
    # export_geo(x,xb,y_top,y_bot,file_name,create_mesh=True,split_lines=True,create_zones=True)
    export_geo(x,xb,y_top,y_bot,file_name,blunt=blunt,create_mesh=True,split_lines=True,create_zones=True)

    plt.show()
