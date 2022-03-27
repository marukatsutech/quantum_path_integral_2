# Quantum Mechanics, Path integral - 2nd
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d


def update_spirals():
    global y_spiral1, z_spiral1, spiral1, y_spiral2, z_spiral2, spiral2
    global y_spiral3, z_spiral3, spiral3, y_spiral4, z_spiral4, spiral4
    global y_spiral_superposed, z_spiral_superposed, spiral_superposed
    y_spiral1 = np.cos(y_theta0_clock1)
    z_spiral1 = np.sin(y_theta0_clock1)
    spiral1.set_xdata(x)
    spiral1.set_ydata(y_spiral1 + offset_spiral1)
    spiral1.set_3d_properties(z_spiral1)
    y_spiral2 = np.cos(y_theta0_clock2)
    z_spiral2 = np.sin(y_theta0_clock2)
    spiral2.set_xdata(x)
    spiral2.set_ydata(y_spiral2 + offset_spiral2)
    spiral2.set_3d_properties(z_spiral2)
    y_spiral3 = np.cos(y_theta0_clock3)
    z_spiral3 = np.sin(y_theta0_clock3)
    spiral3.set_xdata(x)
    spiral3.set_ydata(y_spiral3 + offset_spiral3)
    spiral3.set_3d_properties(z_spiral3)
    y_spiral4 = np.cos(y_theta0_clock4)
    z_spiral4 = np.sin(y_theta0_clock4)
    spiral4.set_xdata(x)
    spiral4.set_ydata(y_spiral4 + offset_spiral4)
    spiral4.set_3d_properties(z_spiral4)
    y_spiral_superposed = (y_spiral1 + y_spiral2 + y_spiral3 + y_spiral4) / 4.
    z_spiral_superposed = (z_spiral1 + z_spiral2 + z_spiral3 + z_spiral4) / 4.
    spiral_superposed .set_xdata(x)
    spiral_superposed .set_ydata(y_spiral_superposed + offset_spiral_superposed)
    spiral_superposed .set_3d_properties(z_spiral_superposed)


def update_curve_theta():
    global z_theta0
    global y_theta0_clock1, curve_theta0_clock1, y_theta0_clock2, curve_theta0_clock2
    global y_theta0_clock3, curve_theta0_clock3, y_theta0_clock4, curve_theta0_clock4
    z_theta0 = x * 0. + z_min
    if function_number == 0:
        y_theta0_clock1 = (x - x0_clock1) / t * 2. * np.pi + theta0_clock1
        y_theta0_clock2 = (x - x0_clock2) / t * 2. * np.pi + theta0_clock2
        y_theta0_clock3 = (x - x0_clock3) / t * 2. * np.pi + theta0_clock3
        y_theta0_clock4 = (x - x0_clock4) / t * 2. * np.pi + theta0_clock4
    elif function_number == 1:
        y_theta0_clock1 = np.abs(x - x0_clock1) * 2. * np.pi / t + theta0_clock1
        y_theta0_clock2 = np.abs(x - x0_clock2) * 2. * np.pi / t + theta0_clock2
        y_theta0_clock3 = np.abs(x - x0_clock3) * 2. * np.pi / t + theta0_clock3
        y_theta0_clock4 = np.abs(x - x0_clock4) * 2. * np.pi / t + theta0_clock4
    else:
        y_theta0_clock1 = (x - x0_clock1) ** 2 * 2. * np.pi / t + theta0_clock1
        y_theta0_clock2 = (x - x0_clock2) ** 2 * 2. * np.pi / t + theta0_clock2
        y_theta0_clock3 = (x - x0_clock3) ** 2 * 2. * np.pi / t + theta0_clock3
        y_theta0_clock4 = (x - x0_clock4) ** 2 * 2. * np.pi / t + theta0_clock4
    update_spirals()
    y_theta0_clock1[y_min > y_theta0_clock1] = np.nan
    y_theta0_clock1[y_max < y_theta0_clock1] = np.nan
    curve_theta0_clock1.set_xdata(x)
    curve_theta0_clock1.set_ydata(y_theta0_clock1)
    curve_theta0_clock1.set_3d_properties(z_theta0)
    y_theta0_clock2[y_min > y_theta0_clock2] = np.nan
    y_theta0_clock2[y_max < y_theta0_clock2] = np.nan
    curve_theta0_clock2.set_xdata(x)
    curve_theta0_clock2.set_ydata(y_theta0_clock2)
    curve_theta0_clock2.set_3d_properties(z_theta0)
    y_theta0_clock3[y_min > y_theta0_clock3] = np.nan
    y_theta0_clock3[y_max < y_theta0_clock3] = np.nan
    curve_theta0_clock3.set_xdata(x)
    curve_theta0_clock3.set_ydata(y_theta0_clock3)
    curve_theta0_clock3.set_3d_properties(z_theta0)
    y_theta0_clock4[y_min > y_theta0_clock4] = np.nan
    y_theta0_clock4[y_max < y_theta0_clock4] = np.nan
    curve_theta0_clock4.set_xdata(x)
    curve_theta0_clock4.set_ydata(y_theta0_clock4)
    curve_theta0_clock4.set_3d_properties(z_theta0)


def update_hands():
    global y_hand1, z_hand1, hand_clock1, y_hand2, z_hand2, hand_clock2
    global y_hand3, z_hand3, hand_clock3, y_hand4, z_hand4, hand_clock4
    y_hand1 = np.cos(theta0_clock1)
    z_hand1 = np.sin(theta0_clock1)
    hand_clock1.remove()
    hand_clock1 = ax1.add_line(
        art3d.Line3D([x0_clock1, x0_clock1], [0, y_hand1], [0, z_hand1], color='red', linewidth=1))
    y_hand2 = np.cos(theta0_clock2)
    z_hand2 = np.sin(theta0_clock2)
    hand_clock2.remove()
    hand_clock2 = ax1.add_line(
        art3d.Line3D([x0_clock2, x0_clock2], [0, y_hand2], [0, z_hand2], color='orange', linewidth=1))
    y_hand3 = np.cos(theta0_clock3)
    z_hand3 = np.sin(theta0_clock3)
    hand_clock3.remove()
    hand_clock3 = ax1.add_line(
        art3d.Line3D([x0_clock3, x0_clock3], [0, y_hand3], [0, z_hand3], color='green', linewidth=1))
    y_hand4 = np.cos(theta0_clock4)
    z_hand4 = np.sin(theta0_clock4)
    hand_clock4.remove()
    hand_clock4 = ax1.add_line(
        art3d.Line3D([x0_clock4, x0_clock4], [0, y_hand4], [0, z_hand4], color='blue', linewidth=1))


def update_clock():
    global c1, c2, c3, c4, clock1, clock2, clock3, clock4
    global ad_line_clock1, ad_line_clock2, ad_line_clock3, ad_line_clock4
    clock1.remove()
    c1 = Circle((0, 0), 1, ec='red', fill=False)
    clock1 = ax1.add_patch(c1)
    art3d.pathpatch_2d_to_3d(c1, z=x0_clock1, zdir="x")
    clock2.remove()
    c2 = Circle((0, 0), 1, ec='orange', fill=False)
    clock2 = ax1.add_patch(c2)
    art3d.pathpatch_2d_to_3d(c2, z=x0_clock2, zdir="x")
    clock3.remove()
    c3 = Circle((0, 0), 1, ec='green', fill=False)
    clock3 = ax1.add_patch(c3)
    art3d.pathpatch_2d_to_3d(c3, z=x0_clock3, zdir="x")
    clock4.remove()
    c4 = Circle((0, 0), 1, ec='blue', fill=False)
    clock4 = ax1.add_patch(c4)
    art3d.pathpatch_2d_to_3d(c4, z=x0_clock4, zdir="x")
    ad_line_clock1.remove()
    ad_line_clock1 = ax1.add_line(art3d.Line3D([x0_clock1, x0_clock1], [0, 0], [-1, z_min], color='red',
                                               ls='--', linewidth=1))
    ad_line_clock2.remove()
    ad_line_clock2 = ax1.add_line(art3d.Line3D([x0_clock2, x0_clock2], [0, 0], [-1, z_min], color='orange',
                                               ls='--', linewidth=1))
    ad_line_clock3.remove()
    ad_line_clock3 = ax1.add_line(art3d.Line3D([x0_clock3, x0_clock3], [0, 0], [-1, z_min], color='green',
                                               ls='--', linewidth=1))
    ad_line_clock4.remove()
    ad_line_clock4 = ax1.add_line(art3d.Line3D([x0_clock4, x0_clock4], [0, 0], [-1, z_min], color='blue',
                                               ls='--', linewidth=1))


def change_theta0_clock4(th0):
    global theta0_clock4
    theta0_clock4 = th0
    update_hands()
    update_curve_theta()


def change_theta0_clock3(th0):
    global theta0_clock3
    theta0_clock3 = th0
    update_hands()
    update_curve_theta()


def change_theta0_clock2(th0):
    global theta0_clock2
    theta0_clock2 = th0
    update_hands()
    update_curve_theta()


def change_theta0_clock1(th0):
    global theta0_clock1
    theta0_clock1 = th0
    update_hands()
    update_curve_theta()


def change_dx(delta_x):
    global dx, x0_clock1, x0_clock2, x0_clock3, x0_clock4
    dx = delta_x
    x0_clock1 = -dx / 2.
    x0_clock2 = -dx / 6.
    x0_clock3 = dx / 6.
    x0_clock4 = dx / 2.
    update_clock()
    update_hands()
    update_curve_theta()


def change_t(tm):
    global t
    t = tm
    update_curve_theta()


def change_function(fn):
    global function_number
    function_number = fn
    if function_number == 0:
        ax1.set_title('Quantum Mechanics, Path integral w = x/t')
    elif function_number == 1:
        ax1.set_title('Quantum Mechanics, Path integral w = abs(x/t)')
    else:
        ax1.set_title('Quantum Mechanics, Path integral w = x**2/t')
    update_curve_theta()


def switch():
    global is_running
    if is_running:
        is_running = False
    else:
        is_running = True


def update(f):
    global cnt
    # tx_step.set_text("Step=" + str(cnt))
    if is_running:
        cnt += 1


# Global variables
num_of_points = 4000

x_min = -10.
x_max = 10.
y_min = -4.
y_max = 20.
z_min = -6.
z_max = 6.

is_running = False
cnt = 0

dx_init = 2.
dx = dx_init
x0_clock1 = -dx / 2.
x0_clock2 = -dx / 6.
x0_clock3 = dx / 6.
x0_clock4 = dx / 2.

theta0_clock1_init = 0.
theta0_clock1 = theta0_clock1_init
theta0_clock2_init = 0.
theta0_clock2 = theta0_clock2_init
theta0_clock3_init = 0.
theta0_clock3 = theta0_clock3_init
theta0_clock4_init = 0.
theta0_clock4 = theta0_clock4_init

t_init = 1.
t = t_init

num_of_points = 10000

function_number = 0

offset_spiral1 = 3.
offset_spiral2 = 6.
offset_spiral3 = 9.
offset_spiral4 = 12.
offset_spiral_superposed = 16.

inc_spinbox = 0.1
inc_spinbox_t = 0.1

# Generate figure and axes
fig = Figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.set_box_aspect((2, 2, 1))
ax1.grid()
ax1.set_title('Quantum Mechanics, Path integral w = x/t')
ax1.set_xlabel('x')
ax1.set_ylabel('Phi (imaginary) or theta of clocks')
ax1.set_zlabel('Phi (Real)')
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_zlim(z_min, z_max)

# Generate items
# tx_step = ax1.text(x_min, y_max, z_max * 1.6, "Step=" + str(0))

c1 = Circle((0, 0), 1, ec='red', fill=False)
clock1 = ax1.add_patch(c1)
art3d.pathpatch_2d_to_3d(c1, z=x0_clock1, zdir="x")
c2 = Circle((0, 0), 1, ec='orange', fill=False)
clock2 = ax1.add_patch(c2)
art3d.pathpatch_2d_to_3d(c2, z=x0_clock2, zdir="x")
c3 = Circle((0, 0), 1, ec='green', fill=False)
clock3 = ax1.add_patch(c3)
art3d.pathpatch_2d_to_3d(c3, z=x0_clock3, zdir="x")
c4 = Circle((0, 0), 1, ec='blue', fill=False)
clock4 = ax1.add_patch(c4)
art3d.pathpatch_2d_to_3d(c4, z=x0_clock4, zdir="x")

y_hand1 = np.cos(theta0_clock1)
z_hand1 = np.sin(theta0_clock1)
hand_clock1 = ax1.add_line(art3d.Line3D([x0_clock1, x0_clock1], [0, y_hand1], [0, z_hand1], color='red', linewidth=1))
ad_line_clock1 = ax1.add_line(art3d.Line3D([x0_clock1, x0_clock1], [0, 0], [-1, z_min], color='red',
                                           ls='--', linewidth=1))
y_hand2 = np.cos(theta0_clock2)
z_hand2 = np.sin(theta0_clock2)
hand_clock2 = ax1.add_line(art3d.Line3D([x0_clock2, x0_clock2], [0, y_hand2], [0, z_hand2], color='orange', linewidth=1))
ad_line_clock2 = ax1.add_line(art3d.Line3D([x0_clock2, x0_clock2], [0, 0], [-1, z_min], color='orange',
                                           ls='--', linewidth=1))
y_hand3 = np.cos(theta0_clock3)
z_hand3 = np.sin(theta0_clock3)
hand_clock3 = ax1.add_line(art3d.Line3D([x0_clock3, x0_clock3], [0, y_hand3], [0, z_hand3], color='green', linewidth=1))
ad_line_clock3 = ax1.add_line(art3d.Line3D([x0_clock3, x0_clock3], [0, 0], [-1, z_min], color='green',
                                           ls='--', linewidth=1))
y_hand4 = np.cos(theta0_clock4)
z_hand4 = np.sin(theta0_clock4)
hand_clock4 = ax1.add_line(art3d.Line3D([x0_clock4, x0_clock4], [0, y_hand4], [0, z_hand4], color='blue', linewidth=1))
ad_line_clock4 = ax1.add_line(art3d.Line3D([x0_clock4, x0_clock4], [0, 0], [-1, z_min], color='blue',
                                           ls='--', linewidth=1))

x = np.linspace(x_min, x_max, num_of_points)

y_theta0_clock1 = (x - x0_clock1) * 2. * np.pi / t + theta0_clock1
y_theta0_clock2 = (x - x0_clock2) * 2. * np.pi / t + theta0_clock2
y_theta0_clock3 = (x - x0_clock3) * 2. * np.pi / t + theta0_clock3
y_theta0_clock4 = (x - x0_clock4) * 2. * np.pi / t + theta0_clock4
z_theta0 = x * 0. + z_min

y_spiral1 = np.cos(y_theta0_clock1)
z_spiral1 = np.sin(y_theta0_clock1)
spiral1, = ax1.plot(x, y_spiral1 + offset_spiral1, z_spiral1, color='red', linewidth=1, label='Clock1, Phi1')
y_spiral2 = np.cos(y_theta0_clock2)
z_spiral2 = np.sin(y_theta0_clock2)
spiral2, = ax1.plot(x, y_spiral2 + offset_spiral2, z_spiral2, color='orange', linewidth=1, label='Clock2, Phi2')
y_spiral3 = np.cos(y_theta0_clock3)
z_spiral3 = np.sin(y_theta0_clock3)
spiral3, = ax1.plot(x, y_spiral3 + offset_spiral3, z_spiral3, color='green', linewidth=1, label='Clock3, Phi3')
y_spiral4 = np.cos(y_theta0_clock4)
z_spiral4 = np.sin(y_theta0_clock4)
spiral4, = ax1.plot(x, y_spiral4 + offset_spiral4, z_spiral4, color='blue', linewidth=1, label='Clock4, Phi4')
y_spiral_superposed = (y_spiral1 + y_spiral2 + y_spiral3 + y_spiral4) / 4.
z_spiral_superposed = (z_spiral1 + z_spiral2 + z_spiral3 + z_spiral4) / 4.
spiral_superposed, = ax1.plot(x, y_spiral_superposed + offset_spiral_superposed, z_spiral_superposed,
                              color='black', linewidth=1, label='Superposed phi (ave.)')

y_theta0_clock1[y_min > y_theta0_clock1] = np.nan
y_theta0_clock1[y_max < y_theta0_clock1] = np.nan
curve_theta0_clock1, = ax1.plot(x, y_theta0_clock1, z_theta0, color='red', ls='--', linewidth=1,
                                label='Theta of clock1')
y_theta0_clock2[y_min > y_theta0_clock2] = np.nan
y_theta0_clock2[y_max < y_theta0_clock2] = np.nan
curve_theta0_clock2, = ax1.plot(x, y_theta0_clock2, z_theta0, color='orange', ls='--', linewidth=1,
                                label='Theta of clock2')
y_theta0_clock3[y_min > y_theta0_clock3] = np.nan
y_theta0_clock3[y_max < y_theta0_clock3] = np.nan
curve_theta0_clock3, = ax1.plot(x, y_theta0_clock3, z_theta0, color='green', ls='--', linewidth=1,
                                label='Theta of clock2')
y_theta0_clock4[y_min > y_theta0_clock4] = np.nan
y_theta0_clock4[y_max < y_theta0_clock4] = np.nan
curve_theta0_clock4, = ax1.plot(x, y_theta0_clock4, z_theta0, color='blue', ls='--', linewidth=1,
                                label='Theta of clock4')

ax1.legend()

# Embed in Tkinter
root = tk.Tk()
root.title("Quantum Mechanics, Path integral")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

frm_fn = ttk.Labelframe(root, relief="ridge", text="Function", labelanchor="n", width=100)
frm_fn.pack(side='left')
var_function = tk.IntVar(value=0)
rdb0_function = tk.Radiobutton(frm_fn, text="w=x/t", command=lambda: change_function(var_function.get()),
                               variable=var_function, value=0)
rdb1_function = tk.Radiobutton(frm_fn, text="w=abs(x/t)", command=lambda: change_function(var_function.get()),
                               variable=var_function, value=1)
rdb2_function = tk.Radiobutton(frm_fn, text="w=x**2/t", command=lambda: change_function(var_function.get()),
                               variable=var_function, value=2)
rdb0_function.pack(anchor=tk.W)
rdb1_function.pack(anchor=tk.W)
rdb2_function.pack(anchor=tk.W)

frm_clk = ttk.Labelframe(root, relief="ridge", text="Theta", labelanchor="n", width=100)
frm_clk .pack(side='left')
lbl_clock1 = tk.Label(frm_clk, text="Clock1:")
lbl_clock1.pack(side='left')
var_clock1 = tk.StringVar(root)  # variable for spinbox-value
var_clock1.set(theta0_clock1_init)  # Initial value
spn_clock1 = tk.Spinbox(
    frm_clk, textvariable=var_clock1, format="%.2f", from_=-10., to=10., increment=inc_spinbox,
    command=lambda: change_theta0_clock1(float(var_clock1.get())), width=4
    )
spn_clock1.pack(side='left')
lbl_clock2 = tk.Label(frm_clk, text="Clock2:")
lbl_clock2.pack(side='left')
var_clock2 = tk.StringVar(root)  # variable for spinbox-value
var_clock2.set(theta0_clock2_init)  # Initial value
spn_clock2 = tk.Spinbox(
    frm_clk, textvariable=var_clock2, format="%.2f", from_=-10., to=10., increment=inc_spinbox,
    command=lambda: change_theta0_clock2(float(var_clock2.get())), width=4
    )
spn_clock2.pack(side='left')
lbl_clock3 = tk.Label(frm_clk, text="Clock3:")
lbl_clock3.pack(side='left')
var_clock3 = tk.StringVar(root)  # variable for spinbox-value
var_clock3.set(theta0_clock3_init)  # Initial value
spn_clock3 = tk.Spinbox(
    frm_clk, textvariable=var_clock3, format="%.2f", from_=-10., to=10., increment=inc_spinbox,
    command=lambda: change_theta0_clock3(float(var_clock3.get())), width=4
    )
spn_clock3.pack(side='left')
lbl_clock4 = tk.Label(frm_clk, text="Clock4:")
lbl_clock4.pack(side='left')
var_clock4 = tk.StringVar(root)  # variable for spinbox-value
var_clock4.set(theta0_clock4_init)  # Initial value
spn_clock4 = tk.Spinbox(
    frm_clk, textvariable=var_clock4, format="%.2f", from_=-10., to=10., increment=inc_spinbox,
    command=lambda: change_theta0_clock4(float(var_clock4.get())), width=4
    )
spn_clock4.pack(side='left')

lbl_dx = tk.Label(root, text="dx:")
lbl_dx.pack(side='left')
var_dx = tk.StringVar(root)  # variable for spinbox-value
var_dx.set(dx_init)  # Initial value
spn_dx = tk.Spinbox(
    root, textvariable=var_dx, format="%.2f", from_=inc_spinbox, to=10., increment=inc_spinbox,
    command=lambda: change_dx(float(var_dx.get())), width=4
    )
spn_dx.pack(side='left')

lbl_t = tk.Label(root, text="t:")
lbl_t.pack(side='left')
var_t = tk.StringVar(root)  # variable for spinbox-value
var_t.set(t_init)  # Initial value
spn_t = tk.Spinbox(
    root, textvariable=var_t, format="%.2f", from_=inc_spinbox_t, to=100., increment=inc_spinbox_t,
    command=lambda: change_t(float(var_t.get())), width=4
    )
spn_t.pack(side='left')

# btn_play = tk.Button(root, text="Play/Pause", command=switch)
# btn_play.pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=200)
root.mainloop()
