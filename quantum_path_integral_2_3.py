# Quantum Mechanics, Path integral - 2_3

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from mpl_toolkits.mplot3d import proj3d


def change_function(value):
    global flag_function
    flag_function = int(value)
    update_qtm1()
    update_path()


def step():
    global cnt_step
    pass
    cnt_step += 1


def update_path():
    global qtm1
    global plt_path_list
    global y1_1
    for k in range(1, 10):
        qtm1a_buffer = qtm1 * 0.
        qtm1b_buffer = qtm1 * 0.
        # Quantum A
        for i in range(len(x_a)):
            dx_squared = (x_a[i] - x) ** 2
            if t != 0.:
                theta = np.mod(np.deg2rad(mass * dx_squared / (2. * h * (k * np.power(10., te)))), (2. * np.pi))
            else:
                theta = 0.
                print('Error: divided by zero')
            rot = np.sin(theta) + np.cos(theta) * 1j
            qtm1a_buffer = qtm1a_buffer + qtm0a[i] * rot
        # Quantum B
        for j in range(len(x_b)):
            dx_squared = (x_b[j] - x) ** 2
            if t != 0.:
                theta = np.mod(np.deg2rad(mass * dx_squared / (2. * h * (k * np.power(10., te)))), (2. * np.pi))
            else:
                theta = 0.
                print('Error: divided by zero')
            rot = np.sin(theta) + np.cos(theta) * 1j
            qtm1b_buffer = qtm1b_buffer + qtm0b[j] * rot
        if flag_function == 0:
            qtm1_ = qtm1a_buffer / np.sum(np.abs(qtm1a_buffer)) + qtm1b_buffer / np.sum(np.abs(qtm1b_buffer))
        else:
            qtm1_ = qtm1a_buffer / np.sum(np.abs(qtm1a_buffer)) - qtm1b_buffer / np.sum(np.abs(qtm1b_buffer))
        # y1_ = qtm1_.imag
        z1_ = qtm1_.real
        y1_1_ = x1 * 0. + k
        plt_path_list[k - 1].set_xdata(x)
        plt_path_list[k - 1].set_ydata(y1_1_)
        plt_path_list[k - 1].set_3d_properties(z1_)


def update_qtm1():
    global qtm1
    global y1, z1, y1_offset
    global plt_qtm1, plt_qtm1_path_real
    qtm1a_buffer = qtm1 * 0.
    qtm1b_buffer = qtm1 * 0.
    # Quantum A
    for i in range(len(x_a)):
        dx_squared = (x_a[i] - x) ** 2
        if t != 0.:
            theta = np.mod(np.deg2rad(mass * dx_squared / (2. * h * t)), (2. * np.pi))
        else:
            theta = 0.
            print('Error: divided by zero')
        rot = np.sin(theta) + np.cos(theta) * 1j
        qtm1a_buffer = qtm1a_buffer + qtm0a[i] * rot
    # Quantum B
    for j in range(len(x_b)):
        dx_squared = (x_b[j] - x) ** 2
        if t != 0.:
            theta = np.mod(np.deg2rad(mass * dx_squared / (2. * h * t)), (2. * np.pi))
        else:
            theta = 0.
            print('Error: divided by zero')
        rot = np.sin(theta) + np.cos(theta) * 1j
        qtm1b_buffer = qtm1b_buffer + qtm0b[j] * rot
    if flag_function == 0:
        qtm1 = qtm1a_buffer / np.sum(np.abs(qtm1a_buffer)) + qtm1b_buffer / np.sum(np.abs(qtm1b_buffer))
    else:
        qtm1 = qtm1a_buffer / np.sum(np.abs(qtm1a_buffer)) - qtm1b_buffer / np.sum(np.abs(qtm1b_buffer))
    y1 = qtm1.imag
    z1 = qtm1.real
    y1_offset = y1 + offset1
    plt_qtm1.set_xdata(x)
    plt_qtm1.set_ydata(y1_offset)
    plt_qtm1.set_3d_properties(z1)

    y1_path = x * 0. + tm
    plt_qtm1_path_real.set_xdata(x)
    plt_qtm1_path_real.set_ydata(y1_path)
    plt_qtm1_path_real.set_3d_properties(z1)


def update_qtm0a():
    global gaussian0a, gaussian0a_std
    global z0a, y0a, y0a, qtm0a
    global ya, za, ya_offset
    global plt_qtm0a
    global x_a
    global ya_zero, plt_qtm0a_real
    x_a = np.arange(mu0a - width_a / 2., mu0a + width_a / 2.)
    gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x_a - mu0a) ** 2 / (2 * sigma0a ** 2))
    gaussian0a_std = gaussian0a / np.sum(np.abs(gaussian0a))
    z0a = np.sin(k0a * x_a) * gaussian0a_std
    y0a = np.cos(k0a * x_a) * 1j * gaussian0a_std
    qtm0a = z0a + y0a
    ya = qtm0a.imag
    za = qtm0a.real
    ya_offset = ya + offset
    plt_qtm0a.set_xdata(x_a)
    plt_qtm0a.set_ydata(ya_offset)
    plt_qtm0a.set_3d_properties(za)

    ya_zero = x_a * 0. + offset
    plt_qtm0a_real.set_xdata(x_a)
    plt_qtm0a_real.set_ydata(ya_zero)
    plt_qtm0a_real.set_3d_properties(za)


def update_qtm0b():
    global gaussian0b, gaussian0b_std
    global z0b, y0b, y0b, qtm0b
    global yb, zb, yb_offset
    global plt_qtm0b
    global x_b
    global yb_zero, plt_qtm0b_real
    x_b = np.arange(mu0b - width_b / 2., mu0b + width_b / 2.)
    gaussian0b = 1 / (np.sqrt(2 * np.pi) * sigma0b) * np.exp(- (x_b - mu0b) ** 2 / (2 * sigma0b ** 2))
    gaussian0b_std = gaussian0b / np.sum(np.abs(gaussian0b))
    z0b = np.sin(k0b * x_b) * gaussian0b_std
    y0b = np.cos(k0b * x_b) * 1j * gaussian0b_std
    qtm0b = z0b + y0b
    yb = qtm0b.imag
    zb = qtm0b.real
    yb_offset = yb + offset
    plt_qtm0b.set_xdata(x_b)
    plt_qtm0b.set_ydata(yb_offset)
    plt_qtm0b.set_3d_properties(zb)

    yb_zero = x_a * 0. + offset
    plt_qtm0b_real.set_xdata(x_b)
    plt_qtm0b_real.set_ydata(yb_zero)
    plt_qtm0b_real.set_3d_properties(zb)


def set_k0a(value):
    global k0a
    k0a = float(value)
    update_qtm0a()
    update_qtm1()
    update_path()


def set_x0a(value):
    global mu0a
    mu0a = float(value)
    update_qtm0a()
    update_qtm1()
    update_path()


def set_sigma0a(value):
    global sigma0a
    sigma0a = float(value)
    update_qtm0a()
    update_qtm1()
    update_path()


def set_k0b(value):
    global k0b
    k0b = float(value)
    update_qtm0b()
    update_qtm1()
    update_path()


def set_x0b(value):
    global mu0b
    mu0b = float(value)
    update_qtm0b()
    update_qtm1()
    update_path()


def set_sigma0b(value):
    global sigma0b
    sigma0b = float(value)
    update_qtm0b()
    update_qtm1()
    update_path()


def set_tm(value):
    global tm, t
    tm = float(value)
    t = tm * np.power(10., te)
    update_qtm1()
    update_path()


def set_te(value):
    global te, t
    te = float(value)
    t = tm * np.power(10., te)
    update_qtm1()
    update_path()


def function_selected(event):
    pass


# Animation control
def reset():
    global is_play, cnt
    global cnt_step
    is_play = False
    # cnt = 0
    cnt_step = 1
    update_qtm1()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    # global cnt
    global txt_step
    txt_step.set_text("Delta t:" + str(tm) + "E" + str(te))
    if is_play:
        # cnt += 1
        pass


# Global variables
# Animation control
cnt = 0
is_play = False
cnt_step = 1

# Data structure
range_x = 3000
range_y = 10
x = np.arange(- range_x, range_x)
y = np.arange(1, range_y)

range_quantum_yz = 0.02

flag_function = 0

# Quantum parameter
mass = 9.1093837015 * 1.0E-31
tm, te = 1., 5
t = tm * np.power(10., te)
h = 6.62607015 * 1.0E-34

# Quantum
# Quantum 0a
k0a = 0.2
sigma0a = 40.
mu0a = - range_x / 2.
width_a = 150
x_a = np.arange(mu0a - width_a / 2., mu0a + width_a / 2.)
gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x_a - mu0a) ** 2 / (2 * sigma0a ** 2))
gaussian0a_std = gaussian0a / np.sum(np.abs(gaussian0a))
z0a = np.sin(k0a * x_a) * gaussian0a_std
y0a = np.cos(k0a * x_a) * 1j * gaussian0a_std
qtm0a = z0a + y0a

# Quantum 0b
# k0b = - 12.
k0b = - 0.2
sigma0b = 40.
mu0b = range_x / 2.
width_b = 150
x_b = np.arange(mu0b - width_b / 2., mu0b + width_b / 2.)
gaussian0b = 1 / (np.sqrt(2 * np.pi) * sigma0b) * np.exp(- (x_b - mu0b) ** 2 / (2 * sigma0b ** 2))
gaussian0b_std = gaussian0b / np.sum(np.abs(gaussian0b))
z0b = np.sin(k0b * x_b) * gaussian0b_std
y0b = np.cos(k0b * x_b) * 1j * gaussian0b_std
qtm0b = z0b + y0b

# Quantum 1 (Initial)
z1 = x * 0.
y1 = x * 0. * 1j
qtm1 = z1 + y1


# Generate figure and axes
title_tk = "Quantum Mechanics, Path integral"
title_ax0 = "Quantum"
title_ax1 = "Quantum (real part)"

x_min0 = - range_x
x_max0 = range_x
y_min0 = - range_quantum_yz * 2.
y_max0 = range_quantum_yz * 2.
z_min0 = - range_quantum_yz
z_max0 = range_quantum_yz

x_min1 = - range_x
x_max1 = range_x
y_min1 = 0.
y_max1 = range_y
z_min1 = - range_quantum_yz
z_max1 = range_quantum_yz

fig = Figure()
ax0 = fig.add_subplot(121, projection='3d')
ax0.set_box_aspect((6, 4, 2))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel('x')
ax0.set_ylabel('y')
ax0.set_zlabel('z')
ax0.set_xlim(x_min0, x_max0)
ax0.set_ylim(y_min0, y_max0)
ax0.set_zlim(z_min0, z_max0)

ax1 = fig.add_subplot(122, projection='3d')
ax1.set_box_aspect((4, 4, 1))
ax1.grid()
ax1.set_title(title_ax1)
ax1.set_xlabel('x')
ax1.set_ylabel('Delta t' + "(E" + str(te) + ")")
ax1.set_zlabel('z')
ax1.set_xlim(x_min1, x_max1)
ax1.set_ylim(y_min1, y_max1)
ax1.set_zlim(z_min1, z_max1)

# Text items
txt_step = ax0.text2D(x_min0, y_max0, "Delta t * " + str(cnt_step))
xz, yz, _ = proj3d.proj_transform(x_min0, y_max0, z_max0, ax0.get_proj())
txt_step.set_position((xz, yz))

# Plot items
# Quantum 0a
offset = - range_quantum_yz
xa = x_a
ya = qtm0a.imag
za = qtm0a.real
ya_offset = ya + offset
plt_qtm0a, = ax0.plot(xa, ya_offset, za, color='blue', linewidth=0.5, label='Quantum A(t0)')

ya_zero = x_a * 0. + offset
plt_qtm0a_real, = ax1.plot(xa, ya_zero, za, color='blue', linewidth=0.5, label='Quantum A(t0)')

# Quantum 0b
offset = - range_quantum_yz
xb = x_b
yb = qtm0b.imag
zb = qtm0b.real
yb_offset = yb + offset
plt_qtm0b, = ax0.plot(xb, yb_offset, zb, color='red', linewidth=0.5, label='Quantum B(t0)')

yb_zero = x_b * 0. + offset
plt_qtm0b_real, = ax1.plot(xb, yb_zero, zb, color='red', linewidth=0.5, label='Quantum B(t0)')

# Quantum 1
offset1 = range_quantum_yz
x1 = x
y1 = qtm1.imag
z1 = qtm1.real
y1_offset = y1 + offset1
plt_qtm1, = ax0.plot(x1, y1_offset, z1, color='green', linewidth=0.5, label='Quantum(delta t)')

y1_1 = x1 * 0. + tm
plt_qtm1_path_real, = ax1.plot(x1, y1_1, z1, color='green', linewidth=0.5, label='Quantum(delta t)')
update_qtm1()

# Path
plt_path_list = []
for i_ in range(1, 10):
    y1_1 = x1 * 0. + 5 * i_
    plt_qtm1_path_real_, = ax1.plot(x1, y1_1, z1, color='gray', linewidth=0.5, label='Quantum(delta t)')
    plt_path_list.append(plt_qtm1_path_real_)
update_path()

# Legend
ax0.legend(loc='lower right', fontsize=6)

# Embed in Tkinter
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Animation
'''
frm_anim = ttk.Labelframe(root, relief='ridge', text='Animation', labelanchor='n')
frm_anim.pack(side='left', fill=tk.Y)
# btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
# btn_play.pack(side='left')
btn_step = tk.Button(frm_anim, text="Step", command=step)
btn_step.pack(side='left')
btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
btn_reset.pack(side='left')
'''
# Quantum A
frm_qtm_a = ttk.Labelframe(root, relief="ridge", text="Quantum A", labelanchor="n")
frm_qtm_a.pack(side='left', fill=tk.Y)

lbl_k0a = tk.Label(frm_qtm_a, text="k:")
lbl_k0a.pack(side='left')
var_k0a = tk.StringVar(root)
var_k0a.set(str(k0a))
spn_k0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_k0a, format="%.2f", from_=-2., to=2., increment=0.01,
    command=lambda: set_k0a(var_k0a.get()), width=6
    )
spn_k0a.pack(side='left')

lbl_x0a = tk.Label(frm_qtm_a, text="x:")
lbl_x0a.pack(side='left')
var_x0a = tk.StringVar(root)
var_x0a.set(str(mu0a))
spn_x0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_x0a, format="%.1f", from_=x_min0, to=x_max0, increment=100,
    command=lambda: set_x0a(var_x0a.get()), width=6
    )
spn_x0a.pack(side='left')

lbl_s0a = tk.Label(frm_qtm_a, text="sigma:")
lbl_s0a.pack(side='left')
var_s0a = tk.StringVar(root)
var_s0a.set(str(sigma0a))
spn_s0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_s0a, format="%.1f", from_=-600, to=600, increment=10,
    command=lambda: set_sigma0a(var_s0a.get()), width=6
    )
spn_s0a.pack(side='left')

# Quantum b
frm_qtm_b = ttk.Labelframe(root, relief="ridge", text="Quantum b", labelanchor="n")
frm_qtm_b.pack(side='left', fill=tk.Y)

lbl_k0b = tk.Label(frm_qtm_b, text="k:")
lbl_k0b.pack(side='left')
var_k0b = tk.StringVar(root)
var_k0b.set(str(k0b))
spn_k0b = tk.Spinbox(
    frm_qtm_b, textvariable=var_k0b, format="%.2f", from_=-2., to=2., increment=0.01,
    command=lambda: set_k0b(var_k0b.get()), width=6
    )
spn_k0b.pack(side='left')

lbl_x0b = tk.Label(frm_qtm_b, text="x:")
lbl_x0b.pack(side='left')
var_x0b = tk.StringVar(root)
var_x0b.set(str(mu0b))
spn_x0b = tk.Spinbox(
    frm_qtm_b, textvariable=var_x0b, format="%.1f", from_=x_min0, to=x_max0, increment=100,
    command=lambda: set_x0b(var_x0b.get()), width=6
    )
spn_x0b.pack(side='left')

lbl_s0b = tk.Label(frm_qtm_b, text="sigma:")
lbl_s0b.pack(side='left')
var_s0b = tk.StringVar(root)
var_s0b.set(str(sigma0a))
spn_s0b = tk.Spinbox(
    frm_qtm_b, textvariable=var_s0b, format="%.1f", from_=-600, to=600, increment=10,
    command=lambda: set_sigma0b(var_s0b.get()), width=6
    )
spn_s0b.pack(side='left')

# Delta t
frm_t = ttk.Labelframe(root, relief="ridge", text="Delta t", labelanchor="n")
frm_t.pack(side='left', fill=tk.Y)

lbl_tm = tk.Label(frm_t, text="Mantissa:")
lbl_tm.pack(side='left')
var_tm = tk.StringVar(root)
var_tm.set(str(tm))
spn_tm = tk.Spinbox(
    frm_t, textvariable=var_tm, format="%.1f", from_=-100., to=100., increment=0.1,
    command=lambda: set_tm(var_tm.get()), width=6
    )
spn_tm.pack(side='left')

lbl_te = tk.Label(frm_t, text="Exponent:")
lbl_te.pack(side='left')
var_te = tk.StringVar(root)
var_te.set(str(te))
spn_te = tk.Spinbox(
    frm_t, textvariable=var_te, format="%.1f", from_=-10., to=10., increment=0.1,
    command=lambda: set_te(var_te.get()), width=6
    )
spn_te.pack(side='left')

frm_fn = ttk.Labelframe(root, relief="ridge", text="Function", labelanchor="n", width=100)
frm_fn.pack(side='left')
var_function = tk.IntVar(value=0)
rdb0_function = tk.Radiobutton(frm_fn, text="Qtm dt=QtmA+QtmB", command=lambda: change_function(var_function.get()),
                               variable=var_function, value=0)
rdb1_function = tk.Radiobutton(frm_fn, text="Qtm dt=QtmA-QtmB", command=lambda: change_function(var_function.get()),
                               variable=var_function, value=1)
rdb0_function.pack(anchor=tk.W)
rdb1_function.pack(anchor=tk.W)

# main loop
anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
root.mainloop()

