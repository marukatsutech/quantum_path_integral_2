# Quantum Mechanics, Path integral - 2_1

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk


def step():
    global x1, y1, z1, y1_offset
    global qtm1, qtm1_buffer
    global plt_qtm1
    qtm1_buffer = np.dot(qtm1, phase_convert_array)
    qtm1 = qtm1_buffer / np.sum(np.abs(qtm1_buffer))
    x1 = x
    y1 = qtm1.imag
    z1 = qtm1.real
    y1_offset = y1 + offset
    plt_qtm1.set_xdata(x1)
    plt_qtm1.set_ydata(y1_offset)
    plt_qtm1.set_3d_properties(z1)


def update_qtm1():
    global x1, y1, z1, y1_offset
    global qtm1, qtm1_buffer
    global plt_qtm1
    qtm1_buffer = np.dot(qtm0a, phase_convert_array)
    qtm1 = qtm1_buffer / np.sum(np.abs(qtm1_buffer))
    x1 = x
    y1 = qtm1.imag
    z1 = qtm1.real
    y1_offset = y1 + offset
    plt_qtm1.set_xdata(x1)
    plt_qtm1.set_ydata(y1_offset)
    plt_qtm1.set_3d_properties(z1)


def update_qtm0a():
    global gaussian0a, gaussian0a_std
    global z0a, y0a, y0a, qtm0a
    global ya, za, ya_offset
    global plt_qtm0a
    gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x - mu0a) ** 2 / (2 * sigma0a ** 2))
    gaussian0a_std = gaussian0a / np.sum(np.abs(gaussian0a))
    z0a = np.sin(k0a * x * k_x_ratio) * gaussian0a_std
    y0a = np.cos(k0a * x * k_x_ratio) * 1j * gaussian0a_std
    qtm0a = z0a + y0a
    ya = qtm0a.imag
    za = qtm0a.real
    ya_offset = ya - offset
    plt_qtm0a.set_xdata(x)
    plt_qtm0a.set_ydata(ya_offset)
    plt_qtm0a.set_3d_properties(za)


def set_k0a(value):
    global k0a
    k0a = float(value)
    # function_selected(True)
    update_qtm0a()
    update_qtm1()


def set_x0a(value):
    global mu0a
    mu0a = float(value)
    # function_selected(True)
    update_qtm0a()
    update_qtm1()


def set_sigma0a(value):
    global sigma0a
    sigma0a = float(value)
    # function_selected(True)
    update_qtm0a()
    update_qtm1()


def set_tm(value):
    global tm, t
    tm = float(value)
    t = tm * np.power(10., te)
    function_selected(True)
    update_qtm1()


def set_te(value):
    global te, t
    te = float(value)
    t = tm * np.power(10., te)
    function_selected(True)
    update_qtm1()


def function_selected(event):
    global combo_func
    global zz1, zz1_min, zz1_max
    global z_min1, z_max1
    global ax1
    global plt1
    global theta_conv, phase_convert_array
    # print(combo_func.get())
    zz1_conv = (xx - yy) ** 2
    if t != 0.:
        theta_conv = np.mod(np.deg2rad(mass * zz1_conv / (2. * h * t)), (2. * np.pi))
        phase_convert_array = np.sin(theta_conv) + np.cos(theta_conv) * 1j
    else:
        theta_conv = 0. * zz1_conv
        phase_convert_array = np.sin(theta_conv) + np.cos(theta_conv) * 1j
    if combo_func.get() == option_func[0]:
        ax1.set_title(title_ax1 + "(x - y)**2")
        zz1 = zz1_conv
    elif combo_func.get() == option_func[1]:
        if t != 0.:
            ax1.set_title(title_ax1 + "Rotation(real part):cos(mass*(x - y)**2/(2*h*t))")
            theta = np.mod(np.deg2rad(mass * zz1_conv / (2. * h * t)), (2. * np.pi))
            zz1 = np.cos(theta)
        else:
            ax1.set_title(title_ax1 + "Not available(t=0)")
            zz1 = zz1_conv * 0.
    else:
        pass
    zz1_min = np.min(zz1)
    zz1_max = np.max(zz1)
    if zz1_min == zz1_max:
        zz1_min -= 1.
        zz1_max += 1.
    if combo_func.get() == option_func[0]:
        z_min1 = np.round(zz1_min)
        z_max1 = np.round(zz1_max)
    else:
        z_min1 = np.round(zz1_min) * 10.
        z_max1 = np.round(zz1_max) * 10.
    ax1.set_zlim(z_min1, z_max1)
    plt1.remove()
    plt1 = ax1.plot_wireframe(xx, yy, zz1, linewidth=1, rstride=200, cstride=200)


# Animation control
def reset():
    global is_play, cnt
    is_play = False
    # cnt = 0
    update_qtm1()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt
    # global txt_step
    # txt_step.set_text("Step=" + str(cnt))
    if is_play:
        cnt += 1


# Global variables
# Animation control
cnt = 0
is_play = False

# Data structure
range_xy = 3000
x = np.arange(- range_xy, range_xy)
y = np.arange(- range_xy, range_xy)
xx, yy = np.meshgrid(x, y)
zz1 = (xx - yy) ** 2

zz1_min = np.min(zz1)
zz1_max = np.max(zz1)

range_quantum_yz = 0.002

# Quantum parameter
mass = 9.1093837015 * 1.0E-31
tm, te = 1., 6
t = tm * np.power(10., te)
h = 6.62607015 * 1.0E-34

# Convert array
theta_conv = np.mod(np.deg2rad(mass * zz1 / (2. * h * t)), (2. * np.pi))
phase_convert_array = np.sin(theta_conv) + np.cos(theta_conv) * 1j

# Quantum 0a
k0a = 12.
sigma0a = 200.
mu0a = 0.
gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x - mu0a) ** 2 / (2 * sigma0a ** 2))
gaussian0a_std = gaussian0a / np.sum(np.abs(gaussian0a))
k_x_ratio = 1 / 360.
z0a = np.sin(k0a * x * k_x_ratio) * gaussian0a_std
y0a = np.cos(k0a * x * k_x_ratio) * 1j * gaussian0a_std
qtm0a = z0a + y0a

# Quantum 1 (Initial)
qtm1_buffer = np.dot(qtm0a, phase_convert_array)
qtm1 = qtm1_buffer / np.sum(np.abs(qtm1_buffer))

# Generate figure and axes
title_tk = "Quantum Mechanics, Path integral"
title_ax0 = "Quantum\n"
title_ax1 = "Conversion matrix monitor:\n"

x_min0 = - range_xy
x_max0 = range_xy
y_min0 = - range_quantum_yz * 2.
y_max0 = range_quantum_yz * 2.
z_min0 = - range_quantum_yz
z_max0 = range_quantum_yz

x_min1 = - range_xy
x_max1 = range_xy
y_min1 = x_min1
y_max1 = x_max1
z_min1 = np.round(zz1_min)
z_max1 = np.round(zz1_max)

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
ax1.set_box_aspect((4, 4, 4))
ax1.grid()
ax1.set_title(title_ax1 + "(x0-x1) ** 2")
ax1.set_xlabel('x0')
ax1.set_ylabel('x1')
ax1.set_zlabel('z')
ax1.set_xlim(x_min1, x_max1)
ax1.set_ylim(y_min1, y_max1)
ax1.set_zlim(z_min1, z_max1)

# Text items
'''
txt_step = ax0.text2D(x_min0, y_max0, "Step=" + str(0))
xz, yz, _ = proj3d.proj_transform(x_min0, y_max0, z_max0, ax0.get_proj())
txt_step.set_position((xz, yz))
'''

# Plot items
# Convert array
plt1 = ax1.plot_wireframe(xx, yy, zz1, linewidth=1, rstride=200, cstride=200)

# Quantum 0a
offset = - range_quantum_yz
xa = x
ya = qtm0a.imag
za = qtm0a.real
ya_offset = ya + offset
plt_qtm0a, = ax0.plot(xa, ya_offset, za, color='blue', linewidth=0.5, label='Quantum(t0)')

# Quantum 1
offset = range_quantum_yz
x1 = x
y1 = qtm1.imag
z1 = qtm1.real
y1_offset = y1 + offset
plt_qtm1, = ax0.plot(x1, y1_offset, z1, color='green', linewidth=0.5, label='Quantum(t1)')

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
frm_anim = ttk.Labelframe(root, relief='ridge', text='Animation', labelanchor='n')
frm_anim.pack(side='left', fill=tk.Y)
# btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
# btn_play.pack(side='left')
btn_step = tk.Button(frm_anim, text="Step", command=step)
btn_step.pack(side='left')
btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
btn_reset.pack(side='left')

# Quantum A
frm_qtm_a = ttk.Labelframe(root, relief="ridge", text="Quantum A", labelanchor="n")
frm_qtm_a.pack(side='left', fill=tk.Y)

lbl_k0a = tk.Label(frm_qtm_a, text="k:")
lbl_k0a.pack(side='left')
var_k0a = tk.StringVar(root)
var_k0a.set(str(k0a))
spn_k0a = tk.Spinbox(
    frm_qtm_a, textvariable=var_k0a, format="%.1f", from_=-20., to=20., increment=0.5,
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
    frm_qtm_a, textvariable=var_s0a, format="%.1f", from_=-400, to=400, increment=10,
    command=lambda: set_sigma0a(var_s0a.get()), width=6
    )
spn_s0a.pack(side='left')

# Convert matrix monitor option
frm_func = ttk.Labelframe(root, relief='ridge', text='Conversion matrix monitor', labelanchor='n')
frm_func.pack(side='left', fill=tk.Y)
option_func = ["(x0-x1)**2", "(x0-x1)**2->rotation"]
variable = tk.StringVar()
combo_func = ttk.Combobox(frm_func, values=option_func, textvariable=variable, width=30)
combo_func.set(option_func[0])
combo_func.bind("<<ComboboxSelected>>", function_selected)
combo_func.pack()

# Delta t
frm_t = ttk.Labelframe(root, relief="ridge", text="Delta t", labelanchor="n")
frm_t.pack(side='left', fill=tk.Y)

lbl_tm = tk.Label(frm_t, text="Mantissa:")
lbl_tm.pack(side='left')
var_tm = tk.StringVar(root)
var_tm.set(str(tm))
spn_tm = tk.Spinbox(
    frm_t, textvariable=var_tm, format="%.1f", from_=-9., to=9., increment=0.1,
    command=lambda: set_tm(var_tm.get()), width=6
    )
spn_tm.pack(side='left')

lbl_te = tk.Label(frm_t, text="Exponent:")
lbl_te.pack(side='left')
var_te = tk.StringVar(root)
var_te.set(str(te))
spn_te = tk.Spinbox(
    frm_t, textvariable=var_te, format="%.1f", from_=-10., to=10., increment=1.,
    command=lambda: set_te(var_te.get()), width=6
    )
spn_te.pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
root.mainloop()

