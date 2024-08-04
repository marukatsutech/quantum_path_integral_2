# Quantum Mechanics, Path integral in potential

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from mpl_toolkits.mplot3d import proj3d


def step():
    global cnt
    apply_path_integral()
    cnt += 1


def apply_path_integral():
    global qtm0a
    global ya, za, plt_qtm0a
    qtm0a_buffer = qtm0a * 0.
    # Quantum A
    for i in range(len(x)):
        dx_squared = (x[i] - x) ** 2
        if t != 0.:
            theta = np.mod(2. * np.pi * mass * dx_squared / (2. * h * t), (2. * np.pi))
        else:
            theta = 0.
            print('Error: divided by zero')
        rot = np.sin(theta) + np.cos(theta) * 1j
        qtm0a_buffer = qtm0a_buffer + qtm0a[i] * rot
    qtm0a_buffer = qtm0a_buffer * rot_potential
    qtm0a = qtm0a_buffer / np.sum(np.abs(qtm0a_buffer))
    ya = qtm0a.imag
    za = qtm0a.real
    plt_qtm0a.set_xdata(x)
    plt_qtm0a.set_ydata(ya)
    plt_qtm0a.set_3d_properties(za)


def update_qtm0a():
    global gaussian0a, amplitude0a
    global z0a, y0a, y0a, qtm0a
    global ya, za, plt_qtm0a
    gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x - mu0a) ** 2 / (2 * sigma0a ** 2))
    amplitude0a = gaussian0a / np.sum(np.abs(gaussian0a))
    z0a = np.sin(k0a * x) * amplitude0a
    y0a = np.cos(k0a * x) * 1j * amplitude0a
    qtm0a = z0a + y0a
    ya = qtm0a.imag
    za = qtm0a.real
    plt_qtm0a.set_xdata(x)
    plt_qtm0a.set_ydata(ya)
    plt_qtm0a.set_3d_properties(za)


def set_k0a(value):
    global k0a
    k0a = float(value)
    update_qtm0a()


def set_x0a(value):
    global mu0a
    mu0a = float(value)
    update_qtm0a()


def set_sigma0a(value):
    global sigma0a
    sigma0a = float(value)
    update_qtm0a()


def set_tm(value):
    global tm, t
    tm = float(value)
    t = tm * np.power(10., te)


def set_te(value):
    global te, t
    te = float(value)
    t = tm * np.power(10., te)


def function_selected(event):
    global theta_potential, rot_potential, plt_potential
    if combo_func.get() == option_func[0]:
        ax1.set_title(title_ax1 + option_func[0])
        theta_potential = 0. * x
    elif combo_func.get() == option_func[1]:
        ax1.set_title(title_ax1 + option_func[1])
        theta_potential = + 0.0001 * (x ** 2)
    elif combo_func.get() == option_func[2]:
        ax1.set_title(title_ax1 + option_func[2])
        theta_potential = - 0.0001 * (x ** 2)
    elif combo_func.get() == option_func[3]:
        ax1.set_title(title_ax1 + option_func[3])
        theta_potential = + 0.1 * x
    elif combo_func.get() == option_func[4]:
        ax1.set_title(title_ax1 + option_func[4])
        theta_potential = - 0.1 * x
    elif combo_func.get() == option_func[5]:
        ax1.set_title(title_ax1 + option_func[5])
        sigma = 500.
        theta_potential = + 50000. * 1. / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (x - 0.) ** 2 / (2 * sigma ** 2))
    else:
        ax1.set_title(title_ax1 + option_func[6])
        sigma = 500.
        theta_potential = - 50000. * 1. / (np.sqrt(2 * np.pi) * sigma) * np.exp(- (x - 0.) ** 2 / (2 * sigma ** 2))
    rot_potential = np.sin(theta_potential) + np.cos(theta_potential) * 1j
    plt_potential.set_ydata(theta_potential)


# Animation control
def reset():
    global is_play, cnt
    global cnt_step
    global k0a, sigma0a, mu0a
    is_play = False
    cnt = 0
    k0a = k0a_init
    sigma0a = sigma0a_init
    mu0a = mu0a_init
    update_qtm0a()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt
    txt_step.set_text("Delta t * " + str(cnt))
    if is_play:
        cnt += 1
        step()


# Global variables
# Animation control
cnt = 0
is_play = False
cnt_step = 1

# Data structure
range_x = 1000
x = np.arange(- range_x, range_x)

range_quantum_yz = 0.02
range_potential = 100

# Quantum parameter
mass = 9.1093837015 * 1.0E-31   # Unit: kg
tm, te = 1., 7
t = tm * np.power(10., te)  # Unit: s
h = 6.62607015 * 1.0E-34    # Unit: m2 kg / s

# Quantum
# Quantum 0a
k0a_init = 0.0
sigma0a_init = 20.
mu0a_init = - range_x / 2.

k0a = k0a_init
sigma0a = sigma0a_init
mu0a = mu0a_init
gaussian0a = 1 / (np.sqrt(2 * np.pi) * sigma0a) * np.exp(- (x - mu0a) ** 2 / (2 * sigma0a ** 2))
amplitude0a = gaussian0a / np.sum(np.abs(gaussian0a))
z0a = np.sin(k0a * x) * amplitude0a
y0a = np.cos(k0a * x) * 1j * amplitude0a
qtm0a = z0a + y0a

# Potential
theta_potential = 0. * x
rot_potential = np.sin(theta_potential) + np.cos(theta_potential) * 1j

# Generate figure and axes
title_tk = "Quantum Mechanics, Path integral in potential"
title_ax0 = "Quantum"
title_ax1 = "Potential:"

x_min0 = - range_x
x_max0 = range_x
y_min0 = - range_quantum_yz
y_max0 = range_quantum_yz
z_min0 = - range_quantum_yz
z_max0 = range_quantum_yz

x_min1 = - range_x
x_max1 = range_x
y_min1 = - range_potential
y_max1 = range_potential

fig = Figure()
ax0 = fig.add_subplot(121, projection='3d')
ax0.set_box_aspect((6, 2, 2))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel('x')
ax0.set_ylabel('y')
ax0.set_zlabel('z')
ax0.set_xlim(x_min0, x_max0)
ax0.set_ylim(y_min0, y_max0)
ax0.set_zlim(z_min0, z_max0)

ax1 = fig.add_subplot(122)
ax1.grid()
ax1.set_title(title_ax1)
ax1.set_xlabel('x')
ax1.set_ylabel('Additional phase')
ax1.set_xlim(x_min1, x_max1)
ax1.set_ylim(y_min1, y_max1)

# Text items
txt_step = ax0.text2D(x_min0, y_max0, "Delta t * " + str(cnt))
xz, yz, _ = proj3d.proj_transform(x_min0, y_max0, z_max0, ax0.get_proj())
txt_step.set_position((xz, yz))

# Plot items
# Quantum 0a
ya = qtm0a.imag
za = qtm0a.real
plt_qtm0a, = ax0.plot(x, ya, za, color='blue', linewidth=0.5, label='Quantum A')

# Potential
plt_potential, = ax1.plot(x, theta_potential, label='Potential')


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
btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
btn_play.pack(side='left')
btn_step = tk.Button(frm_anim, text="Step", command=step)
btn_step.pack(side='left')
btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
btn_reset.pack(side='left')

# Quantum A
frm_qtm_a = ttk.Labelframe(root, relief="ridge", text="Quantum A (initial)", labelanchor="n")
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

# Potential function option
frm_func = ttk.Labelframe(root, relief='ridge', text='Potential option', labelanchor='n')
frm_func.pack(side='left', fill=tk.Y)
option_func = ["None", "+0.0001*(x**2)", "-0.0001*(x**2)", "+0.1*x", "-0.1*x", "+30000*Gaussian(sigma:500)",
               "-30000*Gaussian(sigma:500)"]
variable = tk.StringVar()
combo_func = ttk.Combobox(frm_func, values=option_func, textvariable=variable, width=30)
combo_func.set(option_func[0])
combo_func.bind("<<ComboboxSelected>>", function_selected)
combo_func.pack()

# main loop
anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
root.mainloop()

