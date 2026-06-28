""" Quantum Mechanics, Path integral - 3rd """
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk

""" Global variables """

""" Animation control """
is_play = True

""" Other parameters """
h_base = 6.626   # Planck constant 6.626 * 10**-34 kg*m**2 / s
h_power = -34
h_bar = h_base * 10. ** h_power / (2. * np.pi)

mass_base = 50  # electron 9.10938 10**-31 kg
mass_power = -31
mass = mass_base * 10. ** mass_power

delta_x_base = 1.
delta_x_power = -1.
delta_x = delta_x_base * 10. ** delta_x_power

delta_t_base = 1.
delta_t_power = 1.
delta_t = delta_t_base * 10. ** delta_t_power

k = 0.
exponent = 2
num = 30

""" Create figure and axes """
title_ax0 = "Phases of action"
title_ax1 = "Time evolution"
title_tk = "Path integral"

x_min = -0.8
x_max = 0.8
y_min = -4.
y_max = 60.
z_min = -3.
z_max = 3.

fig = Figure()
# fig = Figure(facecolor='black')
ax0 = fig.add_subplot(121, projection='3d')
ax0.set_box_aspect((12, 12, 1))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel("x")
ax0.set_ylabel("Actions in Delta x")
ax0.set_zlabel("Imaginary part")
ax0.set_xlim(x_min, x_max)
ax0.set_ylim(y_min, y_max)
ax0.set_zlim(z_min, z_max)

# ax0.set_facecolor('black')
# ax0.axis('off')

ax1 = fig.add_subplot(122, projection='3d')
ax1.set_box_aspect((12, 12, 1))
ax1.grid()
ax1.set_title(title_ax1, color="black")
ax1.set_xlabel("x")
ax1.set_ylabel("t")
ax1.set_zlabel("Imaginary part")
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_zlim(z_min, z_max)

# ax1.set_facecolor('black')
# ax1.axis('off')

""" Embed in Tkinter """
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()


def on_move(event):
    if event.inaxes == ax0:
        ax1.view_init(elev=ax0.elev, azim=ax0.azim)
    elif event.inaxes == ax1:
        ax1.view_init(elev=ax1.elev, azim=ax1.azim)
    fig.canvas.draw_idle()


fig.canvas.mpl_connect('motion_notify_event', on_move)

""" Global objects of Tkinter """

""" Classes and functions """


class ActionSuperposed:
    def __init__(self, ax=None, dx=None, mass=None, dt=None, phase=None, x_start=None, x_end=None, label=None,
                 color=None):
        self.ax = ax
        self.dx = dx
        self.dt = dt
        self.mass = mass
        self.phase = phase
        self.x_start = x_start
        self.x_end = x_end
        self.label = label
        self.color = color
        self.radius = 0.8
        self.pitch = 0.001
        self.k = 0.
        self.exponent = 2

        # --- Helix (Superposed Propagator) ---
        self.action_phases = []
        for i in range(num):
            self.plt_helix, = self.ax.plot([], [], [], lw=0.5, color=self.color, ls='-', alpha=1.)
            self.action_phases.append(self.plt_helix)

        self.update_diagrams()

    def update_diagrams(self):
        # --- Update Helix (Propagator) ---
        x_position = np.linspace(self.x_start, self.x_end, 4000)

        # Pre-calculate the position of each slit/path
        x_start = - self.dx / 2.
        x_step = self.dx / (num - 1) if num > 1 else 0
        x_paths = x_start + np.arange(num) * x_step

        for i in range(num):
            helix_x = x_position
            sum_cos = np.zeros_like(x_position)
            sum_sin = np.zeros_like(x_position)

            # Superposition of phases from each path (numerical path integration)
            for j in range(num):
                x = x_paths[j]
                current_phase = self.k * j * 2. * np.pi / (num - 1) if num > 1 else 0

                # Calculate phase theta based on Action S
                # Action S = m * (x - x0)^exponent / (2 * dt)
                # Phase theta = S / h_bar
                theta = self.mass * (x_position - x) ** self.exponent / (2 * h_bar * self.dt * (i + 1)) + current_phase

                sum_cos += np.cos(theta)
                sum_sin += np.sin(theta)

            # --- Normalization and Offset Processing ---
            # Average the superposed waves and apply the y-axis offset (i * 2)
            adjust = np.sqrt(1 / (i + 1)) * 4.
            helix_y = (i * 2) + (self.radius * sum_cos / num) * adjust
            helix_z = (self.radius * sum_sin / num) * adjust

            self.action_phases[i].set_data_3d(helix_x, helix_y, helix_z)

    def set_mass(self, value):
        self.mass = value
        self.update_diagrams()

    def set_dx(self, value):
        self.dx = value
        self.update_diagrams()

    def set_dt(self, value):
        self.dt = value
        self.update_diagrams()

    def set_k(self, value):
        self.k = value
        self.update_diagrams()

    def set_exponent(self, value):
        self.exponent = value
        self.update_diagrams()


class ActionPhase:
    def __init__(self, ax=None, x=None, y=None, mass=None, dt=None, phase=None, x_start=None, x_end=None, label=None, color=None):
        self.ax = ax
        self.x = x
        self.y = y
        self.mass = mass
        self.dt = dt
        self.phase = phase
        self.origin_circle = np.array([self.x, self.y, 0.])  # Center of phase circle
        self.x_start = x_start
        self.x_end = x_end
        self.label = label
        self.color = color
        self.radius = 0.8
        self.pitch = 0.001

        self.exponent = 2

        # --- Circle (Phase Circle)---
        self.plt_circle, = self.ax.plot([], [], [], lw=1, color=self.color, ls='-', alpha=1.)

        # --- Quiver (Phase Arrow)---
        self.quiver_obj = None

        # --- Helix (Propagator)---
        self.plt_helix, = self.ax.plot([], [], [], lw=0.5, color=self.color, ls='-', alpha=1.)

        self.update_diagrams()

    def update_diagrams(self):
        # --- Update Phase Circle ---

        theta = np.linspace(0, 2 * np.pi, 64)

        # Center of circle
        cx = self.x
        cy = self.y
        cz = 0.0

        # Circle in y-z plane
        circle_x = np.full_like(theta, cx)
        circle_y = cy + self.radius * np.cos(theta)
        circle_z = cz + self.radius * np.sin(theta)

        self.plt_circle.set_data_3d(circle_x, circle_y, circle_z)

        # --- Update Quiver (Phase Arrow) ---
        if self.quiver_obj:
            self.quiver_obj.remove()

        self.quiver_obj = self.ax.quiver(
            self.origin_circle[0], self.origin_circle[1], self.origin_circle[2],
            0., np.cos(self.phase), np.sin(self.phase),
            length=self.radius, color=self.color, linewidth=1.5,
            arrow_length_ratio=0.1, normalize=False, ls='-'
        )

        # --- Update Helix (Propagator) ---
        x_position = np.linspace(self.x_start, self.x_end, 4000)

        # Action S = m (x - x0)^2 / (2 Δt)
        # Phase θ = m (x - x0)^2 / (2 ħ Δt)
        theta = self.mass * (x_position - self.x) ** self.exponent / (2 * h_bar * self.dt) + self.phase

        helix_x = x_position
        helix_y = self.y + self.radius * np.cos(theta)
        helix_z = self.radius * np.sin(theta)

        # 描画
        self.plt_helix.set_data_3d(helix_x, helix_y, helix_z)

    def set_mass(self, value):
        self.mass = value
        self.update_diagrams()

    def set_x(self, value):
        self.origin_circle[0] = value
        self.x = value
        self.update_diagrams()

    def set_dt(self, value):
        self.dt = value
        self.update_diagrams()

    def set_phase(self, value):
        self.phase = value
        self.update_diagrams()

    def set_exponent(self, value):
        self.exponent = value
        self.update_diagrams()


def set_mass_base(value):
    global mass, mass_base
    mass_base = value
    mass = mass_base * 10. ** mass_power
    for obj in action_phases:
        obj.set_mass(mass)
        obj.update_diagrams()

    action_superposed.set_mass(mass)
    action_superposed.update_diagrams()


def set_mass_power(value):
    global mass, mass_power
    mass_power = value
    mass = mass_base * 10. ** mass_power
    for obj in action_phases:
        obj.set_mass(mass)
        obj.update_diagrams()

    action_superposed.set_mass(mass)
    action_superposed.update_diagrams()


def set_delta_x_base(value):
    global delta_x, delta_x_base
    delta_x_base = value
    delta_x = delta_x_base * 10. ** delta_x_power
    i = 0
    for obj in action_phases:
        x_start = - delta_x / 2.
        x_step = delta_x / (num - 1)
        x = x_start + i * x_step
        obj.set_x(x)
        i += 1
        obj.update_diagrams()

    update_delta_x_guides()

    action_superposed.set_dx(delta_x)
    action_superposed.update_diagrams()


def set_delta_x_power(value):
    global delta_x, delta_x_power
    delta_x_power = value
    delta_x = delta_x_base * 10. ** delta_x_power
    i = 0
    for obj in action_phases:
        x_start = - delta_x / 2.
        x_step = delta_x / (num - 1)
        x = x_start + i * x_step
        obj.set_x(x)
        i += 1
        obj.update_diagrams()

    update_delta_x_guides()

    action_superposed.set_dx(delta_x)
    action_superposed.update_diagrams()


def set_delta_t_base(value):
    global delta_t, delta_t_base
    delta_t_base = value
    delta_t = delta_t_base * 10. ** delta_t_power
    i = 0
    for obj in action_phases:
        obj.set_dt(delta_t)
        i += 1
        obj.update_diagrams()

    update_diagrams()

    action_superposed.set_dt(delta_t)
    action_superposed.update_diagrams()


def set_delta_t_power(value):
    global delta_t, delta_t_power
    delta_t_power = value
    delta_t = delta_t_base * 10. ** delta_t_power
    i = 0
    for obj in action_phases:
        obj.set_dt(delta_t)
        i += 1
        obj.update_diagrams()

    update_diagrams()

    action_superposed.set_dt(delta_t)
    action_superposed.update_diagrams()


def set_k(value):
    global k
    k = value
    i = 0
    for obj in action_phases:
        phase = k * i * 2. * np.pi / (num - 1)
        obj.set_phase(phase)
        i += 1
        obj.update_diagrams()

    action_superposed.set_k(k)
    action_superposed.update_diagrams()


def set_exponent(value):
    global exponent
    exponent = value
    i = 0
    for obj in action_phases:
        obj.set_exponent(exponent)
        i += 1
        obj.update_diagrams()

    action_superposed.set_exponent(exponent)
    action_superposed.update_diagrams()


def update_delta_x_guides():
    x_left = -delta_x / 2.
    x_right = delta_x / 2.

    y_vals = np.array([y_min, y_max])
    z_vals = np.array([0., 0.])

    guide_left.set_data_3d(
        np.full_like(y_vals, x_left), y_vals, z_vals
    )
    guide_right.set_data_3d(
        np.full_like(y_vals, x_right), y_vals, z_vals
    )


def create_parameter_setter():
    # Mass
    frm_mass = ttk.Labelframe(root, relief="ridge", text="Mass", labelanchor="n", width=100)
    frm_mass.pack(side='left')

    tk.Label(frm_mass, text="Base").pack(side='left')
    var_mass_base = tk.StringVar(root)  # variable for spinbox-value
    var_mass_base.set(str(mass_base))  # Initial value

    spn_mass_base = tk.Spinbox(
        frm_mass, textvariable=var_mass_base, format="%.1f", from_=-100., to=100., increment=1.,
        command=lambda: set_mass_base(float(var_mass_base.get())), width=8
    )
    spn_mass_base.pack(side='left')

    tk.Label(frm_mass, text="Power").pack(side='left')
    var_mass_power = tk.StringVar(root)  # variable for spinbox-value
    var_mass_power.set(str(mass_power))  # Initial value

    spn_mass_power = tk.Spinbox(
        frm_mass, textvariable=var_mass_power, format="%.1f", from_=-10., to=10., increment=0.1,
        command=lambda: set_mass_power(float(var_mass_power.get())), width=6
    )
    spn_mass_power.pack(side='left')

    # Delta x
    frm_delta_x = ttk.Labelframe(root, relief="ridge", text="Delta x", labelanchor="n", width=100)
    frm_delta_x.pack(side='left')

    tk.Label(frm_delta_x, text="Base").pack(side='left')
    var_delta_x_base = tk.StringVar(root)  # variable for spinbox-value
    var_delta_x_base.set(str(delta_x_base))  # Initial value

    spn_delta_x_base = tk.Spinbox(
        frm_delta_x, textvariable=var_delta_x_base, format="%.1f", from_=-10., to=10., increment=0.1,
        command=lambda: set_delta_x_base(float(var_delta_x_base.get())), width=6
    )
    spn_delta_x_base.pack(side='left')

    tk.Label(frm_delta_x, text="Power").pack(side='left')
    var_delta_x_power = tk.StringVar(root)  # variable for spinbox-value
    var_delta_x_power.set(str(delta_x_power))  # Initial value

    spn_delta_x_power = tk.Spinbox(
        frm_delta_x, textvariable=var_delta_x_power, format="%.1f", from_=-10., to=10., increment=0.1,
        command=lambda: set_delta_x_power(float(var_delta_x_power.get())), width=6
    )
    spn_delta_x_power.pack(side='left')

    # k (Wave Number)
    frm_k = ttk.Labelframe(root, relief="ridge", text="k(Wave number)", labelanchor="n", width=100)
    frm_k.pack(side='left')

    tk.Label(frm_k, text="k").pack(side='left')
    var_k = tk.StringVar(root)  # variable for spinbox-value
    var_k.set(str(k))  # Initial value

    spn_k = tk.Spinbox(
        frm_k, textvariable=var_k, format="%.1f", from_=-10., to=10., increment=0.5,
        command=lambda: set_k(float(var_k.get())), width=6
    )
    spn_k.pack(side='left')

    # Delta t
    frm_delta_t = ttk.Labelframe(root, relief="ridge", text="Delta t", labelanchor="n", width=100)
    frm_delta_t.pack(side='left')

    tk.Label(frm_delta_t, text="Base").pack(side='left')
    var_delta_t_base = tk.StringVar(root)  # variable for spinbox-value
    var_delta_t_base.set(str(delta_t_base))  # Initial value

    spn_delta_t_base = tk.Spinbox(
        frm_delta_t, textvariable=var_delta_t_base, format="%.1f", from_=-10., to=10., increment=0.1,
        command=lambda: set_delta_t_base(float(var_delta_t_base.get())), width=6
    )
    spn_delta_t_base.pack(side='left')

    tk.Label(frm_delta_t, text="Power").pack(side='left')
    var_delta_t_power = tk.StringVar(root)  # variable for spinbox-value
    var_delta_t_power.set(str(delta_t_power))  # Initial value

    spn_delta_t_power = tk.Spinbox(
        frm_delta_t, textvariable=var_delta_t_power, format="%.1f", from_=-10., to=10., increment=0.1,
        command=lambda: set_delta_t_power(float(var_delta_t_power.get())), width=6
    )
    spn_delta_t_power.pack(side='left')

    # Exponent
    frm_exponent = ttk.Labelframe(root, relief="ridge", text="Exponent", labelanchor="n", width=100)
    frm_exponent.pack(side='left')

    tk.Label(frm_exponent, text="").pack(side='left')
    var_exponent = tk.StringVar(root)  # variable for spinbox-value
    var_exponent.set(str(exponent))  # Initial value

    spn_exponent = tk.Spinbox(
        frm_exponent, textvariable=var_exponent, format="%.1f", from_=1, to=2, increment=1,
        command=lambda: set_exponent(float(var_exponent.get())), width=4
    )
    spn_exponent.pack(side='left')


def draw_static_diagrams():
    pass


def update_diagrams():
    for obj in action_phases:
        obj.update_diagrams()

    action_superposed.update_diagrams()


def update(f):
    if is_play:
        pass
        # update_diagrams()


""" main loop """
if __name__ == "__main__":
    draw_static_diagrams()
    create_parameter_setter()

    # Δx Guide line
    guide_left, = ax0.plot([], [], [], color="blue", lw=1)
    guide_right, = ax0.plot([], [], [], color="blue", lw=1)
    update_delta_x_guides()

    action_phases = []
    for i_ in range(num):
        x_start = - delta_x / 2.
        x_step = delta_x / (num - 1)
        x = x_start + i_ * x_step
        action_phase = ActionPhase(ax0, x, i_ * 2., mass, delta_t, 0., x_min, x_max, "", "red")
        action_phases.append(action_phase)

    action_superposed = ActionSuperposed(ax1, delta_x, mass, delta_t, 0., x_min, x_max, "", "red")

    anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
    root.mainloop()