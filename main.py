import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.widgets import TextBox

from bisect import bisect

m = 1
mm = 1e-3
um = 1e-6
nm = 1e-9



def length_value_unit(length_in_meters):
    factor = 1
    if length_in_meters >= factor*m:
        return((length_in_meters, 'm'))
    if length_in_meters >= factor*mm and length_in_meters < factor*m:
        return((length_in_meters / mm, 'mm'))
    if length_in_meters >= factor*um and length_in_meters < factor*mm:
        return((length_in_meters / um, 'um'))
    if length_in_meters < factor*um:
        return(length_in_meters / nm, 'nm')

def parameters_text(z):
    z0, w0, zr, NA, waist, curvature = beam_parameters.properties(z)
    z, z_unit = length_value_unit(z)
    z0, z0_unit = length_value_unit(z0)
    w0, w0_unit = length_value_unit(w0)
    zr, zr_unit = length_value_unit(zr)
    waist, waist_unit = length_value_unit(waist)
    curvature, curvature_unit = curvature * m, 'm-1'
    string = "z = {:.5g} {}\n\nz_center = {:.3g} {}\n\nw0 = {:.3g} {}\n\nzr = {:.3g} {}\n\nNA = {:.3g}\n\nWaist = {:.3g} {}\n\nCurvature = {:.3g} {}".format(z, z_unit, z0, z0_unit, w0, w0_unit, zr, zr_unit, NA, waist, waist_unit, curvature, curvature_unit)
    return(string)

def w(z, w0, z0, wavelength):
    zr = np.pi*w0**2/wavelength
    return(w0*np.sqrt(1+((z-z0)/zr)**2))



class BeamPropagation:
    def __init__(self, propagation_distance, beam_waist, beam_center_position, wavelength, lens_focal_lengths = [], lens_positions = []):
        self.L = propagation_distance
        self.w0 = beam_waist
        self.z0 = beam_center_position
        self.wavelength = wavelength

        self.update_lenses(lens_focal_lengths, lens_positions)


    
    def update_lenses(self, lens_focal_lengths, lens_positions):

        self.lens_z, self.lens_f = zip(*sorted(zip(lens_positions, lens_focal_lengths)))

        if len(self.lens_f) != len(self.lens_z):
            raise Exception("The lengths of the lens positions list and the lens focal lengths list must be equal")

        lens_delta_z = [self.lens_z[0] if i == 0 else self.lens_z[i] - self.lens_z[i-1] for i in range(len(self.lens_z))]
        ABCD_matrices = [np.dot(np.array([[1, 0], [-1/self.lens_f[i], 1]]), np.array([[1, lens_delta_z[i]], [0, 1]])) for i in range(len(self.lens_z))]

        zr = np.pi*self.w0**2/self.wavelength
        Q = [-self.z0 + 1j*zr]
        self.beam_positions = [self.z0]
        self.beam_waists = [self.w0]

        for i in range(len(self.lens_z)):
            q_in = Q[-1]
            ABCD_matrix = ABCD_matrices[i]
            A, B, C, D = ABCD_matrix[0, 0], ABCD_matrix[0, 1], ABCD_matrix[1, 0], ABCD_matrix[1, 1]
            q_out = (q_in*A + B)/(q_in*C + D)
            Q.append(q_out)
            self.beam_positions.append(self.lens_z[i] - np.real(q_out))
            self.beam_waists.append(np.sqrt(self.wavelength*np.imag(q_out)/np.pi))

    def properties(self, z):
        i = bisect(self.lens_z, z)
        z0 = self.beam_positions[i]
        w0 = self.beam_waists[i]
        zr = np.pi*w0**2/self.wavelength
        NA = w0/zr
        waist = w(z, w0, z0, self.wavelength)
        curvature = (z-z0)/((z-z0)**2 + zr**2)
        return(z0, w0, zr, NA, waist, curvature)

    def __call__(self, z):
        i = bisect(self.lens_z, z)
        return(w(z, self.beam_waists[i], self.beam_positions[i], self.wavelength))



#########################################################
#################### Parameters #########################

L = 3 * m             # Simulation range
w0 = 500 * um            # Initial beam waist
z0 = 0 * mm             # Initial beam position
wavelength = 810 * nm   # Wavelength

lens_f = [150*mm, 300*mm, 600*mm]         # List of focal lengths of your lenses
lens_z = [150*mm, 600*mm, 1500*mm]       # List of positions of your lenses

nb_points = 1024                # Number of points for graph

#########################################################


# Compute and plot

beam_parameters = BeamPropagation(L, w0, z0, wavelength, lens_f, lens_z)
beam_parameters_vect = np.vectorize(beam_parameters)


nb_of_lenses = len(lens_f)




Z = np.linspace(0, L, nb_points)
radius_array = beam_parameters_vect(Z)
max_radius = np.max(radius_array)

y_unit = mm

fig, ax = plt.subplots()
plt.subplots_adjust(bottom = 0.05*nb_of_lenses + 0.2)
plt.subplots_adjust(right = 0.7)
ax.set_xlim([0, L / m])
ax.set_ylim([-1.2*max_radius / y_unit, 1.2*max_radius / y_unit])
ax.set_facecolor([255/255, 255/255, 255/255])

radius_plus, = ax.plot(Z / m, radius_array / y_unit, color = 'r')
radius_minus, = ax.plot(Z / m, -radius_array / y_unit, color = 'r')
fill = ax.fill_between(Z, -radius_array / y_unit, radius_array / y_unit, color = [255/255, 32/255, 0], alpha = 0.5)

ax_limits = ax.get_position()
text = fig.text(0.72, (ax_limits.y1 + ax_limits.y0)/2, parameters_text(0), horizontalalignment="left", verticalalignment="center")

measurement_bar = ax.axvline(0, color = [120/255, 198/255, 121/255], alpha = 0.8)


cursor_position = 0

def onclick(event):
    global cursor_position
    if event.inaxes is ax and event.button==1:
        if type(event.xdata) == np.float64:
            cursor_position = event.xdata
            text.set_text(parameters_text(cursor_position))
            measurement_bar.set_xdata([cursor_position])
            fig.canvas.draw()
            


cid = fig.canvas.mpl_connect('motion_notify_event', onclick)
cid = fig.canvas.mpl_connect('button_press_event', onclick)



test = []
lens_slider = []
lens_text = []
text_box = []


def update_function(index):
    def update_lens(val):
        global max_radius
        global fill
        global cursor_position

        lens_position = val

        lens_z[index] = lens_position
        beam_parameters.update_lenses(lens_f, lens_z)
        radius_array = beam_parameters_vect(Z)
        if np.max(radius_array) > max_radius:
            max_radius = np.max(radius_array)
            ax.set_ylim([-1.2*max_radius / y_unit, 1.2*max_radius / y_unit])
        radius_plus.set_ydata(radius_array / y_unit)
        radius_minus.set_ydata(-radius_array / y_unit)

        text.set_text(parameters_text(cursor_position))

        fill.remove()
        fill = ax.fill_between(Z, -radius_array / y_unit, radius_array / y_unit, color = [255/255, 32/255, 0], alpha = 0.5)

        text_box[index].set_val('{:.4g}'.format(val))

        arrows_up.set_offsets([list(i) for i in zip(*[lens_z, [0]*nb_of_lenses])])
        arrows_up.set_UVC([0]*nb_of_lenses, [max_radius / y_unit]*nb_of_lenses)
        arrows_down.set_offsets([list(i) for i in zip(*[lens_z, [0]*nb_of_lenses])])
        arrows_down.set_UVC([0]*nb_of_lenses, [-max_radius / y_unit]*nb_of_lenses)

        for i in range(len(lens_z)):
            lens_text[i].set_position((lens_z[i], max_radius / y_unit))
        

    return(update_lens)

def submit_function(index):
    def submit(text):
        try:
            new_val = float(text)
        except:
            return
        if new_val < 0:
            new_val = 0
        if new_val > L:
            new_val = L
        lens_slider[index].set_val(new_val)
    return(submit)

arrows_up = ax.quiver(lens_z, [0]*nb_of_lenses, [0]*nb_of_lenses, [max_radius / y_unit]*nb_of_lenses, scale_units='xy', scale=1)
arrows_down = ax.quiver(lens_z, [0]*nb_of_lenses, [0]*nb_of_lenses, [-max_radius / y_unit]*nb_of_lenses, scale_units='xy', scale=1)


for index, lens_position in enumerate(lens_z):
    ax_lens_slider = plt.axes([0.125, 0.05*(nb_of_lenses - index), 0.575, 0.03], facecolor = 'white')
    lens_slider.append(Slider(ax_lens_slider, "Lens {}".format(index+1), 0, L, valinit = lens_position))
    lens_slider[-1].valtext.set_visible(False)
    lens_slider[-1].on_changed(update_function(index))

    axbox = fig.add_axes([0.71, 0.05*(nb_of_lenses - index), 0.1, 0.04])
    text_box.append(TextBox(axbox, '', initial='{:.4g}'.format(lens_position), textalignment="left"))
    text_box[-1].on_submit(submit_function(index))

    fig.text(0.82, 0.05*(nb_of_lenses - index) + 0.025,  'm', verticalalignment="center")

    lens_text.append(ax.text(lens_position, max_radius / y_unit, "f{} = {:.0f} mm".format(index+1, lens_f[index] / y_unit), horizontalalignment='center'))



ax.set_xlabel('z-axis position (m)')
ax.set_ylabel('Distance from optical axis (mm)')

plt.show()