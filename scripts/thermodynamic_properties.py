import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)



fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)


fig2, ax = plt.subplots()

#label axes

ax1.set_title("Translational KE")
ax2.set_title("Rotational KE")
ax3.set_title("Vibrational KE")
ax4.set_title("Potential KE")

ax.set_title("Total E")

#axes limits
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0, top=10E-20)

ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0, top=10E-20)

ax3.set_xlim(left=0)
ax3.set_ylim(bottom=0, top=10E-20)

ax4.set_xlim(left=0)
ax4.set_ylim(bottom=0, top=10E-20)

ax.set_xlim(left=0)
ax.set_ylim(bottom=0, top=10E-20)

#calculate energies
#translational

def kinetic_trans(path, i):

    data = pd.read_csv(path, sep=',',header=None, skiprows=1)

    u_1 = data.loc[i, 6]
    u_2 = data.loc[i, 7]

    v_1 = data.loc[i, 8]
    v_2 = data.loc[i, 9]

    w_1 = data.loc[i, 10]
    w_2 = data.loc[i, 11]

    m_1 = data.loc[i, 12]
    m_2 = data.loc[i, 13]

    m = m_1

    ke_t = m/4*((u_1 + u_2)**2 + (v_1 + v_2)**2 + (w_1 + w_2)**2)

    return ke_t

#vibrational

def kinetic_vibr(path, i):

    data = pd.read_csv(path, sep=',',header=None, skiprows=1)

    x_1 = data.loc[i, 0]
    x_2 = data.loc[i, 1]

    y_1 = data.loc[i, 2]
    y_2 = data.loc[i, 3]

    z_1 = data.loc[i, 4]
    z_2 = data.loc[i, 5]

    u_1 = data.loc[i, 6]
    u_2 = data.loc[i, 7]

    v_1 = data.loc[i, 8]
    v_2 = data.loc[i, 9]

    w_1 = data.loc[i, 10]
    w_2 = data.loc[i, 11]

    m_1 = data.loc[i, 12]
    m_2 = data.loc[i, 13]

    m = m_1

    num = (u_2 - u_1)/2*(x_2 - x_1) + (v_2 - v_1)/2*(y_2 - y_1) + (w_2 - w_1)/2*(z_2 - z_1)
    den = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)

    ke_vibr = m*((num/den)**2)

    #ke_vibr = m/4*((u_2 - u_1)**2 + (v_2 - v_1)**2 + (w_2 - w_1)**2)

    return ke_vibr

#potential

def potential(path, i):

    data = pd.read_csv(path, sep=',',header=None, skiprows=1)

    x_1 = data.loc[i, 0]
    x_2 = data.loc[i, 1]

    y_1 = data.loc[i, 2]
    y_2 = data.loc[i, 3]

    z_1 = data.loc[i, 4]
    z_2 = data.loc[i, 5]

    k = data.loc[i, 14]
    l_0 = data.loc[i, 15]

    pe = 0.5*k*(math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2) - l_0)**2

    return pe

#rotational

def kinetic_rot(path, i):

    data = pd.read_csv(path, sep=',',header=None, skiprows=1)

    x_1 = data.loc[i, 0]
    x_2 = data.loc[i, 1]

    y_1 = data.loc[i, 2]
    y_2 = data.loc[i, 3]

    z_1 = data.loc[i, 4]
    z_2 = data.loc[i, 5]

    u_1 = data.loc[i, 6]
    u_2 = data.loc[i, 7]

    v_1 = data.loc[i, 8]
    v_2 = data.loc[i, 9]

    w_1 = data.loc[i, 10]
    w_2 = data.loc[i, 11]

    m_1 = data.loc[i, 12]
    m_2 = data.loc[i, 13]

    m = m_1

    l = den = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)

    r_2 = 1/4*((x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2)

    I = 0.5*m*(l**2)

    #w_x = 4/(l**2)*((w_2 - w_1)*(y_2 - y_1) - (v_2 - v_1)*(z_2 - z_1))
    #w_y = 4/(l**2)*((u_2 - u_1)*(z_2 - z_1) - (w_2 - w_1)*(x_2 - x_1))
    #w_z = 4/(l**2)*((v_2 - v_1)*(x_2 - x_1) - (u_2 - u_1)*(y_2 - y_1))

    w_x = 1/r_2/4*((w_2 - w_1)*(y_2 - y_1) - (v_2 - v_1)*(z_2 - z_1))
    w_y = 1/r_2/4*((u_2 - u_1)*(z_2 - z_1) - (w_2 - w_1)*(x_2 - x_1))
    w_z = 1/r_2/4*((v_2 - v_1)*(x_2 - x_1) - (u_2 - u_1)*(y_2 - y_1))

    ke_rot = 0.5*I*((w_x**2) + (w_y**2) + (w_z**2))

    return ke_rot

#total mechanical energy

def total(path, i):

    return kinetic_rot(path, i) + potential(path, i) + kinetic_vibr(path, i) + kinetic_trans(path, i)


def kinetic_atom(path, i):

    data = pd.read_csv(path, sep=',',header=None, skiprows=1)

    u_1 = data.loc[i, 6]
    u_2 = data.loc[i, 7]

    v_1 = data.loc[i, 8]
    v_2 = data.loc[i, 9]

    w_1 = data.loc[i, 10]
    w_2 = data.loc[i, 11]

    m_1 = data.loc[i, 12]

    m = m_1

    ke_1 = 0.5*m*(u_1**2 + v_1**2 + w_1**2)
    ke_2 = 0.5*m*(u_2**2 + v_2**2 + w_2**2)

    return ke_1 + ke_2

#animation


def iterate(i):

    path = 'results/frame_' + str(i) + '.csv'

""" Test for i = 1 molecule """

def test_mol(frame, i):

    path = 'results/frame_' + str(frame) + '.csv'

    k_trans = kinetic_trans(path, i)
    k_rot = kinetic_rot(path, i)
    k_vibr = kinetic_vibr(path, i)

    p = potential(path, i)

    k_atom = kinetic_atom(path, i)

    k_total = k_trans +  k_vibr + k_rot

    print("Translational KE of " + str(i) + " = " + str(k_trans))
    print("Rotational KE of    " + str(i) + " = " + str(k_rot))
    print("Vibrational KE of   " + str(i) + " = " + str(k_vibr))

    print()

    print("Potential E of      " + str(i) + " = " + str(p))

    print()

    print("Total KE of         " + str(i) + " = " + str(k_total))

    print()
    print()

    print("Atomic KE of        " + str(i) + " = " + str(k_atom))

def average_energies_c(frame):

    path = 'results/frame_' + str(frame) + '.csv'

    data = pd.read_csv(path, sep=',',header=None, skiprows=1)

    rows = len(data)

    #define variables
    k_trans = 0
    k_rot = 0
    k_vibr = 0

    p = 0

    e_total = 0

    k_atom = 0

    #traverse list of molecules and add contribution
    for i in range(0, rows):

        k_trans += data.loc[i, 16]
        k_rot += data.loc[i, 17]
        k_vibr += data.loc[i, 18]

        p += data.loc[i, 19]

        e_total += data.loc[i, 20]

        #k_atom += kinetic_atom(path, i)

    #normalise by number of molecules
    k_trans = k_trans/rows
    k_rot = k_rot/rows
    k_vibr = k_vibr/rows

    p = p/rows

    e_total = e_total/rows

    #k_atom  = k_atom/rows

    return [k_trans, k_rot, k_vibr, p, e_total]


def average_energies(frame):

    path = 'results/frame_' + str(frame) + '.csv'

    #define variables
    k_trans = 0
    k_rot = 0
    k_vibr = 0

    p = 0

    e_total = 0

    k_atom = 0

    #traverse list of molecules and add contribution
    for i in range(0, 343):

        k_trans += kinetic_trans(path, i)
        k_rot += kinetic_rot(path, i)
        k_vibr += kinetic_vibr(path, i)

        p += potential(path, i)

        e_total += total(path, i)

        k_atom += kinetic_atom(path, i)

    #normalise by number of molecules
    k_trans = k_trans/343
    k_rot = k_rot/343
    k_vibr = k_vibr/343

    p = p/343

    e_total = e_total/343

    k_atom  = k_atom/343

    """
    #print results
    print("Average Translational KE of frame " + str(frame) + " = " + str(k_trans))
    print("Average Rotational KE of frame    " + str(frame) + " = " + str(k_rot))
    print("Average Vibrational KE of frame   " + str(frame) + " = " + str(k_vibr))

    print()

    print("Average Potential E of frame      " + str(frame) + " = " + str(p))

    print()

    print("Average Total E of frame          " + str(frame) + " = " + str(e_total))

    print()

    print("Average Total KE of frame         " + str(frame) + " = " + str(k_atom))

    """
    return [k_trans, k_rot, k_vibr, p, e_total]

#test_mol(14800, 34)
def average(lst):
    return sum(lst) / len(lst)

def plot_energies(start, end, step):

    #define energy vectors
    k_trans = []
    k_rot = []
    k_vibr = []
    p = []

    e_total = []

    time = range(start, end, step)

    ax1.set_xlim(left=0, right=end)
    ax1.set_ylim(bottom=0, top=10E-20)

    ax2.set_xlim(left=0, right=end)
    ax2.set_ylim(bottom=0, top=10E-20)

    ax3.set_xlim(left=0, right=end)
    ax3.set_ylim(bottom=0, top=10E-20)

    ax4.set_xlim(left=0, right=end)
    ax4.set_ylim(bottom=0, top=10E-20)

    ax.set_xlim(left=0, right=end)
    ax.set_ylim(bottom=0, top=10E-20)


    for frame in range(start, end, step):

        path = 'results/frame_' + str(frame) + '.csv'

        print("Time step = " + str(frame*1E-15))

        energies_step = average_energies_c(frame)

        k_trans.append(energies_step[0])
        k_rot.append(energies_step[1])
        k_vibr.append(energies_step[2])
        p.append(energies_step[3])
        e_total.append(energies_step[4])

    #best fit lines
    #m_t, b_t = np.polyfit(time[11000:], k_trans[11000:], 1)
    #m_r, b_r = np.polyfit(time[11000:], k_rot[11000:], 1)
    #m_v, b_v = np.polyfit(time[11000:], k_vibr[11000:], 1)
    #m_p, b_p = np.polyfit(time[11000:], p[11000:], 1)
    #m_et, b_et = np.polyfit(time[11000:], e_total[11000:], 1)


    ax1.plot(time, k_trans)
    #ax1.plot(time, k_rot, 'o')
    ax2.plot(time, [sum(x) for x in zip(k_rot, k_trans)], 'k')
    #ax1.plot(m_t, b_t, 'o')

    ax2.plot(time, k_rot)
    #ax2.plot(m_r, b_r, 'o')

    ax3.plot(time, k_vibr)
    #ax3.plot(m_v, b_v, 'o')

    ax4.plot(time, p)
    #ax4.plot(m_p, b_p, 'o')

    ax.plot(time, e_total)
    #ax.plot(m_et, b_et, 'o')

    #print across frame averages
    print("Translational Kinetic Energy: " + str(average(k_trans[5000:])))
    print("Rotational Kinetic Energy   : " + str(average(k_rot[5000:])))
    print("Vibrational Kinetic Energy  : " + str(average(k_vibr[5000:])))
    print("Potential Energy            : " + str(average(p[5000:])))

#average_energies(34640)

plot_energies(0, 1000000, 40)
plt.show()
