import streamlit as st
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

import mpmath as mp
mp.dps = 30  # Set decimal places for high precision calculations

#from matplotlib import rc
#rc('text', usetex=True)
plt.rcParams.update({'xtick.labelsize': 12, 
                     'ytick.labelsize': 12,
                     'axes.labelsize': 16})#,
                     #'font.family': 'cmr10'})

def P_3(theta):
    return 5 * mp.cos(theta)**2 - 1

def P_5(theta):
    return 21 * mp.cos(theta)**4 - 14 * mp.cos(theta)**2 + 1

def omega_tac(r, theta, rd1, rd3, w1, w3, deltaOmega3, deltaOmega5, B1, B3, B5, C, Omega_c):
    B = B1 + B3 * P_3(theta) + B5 * P_5(theta)
    deltaOmega = deltaOmega3 * P_3(theta) + deltaOmega5 * P_5(theta)
    rt = rd1 + rd3 * P_3(theta)
    w = w1 + w3 * P_3(theta)
    exp = (1 + mp.exp((rt - r) / w))
    if r <= 0.70:
        return Omega_c + deltaOmega / exp
    elif 0.70 < r <= 0.95:
        return Omega_c + B * (r - 0.7) + deltaOmega / exp
    else:
        return Omega_c + 0.25 * B - C * (r - 0.95) + deltaOmega / exp

@st.cache_data
def get_grid():
    r_values = np.linspace(0.5, 1, 100)
    theta_values = np.linspace(0, np.pi, 100)
    selected_thetas = np.deg2rad([0, 15, 30, 45, 60])
    selected_rs = [0.7, 0.75, 0.8, 0.85, 0.98]
    return r_values, theta_values, selected_thetas, selected_rs

def plot_interactive(rd1, rd3, w1, w3, deltaOmega3, deltaOmega5, B1, B3, B5, C, Omega_c):
    r_values, theta_values, selected_thetas, selected_rs = get_grid()
        
    fig, axs = plt.subplots(3, 1, figsize=(10, 10))
    
    for theta in selected_thetas:
        omega_r = [omega_tac(r, theta, rd1, rd3, w1, w3, deltaOmega3, deltaOmega5, B1, B3, B5, C, Omega_c) for r in r_values]
        axs[0].plot(r_values, omega_r, lw=3, label=f'${int(round(np.rad2deg(theta), 0))}^\circ$')

    rt = rd1 + rd3 * P_3(theta)
    axs[0].axvline(rt, ls='--', c='k', lw=1.5, label=r'$r_t$')
    axs[0].set_xlim([0.5, 1])
    axs[0].set_xlabel(r'$r/\rm{R}_\odot$')
    axs[0].set_ylabel(r'$\Omega$ [nHz]')
    axs[0].legend(fontsize=14)
    
    for r in selected_rs:
        omega_theta = [omega_tac(r, theta, rd1, rd3, w1, w3, deltaOmega3, deltaOmega5, B1, B3, B5, C, Omega_c) for theta in theta_values]
        axs[1].plot(np.rad2deg(theta_values), omega_theta, lw=3, label=f'${r}' + r'~\rm{R}_\odot$')
    axs[1].set_xlim([0, 180])
    axs[1].set_xlabel(r'co-latitude $\theta$ [$^{\circ}$]')
    axs[1].set_ylabel(r'$\Omega$ [nHz]')
    axs[1].legend(fontsize=14)
    
    R, Theta = np.meshgrid(r_values, np.rad2deg(theta_values))
    Omega = np.vectorize(lambda r, theta: float(omega_tac(r, np.deg2rad(theta), rd1, rd3, w1, w3, deltaOmega3, deltaOmega5, B1, B3, B5, C, Omega_c)))(R, Theta)
    pcm = axs[2].pcolormesh(R, Theta, Omega, shading='auto', cmap='coolwarm')
    plt.colorbar(pcm, ax=axs[2], label=r'$\Omega$ [nHz]')
    axs[2].axvline(rt, ls='--', c='k', lw=1.5, label=r'$r_t$')
    axs[2].set_xlabel(r'$r/\rm{R}_\odot$')
    axs[2].set_ylabel(r'co-latitude $\theta$ [$^{\circ}$]')
    axs[2].legend(loc='upper left', fontsize=14)
    
    plt.tight_layout()
    st.pyplot(fig)


st.title("Solar Rotation Profile")

rd1 = st.sidebar.slider(r"$r_{d1}$", min_value=0.6965 - 0.0026, max_value=0.7034 + 0.0024, value=0.6978, step=0.0001)
rd3 = st.sidebar.slider(r"$r_{d3}$", min_value=0.0029 - 0.0010, max_value=0.0137 + 0.0019, value=0.0033, step=0.0001)
w1 = st.sidebar.slider(r"$w_1$", min_value=0.0063 - 0.0032, max_value=0.0093 + 0.0020, value=0.0071, step=0.0001)
w3 = st.sidebar.slider(r"$w_3$", min_value=0.0048 - 0.0017, max_value=0.0112 + 0.0033, value=0.0052, step=0.0001)
deltaOmega3 = st.sidebar.slider(r"$\delta\Omega_3$", min_value=-22.75 - 0.28, max_value=-21.20 + 0.24, value=-22.28, step=0.01)
deltaOmega5 = st.sidebar.slider(r"$\delta\Omega_5$", min_value=-3.74 - 0.11, max_value=-2.42 + 0.12, value=-3.52, step=0.01)
B1 = st.sidebar.slider(r"$B_1$", min_value=-50, max_value=50, value=0, step=1)
B3 = st.sidebar.slider(r"$B_3$", min_value=-50, max_value=50, value=0, step=1)
B5 = st.sidebar.slider(r"$B_5$", min_value=-50, max_value=50, value=0, step=1)
C = st.sidebar.slider(r"$C$", min_value=-50, max_value=50, value=0, step=1)
Omega_c = st.sidebar.slider(r"$\Omega_c$", min_value=360.0, max_value=480.0, value=435.0, step=0.1)

plot_interactive(rd1, rd3, w1, w3, deltaOmega3, deltaOmega5, B1, B3, B5, C, Omega_c)

st.markdown(
    "See [Antia & Basu 2011](https://ui.adsabs.harvard.edu/abs/2011ApJ...735L..45A/abstract) for more information.",
    unsafe_allow_html=True
)
