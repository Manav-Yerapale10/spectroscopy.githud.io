import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

st.set_page_config(page_title="2D NMR HSQC Simulator", layout="centered")
st.title("🧬 2D NMR HSQC Simulator for Protein-Ligand Interaction")

st.markdown("""
This tool allows you to input 1H and 15N chemical shifts and simulate an HSQC spectrum
**before and after ligand binding**. Ideal for protein-ligand interaction visualization.
""")

# Input fields
num_peaks = st.number_input("Number of peaks", min_value=1, max_value=30, value=5)
h_shifts = st.text_input("1H Chemical Shifts (comma-separated)", "7.1, 7.5, 8.0, 8.2, 8.5")
n_shifts = st.text_input("15N Chemical Shifts (comma-separated)", "105.3, 110.2, 115.7, 118.1, 120.3")
delta_h = st.text_input("Δ1H Shifts (optional, comma-separated)", "0.1, -0.05, 0.15, 0.0, -0.1")
delta_n = st.text_input("Δ15N Shifts (optional, comma-separated)", "0.4, 0.6, -0.3, 0.2, -0.1")

# Parse input
try:
    h_shifts = [float(x.strip()) for x in h_shifts.split(",")][:num_peaks]
    n_shifts = [float(x.strip()) for x in n_shifts.split(",")][:num_peaks]
    delta_h = [float(x.strip()) for x in delta_h.split(",")][:num_peaks]
    delta_n = [float(x.strip()) for x in delta_n.split(",")][:num_peaks]
except Exception as e:
    st.error("⚠️ Please check your inputs. All fields must contain comma-separated numbers.")
    st.stop()

# Define simulation
N = 128
h_axis = np.linspace(6, 10, N)
n_axis = np.linspace(90, 130, N)
H, N_grid = np.meshgrid(h_axis, n_axis)

def generate_spectrum(hs, ns):
    spectrum = np.zeros_like(H)
    for h, n in zip(hs, ns):
        peak = np.exp(-((H - h)**2 + (N_grid - n)**2) / 0.01)
        spectrum += peak
    return gaussian_filter(spectrum, sigma=1)

# Simulate spectra
spectrum_free = generate_spectrum(h_shifts, n_shifts)
h_bound = [h + dh for h, dh in zip(h_shifts, delta_h)]
n_bound = [n + dn for n, dn in zip(n_shifts, delta_n)]
spectrum_bound = generate_spectrum(h_bound, n_bound)

# Plot spectra
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
axs[0].imshow(spectrum_free, extent=[h_axis[0], h_axis[-1], n_axis[0], n_axis[-1]], origin='lower', cmap='viridis')
axs[0].set_title("Free Protein")
axs[0].set_xlabel("1H (ppm)")
axs[0].set_ylabel("15N (ppm)")

axs[1].imshow(spectrum_bound, extent=[h_axis[0], h_axis[-1], n_axis[0], n_axis[-1]], origin='lower', cmap='plasma')
axs[1].set_title("Protein + Ligand Bound")
axs[1].set_xlabel("1H (ppm)")
axs[1].set_ylabel("15N (ppm)")

st.pyplot(fig)

st.markdown("""
*Tip: You can simulate ligand binding by changing the Δ shifts.*
""")
