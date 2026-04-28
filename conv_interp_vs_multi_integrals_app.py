import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Box1DKernel, CustomKernel, convolve
from astropy.modeling.models import Gaussian1D
from scipy.integrate import quad
from scipy.special import erf

st.set_page_config(page_title="PSF with IPSV Test", layout="wide")

st.title("PSF with IPSV Test")

st.sidebar.header("Settings")

# User chooses the mode
mode = st.sidebar.radio(
    "Kernel Type", options=["Top-Hat (width = 1 pixel)", "Gaussian IPSV"]
)
use_ipsv = mode == "Gaussian IPSV"

# General settings
st.sidebar.subheader("General Parameters")
sigma = st.sidebar.number_input("sigma (PSF)", value=0.5, step=0.1)
sample_step = 1.0  # Fixed as requested
dx = st.sidebar.number_input("dx", value=0.0001, step=0.0001, format="%.4f")
x_min = -3.5
x_max = 3.5

# Mode specific settings
st.sidebar.subheader("Kernel Parameters")
if use_ipsv:
    ipsv_sigma = st.sidebar.number_input("ipsv_sigma", value=0.8, step=0.1)
    ipsv_offset = st.sidebar.number_input("ipsv_offset", value=0.03, step=0.01)
else:
    tophat_width = 1.0

centroid_option = st.sidebar.radio(
    "Centroid", options=["Random (-0.5 to 0.5)", "Fixed (0.0)"]
)
if centroid_option == "Random (-0.5 to 0.5)":
    rng = np.random.default_rng()
    centroid = rng.uniform(-0.5, 0.5)
else:
    centroid = 0.0

st.write(f"**PSF centroid:** {centroid:.3f}")

# Main calculation logic
x = np.arange(x_min, x_max + dx, dx)
psf_model = Gaussian1D(amplitude=1.0, mean=centroid, stddev=sigma)
psf = psf_model(x)
half_pix = sample_step / 2.0

if use_ipsv:
    ipsv_pixel_norm = (
        ipsv_sigma
        * np.sqrt(np.pi / 2.0)
        * (
            erf((half_pix - ipsv_offset) / (ipsv_sigma * np.sqrt(2)))
            + erf((half_pix + ipsv_offset) / (ipsv_sigma * np.sqrt(2)))
        )
    )
    ipsv_amplitude = 1.0 / ipsv_pixel_norm

    ipsv_sigma_px = ipsv_sigma / dx
    kernel_half = int(6 * ipsv_sigma_px)
    kernel_x = np.arange(-kernel_half, kernel_half + 1) * dx

    kernel_array = ipsv_amplitude * np.exp(
        -0.5 * ((kernel_x + ipsv_offset) / ipsv_sigma) ** 2
    )
    valid_mask = (kernel_x >= -half_pix) & (kernel_x <= half_pix)
    kernel_array = kernel_array * valid_mask

    kernel = CustomKernel(kernel_array)
    kernel_label_short = "Gaussian IPSV"
    kernel_label = f"Gaussian IPSV (σ={ipsv_sigma}, offset={ipsv_offset})"
else:
    kernel_width_px = int(round(tophat_width / dx))
    if kernel_width_px % 2 == 0:
        kernel_width_px += 1
    kernel = Box1DKernel(kernel_width_px)
    kernel_label_short = "Tophat"
    kernel_label = f"Tophat (width={tophat_width})"

if use_ipsv:
    convolved_psf = convolve(psf, kernel, boundary="extend", normalize_kernel=False)
    convolved_psf = convolved_psf * dx
else:
    convolved_psf = convolve(psf, kernel, boundary="extend")

sample_x = np.arange(x_min, x_max + sample_step, sample_step)
sample_y = np.interp(sample_x, x, convolved_psf)

if use_ipsv:

    def pixel_integrand(t, cx):
        ipsv_model = Gaussian1D(
            amplitude=ipsv_amplitude, mean=cx + ipsv_offset, stddev=ipsv_sigma
        )
        return float(psf_model(t)) * float(ipsv_model(t))

else:

    def pixel_integrand(t, cx):
        return float(psf_model(t))


pixel_integrals = np.array(
    [
        quad(pixel_integrand, cx - half_pix, cx + half_pix, args=(cx,))[0]
        for cx in sample_x
    ]
)

fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)

ax.plot(
    x,
    psf,
    color="C0",
    label=f"Original σ={sigma} pix Gaussian PSF at centroid {centroid:.3f} pix",
)

ax.plot(
    x,
    convolved_psf,
    color="C1",
    ls="--",
    label=f"Gaussian PSF convolved with {kernel_label}",
)

ax.scatter(
    sample_x,
    sample_y,
    color="C1",
    s=40,
    zorder=4,
    label=f"sample_y (PSF * {kernel_label_short} interpolated at pixel centers)",
)

ax.scatter(
    sample_x,
    pixel_integrals,
    color="C3",
    marker="x",
    s=60,
    zorder=5,
    linewidths=1.8,
    label=f"pixel_integrals (PSF · {kernel_label_short} integrated over pixel footprint)",
)

for xs, ys in zip(sample_x, sample_y):
    ax.annotate(
        f"{ys:.4f}",
        (xs, ys),
        textcoords="offset points",
        xytext=(0, 6),
        ha="center",
        fontsize=7,
        color="C1",
    )

for i, (xs, pint) in enumerate(zip(sample_x, pixel_integrals)):
    ax.annotate(
        f"ø={pint:.4f}",
        (xs, 0),
        textcoords="offset points",
        xytext=(0, -14),
        ha="center",
        fontsize=7,
        color="C3",
    )

for i, xs in enumerate(sample_x):
    lo, hi = xs - half_pix, xs + half_pix
    mask = (x >= lo) & (x <= hi)
    x_fill = x[mask]
    y_fill = psf_model(x_fill)
    ax.fill_between(
        x_fill,
        0,
        y_fill,
        alpha=0.25,
        color=f"C{i % 9}",
    )
    ax.axvline(lo, color="gray", lw=0.5, ls="-", alpha=0.4)
ax.axvline(sample_x[-1] + half_pix, color="gray", lw=0.5, ls="-", alpha=0.4)

ax.axvline(centroid, color="gray", ls=":", lw=1.2, label=f"Centroid of Gaussian PSF")
ax.set_xlabel("X (pixels)")
ax.set_ylabel("Amplitude")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)

if use_ipsv:
    ax.set_title("Gaussian PSF + Gaussian IPSV", fontsize=12)
else:
    ax.set_title("Gaussian PSF + Tophat", fontsize=12)

st.pyplot(fig)
