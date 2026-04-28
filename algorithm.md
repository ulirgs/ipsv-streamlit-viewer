# PSF Sampling Methods: Tophat vs Gaussian IPSV

This document describes the two numerical integration algorithms used in the tool to calculate the simulated pixel values from a continuous Point Spread Function (PSF). 

Let the continuous 1D PSF be denoted as $\text{PSF}(x)$, and the size of a single pixel be $\Delta x$ (`sample_step` in the code, which is set to 1.0). The center of a given pixel is denoted as $x_c$.

## 1. Tophat (Ideal Pixel) Method

The Tophat method assumes that a pixel has uniform sensitivity across its entire surface and drops sharply to zero outside its boundaries. 

The simulated intensity $I(x_c)$ measured by the pixel centered at $x_c$ is calculated by integrating the continuous $\text{PSF}(x)$ over the physical footprint of the pixel, from $x_c - \frac{\Delta x}{2}$ to $x_c + \frac{\Delta x}{2}$:

$$
I(x_c) = \int_{x_c - \frac{\Delta x}{2}}^{x_c + \frac{\Delta x}{2}} \text{PSF}(x) \, dx
$$

This is the standard, ideal formulation for pixel sampling without any hardware non-idealities.

---

## 2. Gaussian IPSV (Intra-Pixel Sensitivity Variation) Method

In real detectors, the sensitivity of a pixel is often not perfectly uniform; it typically peaks near the center of the pixel and falls off towards the edges. This phenomenon is known as Intra-Pixel Sensitivity Variation (IPSV). 

Here, the IPSV is modeled as a Gaussian function. The simulated intensity $I(x_c)$ is the integral of the PSF *weighted* by the pixel's local sensitivity response:

$$
I(x_c) = \int_{x_c - \frac{\Delta x}{2}}^{x_c + \frac{\Delta x}{2}} \text{PSF}(x) \cdot \text{IPSV}(x; x_c) \, dx
$$

### The IPSV Kernel

The sensitivity profile for the pixel centered at $x_c$ is given by a Gaussian kernel. To account for potential manufacturing asymmetries, the kernel can be offset from the geometric center of the pixel by a value $\delta_{\text{offset}}$. 

$$
\text{IPSV}(x; x_c) = A \cdot \exp\left(-\frac{(x - (x_c + \delta_{\text{offset}}))^2}{2\sigma_{\text{ipsv}}^2}\right)
$$

where:
*   $\sigma_{\text{ipsv}}$ is the standard deviation of the intra-pixel sensitivity.
*   $\delta_{\text{offset}}$ is the centroid offset of the sensitivity peak relative to the geometric center.
*   $A$ is the normalization amplitude.

### Normalization Amplitude ($A$)

To ensure that the total integrated sensitivity of the pixel is conserved (i.e., the "volume" of the pixel's response remains 1), the Gaussian kernel is normalized strictly over the pixel footprint $[-\frac{\Delta x}{2}, \frac{\Delta x}{2}]$, rather than over the infinite real line $(-\infty, \infty)$.

The amplitude $A$ is derived such that:

$$
\int_{x_c - \frac{\Delta x}{2}}^{x_c + \frac{\Delta x}{2}} \text{IPSV}(x; x_c) \, dx = 1
$$

By shifting the coordinate system to the pixel center ($t = x - x_c$), this is equivalent to:

$$
\int_{-\frac{\Delta x}{2}}^{\frac{\Delta x}{2}} A \cdot \exp\left(-\frac{(t - \delta_{\text{offset}})^2}{2\sigma_{\text{ipsv}}^2}\right) dt = 1
$$

Using the error function ($\text{erf}$), the closed-form solution for the normalization constant $A$ is:

$$
A = \frac{1}{\sigma_{\text{ipsv}} \sqrt{\frac{\pi}{2}} \left[ \text{erf}\left( \frac{\frac{\Delta x}{2} - \delta_{\text{offset}}}{\sigma_{\text{ipsv}} \sqrt{2}} \right) + \text{erf}\left( \frac{\frac{\Delta x}{2} + \delta_{\text{offset}}}{\sigma_{\text{ipsv}} \sqrt{2}} \right) \right]}
$$

By integrating the product of the continuous PSF and this normalized IPSV Gaussian inside the pixel bounds, we simulate the realistic charge collection behavior of the detector pixel.
