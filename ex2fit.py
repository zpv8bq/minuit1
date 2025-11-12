import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from matplotlib.backends.backend_pdf import PdfPages

file = r.TFile.Open("experiments.root")
h1 = file.Get("hexp1")
h2 = file.Get("hexp2")


def hist_to_numpy(h):
    x = np.array([h.GetBinCenter(i) for i in range(1, h.GetNbinsX() + 1)])
    y = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])
    err = np.array([h.GetBinError(i) if h.GetBinError(i) > 0 else np.sqrt(max(y[i - 1], 1))
                    for i in range(1, h.GetNbinsX() + 1)])
    return x, y, err


x1, y1, err1 = hist_to_numpy(h1)
x2, y2, err2 = hist_to_numpy(h2)


def chi2(Nsig1, mean, sigma, slope1, offset1, Nsig2, slope2, offset2):
    G1 = Nsig1 * np.exp(-0.5 * ((x1 - mean) / sigma) ** 2)
    G2 = Nsig2 * np.exp(-0.5 * ((x2 - mean) / sigma) ** 2)

    B1 = slope1 * x1 + offset1
    B2 = slope2 * x2 + offset2

    model1 = G1 + B1
    model2 = G2 + B2

    chi2_1 = np.sum(((y1 - model1) / err1) ** 2)
    chi2_2 = np.sum(((y2 - model2) / err2) ** 2)

    return chi2_1 + chi2_2


m = Minuit(chi2, Nsig1=100, mean=50, sigma=5, slope1=-0.1, offset1=10, Nsig2=120, slope2=-0.05, offset2=8)
m.errordef = 1
m.migrad()
m.hesse()

mean = m.values["mean"]
sigma = m.values["sigma"]
emean = m.errors["mean"]
esigma = m.errors["sigma"]
chi2_val = m.fval
ndf = len(x1) + len(x2) - len(m.values)
pval = r.TMath.Prob(chi2_val, ndf)
chi2_r = chi2_val / ndf

print(f"Mean   = {mean:.3f} ± {emean:.3f}")
print(f"Sigma  = {sigma:.3f} ± {esigma:.3f}")
print(f"Red Chi2  = {chi2_r:.2f}")
print(f"p-value    = {pval:.3f}")

with PdfPages("ex2.pdf") as pdf:
    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.errorbar(x1, y1, yerr=err1, fmt='o', color='black', label='Data 1')
    fit1 = (m.values["Nsig1"] * np.exp(-0.5 * ((x1 - mean) / sigma) ** 2) +
            m.values["slope1"] * x1 + m.values["offset1"])
    plt.plot(x1, fit1, 'r-', lw=2, label='Fit')
    plt.xlabel("x")
    plt.ylabel("Counts")
    plt.title("Experiment 1")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.errorbar(x2, y2, yerr=err2, fmt='o', color='black', label='Data 2')
    fit2 = (m.values["Nsig2"] * np.exp(-0.5 * ((x2 - mean) / sigma) ** 2) +
            m.values["slope2"] * x2 + m.values["offset2"])
    plt.plot(x2, fit2, 'r-', lw=2, label='Fit')
    plt.xlabel("x")
    plt.ylabel("Counts")
    plt.title("Experiment 2")
    plt.legend()

    plt.tight_layout()
    pdf.savefig()
    plt.close()

    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')

    summary_text = f"""
    Fit Summary:
    Mean   = {mean:.3f} ± {emean:.3f}
    Sigma  = {sigma:.3f} ± {esigma:.3f}
    reduced Chi2    = {chi2_r:.2f}
    p-value    = {pval:.3f}


    The  p-value is  {pval:.3f} which 
    demonstrates a poor fit.
    The reduced chi2 is close to 2 which also demonstrates a
    poor fit for signal
    """

    ax.text(0.1, 0.9, summary_text, fontsize=12, va='top')
    pdf.savefig()
    plt.close()
