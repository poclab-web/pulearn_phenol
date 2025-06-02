import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path

def Parameter_tuning_vis(p_result, n_result, mcc, n_estimators, descriptor):
    """
    The function to output a line graph of the results of parameter tuning.
    """

    num_u = p_result.index
    
    fig = plt.figure()
    ax = fig.add_subplot()

    p1 = ax.plot(num_u, p_result, label="$\it{TPR}$",
                 linewidth=1, linestyle="solid", markersize=2, marker="o", c="#ff4b00")
    
    p2 = ax.plot(num_u, n_result, label="$\it{TNR}$",
                 linewidth=1, linestyle="dotted", markersize=2, marker="o", c="#005aff")
    
    p3 = ax.plot(num_u, mcc, label="$\it{MCC}$",
                 linewidth=1, linestyle="dashed", markersize=2, marker="o", c="#03af7a")

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([i for i in range(0, 21)])
    ax.set_xlabel("$\it{a}$ ($\it{K}$ = $\it{a}$ * $\it{NP}$)", fontsize = 16)
    ax.set_ylabel("Score", fontsize = 16)
    ax.set_title(f"{descriptor}" + " ($\it{T}$ = " + f"{n_estimators})", fontsize = 18)
    fig.tight_layout() 
    
    plt.savefig(f"Parameter_tuning/Parameter_tuning_{descriptor}_{n_estimators}.png", dpi=300)
    plt.clf()

def dps_distribution(dps, y0_d, y0_t, y1_t, title):
    """
    The function to output a histogram representing the DPS distribution of Unlabeled.
    """
    
    fig, ax = plt.subplots(nrows=2, figsize=(3,4), dpi=160, sharex='col',
                           gridspec_kw={'height_ratios': (1,1)} )
    ax[0].hist(dps, bins=np.arange(0, 1.0+0.1, 0.1), color="#03af7a", edgecolor="k", linewidth=0.2)
    ax[1].hist(dps, bins=np.arange(0, 1.0+0.1, 0.1), color="#03af7a", edgecolor="k", linewidth=0.2)
    fig.subplots_adjust(hspace=0.0)
    ax[1].set_ylim(0, y1_t)
    ax[1].set_yticks([0, 500, 1000])
    ax[0].set_ylim(y0_d, y0_t)
    ax[0].set_yticks([2500, 3000, 3500])
    ax[1].spines['top'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    
    d1 = 0.02 # Amount of X-axis overhang
    d2 = 0.03 # Nyoro Wave Height
    wn = 21   # Number of Nyoro waves (specify odd values)
    pp = (0,d2,0,-d2)
    px = np.linspace(-d1,1+d1,wn)
    py = np.array([1+pp[i%4] for i in range(0,wn)])
    p = Path(list(zip(px,py)), [Path.MOVETO]+[Path.CURVE3]*(wn-1))
    
    line1 = mpatches.PathPatch(p, lw=4, edgecolor='black',
                               facecolor='None', clip_on=False,
                               transform=ax[1].transAxes, zorder=10)
    line2 = mpatches.PathPatch(p,lw=3, edgecolor='white',
                               facecolor='None', clip_on=False,
                               transform=ax[1].transAxes, zorder=10,
                               capstyle='round')
    a = ax[1].add_patch(line1)
    a = ax[1].add_patch(line2)
    ax[0].set_title(f'{title}', fontsize = 18)
    ax[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax[1].set_xlabel('$\it{DPS}$', fontsize = 14)
    ax[0].vlines(0.5, y0_d, y0_t, colors='k', linestyle='dashed', linewidth=0.5)
    ax[1].vlines(0.5, 0, y1_t, colors='k', linestyle='dashed', linewidth=0.5)
    plt.savefig(f"Predict_unlabeled/DPS_distribution_{title}.png", dpi=300, bbox_inches='tight')
    plt.clf()