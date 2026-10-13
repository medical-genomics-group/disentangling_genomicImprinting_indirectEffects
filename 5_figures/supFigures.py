"""
Supplementary figures 
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
import csv
import pickle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from collections import Counter
from matplotlib.patches import Patch
import matplotlib.colors as mc
import colorsys

def set_style(textwidth_pts=418.0,ratio=0.8):
    pts_per_inch = 72.27
    
    plt.rcParams.update({
        # Sizing
        "figure.figsize":      [textwidth_pts / pts_per_inch, textwidth_pts / pts_per_inch * ratio],
        # Fonts
        "font.family":         "serif",      
        "font.size":           14,
        "axes.labelsize":      14,
        "axes.titlesize":      14,
        "xtick.labelsize":     12,
        "ytick.labelsize":     12,
        "legend.fontsize":     12,
        "figure.titlesize":    16,
        # PDF font embedding
        "pdf.fonttype":        42,
        "ps.fonttype":         42,
        # Lines
        "lines.linewidth":     1.0,
        "lines.markersize":    4,
        "axes.linewidth":      0.8,
        # Ticks
        "xtick.direction":     "in",   
        "ytick.direction":     "in",
        "xtick.major.size":    3.5,
        "ytick.major.size":    3.5,
        "xtick.major.width":   0.8,
        "ytick.major.width":   0.8,
        # Legend
        "legend.frameon":      False,  
        "legend.handlelength": 1.5,
    })

def subfig_bbox(subfig, fig):
    bbox     = subfig.get_window_extent()
    fig_bbox = fig.get_window_extent()
    x = bbox.x1 / fig_bbox.width
    y = bbox.y1 / fig_bbox.height
    return (x,y)

######
outdir="POE_mQTLs/outputs/figures/supplementary/"

width=8.3;length=11.7

colors=["#FA6C76", "#EC3F4C", "#CB1422","#9A0510", #reds
        "#6B89BC", "#4d6c9d", "#2e5180", "#073763", #blues
        "#8F54B3","#703096","#56187B","#3D085C", #purples
        "#FFF45E","#FFEE00","#D7C800","#736B00",##FFFE6F","#E1E03F", "#B8B71C", "#8A8905", #yellows
        "whitesmoke","gainsboro","silver","dimgray","black", #greys
        "#70CACA","#165454", #cyan
        "#07485B","#DBA507","#5E7C4D",
        "#3ABC0A","#01787E","#B2D40B"]  #jodie susie and lasso (maybe change)

size_subfig=16#20
size_labels=12
size_axes=9
size_title=14

#show palette
plt.figure(figsize=(10, 2))
for i, c in enumerate(colors):
    plt.bar(i, 1, color=c)

plt.xticks(range(len(colors)), colors, rotation=45, ha="right")
plt.yticks([])
plt.title("Color Palette Preview")
plt.tight_layout()
plt.show()


################
#functions
def plotViolin(data,title,y,labels,colors,save=False,h_thr=False,expected=False,
               axes=None,separators=None,group_labels=False,backColors=colors[16:21],extra_padding=0,ylim=None):
    if axes is None:
       fig, ax = plt.subplots(figsize=(8,6)) #unique plot
    else:
        ax=axes

    ax.margins(x=0)
    
    if separators: #vector of lines to separate different models
        separators.insert(0,0) #add zero at first
        separators.append(len(data)) #add how long is it
        for i in range(0,len(separators)-1):
            ax.axvspan(separators[i]+0.5,separators[i+1]+0.5,
                       color=backColors[i],alpha=0.6)
            ax.axvline(separators[i+1]+0.5,color="black",linestyle="-",linewidth=0.8,alpha=0.5) #line at the end
            #label per group
            if group_labels:
                xm=(separators[i]+1+separators[i+1])/2
                ax.text(x=xm,y=-0.18-extra_padding, #-0.18
                        s=group_labels[i],color=backColors[(len(backColors)-len(group_labels))+i],
                        fontsize=size_labels,ha="center") #,fontweight="bold"
                
    
    parts=ax.violinplot(data,showmedians=True)
    #color each violin
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.75) #change here for color bright
    # Set the color of the median lines
    parts['cmedians'].set_colors(colors) #####brighter?
    parts["cmedians"].set_alpha(1)
    parts["cmins"].set_colors(colors)
    parts["cmaxes"].set_colors(colors)
    parts["cbars"].set_colors(colors)
    
    handles_list=[];labels_list=[]
    if h_thr: #horizontal threshold
        h_thr_line=ax.axhline(h_thr,color="red",linestyle="--",alpha=0.8)
        handles_list.append(h_thr_line)
        labels_list.append("Threshold of variance")
    
    if expected: #add grey dist. next to color ones for expected/true values
        for ind in range(0,len(expected)):
            exp_mean=expected[ind][0]
            exp_std=expected[ind][1]
            if ind==len(expected)-1: 
                label="True values"
                labels_list.append(label)
                h_error=ax.errorbar(ind+1.35,exp_mean,yerr=exp_std,fmt="o",color="grey",capsize=4,label=label)
                handles_list.append(h_error)
            else:
                label=None

            ax.errorbar(ind+1.35,exp_mean,yerr=exp_std,fmt="o",color="grey",capsize=4,label=label)
            
    
    if axes is None: 
        legend = plt.legend(
            loc="upper left",
            bbox_to_anchor=(0.02, 0.98),
            #frameon=True,
            #fancybox=True,
            facecolor="white",   # background color of legend box
            edgecolor="black",   # border color
            fontsize=size_axes,
            labels=labels_list,
            handles=handles_list
        )
        legend.get_frame().set_alpha(0.9)
        legend.get_frame().set_linewidth(1.2)
    
    
    ax.set_xticks(range(1, len(labels)+1))
    ax.set_xticklabels(labels,fontsize=size_labels-2) #size of axes
    ax.set_title(title,y=0.8,fontsize=size_labels+2, pad=10) #size title
    if expected:
        ax.set_xlim(0,len(data)+1)

    if ylim:
        ax.set_ylim(0,ylim)
        ax.set_yticks(np.linspace(0,ylim,3))
    
    ax.set_ylabel(y, fontsize=size_labels+2, labelpad=1) #ylabel size
    
    if save:
        plt.savefig(save)
    return ax,handles_list,labels_list

def add_labels(x,y,axis):
    for i in range(len(x)):
        pos=x[i]
        axis.text(pos,y[i],y[i],ha='center',size=size_axes)  # Aligning text at center

def separateByPattern(list_runs,pattern_list,pattern): #subsamples list based on index of list of patterns by specific patter
    subsampled=[val for val, pat in zip(list_runs, pattern_list) if pat == pattern]    
    return subsampled

def plotMeanMethylation(y_dir,x_dir,x_ind,title="",save=False,axes=None): 
    y=zarr.load(y_dir)
    x=zarr.load(x_dir)

    x_snp=x[:,[x_ind*4,x_ind*4+3]]

    #get index of 0,01,10,2
    h0=[]
    h2=[]
    hm=[] #mother +1 01
    hf=[] #father -1 10

    for j in range(0,len(x_snp)):
        c=x_snp[j,0];i=x_snp[j,1]
        if c==0:
            h0.append(j)
        elif c==2:
            h2.append(j)
        else: #c=1
            if i==-1: #paternal
                hf.append(j)
            elif i==1:
                hm.append(j)
            
    y0=y[h0]
    y2=y[h2]
    ym=y[hm]
    yf=y[hf]

    #averages
    y0_mean=np.mean(y0);y2_mean=np.mean(y2);ym_mean=np.mean(ym);yf_mean=np.mean(yf)
    labels=["0","0|1","1|0","2"]
    
    #errors
    down=[np.std(x)/np.sqrt(np.size(x)) for x in [y0,yf,ym,y2]]
    up=[np.std(x)/np.sqrt(np.size(x)) for x in [y0,yf,ym,y2]]
    error=[down,up]
    #plot bar ###################### add error
    
    if axes is None:
       fig, ax = plt.subplots(figsize=(8,6)) #unique plot
       
    else:
        ax=axes
            
    ax.bar(labels,[y0_mean,yf_mean,ym_mean,y2_mean],color=["red","pink","purple","blue"],alpha=0.5)
    ax.errorbar(labels,[y0_mean,yf_mean,ym_mean,y2_mean],yerr=error,fmt='.',color='Black',alpha=0.5)
    
    if axes is None:
        ax.set_ylabel("Average methylation")
        ax.set_xlabel("Genotype")
    ax.set_title(title,fontsize=16)
    
    if save:
        plt.savefig(save)

def selectedSNPsIndex(beta_JODIE_path):
    beta=zarr.load(beta_JODIE_path)
    beta=np.reshape(beta,[int(len(beta)/4),4])
    selected_SNPs_ind=np.where(beta!=0)[0].tolist()
    return(selected_SNPs_ind)

#####Supplemaentary figures Simulations

##figs1:FDR_SNP 
resultsSim=pd.read_pickle(r"sim_final/out/resultsSim.pkl")

labels=["JODIE","SuSiE","LASSO"]

r2_pearson=[resultsSim[c].tolist() for c in ["r2_Pearson_JODIE","r2_Pearson_SUSIE","r2_Pearson_LASSO"]]
ROC=[resultsSim[c].tolist() for c in ["[TPR,FPR,FDR]SNP_JODIE","[TPR,FPR,FDR]SNP_SUSIE","[TPR,FPR,FDR]SNP_LASSO"]]
FDR=[[x[2] for x in sublist] for sublist in ROC]

set_style(textwidth_pts=418.0,ratio=0.85)
fig=plt.figure() 
plotViolin(FDR, "", "FDR", labels, colors=colors[23:26], axes=fig.subplots(1))

plt.show()
plt.savefig(outdir+"S1_FDR_SNP.pdf",format="pdf")

################################
####figs2:CI_patterns
#covariances CI different patterns
expect=[[0,0],[0,0],[0,0],[0,0],[0.2,0.126],[0.2,0.126],[-0.2,0.126],[-0.2,0.126]]##C+ C- CF CM F+ F- M+ M- []mean,std

cov_CI=[x[0,3] for x in resultsSim["var_JODIE"]]
cov_CI=[separateByPattern(cov_CI,resultsSim["new_notation"],pat) for pat in np.unique(resultsSim["new_notation"]).tolist()]

set_style(textwidth_pts=418.0,ratio=0.85)
plt.figure()
plotViolin(cov_CI,"","Covariance CI",
           np.unique(resultsSim["new_notation"]).tolist(),expected=expect,
           colors=["purple"]*len(np.unique(resultsSim["new_notation"]).tolist()))

plt.xlabel("Simulated signal",fontsize=size_labels)
plt.savefig(outdir+"S2_covCI.pdf")

#######################
#figs4:hist_mQTLs histogram pleiotrophy and polygenic for POE-mQTLs
mQTLs=pd.read_csv("POE_mQTLs/outputs/mQTLs_patternAnnotate.csv")

####
#fig5_dist to probes
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter

def fmt_mb(x, pos):
    return f'{x/1e6:.2g}'

def fmt_10mb(x, pos):
    return f'{x/1e7:.2g}'

set_style(textwidth_pts=418.0,ratio=0.5)

distance=np.abs(mQTLs["Position Probe"]-mQTLs["Position SNP"]).tolist()
cis=[i for i in distance if i<=1e6]
trans=[i for i in distance if i>1e6]
fig,axes=plt.subplots(1,2)

axes[0].hist(cis,50)
axes[0].set_ylabel("Frequency",fontsize=size_axes)
axes[0].set_xlim(0,1e6)
axes[0].tick_params(axis='both',labelsize=size_axes)

axes[1].hist(trans,50)
axes[1].set_ylabel("Frequency",fontsize=size_axes)
axes[1].yaxis.set_major_locator(MaxNLocator(integer=True))
axes[1].tick_params(axis='both',labelsize=size_axes)

plt.draw()

axes[0].xaxis.set_major_formatter(FuncFormatter(fmt_mb))
axes[0].set_xlabel("Distance SNP to probe <1Mb", fontsize=size_axes, labelpad=12)
# place the multiplier to the right of the x-axis
axes[0].annotate('×10⁶', xy=(0.99,0.125), xycoords='axes fraction',
                 xytext=(4,-20), textcoords='offset points',
                 fontsize=size_axes, ha='left', va='top')

axes[1].xaxis.set_major_formatter(FuncFormatter(fmt_10mb))
axes[1].set_xlabel("Distance SNP to probe >1Mb", fontsize=size_axes, labelpad=12)
axes[1].annotate('×10⁷', xy=(0.88,0.125), xycoords='axes fraction',
                 xytext=(4,-20), textcoords='offset points',
                 fontsize=size_axes, ha='left', va='top')

x_a,y_a=[x-0.01 for x in subfig_bbox(axes[0],fig)]
x_b,y_b=[x-0.01 for x in subfig_bbox(axes[1],fig)]

x_label = 0.04

fig.text(x_label, y_a+x_label, "a)", fontsize=size_subfig, ha="center", va="top")
fig.text(0.5+x_label, y_a+x_label, "b)", fontsize=size_subfig, ha="center", va="top")

fig.subplots_adjust(bottom=0.2,top=0.8,wspace=0.4,right=0.95,left=0.1)
fig.savefig(outdir+"s5_distProbe.pdf",
    dpi=300,            
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)


#######################################
### number of probes and snps per chr
probes_dir="screen_full/X_matrix/out_methylation/probes_chr_"
snps_dir="screen_full/X_matrix/out_x/snps_"
probes=[];snps=[]
for chro in range(1,24):
    probes.append(len(pd.read_csv(probes_dir+str(chro)+".csv")))
    snps.append(len(np.loadtxt(snps_dir+str(chro)+".txt",dtype=str)))
    
chrom=list(range(1,24))
chrom=[str(i)for i in chrom]
chrom[22]="X"
print(np.sum(probes))
print(np.sum(snps))

width_bar=0.4
set_style(textwidth_pts=418.0,ratio=0.85)
plt.figure()
plt.bar(range(1,24),probes,width_bar,label="probes",color="darkgray")
plt.bar([x+width_bar for x in range(1,24)],snps,width_bar,label="SNPs",color="black")
plt.xticks([x+width_bar/2 for x in range(1,24)],chrom,fontsize=size_axes)
plt.xlim([0.5,24])
plt.xlabel("chromosome",fontsize=size_title)
plt.legend(loc="best",fontsize=size_labels)
plt.xticks(fontsize=size_axes)

plt.savefig(outdir+"s_noProbesSnps.pdf",
    dpi=300,            
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)
##############
#Susie vs JODIE coefficients for mQTLs
### SUSIE
#just autosomes
snps_filt_autosomes=mQTLs[mQTLs["Chromosome"]!=23]
### compare with filtered ones snps
beta_S=[];pip_S=[];ind_S=[];snps_S=[];probe_S=[];
for index, row in snps_filt_autosomes.iterrows():
    chro=str(row["Chromosome"])
    probe=str(int(row["Index probe"]))
    snp=str(row["RSID"])
    
    path_beta="screen_full/"+chro+"/susie/out/beta_"+probe+".zarr/"
    path_pip="screen_full/"+chro+"/susie/out/pip_"+probe+".zarr/"
    path_cis="screen_full/"+chro+"/cis_probes/"+probe+"_probes_cis.csv"
    
    try:
        #load values from susie
        beta_sus=zarr.load(path_beta)
        pip_sus=zarr.load(path_pip)
    
        beta_sus=np.delete(beta_sus,0)
        beta_sus=np.reshape(beta_sus,(int(np.shape(beta_sus)[0]/4),4))
        pip_sus=np.reshape(pip_sus,(int(np.shape(pip_sus)[0]/4),4))
        cis_snps=np.genfromtxt(path_cis,delimiter="\n",dtype=str)
        
        #get values for specific snp
        ind_snp=np.where(cis_snps==snp)[0][0]
        beta_snp=beta_sus[ind_snp].tolist()
        pip_snp=pip_sus[ind_snp].tolist()
        
        #save for DF
        beta_S.append(beta_snp)
        pip_S.append(pip_snp)
        ind_S.append(ind_snp)
        snps_S.append(snp)
        probe_S.append(probe)
        
    except:
        snps_filt_autosomes=snps_filt_autosomes.drop(index=index) #eliminate row if it didnt run with susie

susie=pd.DataFrame({"RSID":snps_S,"Index SNP":ind_S,"Probe ID":probe_S,"betas":beta_S,"pip":pip_S})
import ast
def getValIndex(df,column,index):
    values=df[column].values.tolist()
    list_ind=[ast.literal_eval(x)[index] for x in values]
    return(list_ind)
                     
##Compare betas
beta_jodie=snps_filt_autosomes["beta I"].values.tolist()
beta_susie=[x[3] for x in susie["betas"].values.tolist()]

set_style(textwidth_pts=418.0,ratio=0.85)
plt.figure()
plt.scatter(beta_jodie,beta_susie,color="darkgreen")

plt.xlabel("JODIE coefficients",fontsize=size_labels)
plt.ylabel("SuSiE coefficients",fontsize=size_labels)
plt.xticks(fontsize=size_axes)
plt.yticks(fontsize=size_axes)

plt.savefig(outdir+"s_susieCoeff.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)
##### TPR and FDR threshold variance JODIE
jodie_thr=pd.read_csv(outdir+"jodie_varThr.csv")
variance_thresholds=jodie_thr["variance thr"].tolist()
TPR_thr=jodie_thr["TPR thr"].tolist()
FDR_thr=jodie_thr["FDR thr"].tolist()

#plot no std
set_style(textwidth_pts=418.0,ratio=0.8)
plt.figure()

plt.axvline(0.05,color="red",linestyle="--")
plt.scatter(variance_thresholds,TPR_thr,color="skyblue")
plt.scatter(variance_thresholds,FDR_thr,color="darkblue")
plt.legend(["Selected threshold","TPR","FDR"],fontsize=size_labels,
           loc="upper right",handletextpad=0.1)
plt.plot(variance_thresholds,TPR_thr,color="skyblue")
plt.plot(variance_thresholds,FDR_thr,color="darkblue")

plt.ylabel("Rate")
plt.xlabel("Threshold of variance")
plt.savefig(outdir+"s_varThr.pdf",format="pdf")


##############################
#### h2 simulations
list_h2=[0.1, 0.2, 0.3, 0.4]

move=0.01
colors=["#07485B","#DBA507","#5E7C4D"]

with open("h2/sim_h2.pkl","rb") as f:
    resultsSim=pickle.load(f)

##################
set_style(textwidth_pts=418.0,ratio=0.8)
fig=plt.figure()
(topfig,bottomfig)=fig.subfigures(2,1)
fig.subplots_adjust(hspace=0.15)

h2_FDR=[];h2_TPR=[];h2_FDR_std=[];h2_TPR_std=[]
for h2 in list_h2:
    filtered_h2=resultsSim[resultsSim["h2"]==h2]
    ROC=[filtered_h2[c].tolist()  for c in ["[TPR,FPR,FDR]parameter_JODIE","[TPR,FPR,FDR]parameter_SUSIE","[TPR,FPR,FDR]parameter_LASSO"]]
    TPR=[[x[0] for x in sublist] for sublist in ROC]
    h2_TPR.append([np.mean(x) for x in TPR]);h2_TPR_std.append([np.std(x) for x in TPR])
    FDR=[[x[2] for x in sublist] for sublist in ROC]
    h2_FDR.append([np.mean(x) for x in FDR]);h2_FDR_std.append([np.std(x) for x in FDR])
    

axs=topfig.subplots()

axs.errorbar([x-move for x in list_h2],[x[0] for x in h2_TPR],[x[0] for x in h2_TPR_std],color=colors[0],fmt="none",marker="o") #jodie
axs.errorbar(list_h2,[x[1] for x in h2_TPR],[x[1] for x in h2_TPR_std],color=colors[1],fmt="none",marker="^") #susie
axs.errorbar([x+move for x in list_h2],[x[2] for x in h2_TPR],[x[2] for x in h2_TPR_std],color=colors[2],fmt="none",marker="s") #lasso

axs.scatter([x-move for x in list_h2],[x[0] for x in h2_TPR],color=colors[0],marker="o",s=40)
axs.scatter(list_h2,[x[1] for x in h2_TPR],color=colors[1],marker="^",s=40)
axs.scatter([x+move for x in list_h2],[x[2] for x in h2_TPR],color=colors[2],marker="s",s=40)

axs.legend(["JODIE","SuSiE","LASSO"],ncol=3,loc=[0.3,1],labelspacing=0.05,
           columnspacing=0.05,handletextpad=0.05)
axs.set_ylabel("TPR",fontsize=size_title)
axs.set_xlabel("h\u00b2",fontsize=size_title)
axs.set_ylim(-0.05,1)

###
axs=bottomfig.subplots()
axs.errorbar([x-move for x in list_h2],[x[0] for x in h2_FDR],[x[0] for x in h2_FDR_std],color=colors[0],fmt="none",marker="o") #jodie
axs.errorbar(list_h2,[x[1] for x in h2_FDR],[x[1] for x in h2_FDR_std],color=colors[1],fmt="none",marker="^") #susie
axs.errorbar([x+move for x in list_h2],[x[2] for x in h2_FDR],[x[2] for x in h2_TPR_std],color=colors[2],fmt="none",marker="s") #lasso

axs.scatter([x-move for x in list_h2],[x[0] for x in h2_FDR],color=colors[0],marker="o",s=40)
axs.scatter(list_h2,[x[1] for x in h2_FDR],color=colors[1],marker="^",s=40)
axs.scatter([x+move for x in list_h2],[x[2] for x in h2_FDR],color=colors[2],marker="s",s=40)

axs.set_ylabel("FDR",fontsize=size_title)
axs.set_xlabel("h\u00b2",fontsize=size_title)
axs.set_ylim(-0.05,1)

plt.draw()
x_a,y_a=[x-0.01 for x in subfig_bbox(topfig,fig)]
x_b,y_b=[x-0.01 for x in subfig_bbox(bottomfig,fig)]

fig.text(x_label, y_a, "a)", fontsize=size_subfig, ha="center", va="top")
fig.text(x_label, y_b-0.05, "b)", fontsize=size_subfig, ha="center", va="top")

topfig.subplots_adjust(right=0.95,left=0.15,top=0.8,bottom=0.1)
bottomfig.subplots_adjust(right=0.95,left=0.15,top=0.8,bottom=0.2)

fig.savefig(outdir+"s_h2_performance.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)

########################
### histogram number of probes with var> thr for each effect
probes=pd.read_csv("POE_mQTLs/outputs/POE_probes.csv") #filtered probes screen 2 POE-probes
var_thr=0.05

effects=["C","M","F","I"]

freq=[]
for effect in effects:
    freq.append(len(probes[probes["var"+effect]>var_thr]))
    
set_style(textwidth_pts=418.0,ratio=0.8)
fig,axes=plt.subplots()

axes.bar(effects,freq)
add_labels(effects,freq,axes)
axes.set_ylabel("Probes with significant variance")
axes.set_xlabel("Effects")

fig.savefig(outdir+"s_histEffects.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)

mat=probes[probes["varM"]>var_thr]
pat=probes[probes["varF"]>var_thr]


### histogram cov CI first screeening
path="linearRegression_test/5run/screen_full"
passed1_dir="POE_mQTLs/outputs/probes_pased1.csv"

passed1=pd.read_csv(passed1_dir) #5434 

plt.figure()
plt.hist(passed1["covCI"].tolist(),bins=100)
plt.xlabel("covariance between C and I")
plt.ylabel("Number of probes")
plt.savefig(outdir+"s_hist_covCI.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)

### violin plot how big are the effects for the dependent probes? per effect
#load data
path="POE_mQTLs/outputs/"
## probes
m_probes=pd.read_csv(path+"maternal_probes.csv") #191
f_probes=pd.read_csv(path+"paternal_probes.csv") #426
poe_probes=pd.read_csv(path+"POE_probes.csv") #281

colors=["#532B6E","#8E0B3B","#004799","#1B7E24",
        "#07485B","#DBA507","#5E7C4D"]

    
plt.figure()
fig,axes=plt.subplots()
plotViolin([m_probes["varM"].tolist(),f_probes["varF"].tolist(),poe_probes["varI"].tolist()],
           "","variances of significant effect",["maternal","paternal","imprinting"],colors[1:4],
           axes=axes)

plt.savefig(outdir+"s_varSignificant.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)
