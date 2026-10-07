#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Structured and clean code to get each figure in the article

Figures

"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from collections import Counter
from matplotlib.patches import Patch
import matplotlib.colors as mc
import colorsys
import zarr
import zipfile
import matplotlib.ticker as ticker

def set_style(textwidth_pts=418.0,ratio=0.8):
    pts_per_inch = 72.27
    
    plt.rcParams.update({
        # Sizing
        "figure.figsize":      [textwidth_pts / pts_per_inch, textwidth_pts / pts_per_inch * ratio],
        # Fonts
        "font.family":         "serif",      # matches LaTeX body font
        "font.size":           14,
        "axes.labelsize":      14,
        "axes.titlesize":      14,
        "xtick.labelsize":     9,
        "ytick.labelsize":     9,
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
        "legend.frameon":      False,  #no box around legend
        "legend.handlelength": 1.5,
    })

######

outdir="POE_mQTLs/outputs/figures/"

colors=["#532B6E","#8E0B3B","#004799","#1B7E24",
        "#07485B","#DBA507","#5E7C4D"]

size_subfig=16
size_labels=12
size_axes=9
size_title=14

width=8.3;length=11.7
left_s=0.15  
right_s=0.97
bottom_s=0.1 
top_s=0.8   
x_label = 0.04  
 
#### function violin plots
def plotViolin(data,title,y,labels,colors,save=False,h_thr=False,expected=False,
               axes=None,separators=None,group_labels=False,backColors=colors[16:21],label_loc=0,ylim=None,
               axes_color=False):
    if axes is None:
       fig, ax = plt.subplots(figsize=(8,6)) #unique plot
    else:
        ax=axes

    ax.margins(x=0)
    if axes_color==False:
        axes_color=backColors
    if separators: #vector of lines to separate different models
        separators.insert(0,0) #add zero at first
        separators.append(len(data)) #add how long is it
        for i in range(0,len(separators)-1):
            ax.axvline(separators[i+1]+0.5,color="black",linestyle="-",linewidth=0.8,alpha=0.5) #line at the end
            #label per group
            if group_labels:
                xm=(separators[i]+1+separators[i+1])/2
                ax.text(x=xm,y=label_loc, 
                        s=group_labels[i],
                        ha="center",fontsize=size_labels) 
    
    parts=ax.violinplot(data,showmedians=True)
    #color each violin
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.75) #change here for color bright
    # Set the color of the median lines
    parts['cmedians'].set_colors(colors) 
    parts["cmedians"].set_alpha(1)
    parts["cmins"].set_colors(colors)
    parts["cmaxes"].set_colors(colors)
    parts["cbars"].set_colors(colors)
    
    handles_list=[];labels_list=[]
    if h_thr: #horizontal threshold
        h_thr_line=ax.axhline(h_thr,color="red",linestyle="--",alpha=0.8)
        handles_list.append(h_thr_line)
        labels_list.append("5% threshold")
    
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
            
    
    if axes is None: #only one plot add title and labels x axis
        legend = plt.legend(
            loc="upper left",
            frameon=True,
            fancybox=True,
            facecolor="white",   # background color of legend box
            edgecolor="black",   # border color
            fontsize=size_axes,
            labels=labels_list,
            handles=handles_list
        )
        legend.get_frame().set_alpha(0.9)
        legend.get_frame().set_linewidth(1.2)
        
    ax.set_xticks(range(1, len(labels)+1))
    ax.set_xticklabels(labels) #size of axes #fontsize=size_labels-2
    ax.set_title(title,y=0.65) #size title #fontsize=size_labels+2, 
    if ylim:
        ax.set_ylim(ylim[0],ylim[1])
        ax.set_yticks(np.linspace(0,ylim[1],3))
    
    ax.set_ylabel(y,labelpad=1,fontsize=size_labels) #ylabel size #fontsize=size_labels+2,
   
    if save:
        plt.savefig(save)
    return ax,handles_list,labels_list

def plotErrorBars(data,label,colors,axes=None,ylims=None,markers=["o","^","s"]):
    ##plot mean + error bars
    if axes is None:
       fig,ax = plt.subplots(figsize=(8,6)) #unique plot
    else:
        ax=axes
    
    for i in range(0,len(data)):
        ax.errorbar(i,np.mean(data[i]),np.std(data[i]),fmt="none",ecolor=colors[i])
        ax.scatter(i,np.mean(data[i]),color=colors[i],marker=markers[i],s=40)
        
    ax.set_xlabel(label,labelpad=-0.05) #,fontsize=size_title
    ax.set_xlim(-0.5,len(data)-0.5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    
    if ylims:
        ax.set_ylim(ylims[0],ylims[1])
        ax.set_yticks(np.linspace(ylims[0],ylims[1],3))

def subfig_bbox(subfig, fig):
    bbox     = subfig.get_window_extent()
    fig_bbox = fig.get_window_extent()
    x = bbox.x1 / fig_bbox.width
    y = bbox.y1 / fig_bbox.height
    return (x,y)

#show palette
plt.figure(figsize=(10, 2))
for i, c in enumerate(colors):
    plt.bar(i, 1, color=c)

plt.xticks(range(len(colors)), colors, rotation=45, ha="right")
plt.yticks([])
plt.title("Color Palette Preview")
plt.tight_layout()
plt.show()


######################################################
############## FIG poe probes
#####################################################
#function get chromosome positions
def getChroPos(df,chr_col,pos_col,pad=20000000):
    chr_pos=pd.DataFrame({chr_col:np.unique(df[chr_col]).tolist(),
                            "start":(df.groupby(chr_col)[pos_col].min()).tolist(),
                            "end":(df.groupby(chr_col)[pos_col].max()+pad).tolist()})
    return(chr_pos)

#function "manhattan"
def plotAllManhattan(df,y_col,y_title,y_lim=False,chrom_pos=None,thr=0.05,chr_col="Chromosome",pos_col="Position",title=False,filtered=None,knowns=None,position=True,save=False,figure=False):
    #get position of the chromosomes
    if chrom_pos is None:
        chrom_pos=getChroPos(df,chr_col,pos_col)
    #add row in 
    chro_add=[];pos_add=[]
    for index,row in chrom_pos.iterrows():
        if row["start"] not in df[pos_col]:
            chro_add.append(row[chr_col])
            pos_add.append(row["start"])
        if row["end"] not in df[pos_col]:
            chro_add.append(row[chr_col])
            pos_add.append(row["end"])
            
    add=pd.DataFrame({chr_col:chro_add,pos_col:pos_add})
    
    df=pd.concat([df,add],ignore_index=True)
   
    #####
    #position
    df=df.sort_values([chr_col,pos_col])
    df.reset_index(inplace=True, drop=True); 
    df['i']=df.index
    
    chrom_offsets = (
        df.groupby(chr_col)[pos_col]
        .max() 
        .cumsum()
        .shift(fill_value=0)
    )

    #offsets to each row
    df['cumulative_pos'] = df[chr_col].map(chrom_offsets)+df[pos_col]
    if position: #plot based on cumulative position
        plot_based="cumulative_pos"
    else: #plot x based on order
        plot_based="i"

    #######   Plot
    #colors
    colors=["green"]*23
    greys=["grey"]*23
    #legend
    legend_elements=[]
    
    ### create or take figure/axes
    if figure==False:
        fig,plot=plt.subplots()
    else:
        plot=figure.subplots();fig=figure

    
    if filtered: #change alpha depending on column given by filtered (boolean)
        #filtered ones (screen 1 passed not screen 2 with MF due to varI<0.05 or not confident value)
        sns.scatterplot(
            data=df[df[filtered]==True], x=plot_based, y=y_col,
            hue=chr_col,palette=greys[:-1],
            s=10,edgecolor="none",
            alpha=0.6,
            legend=False,
            ax=plot
        )  
        legend_elements.append(Line2D([0],[0],marker='o',label='filtered probes',
                          markerfacecolor="grey",alpha=0.6,markersize=8,linestyle="None",markeredgecolor="none"))
    
        #POE probes (passed both analysis)
        filt=df[df[filtered]==False]
        sns.scatterplot(
            data=filt, x=plot_based, y=y_col,marker="D",
            hue=chr_col, palette=[colors[i-1] for i in np.unique(filt[chr_col].tolist())],
            s=12,edgecolor="none",
            legend=False,
            ax=plot
            )
        legend_elements.append(Line2D([0],[0], marker='D', label='POE dependent probes',
                          markerfacecolor=colors[0],markersize=8,linestyle="None",markeredgecolor="none"))
        
    else:
        sns.scatterplot(
            data=df, x=plot_based, y=y_col,
            hue=chr_col, palette=[colors[i-1] for i in np.unique(df[chr_col].tolist())],
            s=8,edgecolor="none",
            legend=False,
            ax=plot
        )  

    plot.grid(False,axis="x")

    #add shading of known imprinted regions
    if knowns is not None:
        area=4500000#make known imprinted genes areas more apperant 
        colors_known=["dodgerblue","purple","blue","green","red"] #for different methods/datasets
        for ind_k,known in enumerate(knowns): #error here when is not list of one takes columns
            known['start_offset']=known['Chromosome'].map(chrom_offsets)+known['start']-area 
            known['end_offset']=known['Chromosome'].map(chrom_offsets)+known['end']+area
          
            for index,row in known.iterrows():
                plt.axvspan(int(row["start_offset"]),int(row["end_offset"]), facecolor=colors_known[ind_k],ymin=0,alpha=0.5,zorder=0)
          
            legend_elements.append(Patch(facecolor=colors_known[ind_k], alpha=0.5, label="Rosenski et al. known iDMR")) 
          
    # Calculate chromosome midpoints
    tags = (df.groupby(chr_col)['cumulative_pos'].max() - 
         df.groupby(chr_col)['cumulative_pos'].min()) / 2 + \
        df.groupby(chr_col)['cumulative_pos'].min()

    # Set axis labels and ticks
    plot.set_xlabel(chr_col.lower(),labelpad=0,fontsize=size_labels) #fontweight="bold" ,fontsize=size_title,
    # x labels
    plot.set_xticks(tags)# position of label
    chro=[str(x) for x in tags.index]
    chro[chro.index("23")]="X"
    plot.set_xticklabels(chro,fontsize=8,rotation=75) #labels #,fontsize=size_axes
    
    plot.set_ylabel(y_title,labelpad=1.2,fontsize=size_labels) #,fontweight="bold" 
    
    
    # Set title
    if title:
        fig.suptitle(title) #,y=0.97,fontsize=size_title

    # x axis limits
    xmin = df['cumulative_pos'].min()
    xmax = df['cumulative_pos'].max()
    padding = (xmax - xmin) * 0.01  # 1% padding
    plot.set_xlim(xmin - padding, xmax + padding)
    
    if thr:
        plot.axhline(thr, ls="--", color="red", label="5% threshold") #draw threshold
        legend_elements.append(Line2D([0], [0], color='red',linestyle='--',label='5% threshold'))
    if y_lim:
        plot.set_ylim([y_lim[0],y_lim[1]])
    if save:
        plt.savefig(save)
        
    return(plot,legend_elements)
        
#load data
probes=pd.read_csv("POE_mQTLs/outputs/poe_probes.csv") #filtered probes screen 2 POE-probes
all_probes=pd.read_csv("POE_mQTLs/outputs/screen2CMFI.csv") #all probes of second screening (filters were not defined yet that's why there is more than screen1 passed)
screen1_passed=pd.read_csv("POE_mQTLs/outputs/probes_pased1.csv") #probes that passed first screening filters

all_probes["Passed1"]=[x in screen1_passed["Probe ID"].tolist() for x in all_probes["Probe ID"].tolist()]
all_probes=all_probes[all_probes["Passed1"]==True]

all_probes["Filtered"]=[x not in probes["Probe ID"].tolist() for x in all_probes["Probe ID"].tolist()]
screens=pd.merge(all_probes,screen1_passed[["Probe ID","varI"]],on="Probe ID")
screens["Diff I"]=screens["varI_x"]-screens["varI_y"]
screens["Diff I abs"]=np.abs(screens["Diff I"])


known_impt_dir="data/knownImp_coordinates.xlsx"
known_imp=pd.read_excel(known_impt_dir,skiprows=[0,1],sheet_name=0)
known_imp["Chromosome"]=[int(x.split("chr")[1]) for x in known_imp["chrom"].tolist()]
known_imp=known_imp[["Chromosome","start","end","ICR name"]]


chrom_pos=getChroPos(all_probes,"Chromosome","Position")

################
### Create plot
set_style(textwidth_pts=418.0,ratio=0.8)
fig=plt.figure() 
(topfig,bottomfig)=fig.subfigures(2,1,height_ratios=[1,1.75])

####################################
#### bottom fig A) and B) ####
top_axs=topfig.subplots(1,2)

#######
# B) violin dist. variances each component POE probes
#######
var_comp=[probes["varC"].tolist(),probes["varM"].tolist(),probes["varF"].tolist(),probes["varI"].tolist()]
plotViolin(var_comp,"","variance",["C","M","F","I"],colors=[colors[i] for i in [0,1,2,3]], 
           axes=top_axs[1],h_thr=0.05,ylim=[0,0.8],label_loc=0.8)

######
# A) violin comparison CI variances on CI and CMFI models 
######
normVal=3.474051424481396
filt=all_probes[(all_probes[["var varC","var varM","var varF","var varI"]]<normVal).all(axis=1)]

both=pd.merge(screen1_passed[["Chromosome","Probe ID","varC","varI"]],
                      filt[["Chromosome","Probe ID","varC","varI"]],
                      on=["Chromosome","Probe ID"])

var_dif=[both["varC_x"].tolist(),both["varI_x"].tolist(),both["varC_y"].tolist(),both["varI_y"].tolist()]
plotViolin(var_dif,"","variance",["C","I","C","I"],colors=[colors[i] for i in [0,3,0,3]],
           axes=top_axs[0],h_thr=0.05,separators=[2],group_labels=["CI model","CMFI model"],
           label_loc=0.85,ylim=[0,0.8])

             
topfig.subplots_adjust(
    bottom=bottom_s,
    top=top_s,
    right=right_s+0.02,
    left=left_s,
    wspace=0.4)

###########################
######## fig C) ########
plot,legend_elements=plotAllManhattan(all_probes,"varI","variance I",filtered="Filtered",
                 knowns=[known_imp],y_lim=[0,0.5],figure=bottomfig)

bottomfig.legend(
    handles=legend_elements,
    loc="upper center",
    bbox_to_anchor=(0.63,0.99),  
    ncol=2,    
    fontsize=9,
    frameon=False
)


bottomfig.subplots_adjust(top=top_s+0.025,bottom=bottom_s+0.03,left=left_s,right=right_s+0.02)

plt.draw()

x_a,y_a=[x-0.01 for x in subfig_bbox(topfig,fig)]
x_b,y_b=[x-0.01 for x in subfig_bbox(bottomfig,fig)]

fig.text(x_label,y_a, "a)", fontsize=size_subfig, ha="center", va="top")
fig.text(0.5+x_label,y_a, "b)", fontsize=size_subfig, ha="center", va="top")
fig.text(x_label, y_b, "c)", fontsize=size_subfig, ha="center", va="top")

plt.show()
fig.savefig(outdir+"figure2.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,
    backend="pdf",
)

######################################################
############## FIG 3 SNPs
#####################################################

################
### Create plot
set_style(textwidth_pts=418.0,ratio=0.8)
fig=plt.figure() 
topfig,bottomfig=fig.subfigures(2,1,height_ratios=[1.5,1])

"""
 3. POE-mQTLs
 a) Scatter plot beta C vs I colors based on classification
 b) Histogram classification of mQTLs
 c) Histogram one POE-probe associated mQTLs 
 d) Histogram one mQTLs associated POE-probes
"""

mQTLs=pd.read_csv("POE_mQTLs/outputs/cluster_df.csv")

##############################
#### fig a)#### Imprinting pattern
unique_patt_noOrd=mQTLs["Selected pattern"].unique().tolist()
unique_patt=[unique_patt_noOrd[i] for i in [1,2,3,0]] #selected order

colors_patt=[colors[i] for i in [1,2,3]]
colors_patt.append("grey")
colors_map=dict(zip(unique_patt,colors_patt))
markers=dict(zip(unique_patt,["^","s","o","x"])) #mother,father,complex,undefined 

top_axs=topfig.subplots(1,2)

for patt in unique_patt:
    if patt=="Undefined":
        subset=mQTLs[mQTLs["Selected pattern"]==patt]
        alpha=0.2
    else:
        subset=mQTLs[mQTLs["Selected pattern"]==patt]
        alpha=0.65
        
    top_axs[0].scatter(subset["beta C"],subset["beta I"],color=colors_map[patt],
                       label=patt,alpha=alpha,marker=markers[patt],s=11)
    top_axs[0].tick_params(axis='both', which='minor', labelsize=size_axes)
    
    
legend_scatter=[plt.scatter([0],[0],color=colors_patt[i],
                            s=40,marker=["^","s","o","x"][i],alpha=0.75) 
                for i in range(0,len(unique_patt))]

top_axs[0].legend(legend_scatter,
           [txt.lower() for txt in unique_patt],
           loc=[0,1],
           ncol=4,
           columnspacing=0.05,
           handletextpad=0.01)

topfig.text(0.08,0.92,"imprinting classification",fontsize=size_labels)
top_axs[0].set_xlabel("direct genetic coefficient",fontsize=size_labels,labelpad=0.05)
top_axs[0].set_ylabel("imprinting coefficient",fontsize=size_labels,labelpad=0.05)
top_axs[0].set_xlim([-1,1])
top_axs[0].set_ylim([-1,1])

topfig.subplots_adjust(top=top_s,bottom=bottom_s,right=right_s,left=left_s,wspace=0.3)

##############################
#### fig b)
#histogram plot different patterns detected
identified=mQTLs.dropna(subset=["Selected pattern","Non imprinting effects"],how="all") 
identified=identified.replace(np.nan,"Direct")

list_patterns=identified["Selected pattern"].tolist() 
list_patterns=[x for x in list_patterns if str(x) != 'nan']
print(np.unique(list_patterns))

order = ['Maternal','Paternal','Complex',"Undefined"]

counts = Counter(list_patterns)
values = [counts[o] for o in order]

bars=top_axs[1].bar(order,values,color=colors_patt,alpha=0.75)

top_axs[1].set_xticks([0,1,2,3],labels=["Mat.","Pat.","Comp.","Und."])
top_axs[1].bar_label(bars,fontsize=9)
top_axs[1].set_ylim(0,1800)
top_axs[1].set_xlabel("imprinting pattern",fontsize=size_labels,labelpad=-0.08)

###########################
#### fig c and d
def add_labels(x,y,axis):
    for i in range(len(x)):
        pos=x[i]
        axis.text(pos,y[i]+0.01,y[i],ha='center',size=9) 

bottom_axs=bottomfig.subplots(1,2)

mQTLs["cluster id"]=[row["Original mQTL"]+"_"+str(row["Cluster"]) for i,row in mQTLs.iterrows()]
clusters=mQTLs[["cluster id","Probe ID"]]
clusters=clusters.drop_duplicates()
same_probe=clusters["Probe ID"].value_counts()

values,bins,bars=bottom_axs[0].hist(same_probe,color="grey")
bottom_axs[0].bar_label(bars,fontsize=9)
bottom_axs[0].set_xlabel("clusters per probe",fontsize=size_labels,labelpad=-0.5) #"mQTLs linked to one probe"
bottom_axs[0].set_ylabel("probes",fontsize=size_labels,labelpad=-1)
bottom_axs[0].set_yticks(np.linspace(0,100,3))
bottom_axs[0].set_ylim(0,100)

###fig d 
same_SNP=clusters["cluster id"].value_counts()

values,bins,bars=bottom_axs[1].hist(same_SNP,color="grey")

bottom_axs[1].bar_label(bars,fontsize=9)
bottom_axs[1].set_xlabel("probes per clusters",fontsize=size_labels,labelpad=-0.5) #Probes linked to one mQTL
bottom_axs[1].set_ylabel("clusters",fontsize=size_labels,labelpad=-1)
bottom_axs[1].set_yticks(np.linspace(0,900,3)) 
bottom_axs[1].set_ylim(0,900)

bottomfig.subplots_adjust(bottom=bottom_s+0.12,top=top_s,wspace=0.3,
                          left=left_s,right=right_s)

plt.draw()

x_a,y_a=[x-0.01 for x in subfig_bbox(topfig,fig)]
x_b,y_b=[x-0.01 for x in subfig_bbox(bottomfig,fig)]

x_label = 0.04  # one global value for all left-side labels

fig.text(x_label,y_a,"a)", fontsize=size_subfig, ha="center", va="top")
fig.text(0.5+x_label,y_a,"b)", fontsize=size_subfig, ha="center", va="top")
fig.text(x_label,y_b,"c)", fontsize=size_subfig, ha="center", va="top")
fig.text(0.5+x_label,y_b,"d)", fontsize=size_subfig, ha="center", va="top")

plt.show()

fig.savefig(outdir+"figure3_cluster.pdf",
    dpi=300,            
    pad_inches=0.02,     
    backend="pdf",
)

######################################################
#####################    FIG 4
####################################################
### Create plot
set_style(textwidth_pts=418.0,ratio=0.75)
fig=plt.figure() 
topfig,bottomfig=fig.subfigures(2,1,height_ratios=[1.2,1]) 


"""
 4. Replication
 a) Overlaps other studies
 b) Violin variances Novel, replicated, exclusive other studies
 c) hist replicated and not with Zeng's methods
"""

def adjustColor(color,amount): 
    try:
        c = mc.cnames[color]
    except:
        c = color
    r, g, b = mc.to_rgb(c)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = max(0, min(1, l * amount)) #1=original 1.5=light 0.5 dark
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return (r, g, b)

def plotImprintingRegions(df,chrom_pos,pos_col,chr_col,source_list,label_list,pad,colors,save=False,figure=False):
    ###chromosome position
    #add row in  0 to fit start and end of chromosome
    chro_add=[];pos_add=[]
    for index,row in chrom_pos.iterrows():
        if row["start"] not in df[pos_col]:
            chro_add.append(row[chr_col])
            pos_add.append(row["start"])
        if row["end"] not in df[pos_col]:
            chro_add.append(row[chr_col])
            pos_add.append(row["end"])
            
    add=pd.DataFrame({chr_col:chro_add,pos_col:pos_add})
    
    df=pd.concat([df,add],ignore_index=True)

    #####
    #position
    df=df.sort_values([chr_col,pos_col])
    df.reset_index(inplace=True, drop=True); 

    chrom_offsets = (
        df.groupby(chr_col)[pos_col]
        .max() 
        .cumsum()
        .shift(fill_value=0)
    )

    #offsets to each chromosome
    df['cumulative_pos'] = df[chr_col].map(chrom_offsets)+df[pos_col]
    #####
    # Calculate chromosome midpoints
    tags = (df.groupby(chr_col)['cumulative_pos'].max() - 
         df.groupby(chr_col)['cumulative_pos'].min()) / 2 + \
        df.groupby(chr_col)['cumulative_pos'].min()

    # x axis limits
    xmin = df['cumulative_pos'].min()
    xmax = df['cumulative_pos'].max()
    padding = (xmax - xmin) * 0.01  # 1% padding

    num=len(source_list)
    if figure==False:
        fig,axs=plt.subplots(num,sharex=True) 

    else:
        axs=figure.subplots(num,sharex=True);fig=figure

    for i,source in zip(range(0,num),source_list): #for each source in df (study_)
        df_s=df.loc[df["Source"]==source].copy()
        df_s['start_offset']=df_s["cumulative_pos"]-pad
        df_s['end_offset']=df_s["cumulative_pos"]+pad
        
        ###colors are different depending on overlap
        if source == "Ours":
            no_overlap_color=colors[0] #red ours and no overlap
        else:
            no_overlap_color=colors[2] #grey not ours no overlap
            
        overlap_color=colors[1] #blue overlap
        
        df_no_overlap=df_s[df_s["Overlap"].isna()]
        df_overlap=df_s[(df_s["Overlap"]=="Position")|(df_s["Overlap"]=="ID")]

       
        for df_so,color_point in zip([df_no_overlap,df_overlap],
                                     [no_overlap_color,overlap_color]):  
            for index,row in df_so.iterrows():
                axs[i].axvspan(
                    int(row["start_offset"]),
                    int(row["end_offset"]),
                    facecolor=color_point,
                    zorder=0,alpha=0.8)
               
        axs[i].set_yticks([])
        axs[i].set_ylabel(label_list[i],fontsize=size_axes,rotation=15,
                          labelpad=1.4,va="center",ha="right")
        
    # Set axis labels and ticks
    axs[num-1].set_xlabel("chromosome",fontsize=size_labels,labelpad=0.5)#,fontweight="bold"
    # x labels
    axs[num-1].set_xticks(tags)# position of label
    chro=[str(int(x)) for x in tags.index]
    chro[chro.index("23")]="X"
    axs[num-1].set_xticklabels(chro,fontsize=9,rotation=75) #labels
    #colors

    plt.setp(axs,xlim=(xmin-padding,xmax+padding))
    fig.subplots_adjust(wspace=0, hspace=0)

    if save:
        fig.savefig(save)

####load data
path_otherStudies="POE_mQTLs/outputs/"
novel=pd.read_csv(path_otherStudies+"novel_studiesComparison.csv")
overlap=pd.read_csv(path_otherStudies+"overlap_studiesComparison.csv")
otherStudies=pd.read_csv(path_otherStudies+"others_studiesComparison.csv")

probes_merged_path="POE_mQTLs/outputs/comparisonStudies.csv"
probes_merged=pd.read_csv(probes_merged_path)

chrom_pos=getChroPos(probes_merged,"Chr","Position")

##############
#ours is a single poiint, give others average distance of same as ours
##############################
### a) overlap other studies
colors_overlap=["#800000","#0066ff","#bfbfbf"]

plotImprintingRegions(probes_merged,chrom_pos,"Position","Chr",
                       source_list=["Ours","Rosenski iDMR","Rosenski p.ASM",
                                    "Cuellar POE-CpG","Cuellar Imprint"],
                       label_list=["Ours","R. et al.\niDMR","R. et al. \np.ASM",
                        "C.P. et al.\nImprint","C.P. et al.\nPOE-CpG"],
                       pad=5e6,colors=colors_overlap,figure=topfig) #5e6

labels=["novel","overlap","other studies"]


legends = [Patch(facecolor=colors_overlap[i], label=labels[i])
           for i in range(len(labels))]

topfig.legend(
    handles=legends,
    loc="upper center",
    bbox_to_anchor=(0.63,0.99),   
    ncol=3,     # all in one row
    fontsize=10,
    frameon=False
)

topfig.subplots_adjust(top=top_s+0.05,bottom=bottom_s+0.05,
                       right=0.98,left=left_s) 

####################################
### b) violin variances: novel, replicated, other studies
bottom_axs=bottomfig.subplots(1,2)

var_replication=[]
for df in [novel,overlap,otherStudies]:
    for var in ["varC","varM","varF","varI"]:
        var_replication.append(df[var].tolist())

colors_overlap[2]="dimgray"
plotViolin(var_replication,"","variance",["C","M","F","I"]*3,colors=[colors[i] for i in [0,1,2,3]]*3,
          h_thr=0.05,separators=[4,8],group_labels=["Novel","Overlap","Others"],
          axes=bottom_axs[0],label_loc=0.52,ylim=[0,0.5],axes_color=colors_overlap) 

###################################
### c)replication Zeng
zeng_rep=pd.read_csv("h2_siblings/2zeng/3analysis/zeng_rep.csv",
                     na_values="",keep_default_na=False)
models=["smp","sibs","sm","sp","null"]
hist_model=[len(zeng_rep[zeng_rep["model"]==x]) for x in models]
labels=["smp","s","sm","sp","null"]
bottom_axs[1].bar(labels,hist_model,
        color=["grey","grey","grey","grey","black"])
bottom_axs[1].set_ylim(0,120)
bottom_axs[1].set_title("selected model",fontsize=size_labels)
bottom_axs[1].set_ylabel("probes",fontsize=size_labels,labelpad=-0.1)

add_labels(range(0,5),hist_model,bottom_axs[1])

bottomfig.subplots_adjust(bottom=bottom_s+0.05,top=top_s,wspace=0.4,
                          left=left_s,right=right_s)

plt.draw()

x_a,y_a=[x-0.01 for x in subfig_bbox(topfig,fig)]
x_b,y_b=[x-0.01 for x in subfig_bbox(bottomfig,fig)]

fig.text(x_label,y_a,"a)", fontsize=size_subfig, ha="center", va="top")
fig.text(x_label,y_b,"b)", fontsize=size_subfig, ha="center", va="top")
fig.text(0.5+x_label,y_b,"c)", fontsize=size_subfig, ha="center", va="top")

plt.show()

fig.savefig(outdir+"figure4.pdf",
    dpi=300,             
    bbox_inches="tight", 
    pad_inches=0.02,     
    backend="pdf",
)

##############################################################
########## Fig 6
##############################################################
from matplotlib.ticker import FormatStrFormatter

#load data
trios="screen_full/trios_cMeth.ped"
trios=np.loadtxt(trios)
ids=trios[:,0].tolist()

######### load files
#load siblings list
sibs_dir="h2_siblings/0data/out_family/siblings.ped"
sibs=np.loadtxt(sibs_dir)

sibs_trios_ind=[i for i in range(0,len(sibs)) if sibs[i,0] in ids] 
sibs_trios=sibs[sibs_trios_ind] #1089

y_dir="screen_full/"
x_dir="screen_full/X_matrix/out_x/"

list_sibs="h2_siblings/0data/out_family/siblings.ped"
list_sibs=np.loadtxt(list_sibs)

trios_ids="screen_full/trios_cMeth.ped"
trios=np.loadtxt(trios_ids)

beta_dir="screen_full/"

probes=pd.read_csv("POE_mQTLs/outputs/POE_probes.csv") #filtered probes screen 2 POE-probes
mQTLs=pd.read_csv("POE_mQTLs/outputs/mQTLs_patternAnnotate.csv")
######################

def loadYX(chro,y_dir,y_ind,x_dir,x_ind):
    y_dir=y_dir+str(chro)+"/out_matrix/"+str(y_ind)+"_y.zarr/"
    y=zarr.load(y_dir)
    #load vcf file
    
    x_dir=x_dir+"x_"+str(chro)+".zarr/"
    x=zarr.load(x_dir)
    x_snp=x[:,x_ind*4:x_ind*4+4]
    
    return(y,x_snp)

def adjustYGenome(chro,beta_dir,x_ind,x_dir,y,components=[0,1,2,3]): 
    zf = zipfile.ZipFile(beta_dir)
    beta=np.genfromtxt(zf.open('mean_beta.csv'), delimiter=',')
    beta=np.delete(beta,(0),axis=0) #intercept value in 0
    beta_flat=beta.flatten()
    
    indexes=[x_ind*4+i for i in components]
    betaj=beta_flat[indexes] 
   
    x_dir=x_dir+"x_"+str(chro)+".zarr/"
    x=zarr.load(x_dir)
    x_j=x[:,indexes] 
    
    ### others residualize
    beta_flat=np.delete(beta_flat,indexes)
    x=np.delete(x,indexes,axis=1)
        
    y_adj=np.matmul(x,beta_flat)
    y_adj=y-y_adj

    yj=np.matmul(x_j,betaj)
    
    return(y_adj,yj)
    
def plotBarPlotsY(y,x_snp,name="",hist=False,axes=None): 
    if axes is None:
       fig, ax = plt.subplots() #unique plot
    else:
        ax=axes

    #get index of 0,01,10,2
    h0=[]
    h2=[]
    hm=[] #mother +1 01
    hf=[] #father -1 10

    for j in range(0,len(x_snp)): #per children
        c=x_snp[j,0];i=x_snp[j,3]
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
    labels=["0","1|0","0|1","2"]
    
    #errors
    down=[np.std(x)/np.sqrt(np.size(x)) for x in [y0,yf,ym,y2]]
    up=[np.std(x)/np.sqrt(np.size(x)) for x in [y0,yf,ym,y2]]
    error=[down,up]
    #plot bar 
    ax.bar(labels,[y0_mean,yf_mean,ym_mean,y2_mean],color=["grey","blue","red","grey"],alpha=0.5)
    ax.errorbar(labels,[y0_mean,yf_mean,ym_mean,y2_mean],yerr=error,fmt='.',color='Black',alpha=0.5)
    ax.set_title(name,fontsize=size_labels)
    ax.set_xlabel("genotype",fontsize=size_axes,labelpad=-1)
    ax.set_ylabel("methylation",fontsize=size_axes,labelpad=-1)
    min_val=np.min([[y0_mean,yf_mean,ym_mean,y2_mean][i]-down[i] for i in
                    range(0,len(up))])
    max_val=np.max([[y0_mean,yf_mean,ym_mean,y2_mean][i]+up[i] for i in
                    range(0,len(up))])
    
    ax.set_yticks([min_val,0,max_val])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

def plotParentals(ym,yf,x_snp,axes=None): 
    if axes is None:
       fig, axs = plt.subplots() #unique plot
    else:
        axs=axes
    #get index of 0,1,2 for maternal and paternal genotypes
    m01=[];m10=[];f10=[];f01=[]
    
    for j in range(0,len(x_snp)):
        i=x_snp[j,3];m=x_snp[j,1];f=x_snp[j,2]#only phenotypes from heterozygous children
        #####parentals
        if i==-1 and m==1: #heterozygous fahter 1|0
            m10.append(j)
        elif i==1 and m==1:
            m01.append(j)
        elif i==-1 and f==1: #heterozygous mother 0|1 
            f10.append(j)
        elif i==1 and f==1:
            f01.append(j)
    
    #get distribution phenotypes
    ym10=ym[m10];ym01=ym[m01];yf10=yf[f10];yf01=yf[f01]
    y_list=[ym10,ym01,yf10,yf01]
    #averages
    y_means=[np.mean(x) for x in y_list]

    #errors
    down=[np.std(x)/np.sqrt(np.size(x)) for x in y_list]
    up=[np.std(x)/np.sqrt(np.size(x)) for x in y_list]
    error=[down,up]
    
    x=[0,1,2,3]
    labels=["1|0","0|1","1|0","0|1"] #parental genotypes
    ####PLOTS 
    ax=axs
    ax.bar(x,y_means,color=["blue","red","blue","red"],alpha=0.5)
    ax.errorbar(x,y_means,yerr=error,fmt='.',color='Black')
    ax.axvline(1.5,color="black")
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels,fontsize=size_axes) 
    ax.set_xlabel("genotype",labelpad=-0.5,fontsize=size_axes)
    
    ax.text(0.5, 1.02, "M=1",transform=ax.get_xaxis_transform(),
            ha="center",fontsize=size_axes)
    ax.text(2.5, 1.02, "F=1",transform=ax.get_xaxis_transform(),
            ha="center",fontsize=size_axes)
   
    ax.set_xlim(-0.5,3.5)
    y_lim=np.max(np.abs(y_means))
    ax.set_ylim(y_lim*-1-0.2,y_lim+0.2)
    ax.set_yticks(np.linspace(y_lim*-1,y_lim,3))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f')) 
    ax.set_ylabel("methylation",fontsize=size_axes,labelpad=-1)
    ax.axhline(0,color="black",linestyle="--",alpha=0.8)
      
def diffSibs_geno(sibs,y,x_snp,ids,name="",axes=None): 
    sibs_neg=[];sibs_0=[];sibs_pos=[]
    
    #get genotype of siblings
    for s1,s2 in zip(sibs[:,0],sibs[:,1]):
        s1=int(s1);s2=int(s2)
        
        if s1 in ids and s2 in ids:
            ind_s1=ids.index(s1)
            ind_s2=ids.index(s2)
            
            ##imprinting sign
            s1_i=x_snp[ind_s1][3]
            s2_i=x_snp[ind_s2][3]
            ##genotype
            s1_c=x_snp[ind_s1][0]
            s2_c=x_snp[ind_s2][0]
            
            if s1_i==-1:
                if s2_i==-1:
                    sibs_0.append(y[ind_s1]-y[ind_s2])
                if s2_c==0:
                    sibs_neg.append(y[ind_s1]-y[ind_s2])
                if s2_c==2:
                    sibs_0.append(y[ind_s1]-y[ind_s2])
                if s2_i==1:
                    sibs_neg.append(y[ind_s1]-y[ind_s2])
            if s1_c==0:
                if s2_i==-1:
                    sibs_pos.append(y[ind_s1]-y[ind_s2])
                if s2_i==0:
                    sibs_0.append(y[ind_s1]-y[ind_s2])
                if s2_i==1:
                    sibs_neg.append(y[ind_s1]-y[ind_s2])
            if s1_c==2:
                sibs_0.append(y[ind_s1]-y[ind_s2])
            if s1_i==1:
                if s2_i==-1:
                    sibs_pos.append(y[ind_s1]-y[ind_s2])
                if s2_c==0:
                    sibs_pos.append(y[ind_s1]-y[ind_s2])
                if s2_c==2:
                    sibs_0.append(y[ind_s1]-y[ind_s2])
                if s2_i==1:
                    sibs_0.append(y[ind_s1]-y[ind_s2])
        else:
            pass
    
    if len(sibs_0)==0 or len(sibs_pos)==0 or len(sibs_neg)==0:
        print("No siblings with given genotypes for SNP")
        return()
    else:
        #plot vilin
        if axes is None:
           fig, ax = plt.subplots() 
        else:
            ax=axes
            
        parts=ax.violinplot([sibs_neg,sibs_0,sibs_pos],showmedians=True)
        
        colors=["blue","grey","red"]
        #color each violin
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_alpha(0.5) 
        # Set the color of the median lines
        parts['cmedians'].set_colors(colors) 
        parts["cmedians"].set_alpha(1)
        parts["cmins"].set_colors(colors)
        parts["cmaxes"].set_colors(colors)
        parts["cbars"].set_colors(colors)
        
        ax.axhline(0,color="black",linestyle="--",alpha=0.8)
        
        labels=["-1","0","+1"]#"All siblings"
        ax.set_xticks(range(1, len(labels)+1))
        ax.set_xticklabels(labels,fontsize=size_axes)
        
        ax.set_ylabel(u"Δ siblings", labelpad=-2,fontsize=size_axes) 
        ax.set_xlabel(u"Δ allele origin",fontsize=size_axes,labelpad=-0.8) 
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f')) 
        
        ax.set_xlim(0.5,3.5)
        ax.set_ylim(-4,4)

################
#### Plot simulations
imp="8_M+"
ind="109_maternal"
beta_imp="sim_final/out_beta_y/beta_"+imp+".zarr/"
beta_ind="sim_final/indirect/out_indirect/beta_"+ind+".zarr/"

#indirect and imprinting
### Example indirect effects
y_dir_ind="sim_final/indirect/out_indirect/y_"+ind+".zarr/"
y_dir_imp="sim_final/out_beta_y/y_"+imp+".zarr/"

x_dir_sim="sim_final/chr22_files/x_22.zarr/"
x_sim=zarr.load(x_dir_sim)
###indirect  
x_ind=5073 
y_ind=zarr.load(y_dir_ind)
x_snp_ind=x_sim[:,x_ind*4:x_ind*4+4]

######################################## FIG 6
###put together both simulations
set_style(textwidth_pts=418.0,ratio=0.9)
fig=plt.figure() 

###
subfigs=fig.subfigures(2,1,hspace=-0.08) #two vertical plots 

###### sim imprinting
###imprinting 
x_imp=4209 
y_imp=zarr.load(y_dir_imp)
x_snp_imp=x_sim[:,x_imp*4:x_imp*4+4]

(leftfig,rightfig)=subfigs[0].subfigures(1,2)  

mainfig,sepfig=leftfig.subfigures(1,2) 

main_axs=mainfig.subplots()
plotBarPlotsY(y_imp,x_snp_imp,axes=main_axs)

sep_axs=sepfig.subplots(2,1)
plotParentals(y_imp,y_imp,x_snp_imp,axes=sep_axs[0])
diffSibs_geno(sibs_trios,y_imp,x_snp_imp,ids,axes=sep_axs[1]) 

mainfig.subplots_adjust(left=0.45,right=1.05,bottom=0.15,top=0.85)
sepfig.subplots_adjust(left=0.35,right=0.95,bottom=0.15,top=0.85,hspace=0.5)

leftfig.suptitle("Paternal imprinting",fontsize=size_labels,y=0.98,x=0.58)
subfigs[0].supylabel("Simulations",fontsize=size_labels,x=0.01)

#########################
##### sim indirect
mainfig,sepfig=rightfig.subfigures(1,2) 

main_axs=mainfig.subplots()
plotBarPlotsY(y_ind,x_snp_ind,axes=main_axs)

sep_axs=sepfig.subplots(2,1)
plotParentals(y_ind,y_ind,x_snp_ind,axes=sep_axs[0])
diffSibs_geno(sibs_trios,y_ind,x_snp_ind,ids,axes=sep_axs[1]) 

mainfig.subplots_adjust(left=0.45,right=1.05,bottom=0.15,top=0.85)
sepfig.subplots_adjust(left=0.35,right=0.95,bottom=0.15,top=0.85,hspace=0.5)

rightfig.suptitle("Indirect maternal effect",fontsize=size_labels,y=0.98,x=0.58)

####################
#### real data
#### indirect
chro=4;y_ind=35417;x_ind=16381;name=""
y,x_snp=loadYX(chro,y_dir,y_ind,x_dir,x_ind)
y_adj,y_snp=adjustYGenome(chro,beta_dir+str(chro)+"/jodie/out/"+str(y_ind)+
                          "/mean_beta.csv.zip",x_ind,x_dir,y,[0,1,2,3]) #all snp

y_adjM,yM=adjustYGenome(chro,beta_dir+str(chro)+"/jodie/out/"+str(y_ind)+"/mean_beta.csv.zip",x_ind,x_dir,y,[1])
y_adjF,yF=adjustYGenome(chro,beta_dir+str(chro)+"/jodie/out/"+str(y_ind)+"/mean_beta.csv.zip",x_ind,x_dir,y,[2])
y_adjI,yI=adjustYGenome(chro,beta_dir+str(chro)+"/jodie/out/"+str(y_ind)+"/mean_beta.csv.zip",x_ind,x_dir,y,[3])


(leftfig,rightfig)=subfigs[1].subfigures(1,2)  #f imp and mat indirect
mainfig,sepfig=rightfig.subfigures(1,2) 

main_axs=mainfig.subplots()
plotBarPlotsY(y_adj,x_snp,name,axes=main_axs)

sep_axs=sepfig.subplots(2,1)
plotParentals(y_adjM,y_adjF, x_snp,axes=sep_axs[0])
diffSibs_geno(sibs_trios,y_adjI,x_snp,ids,axes=sep_axs[1])

mainfig.subplots_adjust(left=0.45,right=1.05,bottom=0.15,top=0.85)
sepfig.subplots_adjust(left=0.35,right=0.95,bottom=0.15,top=0.85,hspace=0.5)


###### example meg3
chro=14;y_ind=22326;x_ind=1631;name=""
snp="rs1004574"

y_dir="screen_full/"
x_dir="LD_fineMap/3matrix/X_matrix/out_x/"

beta_dir="LD_fineMap/4run_regression/out/"

file_mQTLs="LD_fineMap/3.5_1Mb/mQTLs.txt"
file_mQTLs=pd.read_csv(file_mQTLs,sep="\t",header=None)


def loadYX(chro,snp,y_dir,y_ind,x_dir,x_ind):
    y_dir=y_dir+str(chro)+"/out_matrix/"+str(y_ind)+"_y.zarr/"
    y=zarr.load(y_dir)
    
    #load vcf file
    x_dir=x_dir+str(chro)+"_"+snp+".zarr/"
    x=zarr.load(x_dir)
    
    x_snp=x[:,x_ind*4:x_ind*4+4]
    
    return(y,x_snp,x)

def adjustYGenome(chro,beta_dir,x,x_ind,x_dir,y,components=[0,1,2,3]): ## all genome except beta of specific snp an component
    zf = zipfile.ZipFile(beta_dir)
    beta=np.genfromtxt(zf.open('mean_beta.csv'), delimiter=',')
    beta=np.delete(beta,(0),axis=0) #intercept value in 0
    beta_flat=beta.flatten()
    
    indexes=[x_ind*4+i for i in components]
    betaj=beta_flat[indexes] #b^_j
   
    x_j=x[:,indexes] 
    
    ### others residualize
    beta_flat=np.delete(beta_flat,indexes)
    x=np.delete(x,indexes,axis=1)
        
    y_adj=np.matmul(x,beta_flat)
    y_adj=y-y_adj
    
    ##no noise. deterministic
    yj=np.matmul(x_j,betaj)
    
    return(y_adj,yj)
  
y,x_snp,x=loadYX(chro,snp,y_dir,y_ind,x_dir,x_ind)

numQTL=file_mQTLs.index[(file_mQTLs[0]==snp) & (file_mQTLs[3]==y_ind)].tolist()[0]+1
y_adj,y_snp=adjustYGenome(chro,beta_dir+str(numQTL)+"/mean_beta.csv.zip",
                          x,x_ind,x_dir,y,[0,1,2,3]) #all snp

y_adjM,yM=adjustYGenome(chro,beta_dir+str(numQTL)+"/mean_beta.csv.zip",
                        x,x_ind,x_dir,y,[1])
y_adjF,yF=adjustYGenome(chro,beta_dir+str(numQTL)+"/mean_beta.csv.zip",
                        x,x_ind,x_dir,y,[2])
y_adjI,yI=adjustYGenome(chro,beta_dir+str(numQTL)+"/mean_beta.csv.zip",
                        x,x_ind,x_dir,y,[3])

mainfig,sepfig=leftfig.subfigures(1,2) 
main_axs=mainfig.subplots()
plotBarPlotsY(y_adj,x_snp,name,axes=main_axs)


sep_axs=sepfig.subplots(2,1)
plotParentals(y_adjM,y_adjF, x_snp,axes=sep_axs[0])
diffSibs_geno(sibs_trios,y_adjI,x_snp,ids,axes=sep_axs[1])

mainfig.subplots_adjust(left=0.45,right=1.05,bottom=0.15,top=0.85)
sepfig.subplots_adjust(left=0.35,right=0.95,bottom=0.15,top=0.85,hspace=0.5)

subfigs[1].supylabel("Real data",fontsize=size_labels,x=0.01)

plt.draw()


plt.show()

fig.savefig(outdir+"figure6.pdf",
    dpi=300,            
    backend="pdf",
)

