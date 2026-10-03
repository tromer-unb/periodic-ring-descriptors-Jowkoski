from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch
from matplotlib.gridspec import GridSpec

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"; FIG=ROOT/"figures"; FIG.mkdir(exist_ok=True)

C={"navy":"#17324D","blue":"#3B6FB6","teal":"#209486","orange":"#E8913A",
   "coral":"#D95D55","purple":"#7668B5","grey":"#8793A1","light":"#E8EDF2",
   "dark":"#101820","gold":"#D9AA28"}

mpl.rcParams.update({
 "font.family":"DejaVu Sans","font.size":11.0,
 "axes.titlesize":11.5,"axes.titleweight":"bold",
 "axes.labelsize":11.0,"axes.labelweight":"bold",
 "xtick.labelsize":9.8,"ytick.labelsize":9.8,
 "legend.fontsize":9.0,"figure.dpi":160,
 "pdf.fonttype":42,"ps.fonttype":42,
 "axes.linewidth":1.15,"xtick.major.width":1.0,"ytick.major.width":1.0,
})
def clean(ax,grid=False):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out",length=4,width=1.0)
    if grid:
        ax.grid(axis="x",lw=.55,alpha=.20); ax.set_axisbelow(True)
    for q in ax.get_xticklabels()+ax.get_yticklabels(): q.set_fontweight("bold")
def letter(ax,s):
    return None
def save(fig,name):
    fig.savefig(FIG/f"{name}.pdf",bbox_inches="tight",pad_inches=.06)
    fig.savefig(FIG/f"{name}.png",dpi=600,bbox_inches="tight",pad_inches=.06)
    plt.close(fig)
def hbars(ax,names,vals,colors,errs=None,xlabel=None,xmax=None,fmt="{:.3f}"):
    y=np.arange(len(names))
    kw={} if errs is None else {"xerr":errs,"capsize":3}
    bars=ax.barh(y,vals,color=colors,height=.62,**kw)
    ax.set_yticks(y,names); ax.invert_yaxis()
    if xlabel: ax.set_xlabel(xlabel)
    if xmax is None: xmax=max(vals)*1.20
    ax.set_xlim(0,xmax); clean(ax,grid=True)
    for b,v in zip(bars,vals):
        ax.text(v+xmax*.018,b.get_y()+b.get_height()/2,fmt.format(v),
                va="center",ha="left",fontsize=9.4,fontweight="bold",color=C["dark"])
    return bars
# Fig. 1 — descriptor concept
fig=plt.figure(figsize=(7.25,6.35))
gs=GridSpec(2,2,figure=fig,hspace=.34,wspace=.25)
ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
ang=np.arange(6)*np.pi/3
centers=[(-1.45,0),(0,0),(1.45,0),(-.725,1.26),(.725,1.26),(-.725,-1.26),(.725,-1.26)]
for j,(cx,cy) in enumerate(centers):
    pts=np.c_[cx+np.cos(ang),cy+np.sin(ang)]
    ax.add_patch(Polygon(pts,closed=True,facecolor=C["orange"] if j==1 else "none",
                         edgecolor=C["navy"],lw=2.0,alpha=.28 if j==1 else .55))
ax.scatter(np.cos(ang),np.sin(ang),s=42,c=C["navy"],zorder=5)
ax.text(0,-1.72,"periodic face",ha="center",fontweight="bold",color=C["orange"],fontsize=11)
ax.set_title("a   2D: planar ring",loc="left"); ax.set_aspect("equal"); ax.set_xlim(-2.75,2.75); ax.set_ylim(-2.0,2.05); ax.axis("off")

ax=fig.add_subplot(gs[0,1],projection="3d")
t=np.linspace(0,2*np.pi,6,endpoint=False); x=np.cos(t); y=.76*np.sin(t); z=.42*np.sin(2*t)
ax.plot(np.r_[x,x[0]],np.r_[y,y[0]],np.r_[z,z[0]],lw=3.2,c=C["teal"])
ax.scatter(x,y,z,s=38,c=C["navy"],depthshade=False)
ax.plot([x[0],x[1]],[y[0],y[1]],[z[0],z[1]],lw=5.2,c=C["coral"])
ax.set_title("b   3D: periodic shortest ring",loc="left",pad=0); ax.set_axis_off(); ax.view_init(24,-58)

ax=fig.add_subplot(gs[1,0]); letter(ax,"c")
raw=np.array([[-1.1,-.15],[-.4,-.9],[.55,-.72],[1.15,.05],[.42,.95],[-.62,.78]])
raw-=raw.mean(0); rr=np.mean(np.linalg.norm(raw,axis=1)); zc=(raw[:,0]+1j*raw[:,1])/rr
a=.25; wt=rr*(zc+a*a/zc); tr=np.c_[wt.real,wt.imag]
ax.plot(*np.r_[raw,raw[:1]].T,lw=2.6,c=C["navy"],label="raw")
ax.plot(*np.r_[tr,tr[:1]].T,lw=2.6,c=C["orange"],label="Joukowsky")
ax.scatter(raw[:,0],raw[:,1],s=30,c=C["navy"]); ax.scatter(tr[:,0],tr[:,1],s=30,c=C["orange"])
ax.set_title("c   Local-plane shape response",loc="left"); ax.set_aspect("equal"); clean(ax)
ax.set_xticks([]); ax.set_yticks([]); ax.legend(frameon=False,loc="upper left",ncol=2)

ax=fig.add_subplot(gs[1,1]); letter(ax,"d"); ax.axis("off"); ax.set_title("d   Structure-level fingerprint",loc="left")
blocks=[("TOPOLOGY","ring sizes\ncoordination",C["teal"]),
        ("GEOMETRY","area · perimeter\nradius · anisotropy",C["blue"]),
        ("J RESPONSE",r"$A_J, P_J, eta_J$"+"\n"+r"$R_A, R_P$",C["orange"])]
for i,(ttl,txt,col) in enumerate(blocks):
    y0=.68-i*.25
    box=FancyBboxPatch((.08,y0),.84,.20,boxstyle="round,pad=.018,rounding_size=.02",
        transform=ax.transAxes,facecolor=col,edgecolor="none")
    ax.add_patch(box)
    ax.text(.13,y0+.155,ttl,transform=ax.transAxes,color="white",fontsize=10.8,fontweight="bold",va="center")
    ax.text(.54,y0+.045,txt,transform=ax.transAxes,color="white",fontsize=10.0,fontweight="bold",va="center",ha="center")
ax.text(.50,.10,"means  ·  dispersions  ·  ring fractions",transform=ax.transAxes,ha="center",
        fontsize=10.3,fontweight="bold",color=C["dark"])
fig.subplots_adjust(left=.06,right=.98,bottom=.05,top=.96)
save(fig,"fig1_descriptor")
# Fig. 2 — invariance
poly=json.loads((DATA/"polymorph_results.json").read_text())
D=next(q for q in poly if q["name"]=="diamond"); L=next(q for q in poly if q["name"]=="lonsdaleite")
fig=plt.figure(figsize=(7.25,5.35))
gs=GridSpec(2,2,figure=fig,height_ratios=[1,1.1],hspace=.42,wspace=.34)
ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
names=["primitive\ncell","64-atom\ncell"]
y=np.arange(2); h=.22
vals=[np.array([0,36/64]),np.array([0,56/64]),np.array([2,2])]
labs=["cycle basis","minimum basis","periodic shortest"]
cols=[C["grey"],C["purple"],C["teal"]]
for i,(v,lab,col) in enumerate(zip(vals,labs,cols)):
    ax.barh(y+(i-1)*h,v,height=h,color=col,label=lab)
ax.set_yticks(y,names); ax.invert_yaxis(); ax.set_xlabel("rings per atom")
ax.set_title("a   Cell dependence of finite cycle bases",loc="left"); ax.set_xlim(0,2.25); clean(ax,grid=True)
ax.legend(frameon=False,fontsize=8.1,loc="lower right")

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
labels=["diamond\nprimitive","diamond\nconventional","diamond\nsupercell","lonsdaleite","lonsdaleite\n2×2×2"]
vals=np.array([2,2,2,2,2],float)
ax.scatter(np.arange(5),vals,s=72,c=[C["navy"]]*3+[C["orange"]]*2,zorder=4)
ax.plot(np.arange(3),vals[:3],lw=2.1,c=C["navy"]); ax.plot(np.arange(3,5),vals[3:],lw=2.1,c=C["orange"])
ax.axhline(2,color=C["grey"],lw=1,ls="--")
ax.set_xticks(np.arange(5),["D-prim","D-conv","D-super","Lons","Lons-2×"])
ax.set_ylabel("rings per atom"); ax.set_ylim(1.86,2.14); ax.set_title("b   Cell-invariant normalized ring count",loc="left"); clean(ax)

ax=fig.add_subplot(gs[1,:]); letter(ax,"c")
features=[("area_mean","projected area"),("perimeter_mean","perimeter"),
          ("anisotropy_mean","anisotropy"),("j_area_ratio_mean",r"$R_A$")]
ratio=np.array([L[k]/D[k] for k,_ in features])
x=np.arange(len(features)); bars=ax.bar(x,ratio,color=[C["blue"],C["blue"],C["orange"],C["orange"]],width=.62)
ax.axhline(1,color=C["grey"],lw=1.1,ls="--")
ax.set_xticks(x,[lab for _,lab in features]); ax.set_ylabel("lonsdaleite / diamond")
ax.set_ylim(.92,1.12); ax.set_title("c   Polymorph sensitivity at fixed ring size",loc="left"); clean(ax)
for b,v in zip(bars,ratio):
    ax.text(b.get_x()+b.get_width()/2,v+(.006 if v>=1 else -.014),f"{v:.3f}",
            ha="center",va="bottom" if v>=1 else "top",fontweight="bold",fontsize=10)
fig.subplots_adjust(left=.10,right=.98,bottom=.09,top=.95)
save(fig,"fig2_invariance")
# Fig. 3 — 2D results
m2=pd.read_csv(DATA/"metrics_2d_allotropes.csv",header=[0,1],index_col=0)
a2=pd.read_csv(DATA/"ablation_2d.csv",header=[0,1],index_col=0)
d2=pd.read_csv(DATA/"metrics_2d_defects.csv",header=[0,1],index_col=[0,1])
fig=plt.figure(figsize=(7.25,5.55))
gs=GridSpec(2,2,figure=fig,height_ratios=[1,1.05],hspace=.45,wspace=.42)

ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
order=["mean_baseline","metadata","descriptor","combined"]
names=["mean baseline",r"$N+\rho$","ring descriptor","ring + metadata"]
vals=[m2.loc[k,("mae","mean")] for k in order]; err=[m2.loc[k,("mae","std")] for k in order]
hbars(ax,names,vals,[C["grey"],C["gold"],C["teal"],C["navy"]],err,
      "MAE (native target units)",.265)
ax.set_title("a   120 carbon allotropes",loc="left")

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
order=["topology","topology_j","topology_raw","topology_raw_j"]
names=["topology","topology + J","topology + raw","topology + raw + J"]
vals=[a2.loc[k,("mae","mean")] for k in order]; err=[a2.loc[k,("mae","std")] for k in order]
hbars(ax,names,vals,[C["grey"],C["orange"],C["blue"],C["navy"]],err,"MAE",.165)
ax.set_title("b   2D feature ablation",loc="left")

ax=fig.add_subplot(gs[1,:]); letter(ax,"c")
classes=["all","vacancy","B","N"]; xx=np.arange(len(classes)); w=.22
series=[("count_type_baseline",C["grey"],"count/type"),("jowkoski",C["teal"],"periodic ring"),
        ("matminer",C["navy"],"Matminer")]
for j,(model,col,lab) in enumerate(series):
    vals=[d2.loc[(model,k),("mae","mean")] for k in classes]
    ax.bar(xx+(j-1)*w,vals,width=w,color=col,label=lab)
ax.set_xticks(xx,["all defects","vacancy","B substitution","N substitution"])
ax.set_ylabel("formation-energy MAE (eV)"); ax.set_ylim(0,2.75)
ax.set_title("c   Defective graphene",loc="left"); clean(ax)
ax.legend(frameon=False,ncol=3,loc="upper right")
fig.subplots_adjust(left=.14,right=.98,bottom=.08,top=.95)
save(fig,"fig3_2d_results")

# Fig. 4 — 3D results
old=pd.read_csv(DATA/"metrics_3d_old.csv",header=[0,1],index_col=[0,1])
new=pd.read_csv(DATA/"metrics_3d_periodic.csv",header=[0,1],index_col=[0,1])
soap=pd.read_csv(DATA/"soap_3d_summary.csv",index_col=0)
a3=pd.read_csv(DATA/"ablation_3d.csv",header=[0,1],index_col=0)
fig=plt.figure(figsize=(7.25,5.65))
gs=GridSpec(2,2,figure=fig,height_ratios=[1,1.05],hspace=.45,wspace=.42)

ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
names=["count/type","legacy cycles","periodic rings","SOAP","Matminer"]
vals=[new.loc[("count_type_baseline","all"),("mae","mean")],
      old.loc[("jowkoski","all"),("mae","mean")],
      new.loc[("jowkoski","all"),("mae","mean")],
      soap.loc["mean","mae"],new.loc[("matminer","all"),("mae","mean")]]
hbars(ax,names,vals,[C["grey"],C["coral"],C["teal"],C["purple"],C["navy"]],
      None,"formation-energy MAE (eV)",.86)
ax.set_title("a   3D defect benchmark",loc="left")

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
order=["topology","topology_j","topology_raw","topology_raw_j"]
names=["topology","topology + J","topology + raw","topology + raw + J"]
vals=[a3.loc[k,("mae","mean")] for k in order]; err=[a3.loc[k,("mae","std")] for k in order]
hbars(ax,names,vals,[C["grey"],C["orange"],C["blue"],C["navy"]],err,"MAE (eV)",.46)
ax.set_title("b   3D feature ablation",loc="left")

ax=fig.add_subplot(gs[1,:]); letter(ax,"c")
classes=["vacancy","B","N"]; xx=np.arange(3); w=.22
for j,(model,col,lab) in enumerate(series):
    vals=[new.loc[(model,k),("mae","mean")] for k in classes]
    ax.bar(xx+(j-1)*w,vals,width=w,color=col,label=lab)
ax.set_xticks(xx,["vacancy","B substitution","N substitution"])
ax.set_ylabel("formation-energy MAE (eV)"); ax.set_ylim(0,1.85)
ax.set_title("c   Performance by defect chemistry",loc="left"); clean(ax)
ax.legend(frameon=False,ncol=3,loc="upper right")
fig.subplots_adjust(left=.14,right=.98,bottom=.08,top=.95)
save(fig,"fig4_3d_results")
# Fig. 5 — regime dependence and hybrids
dftb=pd.read_csv(DATA/"metrics_dftb.csv",header=[0,1],index_col=0)
hy=pd.read_csv(DATA/"hybrid_matminer.csv",header=[0,1],index_col=0)
fig=plt.figure(figsize=(7.25,4.35))
gs=GridSpec(1,2,figure=fig,wspace=.38)

ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
def_base=new.loc[("count_type_baseline","all"),("mae","mean")]
dftb_base=dftb.loc["mean_baseline",("mae","mean")]
models=["periodic ring","Matminer","SOAP"]; cols=[C["teal"],C["navy"],C["purple"]]
def_vals=[new.loc[("jowkoski","all"),("mae","mean")]/def_base,
          new.loc[("matminer","all"),("mae","mean")]/def_base,soap.loc["mean","mae"]/def_base]
dftb_vals=[dftb.loc["topology_raw_j",("mae","mean")]/dftb_base,
           dftb.loc["matminer",("mae","mean")]/dftb_base,dftb.loc["soap",("mae","mean")]/dftb_base]
xx=np.arange(2); w=.23
for j,(lab,col) in enumerate(zip(models,cols)):
    ax.bar(xx+(j-1)*w,[def_vals[j],dftb_vals[j]],width=w,color=col,label=lab)
ax.axhline(1,c=C["grey"],ls="--",lw=1.2)
ax.set_xticks(xx,["defect\nreconstruction","smooth\ndistortion"])
ax.set_ylabel("MAE / simple baseline"); ax.set_ylim(0,1.16)
ax.set_title("a   Regime dependence",loc="left"); clean(ax)
ax.legend(frameon=False,ncol=1,loc="upper left")

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
names=["Matminer","+ 10 J","+ ring raw","+ ring raw + J"]
vals=[new.loc[("matminer","all"),("mae","mean")],
      hy.loc["matminer_j",("mae","mean")],hy.loc["matminer_top_raw",("mae","mean")],
      hy.loc["matminer_top_raw_j",("mae","mean")]]
errs=[new.loc[("matminer","all"),("mae","std")],hy.loc["matminer_j",("mae","std")],
      hy.loc["matminer_top_raw",("mae","std")],hy.loc["matminer_top_raw_j",("mae","std")]]
hbars(ax,names,vals,[C["grey"],C["orange"],C["teal"],C["navy"]],errs,
      "3D defect MAE (eV)",.245)
ax.set_xlim(.195,.245)
for t in ax.texts: t.set_visible(False)
for yi,v in enumerate(vals):
    ax.text(v+.0010,yi,f"{v:.3f}",va="center",ha="left",fontsize=9.6,fontweight="bold")
ax.set_title("b   Hybrid complementarity",loc="left")
fig.subplots_adjust(left=.12,right=.98,bottom=.14,top=.92)
save(fig,"fig5_regimes_hybrids")

# Fig. 6 — robustness
hold=pd.read_csv(DATA/"holdout_3d.csv"); hs=pd.read_csv(DATA/"holdout_soap.csv")
cut=pd.read_csv(DATA/"cutoff_sensitivity_3d.csv")
cost=json.loads((DATA/"feature_cost_3d.json").read_text())
fig=plt.figure(figsize=(7.25,5.55))
gs=GridSpec(2,2,figure=fig,height_ratios=[1,1.08],hspace=.45,wspace=.38)

ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
ax.plot(cut.cutoff,100*cut.same_ring_distribution_fraction,marker="o",ms=7,lw=2.4,c=C["teal"])
ax.axvline(1.895,c=C["navy"],ls="--",lw=1.2)
ax.set_xlabel("bond cutoff (Å)"); ax.set_ylabel("unchanged ring sets (%)"); ax.set_ylim(96.7,100.25)
ax.set_title("a   Cutoff robustness",loc="left"); clean(ax)

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
for name,col,label in [("jowkoski",C["teal"],"ring"),("matminer",C["navy"],"Matminer"),("soap",C["purple"],"SOAP")]:
    x=cost[name]["n_features"]; y=cost[name]["mean_seconds"]
    ax.scatter(x,y,s=105,c=col,edgecolor="white",linewidth=.8,zorder=3)
    ax.annotate(label,(x,y),xytext=(6,-15 if name=="jowkoski" else 5),textcoords="offset points",fontweight="bold",fontsize=9)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("descriptor dimension"); ax.set_ylabel("extraction time (s)")
ax.set_title("b   Runtime audit",loc="left"); clean(ax)

ax=fig.add_subplot(gs[1,:]); letter(ax,"c")
classes=["all","vacancy","B","N"]; xx=np.arange(4); w=.18
sets=[("count_type_baseline",C["grey"],"count/type"),("jowkoski",C["teal"],"periodic ring"),
      ("matminer",C["navy"],"Matminer")]
for j,(model,col,lab) in enumerate(sets):
    vals=[hold[(hold.model==model)&(hold.defect_type==k)].mae.iloc[0] for k in classes]
    ax.bar(xx+(j-1.5)*w,vals,width=w,color=col,label=lab)
vals=[hs[hs.defect_type==k].mae.iloc[0] for k in classes]
ax.bar(xx+1.5*w,vals,width=w,color=C["purple"],label="SOAP")
ax.set_xticks(xx,["all","vacancy","B substitution","N substitution"])
ax.set_ylabel("MAE (eV)"); ax.set_ylim(0,4.0); ax.set_title("c   Concentration extrapolation (1–4 → 5 defects)",loc="left")
clean(ax); ax.legend(frameon=False,ncol=4,loc="upper right",fontsize=8.5)
fig.subplots_adjust(left=.11,right=.98,bottom=.08,top=.95)
save(fig,"fig6_robustness")
# Supplementary figures
scan=pd.read_csv(DATA/"scan_a_2d_summary.csv",header=[0,1],index_col=0)
corr=pd.read_csv(DATA/"joukowsky_a_correlations.csv")
fig=plt.figure(figsize=(7.25,4.0))
gs=GridSpec(1,2,figure=fig,wspace=.38)
ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
avec=np.array([float(str(k).split("=")[1]) for k in scan.index])
mae=scan[("mae","mean")].to_numpy()
ax.plot(avec,mae,marker="o",ms=7,lw=2.3,c=C["orange"])
raw_ref=pd.read_csv(DATA/"ablation_2d.csv",header=[0,1],index_col=0).loc["topology_raw",("mae","mean")]
ax.axhline(raw_ref,c=C["navy"],ls="--",lw=1.5,label="raw geometry")
ax.set_xlabel("Joukowsky parameter a"); ax.set_ylabel("MAE"); ax.set_title("a   2D parameter scan",loc="left")
clean(ax); ax.legend(frameon=False)

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
for raw,col,label in [("area_mean",C["blue"],"area"),("perimeter_mean",C["teal"],"perimeter"),
                      ("anisotropy_mean",C["orange"],"anisotropy")]:
    g=corr[corr.raw==raw]; ax.plot(g.a,g.pearson_r,marker="o",ms=6,lw=2,label=label,c=col)
ax.set_xlabel("Joukowsky parameter a"); ax.set_ylabel("Pearson r"); ax.set_ylim(.50,1.02)
ax.set_title("b   Raw–transformed correlation",loc="left"); clean(ax); ax.legend(frameon=False)
fig.subplots_adjust(left=.11,right=.98,bottom=.16,top=.91)
save(fig,"figS1_joukowsky_parameter")

summary=json.loads((DATA/"face_vs_mcb_2d_summary.json").read_text())
fig=plt.figure(figsize=(7.25,4.0))
gs=GridSpec(1,2,figure=fig,wspace=.38)
ax=fig.add_subplot(gs[0,0]); letter(ax,"a")
names=["same ring count","same ring-size multiset"]; vals=[summary["same_ring_count"],summary["same_ring_size_multiset"]]
hbars(ax,names,vals,[C["teal"],C["navy"]],None,"structures",130,fmt="{:.0f}")
for t in ax.texts: t.set_visible(False)
ax.set_title("a   Face tracing matches MCB",loc="left")
ax.text(.04,.08,r"max. mean-feature difference $<4\times10^{-15}$",transform=ax.transAxes,
        fontsize=9.5,fontweight="bold")

ax=fig.add_subplot(gs[0,1]); letter(ax,"b")
ax.plot(cut.cutoff,cut.mean_abs_delta_rings_per_atom,marker="o",ms=7,lw=2.2,c=C["teal"],label="rings / atom")
ax.set_xlabel("bond cutoff (Å)"); ax.set_ylabel(r"mean $|\Delta|$ rings / atom",color=C["teal"])
ax2=ax.twinx()
ax2.plot(cut.cutoff,cut.mean_abs_delta_j_area_ratio_mean,marker="s",ms=6,lw=2.0,c=C["orange"],label=r"$R_A$")
ax2.set_ylabel(r"mean $|\Delta R_A|$",color=C["orange"],fontweight="bold")
ax.set_title("b   Cutoff perturbation magnitude",loc="left"); clean(ax)
for q in ax2.get_yticklabels(): q.set_fontweight("bold")
fig.subplots_adjust(left=.13,right=.87,bottom=.16,top=.91)
save(fig,"figS2_algorithm_cutoff")
print("generated",len(list(FIG.glob("*.pdf"))),"PDF figures and",len(list(FIG.glob("*.png"))),"PNG figures")