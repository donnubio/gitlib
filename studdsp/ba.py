import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import FloatVector
import os
import plotly.express as px
import plotly.graph_objects as go

def ba(y1,y2, bs_nrep=1000, fl_rerr=False, ba_ci_conf=0.95, ci_conf=0.95, true_val_type="mean"):

    file_path = os.path.dirname(os.path.realpath(__file__))
    ro.r['source'](file_path + '/bameth.R')
    fun_bacalc_r = ro.globalenv['ba_boot_from_expr']
    y1_r = FloatVector( y1 )
    y2_r = FloatVector( y2 )
    ba_r = fun_bacalc_r(y1_r, y2_r, bs_nrep, False, fl_rerr, ba_ci_conf, ci_conf, true_val_type) #(y1,y2, bs_nrep, fl_trsfrm=F, fl_rerr=F, ba_ci_conf=0.95, ci_conf=0.95)
    D = np.array( ba_r.rx('BA_D') ).flatten()
    A = np.array( ba_r.rx('BA_A') ).flatten()
    ba_pr_names = list( np.array( ba_r.rx('BA_NMS_BA_PR') ).flatten() )
    ba_pr_ = np.array( ba_r.rx('BA_PR_0') ).flatten()
    ba_pr_L_ = np.array( ba_r.rx('BA_PR_PC_L') ).flatten() #perc
    ba_pr_R_ = np.array( ba_r.rx('BA_PR_PC_R') ).flatten() #perc
    ba_pr = {ba_pr_names[i]: ba_pr_[i] for i in range(len(ba_pr_names))}
    ba_pr_L = {ba_pr_names[i]: ba_pr_L_[i] for i in range(len(ba_pr_names))}
    ba_pr_R = {ba_pr_names[i]: ba_pr_R_[i] for i in range(len(ba_pr_names))}    

    return(D,A, ba_pr, ba_pr_L, ba_pr_R)


def BAPlot(A,D,Bias,LoA_U,LoA_L,
           CIBias=None, CILoA_U=None, CILoA_L=None,
           h=300, w=None,
           xlabel=None, ylabel=None, title=None, fontsize=14, showlegend=False):
    layout = go.Layout(autosize=True, width=w, height=h,  margin=dict(l=10, r=10, t=30, b=10), template='simple_white')
    fig = go.Figure(layout=layout)
    fig.add_scatter(x=A,y=D,name="aa",mode="markers")
    fig.add_hline(y=Bias, line = dict(color='royalblue', width=2, dash='solid'))      
    fig.add_hline(y=LoA_U, line = dict(color='red', width=2, dash='dash')) 
    fig.add_hline(y=LoA_L, line = dict(color='red', width=2, dash='dash')) 
    fig.add_trace( go.Scatter( x=[min(A) + (max(A)-min(A))*0.75] * 3, 
                               y=[Bias,LoA_U,LoA_L],
                               text=['Bias='+str(round(Bias,3)), 'LoAU='+str(round(LoA_U,3)), 'LoAL='+str(round(LoA_L,3))],
                               mode='text',
                               textposition=["bottom right","bottom right","bottom right"])) #top
    fig.update_layout(xaxis_title=xlabel,
                    yaxis_title=ylabel,
                    #legend_title="Legend Title",
                    title=title,
                    font=dict(
                    #    family="Courier New, monospace",
                        size=fontsize,
                    #    color="RebeccaPurple"
                    ),
                    showlegend=showlegend 
                    )

    if CIBias:
        x = [min(A), max(A)]; y_upper = [CIBias[1]]*2; y_lower = [CIBias[0]]*2
        fig.add_scatter( x=x+x[::-1], # x, then x reversed
                        y=y_upper+y_lower[::-1], # upper, then lower reversed
                        fill='toself',
                        fillcolor='rgba(0,100,80,0.1)',
                        line=dict(color='rgba(255,255,255,0)'),
                        hoverinfo="skip",
                        showlegend=False
                        )
    if CILoA_U:
        x = [min(A), max(A)]; y_upper = [CILoA_U[1]]*2; y_lower = [CILoA_U[0]]*2
        fig.add_scatter( x=x+x[::-1], # x, then x reversed
                        y=y_upper+y_lower[::-1], # upper, then lower reversed
                        fill='toself',
                        fillcolor='rgba(0,100,80,0.1)',
                        line=dict(color='rgba(255,255,255,0)'),
                        hoverinfo="skip",
                        showlegend=False
                        )        

    if CILoA_L:
        x = [min(A), max(A)]; y_upper = [CILoA_L[1]]*2; y_lower = [CILoA_L[0]]*2
        fig.add_scatter( x=x+x[::-1], # x, then x reversed
                        y=y_upper+y_lower[::-1], # upper, then lower reversed
                        fill='toself',
                        fillcolor='rgba(0,100,80,0.1)',
                        line=dict(color='rgba(255,255,255,0)'),
                        hoverinfo="skip",
                        showlegend=False
                        )     
    fig.show()
