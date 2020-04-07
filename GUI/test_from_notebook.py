import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import streamlit as st
import model
from post_regression import ci_bootstrap


def example(population, infected, removed,
            infected_time_series, alpha=0.95, beta=0.38):

    n_days = len(infected_time_series)
    days_list = np.linspace(0, n_days-1, n_days)

    P = population
    I = infected_time_series
    n_days = len(I)
    t_d = np.linspace(0, n_days-1, n_days)

    n_S0 = P-infected
    n_I0 = infected
    n_R0 = removed
    n0 = [n_S0, n_I0, n_R0]
    p = [alpha,beta]

    # Create empty model
    SIR = model.SIR()
    SIR.set_params(p, n0)
    # Train model
    SIR.fit(t_d, I, P)
    # Build numerical solution
    I_opt = SIR.solve(n_days-1, n_days).fetch()[:, 2]
    R_opt = SIR.r0 #

    # Get the confidence interval through bootstrap
    par_ci, par_list = ci_bootstrap(SIR, t_d, I, P, n_iter=1000)
    alpha_min = par_ci[0][0]
    alpha_max = par_ci[0][1]
    ci_dict = {"IC 95% for alpha": par_ci[0],
               "IC 95% for beta": par_ci[1],
               "IC 95% for r0": par_ci[2]}

    # Data for main plot
    beta_0 = SIR.p[1]
    SIR_minus = model.SIR().set_params([alpha_min, beta_0], n0)
    SIR_plus = model.SIR().set_params([alpha_max, beta_0], n0)
    I_minus = SIR_minus.solve(n_days-1, n_days).fetch()[:, 2]
    I_plus = SIR_plus.solve(n_days-1, n_days).fetch()[:, 2]

    # Data for subplots
    long_term_days = 90
    tf_long = long_term_days-1
    sol_long = SIR.solve(tf_long, long_term_days).fetch()
    N_S_long = sol_long[:, 1]*population
    N_I_long = sol_long[:, 2]*population
    N_R_long = sol_long[:, 3]*population
    tspan_long = np.linspace(0, tf_long, long_term_days)

    # Show main plot
    layout = go.Layout(paper_bgcolor='rgba(0,0,0,0)',
                       plot_bgcolor='rgba(0,0,0,0)',
                       title_text="Fitting of the SIR model" +
                       "against 15 days of UK data",
                       )
    fig = go.Figure(layout=layout)
    fig.add_trace(go.Scatter(x=t_d,
                             y=I,
                             name='UK Data',
                             ))
    fig.add_trace(go.Scatter(x=t_d,
                             y=I_opt,
                             fill=None,
                             mode='lines',
                             name='SIR',
                             ))
    fig.add_trace(go.Scatter(x=t_d,
                             y=I_minus,
                             fill=None,
                             mode='lines',
                             line_color='rgba(200,0,20,0.3)',
                             showlegend=False,
                             name='ICMinus',
                             ))
    fig.add_trace(go.Scatter(x=t_d,
                             y=I_plus,
                             fill='tonexty',
                             mode='lines',
                             fillcolor='rgba(200,0,20,0.2)',
                             line_color='rgba(255,255,255,0)',
                             name='95% IC',
                             ))
    fig.add_trace(go.Scatter(x=t_d,
                             y=I_plus,
                             fill=None,
                             mode='lines',
                             line_color='rgba(200,0,20,0.3)',
                             showlegend=False,
                             name='ICplus',
                             ))
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='black')
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='black')
    fig.update_xaxes(title_text='Day')
    fig.update_yaxes(title_text='Number of people infected')

    # Show Confidence Intervals
    st.subheader('Model Plot:')
    st.plotly_chart(fig)

    st.subheader('Confidence Intervals:')
    st.write(pd.DataFrame.from_dict(ci_dict))

    # Show subplots: Susceptible, Infected, Recovered
    layout2 = go.Layout(paper_bgcolor='rgba(0,0,0,0)',
                        plot_bgcolor='rgba(0,0,0,0)',
                        )
    fig2 = make_subplots(rows=1, cols=3,
                         subplot_titles=("Susceptible",
                                         "Infected",
                                         "Recovered or removed"),
                         shared_yaxes=True)
    fig2.add_trace(go.Scatter(x=tspan_long,
                              y=N_S_long,
                              ), row=1, col=1)
    fig2.add_trace(go.Scatter(x=tspan_long,
                              y=N_I_long,
                              ), row=1, col=2)
    fig2.add_trace(go.Scatter(x=tspan_long,
                              y=N_R_long,
                              ), row=1, col=3)
    fig2.update_xaxes(title_text="Days", row=1, col=1)
    fig2.update_xaxes(title_text="Days", row=1, col=2)
    fig2.update_xaxes(title_text="Days", row=1, col=3)
    fig2.update_yaxes(title_text="Number of people", row=1, col=1)
    fig2.update_layout(showlegend=False)
    st.subheader('Number of susceptible,infected and recovered' +
                 'in a three months period:')
    st.plotly_chart(fig2)

    return
