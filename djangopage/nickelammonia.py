from djangopage.reactions import *
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import least_squares
import scipy.interpolate, scipy.optimize
from django_plotly_dash import DjangoDash

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = DjangoDash('nickelammonia', add_bootstrap_links=True)
# add_bootstrap_links=True


app.layout = html.Div([
    html.Div([
        html.H1(u'Nickel-Ammonia-H\u2082O system')],
        style={
            'text-align':'center',
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '-10px 0 -10px 0',
            'margin-bottom': '2px',
            'width': '100%',
            'color': '#10328f',
        }
        ),

    html.Div([
        html.H6(u"Total nickel concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='nickel_slider',
            min=0.1,
            max=3.0,
            value=1,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
                ),
            ],
        style={
            'padding': '10px 15px 10px',
            'border': 'thin lightgrey solid',
            'margin-bottom': '3px',
            'backgroundColor': 'rgb(250, 250, 250)',
            'width': '100%',
            }
        ),

    html.Div([
        html.H6(u'Total ammonia concentration (kmolm\u207B\u00B3):'),
        dcc.Slider(
            id='ammonia_dropdown',
            min=0.0,
            max=3.0,
            value=1.5,
            step=0,
            marks={i: str(i) for i in [0, 0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
        ),
    ],
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250,250,250)',
            'padding': '10px 15px 10px',
            'margin-bottom': '3px',
            'width': '100%',
            }
    ),

    html.Div([
        html.Div([
            dcc.Graph(
                # animate=True,
                id='speciation plot',
                )
        ],
        style={
            'width': '48%',
            'margin-right': '1%',
        }
        ),
        html.Div([
            dcc.Graph(
                id='Potential-pH curve',
            ),
        ],
        style={
            'width': '48%',
            'margin-left': '1%',
        }
        )
    ],
    className='row',
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '0 0px 0 15px',
            "margin-left": "0px",
            'margin-right': '0px',
            'margin-bottom': '3px'
            }),


], style={
        'margin-top': '40px',
        # "margin-left": "10px",
        # 'margin-right': '10px',
        'margin-bottom': '5px'
})

@app.callback(
    Output('speciation plot', 'figure'),
    [Input('nickel_slider', 'value'),
    Input('ammonia_dropdown', 'value')])

def potential_graph(ni_total, ammonia_total):
    k1 = 10 ** (-9.25)
    k2 = 10 ** (-18.26)
    k3 = 10 ** (-35.91)
    logk = 11.853
    pH_x = list(np.linspace(0, 14, 200))
    T_ = 298
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.
    def concs(ammonia_total, ni_total, pH_x):
        h = 10 ** (-pH_x)
        if ammonia_total != 0.0:
            def equations(p):
                nh4, nio2 = p
                f = 10 ** (logk - 2 * pH_x)
                ni2pfree = f / (1 + f / ni_total)
                nh3 = k1 * nh4 / h
                nin4p2 = k2 * (nh4 ** 4) / (h ** 2)
                nin6p2 = k3 * (nh4 ** 6) / (h ** 4)
                return (ammonia_total - nh3 - 4 * nin4p2 - 6 * nin6p2 - nh4,
                        ni_total - ni2pfree - nin4p2 - nin6p2 - nio2)
            res = least_squares(equations, (0.1, 0.1), bounds=((0, 0), (ammonia_total, ni_total)),method='dogbox', xtol=1e-12)
            nh4 = res.x[0]
            nio2 = res.x[1]
            f = 10 ** (logk - 2 * pH_x)
            ni2pfree = f / (1 + f / ni_total)
            nh3 = k1 * nh4 / h
            nin4p2 = k2 * (nh4 ** 4) / (h ** 2)
            nin6p2 = k3 * (nh4 ** 6) / (h ** 4)
        elif ammonia_total == 0.0:
            f = 10 ** (logk - 2 * pH_x)
            ni2pfree = f / (1 + f / ni_total)
            nio2 = ni_total - ni2pfree
            nin4p2 = nin6p2 = nh3 = nh4 = 0
        return [nh3, nh4, nin4p2, nin6p2, ni2pfree, nio2]

    nh3plot = []
    nh4plot = []
    nin4p2plot = []
    nin6p2plot = []
    ni2pplot = []
    nio2plot = []

    for pHval in pH_x:
        result = concs(ammonia_total, ni_total, pHval)
        nh3plot.append(result[0])
        nh4plot.append(result[1])
        nin4p2plot.append(result[2])
        nin6p2plot.append(result[3])
        ni2pplot.append(result[4])
        nio2plot.append(result[5])

    if ammonia_total != 0.0:
        datasets = [nh3plot,nh4plot,nin4p2plot,nin6p2plot,ni2pplot,nio2plot]
        name = ['NH<sub>3</sub>', 'NH<sub>4</sub><sup>+</sup>', '[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>', '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>', 'Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        fill = [None, None, None, None, None, None]
        color = ['rgb(90, 0, 100)', 'rgb(40, 130, 80)', 'rgb(245, 137, 22)', 'rgb(63, 63, 191)', 'rgb(191, 63, 63)', 'rgb(15, 15, 15)']
    elif ammonia_total == 0.0:
        datasets = [ni2pplot, nio2plot]
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        fill = [None for i in range(len(name))]
        color = ['rgb(191, 63, 63)', 'rgb(243, 238, 77)']

    data1 = []
    for i, dataset in enumerate(datasets):
        data1.append(go.Scatter(
            x=pH_x,
            y=dataset,
            mode='lines',
            hoverinfo='skip',
            fill=fill[i],
            name=name[i],
            showlegend=True,
            line=dict(
                shape='spline',
                width=2.5,
                color=color[i]
            )
        ))

    layout = go.Layout(
        xaxis={'title': 'pH', 'linecolor': 'grey', 'mirror':True},
        yaxis={'title': 'Concentration (kmolm<sup>-3</sup>)', 'linecolor': 'grey', 'mirror':True},
        # transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'l':10},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgb(240,240,240)',
        autosize=False,
        width=600,
        height=500,
    )

    fig1 = go.Figure(data=data1, layout=layout)
    fig1.update_xaxes(gridcolor='LightPink', range=[0, 14],
                     nticks=20, mirror=True, ticks='outside', showline=True)

    fig1.update_yaxes(gridcolor='LightPink', ticks='outside',
                    range=[0, max(ni_total, ammonia_total)*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1

########################################################################################################################
@app.callback(
    Output('Potential-pH curve', 'figure'),
    [Input('nickel_slider', 'value'),
     Input('ammonia_dropdown', 'value')])

def potential_graph(ni_total, ammonia_total):
    k1 = 10 ** (-9.25)
    k2 = 10 ** (-18.26)
    k3 = 10 ** (-35.91)
    logk = 11.853
    pH_x = list(np.linspace(0, 14, 200))
    T_ = 298
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.
    def concs(ammonia_total, ni_total, pH_x):
        h = 10 ** (-pH_x)
        if ammonia_total != 0.0:
            def equations(p):
                nh4, nio2 = p
                f = 10 ** (logk - 2 * pH_x)
                ni2pfree = f / (1 + f / ni_total)
                nh3 = k1 * nh4 / h
                nin4p2 = k2 * (nh4 ** 4) / (h ** 2)
                nin6p2 = k3 * (nh4 ** 6) / (h ** 4)
                return (ammonia_total - nh3 - 4 * nin4p2 - 6 * nin6p2 - nh4,
                        ni_total - ni2pfree - nin4p2 - nin6p2 - nio2)
            res = least_squares(equations, (0.1, 0.1), bounds=((0, 0), (ammonia_total, ni_total)),method='dogbox', xtol=1e-12)
            nh4 = res.x[0]
            nio2 = res.x[1]
            f = 10 ** (logk - 2 * pH_x)
            ni2pfree = f / (1 + f / ni_total)
            nh3 = k1 * nh4 / h
            nin4p2 = k2 * (nh4 ** 4) / (h ** 2)
            nin6p2 = k3 * (nh4 ** 6) / (h ** 4)
        elif ammonia_total == 0.0:
            f = 10 ** (logk - 2 * pH_x)
            ni2pfree = f / (1 + f / ni_total)
            nio2 = ni_total - ni2pfree
            nin4p2 = nin6p2 = nh3 = nh4 = 0
        return [nh3, nh4, nin4p2, nin6p2, ni2pfree, nio2]

    nh3plot = []
    nh4plot = []
    nin4p2plot = []
    nin6p2plot = []
    ni2pplot = []
    nio2plot = []

    for pHval in pH_x:
        concs_result = concs(ammonia_total, ni_total, pHval)
        nh3plot.append(concs_result[0])
        nh4plot.append(concs_result[1])
        nin4p2plot.append(concs_result[2])
        nin6p2plot.append(concs_result[3])
        ni2pplot.append(concs_result[4])
        nio2plot.append(concs_result[5])
    #####################################################################################################
    ni2pfree = ni_total
    nh3 = nh4 = ammonia_total
    nio2 = max(nio2plot)
    nin4p2 = max(nin4p2plot)
    nin6p2 = max(nin6p2plot)

    def statuschecker(nin4p2plot, nin6p2plot, nio2plot, ammonia_total):
        if ammonia_total != 0.0:
            maxnin4p2 = max(nin4p2plot)
            i = nin4p2plot.index(maxnin4p2)
            maxnin6p2 = max(nin6p2plot)
            j = nin6p2plot.index(maxnin6p2)
            if maxnin6p2 > nio2plot[j] and maxnin4p2 > nio2plot[i]:
                status = 0
            elif maxnin6p2 > nio2plot[j] and maxnin4p2 < nio2plot[i]:
                status = 1
            elif maxnin6p2 < nio2plot[j] and maxnin4p2 > nio2plot[i]:
                status = 2  # (This actually never happens)
            elif maxnin6p2 < nio2plot[j] and maxnin4p2 < nio2plot[i]:
                status = 3
        elif ammonia_total == 0.0:
            status = 3
        return status
    status = statuschecker(nin4p2plot, nin6p2plot, nio2plot, ammonia_total)

    def trace_generator(pH_x,ni2pfree,nio2,nin4p2,nin6p2,nh3,nh4,T_):
        interp1 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, ni2pplot)
        interp2 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, nio2plot)
        interp3 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, nin4p2plot)
        interp4 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, nin6p2plot)
        interp5 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, nh3plot)
        interp6 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, nh4plot)

        def difference1(pH_x):
            return np.abs(interp2(pH_x) - interp1(pH_x))

        def difference2(pH_x):
            return np.abs(interp3(pH_x) - interp2(pH_x))

        def difference3(pH_x):
            return np.abs(interp4(pH_x) - interp3(pH_x))

        def difference4(pH_x):
            return np.abs(interp6(pH_x) - interp5(pH_x))

        def difference5(pH_x):
            return np.abs(interp2(pH_x) - interp4(pH_x))

        def difference6(pH_x):
            return np.abs(interp2(pH_x) - interp3(pH_x))

        aaa = []
        if status == 0:
            aaa.append(scipy.optimize.fsolve(difference1, x0=6.0))
            aaa.append(scipy.optimize.fsolve(difference2, x0=8.0))
            aaa.append(scipy.optimize.fsolve(difference3, x0=9.0))
            aaa.append(scipy.optimize.fsolve(difference4, x0=10))
            aaa.append(scipy.optimize.fsolve(difference5, x0=11))

        if status == 1:
            aaa.append(scipy.optimize.fsolve(difference1, x0=6.0))
            aaa.append(scipy.optimize.fsolve(difference3, x0=9.0))
            aaa.append(scipy.optimize.fsolve(difference4, x0=10))
            aaa.append(scipy.optimize.fsolve(difference5, x0=11))

        if status == 2:
            aaa.append(scipy.optimize.fsolve(difference1, x0=6.0))
            aaa.append(scipy.optimize.fsolve(difference2, x0=8.0))
            aaa.append(scipy.optimize.fsolve(difference6, x0=10))

        if status == 3 or ammonia_total == 0.0:
            aaa.append(scipy.optimize.fsolve(difference1, x0=6.0))

        x_intercept = []
        for xvalues in aaa:
            for xvalue in xvalues:
                x_intercept.append(xvalue)

        x_data = []
        if status != 3 or ammonia_total == 0.0:
            for i, item in enumerate(x_intercept):
                if i == 0:
                    x_data.append(np.linspace(0, item, 5))
                elif i >= 1:
                    x_data.append(np.linspace(x_intercept[i - 1], item, 5))
            finalindex = len(x_intercept) - 1
            x_data.append(np.linspace(x_intercept[finalindex], 14, 5))
        elif status == 3 or ammonia_total == 0.0:
            x_data.append(list(np.linspace(0, x_intercept[0], 5)))
            x_data.append(list(np.linspace(x_intercept[0], 14, 5)))

        new_x_data = []
        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_data.append(xvalue)

        if status == 0:
            y_data_bottom = [[R1(ni2pfree, T_) for i in range(len(x_data[0]))],np.linspace(R1(ni2pfree, T_), R2(x_intercept[1], T_), 5),np.linspace(R2(x_intercept[1], T_), R3(x_intercept[2], nh4, nin4p2, T_), 5),
                             np.linspace(R3(x_intercept[2], nh4, nin4p2, T_), R4(x_intercept[3], nh4, nin6p2, T_), 5),np.linspace(R4(x_intercept[3], nh4, nin6p2, T_), R5(nh3, nin6p2, T_), 5),np.linspace(R5(nh3, nin6p2, T_), R2(14, T_), 5)]
        elif status == 1:
            y_data_bottom = [[R1(ni2pfree, T_) for i in range(len(x_data[0]))],np.linspace(R1(ni2pfree, T_), R2(x_intercept[1], T_), 5), np.linspace(R2(x_intercept[1], T_), R4(x_intercept[2], nh4, nin6p2, T_), 5),
                             np.linspace(R4(x_intercept[2], nh4, nin6p2, T_), R5(nh3, nin6p2, T_), 5),np.linspace(R5(nh3, nin6p2, T_), R2(14, T_), 5)]
        elif status == 2:
            y_data_bottom = [[R1(ni2pfree, T_) for i in range(len(x_data[0]))],np.linspace(R1(ni2pfree, T_), R2(x_intercept[1], T_), 5),
                             np.linspace(R2(x_intercept[1], T_), R3(x_intercept[1], nh4, nin4p2, T_), 5),np.linspace(R3(x_intercept[1], nh4, nin4p2, T_), R2(14, T_), 5)]
        elif status == 3 or ammonia_total == 0.0:
            y_data_bottom = [[R1(ni2pfree, T_) for i in range(len(x_data[0]))],np.linspace(R1(ni2pfree, T_), R2(14, T_), 5)]

        new_y_bottom = []
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        if status == 0:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),np.linspace(T1(ni2pfree, x_intercept[0], T_), T2(x_intercept[1], T_), 5),np.linspace(T2(x_intercept[1], T_), T3(nh4, nin4p2, x_intercept[2], T_), 5),
                          np.linspace(T3(nh4, nin4p2, x_intercept[2], T_), T4(x_intercept[3], nh4, nin6p2, T_), 5),np.linspace(T4(x_intercept[3], nh4, nin6p2, T_), T5(nh3, x_intercept[4], nin6p2, T_), 5),np.linspace(T5(nh3, x_intercept[4], nin6p2, T_), T2(14, T_), 5)]

        elif status == 1:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),np.linspace(T1(ni2pfree, x_intercept[0], T_), T2(x_intercept[1], T_), 5),np.linspace(T2(x_intercept[1], T_), T4(x_intercept[2], nh4, nin6p2, T_), 5),
                          np.linspace(T4(x_intercept[2], nh4, nin6p2, T_), T5(nh3, x_intercept[3], nin6p2, T_), 5),np.linspace(T5(nh3, x_intercept[3], nin6p2, T_), T2(14, T_), 5)]
        elif status == 2:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),np.linspace(T1(ni2pfree, x_intercept[0], T_), T2(x_intercept[1], T_), 5),
                          np.linspace(T2(x_intercept[1], T_), T3(nh4, nin4p2, x_intercept[2]), 5),np.linspace(T3(nh4, nin4p2, x_intercept[2]), T2(14, T_), 5)]
        elif status == 3 or ammonia_total == 0.0:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),np.linspace(T1(ni2pfree, x_intercept[0], T_), T2(14, T_), 5)]

        new_y_top = []
        for yvalues in y_data_top:
            for yvalue in yvalues:
                new_y_top.append(yvalue)

        nio3regionx = list(new_x_data) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        nio3regiony = list(new_y_top) + list(np.linspace(T2(14, T_), 2.4, 5)) + list([2.4 for i in range(0, 5)]) + list(np.linspace(2.4, T1(ni2pfree, 0, T_), 5))
        niregionx = list(new_x_data) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        niregiony = list(new_y_bottom) + list(np.linspace(R2(14, T_), -1, 5)) + list([-1 for i in range(0, 5)]) + list(np.linspace(-1, R1(ni2pfree, T_), 5))
        if status == 0:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, T_), R1(ni2pfree, T_), 5))

            nio2regionx1 = list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(x_data[1]) + list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(reversed(x_data[1]))
            nio2regiony1 = list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(y_data_bottom[1]) + list(np.linspace(R2(x_intercept[1], T_), T2(x_intercept[1], T_), 5)) + list(reversed(y_data_top[1]))

            nin4p2regionx1 = list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(x_data[2]) + list(np.linspace(x_intercept[2], x_intercept[2], 5)) + list(reversed(x_data[2]))
            nin4p2regiony1 = list(np.linspace(R2(x_intercept[1], T_), T2(x_intercept[1], T_), 5)) + list(y_data_bottom[2]) + list(np.linspace(R3(x_intercept[2], nh4, nin4p2, T_), T3(nh4, nin4p2, x_intercept[2], T_), 5)) + list(reversed(y_data_top[2]))

            nin6p2regionx1 = list(np.linspace(x_intercept[2], x_intercept[2], 5)) + list(x_data[3]) + list(x_data[4]) + list(np.linspace(x_intercept[4], x_intercept[4], 5)) + list(reversed(x_data[4])) + list(reversed(x_data[3]))
            nin6p2regiony1 = list(np.linspace(T3(nh4, nin4p2, x_intercept[2], T_), R3(x_intercept[2], nh4, nin4p2, T_), 5)) + list(y_data_bottom[3]) + list(y_data_bottom[4]) + list(np.linspace(R5(nh3, nin6p2, T_), T5(nh3, x_intercept[4], nin6p2, T_), 5)) + list(reversed(y_data_top[4])) + list(reversed(y_data_top[3]))

            nio2regionx2 = list(np.linspace(x_intercept[4], x_intercept[4], 5)) + list(x_data[5]) + list([14 for i in range(0, 5)]) + list(reversed(x_data[5]))
            nio2regiony2 = list(np.linspace(T5(nh3, x_intercept[4], nin6p2, T_), R5(nh3, nin6p2, T_), 5)) + list(y_data_bottom[5]) + list(np.linspace(R2(14, T_), T2(14, T_), 5)) + list(reversed(y_data_top[5]))

            xs = [nip2regionx, nio2regionx1, nin4p2regionx1, nin6p2regionx1, nio2regionx2]
            ys = [nip2regiony, nio2regiony1, nin4p2regiony1, nin6p2regiony1, nio2regiony2]

        elif status == 1:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list( np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, T_), R1(ni2pfree, T_), 5))

            nio2regionx1 = list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(x_data[1]) + list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(reversed(x_data[1]))
            nio2regiony1 = list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(y_data_bottom[1]) + list(np.linspace(R2(x_intercept[1], T_), T2(x_intercept[1], T_), 5)) + list(reversed(y_data_top[1]))

            nin6p2regionx1 = list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(x_data[2]) + list(x_data[3]) + list(np.linspace(x_intercept[3], x_intercept[3], 5)) + list(reversed(x_data[3])) + list(reversed(x_data[2]))
            nin6p2regiony1 = list(np.linspace(T2(x_intercept[1], T_), R2(x_intercept[1], T_), 5)) + list(y_data_bottom[2]) + list(y_data_bottom[3]) + list(np.linspace(R5(nh3, nin6p2, T_), T5(nh3, x_intercept[3], nin6p2, T_), 5)) + list(reversed(y_data_top[3])) + list(reversed(y_data_top[2]))

            nio2regionx2 = list(np.linspace(x_intercept[3], x_intercept[3], 5)) + list(x_data[4]) + list([14 for i in range(0, 5)]) + list(reversed(x_data[4]))
            nio2regiony2 = list(np.linspace(T5(nh3, x_intercept[3], nin6p2, T_), R5(nh3, nin6p2, T_), 5)) + list(y_data_bottom[4]) + list(np.linspace(R2(14, T_), T2(14, T_), 5)) + list(reversed(y_data_top[4]))

            xs = [nip2regionx, nio2regionx1, nin6p2regionx1, nio2regionx2]
            ys = [nip2regiony, nio2regiony1, nin6p2regiony1, nio2regiony2]

        elif status == 2:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, T_), R1(ni2pfree, T_), 5))

            nio2regionx1 = list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(x_data[1]) + list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(reversed(x_data[1]))
            nio2regiony1 = list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(y_data_bottom[1]) + list(np.linspace(R2(x_intercept[1], T_), T2(x_intercept[1], T_), 5)) + list(reversed(y_data_top[1]))

            nin4p2regionx1 = list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(x_data[2]) + list(np.linspace(x_intercept[2], x_intercept[2], 5)) + list(reversed(x_data[2]))
            nin4p2regiony1 = list(np.linspace(R2(x_intercept[1], T_), T2(x_intercept[1], T_), 5)) + list(y_data_bottom[2]) + list(np.linspace(R3(x_intercept[2], nh4, nin4p2, T_), T3(nh4, nin4p2, x_intercept[2], T_), 5)) + list(reversed(y_data_top[2]))

            nio2regionx2 = list(np.linspace(x_intercept[2], x_intercept[2], 5)) + list(x_data[3]) + list([14 for i in range(0, 5)]) + list(reversed(x_data[3]))
            nio2regiony2 = list(np.linspace(R3(x_intercept[2], nh4, nin4p2, T_), T3(nh4, nin4p2, x_intercept[2], T_), 5)) + list(y_data_bottom[3]) + list(np.linspace(R2(14, T_), T2(14, T_), 5)) + list(reversed(y_data_top[3]))

            xs = [nip2regionx, nio2regionx1, nin4p2regionx1, nio2regionx2]
            ys = [nip2regiony, nio2regiony1, nin4p2regiony1, nio2regiony2]

        elif status == 3 or ammonia_total == 0.0:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, T_), R1(ni2pfree, T_), 5))

            nio2regionx = list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(x_data[1]) + list([14 for i in range(0, 5)]) + list(reversed(x_data[1]))
            nio2regiony = list(np.linspace(R2(x_intercept[0], T_), T2(x_intercept[0], T_), 5)) + list(y_data_bottom[1]) + list(np.linspace(R2(14, T_), T2(14, T_), 5)) + list(reversed(y_data_top[1]))
            xs = [nip2regionx, nio2regionx]
            ys = [nip2regiony, nio2regiony]
        return [xs, ys, niregionx, niregiony, nio3regionx, nio3regiony]

    xs = trace_generator(pH_x,ni2pfree,nio2,nin4p2, nin6p2, nh3, nh4, T_)[0]
    ys = trace_generator(pH_x, ni2pfree, nio2, nin4p2, nin6p2, nh3, nh4, T_)[1]
    niregionx = trace_generator(pH_x, ni2pfree, nio2, nin4p2, nin6p2, nh3, nh4, T_)[2]
    niregiony = trace_generator(pH_x, ni2pfree, nio2, nin4p2, nin6p2, nh3, nh4, T_)[3]
    nio3regionx = trace_generator(pH_x, ni2pfree, nio2, nin4p2, nin6p2, nh3, nh4, T_)[4]
    nio3regiony = trace_generator(pH_x, ni2pfree, nio2, nin4p2, nin6p2, nh3, nh4, T_)[5]

    if status == 0:
        name = ['Ni<sup>2+</sup>','Ni(OH)<sub>2</sub>','[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>', '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(30, 154, 216, 0.43)', 'rgba(243, 238, 77, 0.5)', 'rgba(114, 102, 234, 0.63)','rgba(245, 169, 104, 0.81)', 'rgba(243, 238, 77, 0.5)']
    elif status == 1:
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>', '[Ni(NH<sub>3</sub>)<sub>6</sub>]<sup>2+</sup>','Ni(OH)<sub>2</sub>']
        color = ['rgba(101, 18, 183, 0.41)', 'rgba(243, 238, 77, 0.5)', 'rgba(191, 38, 42, 0.43)','rgba(114, 204, 234, 0.63)']
    elif status == 2: # this never really happens
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>', 'Ni(OH)<sub>2</sub>', '[Ni(NH<sub>3</sub>)<sub>4</sub>]<sup>2+</sup>']
    elif status == 3 or ammonia_total==0 :
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(127, 14, 191, 0.3)', 'rgba(30, 205, 40, 0.5)']

    data = []
    for i, xvals in enumerate(xs):
        data.append(go.Scatter(
            x=xvals,
            y=ys[i],
            mode='none',
            fill='toself',
            hoverinfo='skip',
            fillcolor=color[i],
            showlegend=True,
            name=name[i]
        ))

    # add water splitting
    ywater = [W1(pH_x, T_), W2(pH_x, T_)]
    for ys in ywater:
        data.append(go.Scatter(
            x = pH_x,
            y = ys,
            mode='lines',
            hoverinfo='skip',
            showlegend=False,
            line=dict(
                shape='spline',
                color='blue',
                width=2,
                dash='dash'
                )
            ))

    layout = go.Layout(
        xaxis={'title': 'pH', 'linecolor':'grey', 'mirror':True},
        yaxis={'title': "Electrode potential, <i>E</i> vs SHE/ V ", 'linecolor':'grey', 'mirror':True},
        transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'r': 20},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(50,0,0,0)',
        autosize=False,
        width=600,
        height=500,
            )
    fig = go.Figure(data=data, layout=layout)
    extrax = [niregionx, nio3regionx]
    extray = [niregiony, nio3regiony]
    namee = ['Ni', 'Ni(OH)<sub>3</sub>']
    colors1 = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)']
    tracesextra = []
    for i, extraxx in enumerate(extrax):
        tracesextra.append(go.Scatter(
            x=extraxx,
            y=extray[i],
            mode='lines',
            name=namee[i],
            #fillcolor=colors1[i],
            fill='toself',
            showlegend=True,
            hoverinfo='skip'
        ))
    for trace in tracesextra:
        fig.add_traces(trace)

    if ammonia_total != 0:
        fig.add_trace(go.Scatter(
            x=list([9.25 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                    shape='spline',
                    color='purple',
                    width=1.5,
                    dash='dash'
                    )
        ))
        fig.add_annotation(
            x=8.2,
            y=2,
            text="NH<sub>4</sub><sup>+</sup>",
            showarrow=False,
            font=dict(
                    family="Courier New bold",
                    size=20,
                    color="purple"
                    ),)
        fig.add_annotation(
            x=10.3,
            y=1.97,
            text="NH<sub>3</sub>",
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=20,
                color="purple"
            )
        )
    fig.add_annotation(
        x=0.4,
        y=1.37,
        showarrow=False,
        text='O<sub>2</sub>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue',
        )
    )
    fig.add_annotation(
        x=0.55,
        y=1.02,
        showarrow=False,
        text='H<sub>2</sub>O',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.add_annotation(
        x=0.37,
        y=0.1,
        showarrow=False,
        text='H<sup>+</sup>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.add_annotation(
        x=0.4,
        y=-0.2,
        showarrow=False,
        text='H<sub>2</sub>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.update_xaxes(gridcolor='LightPink',
                     nticks=20,  ticks='outside', showline=True,
                    zeroline=True, zerolinewidth=2, zerolinecolor='LightPink')
    fig.update_yaxes(nticks=20, gridcolor='LightPink',  ticks='outside', showline=True,
                    zeroline=True, zerolinewidth=2, zerolinecolor='LightPink')
    fig.update_layout(yaxis=dict(range=[-1,2.4]), xaxis=dict(range=[0,14]),
                        title={
                            'text': "Potential-pH diagram",
                            'y':0.95,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top'
                            })
    return fig