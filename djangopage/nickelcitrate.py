from djangopage.reactions1 import *
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from django_plotly_dash import DjangoDash
from scipy.optimize import least_squares
import scipy.interpolate, scipy.optimize

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = DjangoDash('nickelcitrate', add_bootstrap_links=True)
# , add_bootstrap_links=True
# add_bootstrap_links=True


app.layout = html.Div([
    html.Div([
        html.H1(u'Nickel-Citrate-H\u2082O system')],
        style={
            'text-align':'center',
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '-10px 0 -10px 0',
            'margin-bottom': '2px',
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
            marks={n_activity: str(n_activity) for n_activity in [0.1,0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3]},
                ),
            ],
        style={
            'padding': '10px 15px 10px',
            'border': 'thin lightgrey solid',
            'margin-bottom': '3px',
            'backgroundColor': 'rgb(250, 250, 250)',
            }
        ),

    html.Div([
        html.H6(u'Total citrate concentration (kmolm\u207B\u00B3):'),
        dcc.Slider(
            id='citrate_dropdown',
            min=0,
            max=3.0,
            value=1.5,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0, 0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3]},
                 ),
    ],
        style={
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 15px 10px',
            'margin-bottom': '3px',
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

####################################################################################################################
@app.callback(
    Output('speciation plot', 'figure'),
    [Input('nickel_slider', 'value'),
     Input('citrate_dropdown', 'value')])

def speciation_graph(ni_total, citrate_total):
    k1 = 7.129 * (10 ** 11)
    k2 = 9.905 * (10 ** 2)
    k3 = 5.296 * (10 ** 4)
    k4 = 2.481 * (10 ** 6)
    k5 = 2.512 * (10 ** 5)
    k6 = 1.995 * (10 ** 3)
    k7 = 5.623 * (10 ** 1)
    pH_x = np.linspace(0, 14, 200)
    T_ = 298
    #----------------------------------------------------------------------------------------------
    # begin first function, output all species concentrations. One concentration for each pH value.
    def concs(citrate_total, ni_total, pH_x):
        h = 10 ** (-pH_x)
        if citrate_total != 0:
            def equations(p):
                cit3, nio2 = p
                Hcit = h * cit3 * k4
                H2cit = h * Hcit * k3
                H3cit = h * H2cit * k2
                ni2pfree = (k1 * (h ** 2)) / (1 + ((k1 * (h ** 2)) / ni_total))
                NiH2cit = k7 * ni2pfree * H2cit
                NiHcit = k6 * ni2pfree * Hcit
                Nicit = k5 * ni2pfree * cit3
                return (citrate_total - Hcit - H2cit - H3cit - Nicit - NiHcit - NiH2cit - cit3,
                        ni_total - Nicit - NiHcit - NiH2cit - ni2pfree - nio2)
            res = least_squares(equations, (0.1, 0.1), bounds=((0, 0), (citrate_total, ni_total)), method='dogbox',xtol=1e-12)
            cit3 = res.x[0]
            nio2 = res.x[1]
            ni2pfree = (k1 * (h ** 2)) / (1 + ((k1 * (h ** 2)) / ni_total))
            Hcit = h * cit3 * k4
            H2cit = h * Hcit * k3
            H3cit = h * H2cit * k2
            NiH2cit = k7 * ni2pfree * H2cit
            NiHcit = k6 * ni2pfree * Hcit
            Nicit = k5 * ni2pfree * cit3
        elif citrate_total == 0:
            ni2pfree = (k1 * (h ** 2)) / (1 + ((k1 * (h ** 2)) / ni_total))
            nio2 = ni_total - ni2pfree
            cit3 = Hcit =H2cit =H3cit =NiHcit =Nicit =NiH2cit=0
        return [cit3, nio2, ni2pfree, Hcit, H2cit, H3cit, NiH2cit, NiHcit, Nicit]

    cit3freeplot = []
    nio2freeplot = []
    ni2pfreeplot = []
    Hcitfreeplot = []
    H2citfreeplot = []
    H3citfreeplot = []
    NiH2citfreeplot = []
    NiHcitfreeplot = []
    Nicitfreeplot = []

    for pHval in pH_x:
        concs_results = concs(citrate_total, ni_total, pHval)
        cit3freeplot.append(concs_results[0])
        nio2freeplot.append(concs_results[1])
        ni2pfreeplot.append(concs_results[2])
        Hcitfreeplot.append(concs_results[3])
        H2citfreeplot.append(concs_results[4])
        H3citfreeplot.append(concs_results[5])
        NiH2citfreeplot.append(concs_results[6])
        NiHcitfreeplot.append(concs_results[7])
        Nicitfreeplot.append(concs_results[8])

    if citrate_total != 0.0:
        datasets = [cit3freeplot, nio2freeplot, ni2pfreeplot, Hcitfreeplot, H2citfreeplot,
                    H3citfreeplot, NiH2citfreeplot, NiHcitfreeplot, Nicitfreeplot]
        name = ['Cit<sup>3-</sup>', 'Ni(OH)<sub>2</sub>', 'Ni<sup>2+</sup>', 'Hcit<sup>-</sup>',
                'H<sub>2</sub>cit<sup>2-</sup>', 'H<sub>3</sub>cit', 'NiH<sub>2</sub>cit<sup>+</sup>',
                'NiHCit', 'Nicit<sup>-</sup>']
        fill = [None, None, None, None, None, None, None, None, None]
        color = ['rgb(90, 0, 100)', 'rgb(40, 130, 80)', 'rgb(9, 0, 0)', 'rgb(63, 63, 191)', 'rgb(191, 63, 63)',
                 'rgb(66, 81, 245)', 'rgb(218, 66, 245)', 'rgb(245, 144, 66)', 'rgb(245, 66, 90)']
    elif citrate_total == 0.0:
        datasets = [ni2pfreeplot, nio2freeplot]
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
                    range=[0, max(citrate_total, ni_total)*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1
#######################################################################################################################
@app.callback(
    Output('Potential-pH curve', 'figure'),
    [Input('nickel_slider', 'value'),
     Input('citrate_dropdown', 'value')])


def potential_graph(ni_total, citrate_total):
    k1 = 7.129 * (10 ** 11)
    k2 = 9.905 * (10 ** 2)
    k3 = 5.296 * (10 ** 4)
    k4 = 2.481 * (10 ** 6)
    k5 = 2.512 * (10 ** 5)
    k6 = 1.995 * (10 ** 3)
    k7 = 5.623 * (10 ** 1)
    pH_x = np.linspace(0, 14, 200)
    T_ = 298

    def concs(citrate_total, ni_total, pH_x):
        h = 10 ** (-pH_x)
        if citrate_total != 0:
            def equations(p):
                cit3, nio2 = p
                Hcit = h * cit3 * k4
                H2cit = h * Hcit * k3
                H3cit = h * H2cit * k2
                ni2pfree = (k1 * (h ** 2)) / (1 + ((k1 * (h ** 2)) / ni_total))
                NiH2cit = k7 * ni2pfree * H2cit
                NiHcit = k6 * ni2pfree * Hcit
                Nicit = k5 * ni2pfree * cit3
                return (citrate_total - Hcit - H2cit - H3cit - Nicit - NiHcit - NiH2cit - cit3,
                        ni_total - Nicit - NiHcit - NiH2cit - ni2pfree - nio2)
            res = least_squares(equations, (0.1, 0.1), bounds=((0, 0), (citrate_total, ni_total)), method='dogbox',xtol=1e-12)
            cit3 = res.x[0]
            nio2 = res.x[1]
            ni2pfree = (k1 * (h ** 2)) / (1 + ((k1 * (h ** 2)) / ni_total))
            Hcit = h * cit3 * k4
            H2cit = h * Hcit * k3
            H3cit = h * H2cit * k2
            NiH2cit = k7 * ni2pfree * H2cit
            NiHcit = k6 * ni2pfree * Hcit
            Nicit = k5 * ni2pfree * cit3
        elif citrate_total == 0:
            ni2pfree = (k1 * (h ** 2)) / (1 + ((k1 * (h ** 2)) / ni_total))
            nio2 = ni_total - ni2pfree
            cit3 = Hcit = H2cit = H3cit = NiHcit = Nicit = NiH2cit = 0
        return [cit3, nio2, ni2pfree, Hcit, H2cit, H3cit, NiH2cit, NiHcit, Nicit]

    cit3plot = []
    nio2plot = []
    ni2pplot = []
    Hcitplot = []
    H2citplot = []
    H3citplot = []
    NiH2citplot = []
    NiHcitplot = []
    Nicitplot = []

    for pHval in pH_x:
        concs_result = concs(citrate_total, ni_total, pHval)
        cit3plot.append(concs_result[0])
        nio2plot.append(concs_result[1])
        ni2pplot.append(concs_result[2])
        Hcitplot.append(concs_result[3])
        H2citplot.append(concs_result[4])
        H3citplot.append(concs_result[5])
        NiH2citplot.append(concs_result[6])
        NiHcitplot.append(concs_result[7])
        Nicitplot.append(concs_result[8])
        
    ni2pfree = nio2 = ni_total
    cit3 = Hcit = H2cit = H3cit = citrate_total
    Nicit = max(Nicitplot)
    NiHcit = max(NiHcitplot)
    NiH2cit = max(NiH2citplot)

    def statuschecker(NiH2citplot, NiHcitplot, Nicitplot, ni2pplot, nio2plot, citrate_total):
        if citrate_total != 0.0:
            maxNiH2cit = max(NiH2citplot)
            i = NiH2citplot.index(maxNiH2cit)
            maxNicit = max(Nicitplot)
            j = Nicitplot.index(maxNicit)
            maxNiHcit = max(NiHcitplot)
            g = NiHcitplot.index(maxNiHcit)

            if maxNiH2cit > ni2pplot[i] and maxNiHcit > ni2pplot[g]:
                status = 0
            elif maxNiH2cit < ni2pplot[i] and maxNicit > nio2plot[j]:
                status = 1
            elif maxNiH2cit < ni2pplot[i] and maxNicit < nio2plot[j]:
                status = 2
            elif maxNiH2cit > ni2pplot[i] and maxNiHcit < ni2pplot[g]:
                status = 3
        elif citrate_total == 0.0:
            status = 2
        return status
    status = statuschecker(NiH2citplot, NiHcitplot, Nicitplot, ni2pplot, nio2plot, citrate_total)

    def trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_):
        interp1 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, ni2pplot)
        interp2 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, NiH2citplot)
        interp3 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, H3citplot)
        interp4 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, H2citplot)
        interp5 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, Hcitplot)
        interp6 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, cit3plot)
        interp7 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, NiHcitplot)
        interp8 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, Nicitplot)
        interp9 = scipy.interpolate.InterpolatedUnivariateSpline(pH_x, nio2plot)

        def difference1(pH_x):
            return np.abs(interp2(pH_x) - interp1(pH_x))

        def difference2(pH_x):
            return np.abs(interp4(pH_x) - interp3(pH_x))

        def difference3(pH_x):
            return np.abs(interp7(pH_x) - interp2(pH_x))

        def difference4(pH_x):
            return np.abs(interp8(pH_x) - interp7(pH_x))

        def difference5(pH_x):
            return np.abs(interp5(pH_x) - interp4(pH_x))

        def difference6(pH_x):
            return np.abs(interp6(pH_x) - interp5(pH_x))

        def difference7(pH_x):
            return np.abs(interp9(pH_x) - interp8(pH_x))

        def difference8(pH_x):
            return np.abs(interp8(pH_x) - interp1(pH_x))

        def difference9(pH_x):
            return np.abs(interp9(pH_x) - interp1(pH_x))

        aaa = []
        if status == 0:
            aaa.append(scipy.optimize.fsolve(difference1, x0=1.5))
            aaa.append(scipy.optimize.fsolve(difference2, x0=3.0))
            aaa.append(scipy.optimize.fsolve(difference3, x0=3.0))
            aaa.append(scipy.optimize.fsolve(difference4, x0=4.0))
            aaa.append(scipy.optimize.fsolve(difference5, x0=5.0))
            aaa.append(scipy.optimize.fsolve(difference6, x0=6.0))
            aaa.append(scipy.optimize.fsolve(difference7, x0=9.0))

        elif status == 1:
            aaa.append(scipy.optimize.fsolve(difference8, x0=5.5))
            aaa.append(scipy.optimize.fsolve(difference6, x0=6.3))
            aaa.append(scipy.optimize.fsolve(difference7, x0=8.5))

        elif status == 2:
            aaa.append(scipy.optimize.fsolve(difference9, x0=5.5))

        elif status == 3:
            aaa.append(scipy.optimize.fsolve(difference1, x0=1.5))
            aaa.append(scipy.optimize.fsolve(difference1, x0=2.5))
            aaa.append(scipy.optimize.fsolve(difference8, x0=5.0))
            aaa.append(scipy.optimize.fsolve(difference6, x0=6.0))
            aaa.append(scipy.optimize.fsolve(difference7, x0=9.0))

        x_intercept = []
        for xvalues in aaa:
            for xvalue in xvalues:
                x_intercept.append(xvalue)

        x_data = []
        if status != 2 or citrate_total == 0.0:
            for i, item in enumerate(x_intercept):
                if i == 0:
                    x_data.append(np.linspace(0, item, 5))
                elif i >= 1:
                    x_data.append(np.linspace(x_intercept[i - 1], item, 5))
            finalindex = len(x_intercept) - 1
            x_data.append(np.linspace(x_intercept[finalindex], 14, 5))
        elif status == 2 or citrate_total == 0.0:
            x_data.append(list(np.linspace(0, x_intercept[0], 5)))
            x_data.append(list(np.linspace(x_intercept[0], 14, 5)))

        new_x_data = []
        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_data.append(xvalue)

        if status == 0:
            y_data_bottom = [np.linspace(R1(ni2pfree, T_), R1(ni2pfree, T_), 5),
                             np.linspace(R1(ni2pfree, T_), R2(x_intercept[1], T_, NiH2cit, H3cit), 5),
                             np.linspace(R2(x_intercept[1], T_, NiH2cit, H3cit), R3(x_intercept[2], T_, NiHcit, H3cit),
                                         5),
                             np.linspace(R3(x_intercept[2], T_, NiHcit, H3cit), R4(x_intercept[3], T_, NiHcit, H2cit),
                                         5),
                             np.linspace(R4(x_intercept[3], T_, NiHcit, H2cit), R5(x_intercept[4], T_, Nicit, H2cit),
                                         5),
                             np.linspace(R5(x_intercept[4], T_, Nicit, H2cit), R6(x_intercept[5], T_, Nicit, Hcit), 5),
                             np.linspace(R6(x_intercept[5], T_, Nicit, Hcit), R7(T_, Nicit, cit3), 5),
                             np.linspace(R7(T_, Nicit, cit3), R8(14, T_), 5)]

        elif status == 1:
            y_data_bottom = [np.linspace(R1(ni2pfree, T_), R1(ni2pfree, T_), 5),
                             np.linspace(R1(ni2pfree, T_), R6(x_intercept[1], T_, Nicit, Hcit), 5),
                             np.linspace(R6(x_intercept[1], T_, Nicit, Hcit), R7(T_, Nicit, cit3), 5),
                             np.linspace(R7(T_, Nicit, cit3), R8(14, T_), 5)]

        elif status == 2 or citrate_total == 0.0:
            y_data_bottom = [[R1(ni2pfree, T_) for i in range(len(x_data[0]))],
                             np.linspace(R1(ni2pfree, T_), R8(14, T_), 5)]

        elif status == 3:
            y_data_bottom = [np.linspace(R1(ni2pfree, T_), R1(ni2pfree, T_), 5),
                             np.linspace(R1(ni2pfree, T_), R2(x_intercept[1], T_, NiH2cit, H3cit), 5),
                             np.linspace(R2(x_intercept[1], T_, NiH2cit, H3cit), R1(ni2pfree, T_), 5),
                             np.linspace(R1(ni2pfree, T_), R6(x_intercept[3], T_, Nicit, Hcit), 5),
                             np.linspace(R6(x_intercept[3], T_, Nicit, Hcit), R7(T_, Nicit, cit3), 5),
                             np.linspace(R7(T_, Nicit, cit3), R8(14, T_), 5)]

        new_y_bottom = []
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        if status == 0:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),
                          np.linspace(T1(ni2pfree, x_intercept[0], T_), T2(x_intercept[1], NiH2cit, H3cit, T_), 5),
                          np.linspace(T2(x_intercept[1], NiH2cit, H3cit, T_), T3(x_intercept[2], NiHcit, H3cit, T_), 5),
                          np.linspace(T3(x_intercept[2], NiHcit, H3cit, T_), T4(x_intercept[3], NiHcit, H2cit, T_), 5),
                          np.linspace(T4(x_intercept[3], NiHcit, H2cit, T_), T5(x_intercept[4], Nicit, H2cit, T_), 5),
                          np.linspace(T5(x_intercept[4], Nicit, H2cit, T_), T6(x_intercept[5], Nicit, Hcit, T_), 5),
                          np.linspace(T6(x_intercept[5], Nicit, Hcit, T_), T7(x_intercept[6], Nicit, cit3, T_), 5),
                          np.linspace(T7(x_intercept[6], Nicit, cit3, T_), T8(14, T_), 5)]
        elif status == 1:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),
                          np.linspace(T1(ni2pfree, x_intercept[0], T_), T6(x_intercept[1], Nicit, Hcit, T_), 5),
                          np.linspace(T6(x_intercept[1], Nicit, Hcit, T_), T7(x_intercept[2], Nicit, cit3, T_), 5),
                          np.linspace(T7(x_intercept[2], Nicit, cit3, T_), T8(14, T_), 5)]

        elif status == 2 or citrate_total == 0.0:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),
                          np.linspace(T1(ni2pfree, x_intercept[0], T_), T8(14, T_), 5)]

        elif status == 3:
            y_data_top = [np.linspace(T1(ni2pfree, 0, T_), T1(ni2pfree, x_intercept[0], T_), 5),
                          np.linspace(T1(ni2pfree, x_intercept[0], T_), T2(x_intercept[1], NiH2cit, H3cit, T_), 5),
                          np.linspace(T2(x_intercept[1], NiH2cit, H3cit, T_), T1(ni2pfree, x_intercept[2], T_), 5),
                          np.linspace(T1(ni2pfree, x_intercept[2], T_), T6(x_intercept[3], Nicit, Hcit, T_), 5),
                          np.linspace(T6(x_intercept[3], Nicit, Hcit, T_), T7(x_intercept[4], Nicit, cit3, T_), 5),
                          np.linspace(T7(x_intercept[4], Nicit, cit3, T_), T8(14, T_), 5)]

        new_y_top = []
        for yvalues in y_data_top:
            for yvalue in yvalues:
                new_y_top.append(yvalue)

        nio3regionx = list(new_x_data) + list([14 for i in range(0, 5)]) + list(
            reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        nio3regiony = list(new_y_top) + list(np.linspace(T8(14, T_), 2.6, 5)) + list(
            [2.6 for i in range(0, 5)]) + list(np.linspace(2.6, T1(ni2pfree, 0, 298), 5))

        niregionx = list(new_x_data) + list([14 for i in range(0, 5)]) + list(
            reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        niregiony = list(new_y_bottom) + list(np.linspace(R8(14, T_), -1.8, 5)) + list(
            [-1.8 for i in range(0, 5)]) + list(np.linspace(-1.8, R1(ni2pfree, 298), 5))

        if status == 0.0:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(
                reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[0], T_), 5)) + list(
                reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, 298), R1(ni2pfree, T_), 5))

            NiH2citregionx = list(reversed(x_data[2])) + list(reversed(x_data[1])) + list(
                np.linspace(x_intercept[0], x_intercept[0], 5)) + list(x_data[1]) + list(x_data[2]) + list(
                np.linspace(x_intercept[2], x_intercept[2], 5))
            NiH2citregiony = list(reversed(y_data_bottom[2])) + list(reversed(y_data_bottom[1])) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[0], T_), 5)) + list(y_data_top[1]) + list(
                y_data_top[2]) + list(
                np.linspace(T3(x_intercept[2], NiHcit, H3cit, T_), R3(x_intercept[2], T_, NiHcit, H3cit), 5))

            NiHcitregionx = list(x_data[3]) + list(np.linspace(x_intercept[3], x_intercept[3], 5)) + list(
                reversed(x_data[3])) + list(np.linspace(x_intercept[2], x_intercept[2], 5))
            NiHcitregiony = list(y_data_bottom[3]) + list(
                np.linspace(R4(x_intercept[3], T_, NiHcit, H2cit), T4(x_intercept[3], NiHcit, H2cit, T_), 5)) + list(
                reversed(y_data_top[3])) + list(
                np.linspace(R3(x_intercept[2], T_, NiHcit, H3cit), T3(x_intercept[2], NiHcit, H3cit, T_), 5))

            Nicitregionx = list(x_data[4]) + list(x_data[5]) + list(x_data[6]) + list(
                np.linspace(x_intercept[6], x_intercept[6], 5)) + list(reversed(x_data[6])) + list(
                reversed(x_data[5])) + list(reversed(x_data[4])) + list(np.linspace(x_intercept[3], x_intercept[3], 5))
            Nicitregiony = list(y_data_bottom[4]) + list(y_data_bottom[5]) + list(y_data_bottom[6]) + list(
                np.linspace(R7(T_, Nicit, cit3), T7(x_intercept[6], Nicit, cit3, T_), 5)) + list(
                reversed(y_data_top[6])) + list(reversed(y_data_top[5])) + list(reversed(y_data_top[4])) + list(
                np.linspace(T4(x_intercept[3], NiHcit, H2cit, T_), R4(x_intercept[3], T_, NiHcit, H2cit), 5))

            nio2regionx = list(reversed(x_data[7])) + list(np.linspace(x_intercept[6], x_intercept[6], 5)) + list(
                x_data[7]) + list([14 for i in range(0, 5)])
            nio2regiony = list(reversed(y_data_bottom[7])) + list(
                np.linspace(R7(T_, Nicit, cit3), T7(x_intercept[6], Nicit, cit3, T_), 5)) + list(y_data_top[7]) + list(
                np.linspace(T8(14, T_), R8(14, T_), 5))

            xs = [nip2regionx, NiH2citregionx, NiHcitregionx, Nicitregionx, nio2regionx]
            ys = [nip2regiony, NiH2citregiony, NiHcitregiony, Nicitregiony, nio2regiony]

        elif status == 1.0:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(
                reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[0], T_), 5)) + list(
                reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, 298), R1(ni2pfree, T_), 5))

            Nicitregionx = list(x_data[1]) + list(x_data[2]) + list(
                np.linspace(x_intercept[2], x_intercept[2], 5)) + list(reversed(x_data[2])) + list(
                reversed(x_data[1])) + list(np.linspace(x_intercept[0], x_intercept[0], 5))
            Nicitregiony = list(y_data_bottom[1]) + list(y_data_bottom[2]) + list(
                np.linspace(R7(T_, Nicit, cit3), T7(x_intercept[2], Nicit, cit3, T_), 5)) + list(
                reversed(y_data_top[2])) + list(reversed(y_data_top[1])) + list(
                np.linspace(T1(ni2pfree, x_intercept[0], T_), R1(ni2pfree, T_), 5))

            nio2regionx = list(reversed(x_data[3])) + list(np.linspace(x_intercept[2], x_intercept[2], 5)) + list(
                x_data[3]) + list([14 for i in range(0, 5)])
            nio2regiony = list(reversed(y_data_bottom[3])) + list(
                np.linspace(R7(T_, Nicit, cit3), T7(x_intercept[2], Nicit, cit3, T_), 5)) + list(y_data_top[3]) + list(
                np.linspace(T8(14, T_), R8(14, T_), 5))

            xs = [nip2regionx, Nicitregionx, nio2regionx]
            ys = [nip2regiony, Nicitregiony, nio2regiony]

        if status == 2 or citrate_total == 0.0:
            nip2regionx = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(
                reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony = list(y_data_bottom[0]) + list(
                np.linspace(R8(x_intercept[0], T_), T8(x_intercept[0], T_), 5)) + list(reversed(y_data_top[0])) + list(
                np.linspace(T1(ni2pfree, 0, T_), R1(ni2pfree, T_), 5))

            nio2regionx = list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(x_data[1]) + list(
                [14 for i in range(0, 5)]) + list(reversed(x_data[1]))
            nio2regiony = list(np.linspace(R8(x_intercept[0], T_), T8(x_intercept[0], T_), 5)) + list(
                y_data_bottom[1]) + list(np.linspace(R8(14, T_), T8(14, T_), 5)) + list(reversed(y_data_top[1]))

            xs = [nip2regionx, nio2regionx]
            ys = [nip2regiony, nio2regiony]

        elif status == 3.0:
            nip2regionx1 = list(x_data[0]) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(
                reversed(x_data[0])) + list([0 for i in range(0, 5)])
            nip2regiony1 = list(y_data_bottom[0]) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[0], T_), 5)) + list(
                reversed(y_data_top[0])) + list(np.linspace(T1(ni2pfree, 0, 298), R1(ni2pfree, T_), 5))

            NiH2citregionx = list(reversed(x_data[1])) + list(np.linspace(x_intercept[0], x_intercept[0], 5)) + list(
                x_data[1]) + list(np.linspace(x_intercept[1], x_intercept[1], 5))
            NiH2citregiony = list(reversed(y_data_bottom[1])) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[0], T_), 5)) + list(y_data_top[1]) + list(
                np.linspace(T2(x_intercept[1], NiH2cit, H3cit, T_), R2(x_intercept[1], T_, NiH2cit, H3cit), 5))

            nip2regionx2 = list(reversed(x_data[2])) + list(np.linspace(x_intercept[1], x_intercept[1], 5)) + list(
                x_data[2]) + list(np.linspace(x_intercept[2], x_intercept[2], 5))
            nip2regiony2 = list(reversed(y_data_bottom[2])) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[2], T_), 5)) + list(y_data_top[2]) + list(
                np.linspace(T1(ni2pfree, x_intercept[2], T_), R1(ni2pfree, T_), 5))

            Nicitregionx = list(x_data[3]) + list(x_data[4]) + list(
                np.linspace(x_intercept[4], x_intercept[4], 5)) + list(reversed(x_data[4])) + list(
                reversed(x_data[3])) + list(np.linspace(x_intercept[2], x_intercept[2], 5))
            Nicitregiony = list(y_data_bottom[3]) + list(y_data_bottom[4]) + list(
                np.linspace(R7(T_, Nicit, cit3), T7(x_intercept[4], Nicit, cit3, T_), 5)) + list(
                reversed(y_data_top[4])) + list(reversed(y_data_top[3])) + list(
                np.linspace(R1(ni2pfree, T_), T1(ni2pfree, x_intercept[2], T_), 5))

            nio2regionx = list(reversed(x_data[5])) + list(np.linspace(x_intercept[4], x_intercept[4], 5)) + list(
                x_data[5]) + list([14 for i in range(0, 5)])
            nio2regiony = list(reversed(y_data_bottom[5])) + list(
                np.linspace(R7(T_, Nicit, cit3), T7(x_intercept[4], Nicit, cit3, T_), 5)) + list(y_data_top[5]) + list(
                np.linspace(T8(14, T_), R8(14, T_), 5))

            xs = [nip2regionx1, NiH2citregionx, nip2regionx2, Nicitregionx, nio2regionx]
            ys = [nip2regiony1, NiH2citregiony, nip2regiony2, Nicitregiony, nio2regiony]
        return [xs,ys,niregionx,niregiony,nio3regionx,nio3regiony]
    # end function, return data to add to traces, should already be in correct form to
    # cooperate with dash notation

    xs = trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[0]
    ys = trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[1]
    niregionx = trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[2]
    niregiony = trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[3]
    nio3regionx = trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[4]
    nio3regiony = trace_generator(pH_x, ni2pfree, NiH2cit, NiHcit, Nicit, H3cit, H2cit, Hcit, cit3, T_)[5]

    if status == 0:
        name = ['Ni<sup>2+</sup>', 'NiH<sub>2</sub>cit<sup>+</sup>','NiHCit', 'Nicit<sup>-</sup>','Ni(OH)<sub>2</sub>']
        color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)', 'rgba(114, 102, 234, 0.63)','rgba(114, 204, 234, 0.63)', 'rgba(52, 207, 42, 0.61)']
    elif status == 1:
        name = ['Ni<sup>2+</sup>', 'Nicit<sup>-</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(101, 18, 183, 0.41)', 'rgba(30, 192, 42, 0.43)', 'rgba(243, 238, 77, 0.5)']
    elif status == 2 or citrate_total==0 :
        name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)']
    elif status == 3:
        name = ['Ni<sup>2+</sup>', 'NiH<sub>2</sub>cit<sup>+</sup>', 'Ni<sup>2+</sup>','Nicit<sup>-</sup>', 'Ni(OH)<sub>2</sub>']
        color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)', 'rgba(114, 102, 234, 0.63)','rgba(114, 204, 234, 0.63)', 'rgba(52, 207, 42, 0.61)']

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
            fill='toself',
            showlegend=True,
            hoverinfo='skip'
        ))
    for trace in tracesextra:
        fig.add_traces(trace)

    if citrate_total != 0:
        fig.add_trace(go.Scatter(
            x=list([2.996 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                    shape='spline',
                    color='green',
                    width=1.5,
                    dash='dash'
                    )
        ))
        fig.add_trace(go.Scatter(
            x=list([4.7478 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                shape='spline',
                color='green',
                width=1.5,
                dash='dash'
            )
        ))
        fig.add_trace(go.Scatter(
            x=list([6.3946 for i in range(len(np.linspace(-1, 2.4, 5)))]),
            y=list(np.linspace(-1, 2.4, 5)),
            mode='lines',
            showlegend=False,
            hoverinfo='skip',
            line=dict(
                shape='spline',
                color='green',
                width=1.5,
                dash='dash'
            )
        ))

        fig.add_annotation(
            x=2.1,
            y=2.2,
            text='H<sub>3</sub>cit',
            showarrow=False,
            font=dict(
                    family="Courier New bold",
                    size=16,
                    color="purple"
                    ),)
        fig.add_annotation(
            x=3.8,
            y=2.2,
            text='H<sub>2</sub>cit<sup>2-</sup>',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
                color="purple"
            )
        )
        fig.add_annotation(
            x=5.6,
            y=2.2,
            text='Hcit<sup>2-</sup>',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
                color="purple"
            )
        )
        fig.add_annotation(
            x=7.2,
            y=2.2,
            text='Cit<sup>3-</sup>',
            showarrow=False,
            font=dict(
                family="Courier New bold",
                size=16,
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