from djangopage.reactions2 import *
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from django_plotly_dash import DjangoDash

app1 = DjangoDash('copper', add_bootstrap_links=True)
# app1.css.append_css({
#     'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
# })

app1.layout = html.Div([
    html.Div([
        html.H1(u'Copper-H\u2082O system')],
        style={
            'text-align':'center',
            'border': 'thin lightgrey solid',
            'backgroundColor': '#FEFDEB',
            'padding': '-4px 0 -4px 0',
            'margin-bottom': '2px',
            'color': '#10328f',
        }
        ),

    html.Div([
        html.H6(u"Total copper(I) concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='copper_slider1',
            min=0.1,
            max=2.0,
            value=1.0,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2]},
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
        html.H6(u"Total copper(II) concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='copper_slider2',
            min=0.1,
            max=2.0,
            value=1.0,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
              0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4,1.5,
              1.6, 1.7, 1.8, 1.9, 2]},
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
        html.Div([
            dcc.Graph(
                # animate=True,
                id='copperpure1 speciation',
                )
        ],
        style={
            'width': '48%',
            'margin-left': '1%',
        }
        ),
        html.Div([
            dcc.Graph(id='copperpure potential')
        ],
        style={
            'width': '48%',
            'margin-right': '1%',
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
    html.Div([
        html.Div([
            dcc.Graph(
                # animate=True,
                id='copperpure2 speciation',
            )
        ],
            style={
                'width': '48%',
                'margin-left': '1%',
            }
        ),
        html.Div([
            dcc.Graph(id='copperpure3 speciation')
        ],
            style={
                'width': '48%',
                'margin-right': '1%',
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
    ],

    style = {
        'margin-top': '40px',
        # "margin-left": "10px",
        # 'margin-right': '10px',
        'margin-bottom': '5px'
        }
    )


@app1.callback(
    Output('copperpure potential', 'figure'),
    [Input('copper_slider2', 'value')])

def speciation_graph(Cu_total):
    cu2p = cuo2 = Cu_total
    T_ = 298
    pH_x =np.linspace(-2,16,71)

    def trace_generator(pH_x, cu2p, cuo2, T_):

        def interceptgenerator(pH_x, cu2p, cuo2, T_):

            bottom_eqs = [[R1(cu2p, T_) for i in range(len(pH_x))], R2(pH_x, T_, cu2p), R3(pH_x, T_),
                          R4(pH_x, T_, cuo2)]
            y_vars = []
            # y vars has m, +c of every line in bottom eqs
            for i, eq in enumerate(bottom_eqs):
                y_vars.append(np.polyfit(pH_x, eq, 1))  # 1 means linear equaiton

            As = []  # M
            Bs = []  # C
            for i, var in enumerate(y_vars):
                if i >= 1:
                    As.append(np.array([[-y_vars[i - 1][0], 1], [-y_vars[i][0], 1]]))
                    Bs.append(np.array([y_vars[i - 1][1], y_vars[i][1]]))
            # inters has [x_intercept, y_intercept] of all intercepts
            inters = []
            for i, ms in enumerate(As):
                for j, cs in enumerate(Bs):
                    if i == j:
                        inters.append(np.linalg.inv(As[i]).dot(Bs[j]))
            return inters

        # interceptgenerator returns inters is array of x and y values of intercepts
        inters = interceptgenerator(pH_x, cu2p, cuo2, T_)
        xinters = []
        for item in inters:
            xinters.append(item[0])

        x_data = []
        for i, item in enumerate(xinters):
            if i == 0:
                x_data.append(np.linspace(-2, item, 5))
            elif i >= 1:
                x_data.append(np.linspace(xinters[i - 1], item, 5))
        finalindex = len(xinters) - 1
        x_data.append(np.linspace(xinters[finalindex], 16, 5))

        y_data_bottom = [[R1(cu2p, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_, cu2p), R3(x_data[2], T_),
                         R4(x_data[3], T_, cuo2)]
        new_x_bottom = []
        new_y_bottom = []

        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_bottom.append(xvalue)
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        x_data_verticals_1 = np.linspace(inters[1][0], inters[1][0], 5)
        y_data_verticals_1 = np.linspace(inters[1][1], 2.4, 5)
        x_data_verticals_2 = np.linspace(inters[2][0], inters[2][0], 5)
        y_data_verticals_2 = np.linspace(inters[2][1], 2.4, 5)
        T1_xdata = np.linspace(inters[0][0], 16, 5)
        T1_ydata = np.linspace(inters[0][1], B1(16, 298), 5)
        abx = np.linspace(0, x_data[0], 5)

        cu2pregionx = list(x_data[0]) + list(x_data[1]) + list(x_data_verticals_1) + list(
            np.linspace(-2, inters[1][0], 5)) + list([-2 for i in range(0, 5)])
        cu2pregiony = list(y_data_bottom[0]) + list(y_data_bottom[1]) + list(y_data_verticals_1) + list(
            [2.4 for i in range(0, 5)]) + list(reversed(np.linspace(inters[0][1], 2.4, 5)))

        cuoregionx = list(x_data_verticals_1) + list(x_data[2]) + list(x_data_verticals_2) + list(reversed(x_data[2]))
        cuoregiony = list(y_data_verticals_1) + list(y_data_bottom[2]) + list(y_data_verticals_2) + list(
            [2.4 for i in range(0, 5)])

        cuo2regionx = list(x_data_verticals_2) + list(x_data[3]) + list([16 for i in range(0, 5)]) + list(
            reversed(x_data[3]))
        cuo2regiony = list(y_data_verticals_2) + list(y_data_bottom[3]) + list(
            np.linspace(inters[2][1], 2.4, 5)) + list([2.4 for i in range(0, 5)])

        cu2Oregionx = list(x_data[1]) + list(x_data[2]) + list(x_data[3]) + list(reversed(T1_xdata))
        cu2Oregiony = list(y_data_bottom[1]) + list(y_data_bottom[2]) + list(y_data_bottom[3]) + list(
            reversed(T1_ydata))

        curegionx = list(x_data[0]) + list(T1_xdata) + list([16 for i in range(0, 5)]) + list(
            reversed(np.linspace(-2, 16, 5))) + list(-2 for i in range(0, 5))
        curegiony = list(y_data_bottom[0]) + list(T1_ydata) + list(np.linspace(B1(16, 298), -1, 5)) + list(
            -1 for i in range(0, 5)) + list(np.linspace(-1, inters[0][1], 5))

        xs = [cu2pregionx, cuoregionx, cuo2regionx, cu2Oregionx, curegionx]
        ys = [cu2pregiony, cuoregiony, cuo2regiony, cu2Oregiony, curegiony]
        return [xs, ys]
#------------------------------------------------------------------------------------------------
    xs = trace_generator(pH_x, cu2p, cuo2, T_)[0]
    ys = trace_generator(pH_x, cu2p, cuo2, T_)[1]

    name = ['Cu<sup>2+</sup>', 'CuO', 'CuO<sub>2</sub>', 'Cu<sub>2</sub>O','Cu']
    color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)', 'rgba(252, 177, 101, 0.8)','rgba(7, 117, 189, 0.66)', 'rgba(63, 63, 191, 0.5)']

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

    ywater = [W1(pH_x, T_), W2(pH_x, T_)]
    for ys in ywater:
        data.append(go.Scatter(
            x= pH_x,
            y= ys,
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
        # width=600,
        # height=500,
            )

    fig = go.Figure(data=data, layout=layout)
    fig.add_annotation(
        x=-1.5,
        y= 1.45,
        showarrow=False,
        text='O<sub>2</sub>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue',
        )
    )
    fig.add_annotation(
        x=-1.37,
        y= 1.15,
        showarrow=False,
        text='H<sub>2</sub>O',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.add_annotation(
        x=-1.5,
        y= 0.2,
        showarrow=False,
        text='H<sup>+</sup>',
        font=dict(
            family='Courier New bold',
            size=18,
            color='blue'
        )
    )
    fig.add_annotation(
        x= -1.5,
        y=-0.1,
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
    fig.update_layout(yaxis=dict(range=[-1,1.8]), xaxis=dict(range=[-2,16]),
                        title={
                            'text': "Potential-pH diagram",
                            'y':0.95,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top'
                            })
    return fig

@app1.callback(
    Output('copperpure1 speciation', 'figure'),
    [Input('copper_slider2', 'value')])
def Copperpure1(Cu_total):
    pH_x = np.linspace(0, 16, 70)
    cu2p = 10 ** (6.57 - 2 * pH_x) / (1 + (10 ** (6.57 - 2 * pH_x)) / Cu_total)
    hcuo2 = 10 ** (-19.67 + pH_x)
    cuo2 = 10 ** (-32.81 + 2 * pH_x)
    cuo = Cu_total - cu2p - hcuo2 - cuo2

    datasets = [cu2p, hcuo2, cuo2, cuo]
    name = ['Cu<sup>2+</sup>', 'HCuO<sub>2</sub>', 'CuO<sub>2</sub>', 'CuO']
    fill = [None, None, None,None]
    color = ['rgb(210, 80, 80)', 'rgb(90 ,0, 100)','rgb(40, 130, 80)', 'rgb(235, 154, 14)']

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
        yaxis={'title': 'Concentration (kmolm<sup>-3</sup>)', 'linecolor': 'grey',
            'mirror':True},
        # transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'l':10},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgb(240,240,240)',
        autosize=False,
        # width=600,
        # height=500,
    )
    fig1 = go.Figure(data=data1, layout=layout)
    fig1.update_xaxes(gridcolor='white', range=[0, 16],
                     nticks=20, mirror=True, ticks='outside', showline=True)
    fig1.update_yaxes(gridcolor='white', ticks='outside',
                    range=[0, Cu_total*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot [Cu(II)/CuO]",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1

@app1.callback(
    Output('copperpure2 speciation', 'figure'),
    [Input('copper_slider2', 'value')])
def nickelpure1(Cu_total):
    pH_x = np.linspace(0, 14, 50)
    def concs(pH_x, Cu_total):
        cu2p = 10 ** (8.61 - 2 * pH_x) / (1 + (10 ** (8.61 - 2 * pH_x)) /Cu_total)
        hcuo2 = 10 ** (-17.62 + pH_x)
        cuo2 = 10 ** (-30.76 + 2 * pH_x)
        cuoh2 = Cu_total -cuo2 -cu2p -hcuo2
        return [cu2p, hcuo2, cuo2, cuoh2]

    cu2pfreeplot = []
    hcuo2plot = []
    cuo2plot = []
    cuoh2plot = []
    for pHval in pH_x:
        cu2pfreeplot.append(concs(pHval, Cu_total)[0])
        hcuo2plot.append(concs(pHval, Cu_total)[1])
        cuo2plot.append(concs(pHval, Cu_total)[2])
        cuoh2plot.append(concs(pHval, Cu_total)[3])

    datasets = [cu2pfreeplot, hcuo2plot, cuo2plot, cuoh2plot]
    name = ['Cu<sup>2+</sup>', 'HCuO<sub>2</sub>', 'CuO<sub>2</sub>', 'Cu(OH)<sub>2</sub>']
    fill = [None, None, None,None]
    color = ['rgb(210, 80, 80)', 'rgb(90 ,0, 100)','rgb(40, 130, 80)', 'rgb(235, 154, 14)']

    data2 = []
    for i, dataset in enumerate(datasets):
        data2.append(go.Scatter(
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
        yaxis={'title': 'Concentration (kmolm<sup>-3</sup>)', 'linecolor': 'grey',
            'mirror':True},
        # transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'l':10},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgb(240,240,240)',
        autosize=False,
        # width=600,
        # height=500,
    )
    fig2 = go.Figure(data=data2, layout=layout)
    fig2.update_xaxes(gridcolor='white', range=[0, 14],
                     nticks=20, mirror=True, ticks='outside', showline=True)
    fig2.update_yaxes(gridcolor='white', ticks='outside',
                    range=[0, Cu_total*1.05])
    fig2.update_layout(
        title={
            'text': "Speciation plot [Cu(II)/Cu(OH)<sub>2</sub>]",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig2

@app1.callback(
    Output('copperpure3 speciation', 'figure'),
    [Input('copper_slider1', 'value')])
def nickelpure1(Cu_total):
    pH_x = np.linspace(0, 16, 50)
    def concs(pH_x, Cu_total):
        cup = 10 ** (-1.01 - pH_x)
        cu2o = Cu_total - cup
        return [cup, cu2o]

    cupfreeplot = []
    cu2oplot = []
    for pHval in pH_x:
        cupfreeplot.append(concs(pHval, Cu_total)[0])
        cu2oplot.append(concs(pHval, Cu_total)[1])

    datasets = [cupfreeplot, cu2oplot]
    name = ['Cu<sup>+</sup>', 'Cu<sub>2</sub>O']
    fill = [None, None, None,None]
    color = ['rgb(210, 80, 80)', 'rgb(90 ,0, 100)']

    data3 = []
    for i, dataset in enumerate(datasets):
        data3.append(go.Scatter(
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
        yaxis={'title': 'Concentration (kmolm<sup>-3</sup>)', 'linecolor': 'grey',
            'mirror':True},
        # transition = {'duration': 1200},
        font=dict(family='Courier Sans', color='grey'),
        margin={'t': 50, 'l':10},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgb(240,240,240)',
        autosize=False,
        # width=600,
        # height=500,
    )
    fig3 = go.Figure(data=data3, layout=layout)
    fig3.update_xaxes(gridcolor='white', range=[0, 16],
                     nticks=20, mirror=True, ticks='outside', showline=True)
    fig3.update_yaxes(gridcolor='white', ticks='outside',
                    range=[0, Cu_total*1.05])
    fig3.update_layout(
        title={
            'text': "Speciation plot [Cu(I)/Cu<sub>2</sub>O]",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig3
