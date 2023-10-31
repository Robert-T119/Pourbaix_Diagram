from djangopage.reactions import *
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import fsolve
from django_plotly_dash import DjangoDash

app1 = DjangoDash('nickelpure', add_bootstrap_links=True)
# app1.css.append_css({
#     'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
# })

app1.layout = html.Div([
    html.Div([
        html.H1(u'Nickel-H\u2082O system')],
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
        html.H6(u"Total nickel concentration (kmolm\u207B\u00B3):"),
        dcc.Slider(
            id='nickel_slider',
            min=0.1,
            max=2.5,
            value=1.1,
            step=0,
            marks={n_activity: str(n_activity) for n_activity in [0.1, 0.2, 0.3,
                0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5]},
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
                #animate=True,
                id='nickelpure speciation',
                )
        ],
        style={
            'width': '48%',
            'margin-left': '1%',
        }
        ),
        html.Div([
            dcc.Graph(id='nickelpure potential')
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
    Output('nickelpure potential', 'figure'), ##figure????
    [Input('nickel_slider', 'value')])


def nickelpure(nitotal):
    pH_x = np.linspace(0, 14, 50)
    def trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_):

        def interceptgenerator(pH_x, nip2, nin4p2, nin6p2, n, nhp, T_):

            bottom_eqs = [[R1(nip2, T_) for i in range(len(pH_x))], R2(pH_x, T_)]
            y_vars = []
            # y vars has m, +c of every line in bottom eqs
            for i, eq in enumerate(bottom_eqs):
                y_vars.append(np.polyfit(pH_x, eq, 1)) #1 means linear equaiton

            As = [] #M
            Bs = [] #C
            for i, var in enumerate(y_vars):
                if i >= 1:
                    As.append(np.array([[-y_vars[i-1][0], 1], [-y_vars[i][0], 1]]))
                    Bs.append(np.array([y_vars[i-1][1], y_vars[i][1]]))
            # inters has [x_intercept, y_intercept] of all intercepts
            inters = []
            for i, ms in enumerate(As):
                for j, cs in enumerate(Bs):
                    if i == j:
                        inters.append(np.linalg.inv(As[i]).dot(Bs[j]))
            return inters
        # interceptgenerator returns inters is array of x and y values of intercepts
        #------------------------------------------------------------------------------------
        # inters = interceptgenerator(pH_x, nip2, nin4p2, nin6p2, n, nhp, T_, status)[0]
        inters = interceptgenerator(pH_x, nip2, nin4p2, nin6p2, n, nhp, T_)
        xinters = []
        for item in inters:
            xinters.append(item[0])
        x_data = []

        x_data.append(list(np.linspace(0, xinters[0], 5)))
        x_data.append(list(np.linspace(xinters[0], 14, 5)))

        y_data_bottom = [[R1(nip2, T_) for i in range(len(x_data[0]))], R2(x_data[1], T_)]
        new_x_bottom = []
        new_y_bottom = []

        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_bottom.append(xvalue)
        for yvalues in y_data_bottom:
            for yvalue in yvalues:
                new_y_bottom.append(yvalue)

        y_data_top = [T1(nip2, x_data[0], T_), T2(x_data[1], T_)]
        new_x_top = []
        new_y_top = []
        for xvalues in x_data:
            for xvalue in xvalues:
                new_x_top.append(xvalue)
        for yvalues in y_data_top:
            for yvalue in yvalues:
                new_y_top.append(yvalue)

        y_interps = [T1(nip2, inters[0][0], T_)]
        vys = []
        for i, val in enumerate(inters):
            vys.append(list(np.linspace(inters[i][1], y_interps[i], 5)))
        x_data_verticals = []
        for i in range(len(inters)):
            x_data_verticals.append([inters[i][0] for j in range(len(vys[i]))])
        new_x_vert = []
        new_y_vert = []
        for xvalues in x_data_verticals:
            for xvalue in xvalues:
                new_x_vert.append(xvalue)
        for yvalues in vys:
            for yvalue in yvalues:
                new_y_vert.append(yvalue)

        nip2regionx = list(x_data[0]) + list(x_data_verticals[0]) + list(reversed(x_data[0]))
        nip2regiony = list(y_data_bottom[0]) + list(vys[0]) + list(reversed(y_data_top[0]))
        nio3regiony = list(new_y_top) + list(np.linspace(T2(14, T_), 2.4, 5)) + list([2.4 for i in range(0, 5)]) + list(np.linspace(T1(nitotal, 0, 298), 2.4, 5))
        nio3regionx = list(new_x_bottom) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])
        niregiony = list(new_y_bottom) + list(np.linspace(R2(14, T_), -1, 5)) + list([-1 for i in range(0, 5)]) + list(np.linspace(-1, R1(nitotal, 298), 5))
        niregionx = list(new_x_bottom) + list([14 for i in range(0, 5)]) + list(reversed(np.linspace(0, 14, 5))) + list([0 for i in range(0, 5)])

        nio2regionx =  list(reversed(x_data[1])) + list(x_data_verticals[0]) + list(x_data[1])
        nio2regiony = list(reversed(y_data_bottom[1]))  + list(vys[0]) + list(y_data_top[1])
        xs = [nip2regionx, nio2regionx]
        ys = [nip2regiony, nio2regiony]

        return [xs, ys, niregionx, niregiony, nio3regionx, nio3regiony]
    # end function, return data to add to traces, should already be in correct form to
    # cooperate with dash notation
    #------------------------------------------------------------------------------------------------
    nip2 = nip2ppt = nitotal
    nhp = n = ntotal = nin4p2 = nin6p2 = 0
    T_ = 298

    xs = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[0]
    ys = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[1]
    niregionx = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[2]
    niregiony = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[3]
    nio3regionx = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[4]
    nio3regiony = trace_generator(pH_x, nip2, nip2ppt, nin4p2, nin6p2, n, nhp, T_)[5]
    name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
    color = ['rgba(191, 63, 63, 0.5)', 'rgba(243, 238, 77, 0.5)']

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
        # width=600,
        # height=500,
            )

    fig = go.Figure(data=data, layout=layout)
    extrax = [niregionx, nio3regionx]
    extray = [niregiony, nio3regiony]
    namee = ['Ni', 'Ni(OH)<sub>3</sub>']
    colors1 = ['rgba(127, 63, 191, 0.5)', 'rgba(30, 205, 40, 0.5)']
    tracesextra = []
    for i, extraxx in enumerate(extrax):
        tracesextra.append(go.Scatter(
            x = extraxx,
            y = extray[i],
            mode='lines',
            name=namee[i],
            fill='toself',
            showlegend=True,
            hoverinfo='skip'
        ))
    for trace in tracesextra:
        fig.add_traces(trace)

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
        x=0.35,
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

@app1.callback(
    Output('nickelpure speciation', 'figure'),
    [Input('nickel_slider', 'value')])
def nickelpure1(nitotal):
    def concs(pH_x, nitotal):
        logk = 11.853
        f = 10**(logk - 2*pH_x)
        nip2free = f/(1+f/nitotal)
        nio2 = nitotal - nip2free
        return [nip2free, nio2]

    pH_x = list(np.linspace(0, 14, 50))
    nip2freeplot = []
    nio2plot = []
    for pHval in pH_x:
        nip2freeplot.append(concs(pHval, nitotal)[0])
        nio2plot.append(concs(pHval, nitotal)[1])

    datasets = [nip2freeplot, nio2plot]
    name = ['Ni<sup>2+</sup>', 'Ni(OH)<sub>2</sub>']
    fill = [None for i in range(len(name))]
    color = ['rgb(210, 80, 80)', 'rgb(235, 154, 14)']

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
    fig1.update_xaxes(gridcolor='white', range=[0, 14],
                     nticks=20, mirror=True, ticks='outside', showline=True)
    fig1.update_yaxes(gridcolor='white', ticks='outside',
                    range=[0, nitotal*1.05])
    fig1.update_layout(
        title={
            'text': "Speciation plot",
            'y': 0.95,
            'x': 0.45,
            'xanchor': 'center',
            'yanchor': 'top'
            })
    return fig1
