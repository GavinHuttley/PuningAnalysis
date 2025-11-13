from cogent3 import get_app
import numpy

load_json_app = get_app("load_json")

def update_figure_format(fig):
    fig.update_layout(
        template='plotly_white',
        width=1200,  
        height=500,
        margin=dict(l=20, r=20, t=50, b=20),
        autosize=True,
        title_font=dict(size=25, family='CMU Serif', color='black'),
        legend=dict(
        font=dict(size=14))
    )

    fig.update_xaxes(
        title_font=dict(size=30, family='CMU Serif', color='black'), 
        tickfont=dict(size=20),
        gridcolor='lightgrey'
)
    fig.update_yaxes(
        title_font=dict(size=30, family='CMU Serif', color='black'),
        tickfont=dict(size=20),
        gridcolor='lightgrey'
    )

    return fig

