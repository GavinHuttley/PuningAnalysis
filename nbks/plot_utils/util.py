from cogent3 import get_app

load_json_app = get_app("load_json")


def update_figure_format(fig):
    fig.update_layout(
        template="plotly_white",
        width=1200,
        height=500,
        margin={"l": 20, "r": 20, "t": 50, "b": 20},
        autosize=True,
        title_font={"size": 25, "family": "Times New Roman", "color": "black"},
        legend={"font": {"size": 20}},
    )

    fig.update_xaxes(
        title_font={"size": 20, "family": "Times New Roman", "color": "black"},
        tickfont={"size": 20},
        gridcolor="lightgrey",
    )
    fig.update_yaxes(
        title_font={"size": 20, "family": "Times New Roman", "color": "black"},
        tickfont={"size": 20},
        gridcolor="lightgrey",
    )

    return fig
