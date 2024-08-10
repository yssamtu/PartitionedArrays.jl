import matplotlib.pyplot as plt


def draw_legend(ax, figsize, save_file_name: str):
    fig_legend, ax_legend = plt.subplots(figsize=(figsize[0] * 0.5, figsize[1] * 0.5))
    handles, labels = ax.get_legend_handles_labels()
    legend = fig_legend.legend(handles, labels, loc="center", ncol=2, frameon=False)
    ax_legend.axis("off")
    fig_legend.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig_legend.dpi_scale_trans.inverted())
    fig_legend.set_size_inches(bbox.width, bbox.height)
    fig_legend.savefig(f"{save_file_name}_legend.pdf")
    plt.close(fig_legend)
