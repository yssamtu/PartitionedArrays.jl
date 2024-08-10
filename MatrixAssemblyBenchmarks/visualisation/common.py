import matplotlib.pyplot as plt

def change_func_name(data: dict):
    func_map = {
        "psparse": "psparse",
        "petsc_setvalues": "petsc_setvalues",
        "petsc_coo": "petsc_coo",
        "assemble_matrix_no_compressed_snd_and_with_int_vector_cache": "new_map",
        "assemble_matrix_no_compressed_snd_and_with_tuple_vector_cache": "new_pair",
        "assemble_matrix_no_compressed_snd_and_with_auto_cache": "new_auto",
        "assemble_matrix_with_compressed_snd_and_with_int_vector_cache": "new_nodup_map",
        "assemble_matrix_with_compressed_snd_and_with_tuple_vector_cache": "new_nodup_pair",
        "assemble_matrix_with_compressed_snd_and_with_auto_cache": "new_nodup_auto"
    }
    return {func_map[key]: value for key, value in data.items()}

def draw_legend(ax, figsize, save_file_name: str):
    fig_legend, ax_legend = plt.subplots(figsize=(figsize[0] * 0.5, figsize[1] * 0.5))
    handles, labels = ax.get_legend_handles_labels()
    legend = fig_legend.legend(handles, labels, loc="center", ncol=3, frameon=False)
    ax_legend.axis("off")
    fig_legend.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig_legend.dpi_scale_trans.inverted())
    fig_legend.set_size_inches(bbox.width, bbox.height)
    fig_legend.savefig(f"{save_file_name}_legend.pdf")
    plt.close(fig_legend)
