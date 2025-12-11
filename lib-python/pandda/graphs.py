import giant.logs as lg

logger = lg.getLogger(__name__)

try:
    import matplotlib as mpl

    mpl.use("agg")
    mpl.interactive(False)
    from matplotlib import pyplot as plt

    plt.switch_backend("agg")  # yes I know this done twice -- for safety!
    plt.interactive(0)
except (AttributeError, ImportError) as e:
    # Handle NumPy 2.x compatibility issues with matplotlib
    import sys

    error_msg = str(e)
    if "_ARRAY_API" in error_msg or "numpy" in error_msg.lower():
        print(
            "\n" + "=" * 80 + "\n"
            "ERROR: NumPy 2.x compatibility issue with matplotlib!\n"
            "=" * 80 + "\n"
            "Matplotlib in CCP4 Python was compiled against NumPy 1.x but NumPy 2.x is installed.\n\n"
            "SOLUTION: Downgrade NumPy in the CCP4 Python environment:\n"
            "  /Applications/ccp4-9/bin/ccp4-python -m pip install 'numpy<2'\n\n"
            "Or upgrade matplotlib to a version compatible with NumPy 2.x:\n"
            "  /Applications/ccp4-9/bin/ccp4-python -m pip install --upgrade matplotlib\n\n"
            "Original error: " + error_msg + "\n"
            "=" * 80 + "\n",
            file=sys.stderr,
            flush=True,
        )
        raise
    else:
        # Other import errors - try to continue
        logger.warning(f"Error importing matplotlib: {e}")
        try:
            import matplotlib as mpl
            from matplotlib import pyplot as plt
        except Exception as e2:
            logger.error(f"Failed to import matplotlib: {e2}")
            raise
except Exception as e:
    # Catch-all for other errors
    logger.warning(f"Error setting up matplotlib: {e}")
    try:
        import matplotlib as mpl
        from matplotlib import pyplot as plt
    except Exception as e2:
            logger.error(f"Failed to import matplotlib: {e2}")
            raise

import collections  # noqa: E402
import json  # noqa: E402
import pathlib as pl  # noqa: E402

import numpy as np  # noqa: E402


def failure_graph(title):
    fig = plt.figure()
    plt.title("Failed to make {}".format(title))
    return fig


class PlotConfig(object):

    def __init__(self):

        self.colour_map_name = None
        self.plot_style = None
        self.font_family = None
        self.font_name = None

    def configure(
        self,
        colour_map_name=None,
        plot_style=None,
        font_family=None,
        font_name=None,
    ):

        if colour_map_name is not None:
            self.colour_map_name = colour_map_name

        if plot_style is not None:
            self.set_plot_style(
                plot_style=plot_style,
            )

        if font_family is not None:
            self.set_font_family(
                font_family=font_family,
            )

        if font_name is not None:
            self.set_font(
                font_name=font_name,
                font_family=font_family,
            )

    def set_plot_style(
        self,
        plot_style,
    ):

        if plot_style == "xkcd":

            try:

                plt.xkcd()

                logger('Theme set to "xkcd".')

            except Exception as e:

                logger.warning('Failed to set plot style "xkcd".')

            self.plot_style = plot_style

            return None

        try:

            logger('Setting plot_style to "{}"'.format(plot_style))

            plt.style.use(plot_style)

            self.plot_style = plot_style

        except Exception as e:

            logger.warning(
                "Failed to set plot style to {}.\n\t{}".format(plot_style, str(e))
            )

    def set_font_family(
        self,
        font_family,
    ):

        try:

            logger('Setting font_family to "{}"'.format(font_family))

            plt.rc("font", family=font_family)

            self.font_family = font_family

        except Exception as e:

            logger.warning(
                "Failed to set font family to {}.\n\t{}".format(font_family, str(e))
            )

    def set_font(
        self,
        font_name,
        font_family,
    ):

        try:

            if font_family is None:
                logger.warning(
                    "Cannot set font: must provide a font family in order to set font."
                )
                return None

            font_family_str = "font." + str(font_family)

            if font_family_str not in list(plt.rcParams.keys()):
                logger.warning(
                    'Cannot set font: invalid font family provided "{}".'.format(
                        font_family
                    )
                )
                return None

            family_fonts = plt.rcParams[font_family_str]

            if font_name not in family_fonts:
                logger.warning(
                    'font "{}" does not exist in font family "{}". Setting the font may not work. (valid options: {})'.format(
                        font_name,
                        font_family,
                        ", ".join(family_fonts),
                    )
                )

            logger('Setting font_name to "{}"'.format(font_name))

            plt.rcParams[font_family_str].insert(0, font_name)

            self.font_name = font_name

        except Exception as e:

            logger.warning("Failed to set font to {}.\n\t{}".format(font_name, str(e)))

    def get_colour_map(self):
        return mpl.cm.get_cmap(self.colour_map_name)

    def get_colours(self, n):
        cm = self.get_colour_map()
        return cm(np.linspace(0.0, 1.0, n))


config = PlotConfig()


def configure(
    colour_map_name="rainbow",
    plot_style="ggplot",
    font_family="monospace",
    font_name=None,
):

    global config

    config.configure(
        colour_map_name=colour_map_name,
        plot_style=plot_style,
        font_family=font_family,
        font_name=font_name,
    )

    return config


#####


class JsonPlotter(object):
    """Format graph data into json for html plotting"""

    def __init__(self):
        self.data = {}

    def __str__(self):
        return self.as_json()

    def as_json(self):
        return json.dumps(self.data)

    def as_javascript(self):
        return ""


class PanddaPlotter(object):

    output_key = None

    def __init__(self, output_path):

        self.output_path = output_path

    def __call__(self, *args, **kwargs):

        fig = self.plot(*args, **kwargs)

        # try:
        #     fig = self.plot(*args, **kwargs)
        # except:
        #     fig = failure_graph(title='error')

        filename = self.get_path()

        self.save(
            fig=fig,
            filename=filename,
        )

        return {self.output_key: filename}

    def setup(self, nrows=1, ncols=1, **kwargs):

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, **kwargs)

        return fig, axes

    def save(self, fig, filename):

        fig.tight_layout()
        fig.savefig(filename)
        plt.close(fig)

    def get_path(self):

        p = pl.Path(self.output_path)

        if not p.parent.exists():
            p.parent.mkdir(parents=True)

        return str(p)

    def json(self):
        return self.json_plotter.as_json()


class PanddaDatasetPlotter(PanddaPlotter):

    output_key = None

    def __init__(self, output_path_template):

        self.output_path_template = str(output_path_template)
        self.json_plotter = JsonPlotter()

        assert "{label}" in output_path_template  # remove this

    def __call__(self, datasets, *args, **kwargs):

        output_files = collections.OrderedDict()

        for dkey, dataset in sorted(datasets.items()):

            fig = self.plot(dataset=dataset, dataset_label=dkey, *args, **kwargs)

            # try:
            #     fig = self.plot(
            #         dataset = dataset,
            #         dataset_label = dkey,
            #         *args, **kwargs
            #         )
            # except:
            #     fig = failure_graph(title='error')

            filename = self.get_path(
                dataset_label=dkey,
            )

            self.save(
                fig=fig,
                filename=filename,
            )

            output_files[dkey] = filename

        return {self.output_key: output_files}

    def get_label(self, dataset=None):

        if dataset is None:
            return None

        return dataset.tag

    def get_path(self, dataset_label=None, **kwargs):

        p = pl.Path(self.output_path_template.format(label=dataset_label, **kwargs))

        if not p.parent.exists():
            p.parent.mkdir(parents=True)

        return str(p)


class PanddaMultiDatasetPlotter(PanddaPlotter):

    output_key = None

    def __init__(self, output_path):

        self.output_path = str(output_path)
        self.json_plotter = JsonPlotter()

    def __call__(self, datasets, *args, **kwargs):

        dataset_keys = sorted(datasets.keys())
        dataset_list = [datasets[k] for k in dataset_keys]

        fig = self.plot(
            dataset_labels=dataset_keys, dataset_list=dataset_list, *args, **kwargs
        )

        # try:
        #     fig = self.plot(
        #         dataset_labels = dataset_keys,
        #         dataset_list = dataset_list,
        #         *args, **kwargs
        #         )
        # except:
        #     fig = failure_graph(title='error')

        filename = self.get_path()

        self.save(
            fig=fig,
            filename=filename,
        )

        return {self.output_key: filename}
