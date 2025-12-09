### TMP HACK ############################################
# import giant.logs as lg
# logger = lg.getLogger(__name__)
# lg.setup_logging_basic(None)
###
import os, sys, copy
import logging as lg
logger = lg.getLogger()
logger.setLevel(lg.DEBUG)
ch = lg.StreamHandler(stream=sys.stdout)
ch.setFormatter(lg.Formatter(fmt='%(message)s'))
ch.setLevel(lg.DEBUG)
logger.addHandler(ch)
### TMP HACK ############################################

import os, json, collections
import pathlib as pl
import pandas as pd 
import numpy as np

### TMP HACK ############################################
# from pandda import (
#     resources,
#     )
### TMP HACK ############################################

### TMP HACK ############################################
# from pandda.inspect.exceptions import (
#     MissingFile,
#     )
###
class MissingFile(Exception):
    pass
### TMP HACK ############################################

### TMP HACK ############################################
# from pandda.inspect.parser import (
#     parser
#     )
###
import argparse

parser = argparse.ArgumentParser(
    prog = "pandda.inspect",
    )

parser.add_argument(
    "-m", 
    "--mode", 
    help = "mode to run in", 
    choices = [
        "events", "datasets",
        ],
    default = "events",
    type = str,
    )

parser.add_argument(
    "-p",
    "--pandda_directory",
    help = "input/output pandda directory", 
    default = ".",
    type = str,
    )

# parser.parse_args()
### TMP HACK ############################################


### TMP HACK ############################################
# IMG_DIR = pl.Path(os.path.realpath(resources.__path__[0]))
# LOGO_PATH = resources.SMALL_LOGO_PATH
###
IMG_DIR = pl.Path(sys.argv[-1]).parent / '../resources'
print(IMG_DIR)
LOGO_PATH = str(IMG_DIR / 'pandda-logo-small.png')
### TMP HACK ############################################

# ANGLES FOR COLOURS
COLOUR_YELLOW = 20.0
COLOUR_GREEN = 80.0
COLOUR_BLUE = 200.0
COLOUR_MAGENTA = 260.0

# COLOUR FOR THE PROTEIN
MOL_COLOUR = COLOUR_GREEN
# COLOUR FOR THE LIGANDS
LIG_COLOUR = COLOUR_MAGENTA

# Ordering the dict is IMPORTANT
pandda_column_labels = collections.OrderedDict(
    [
        (
            'dataset', {
                'F' : 'FDATASET', 
                'PHI': 'PHDATASET', 
                'diff': 0,
                },
            ),
        (
            'ground', {
                'F' : 'FGROUND', 
                'PHI': 'PHGROUND', 
                'diff': 0,
                },
            ),
        (
            'zvalues', {
                'F' : 'FZVALUES', 
                'PHI': 'PHZVALUES', 
                'diff': 1,
                },
            ),
        (
            'event', {
                'F' : 'FEVENT', 
                'PHI': 'PHEVENT', 
                'diff': 0,
                },
            ),
        ]
    )

# =========================================================================
# GTK FUNCTIONS
# =========================================================================

def catchup(block=False):
    while gtk.events_pending():
        gtk.main_iteration(block)

def nonmodal_msg(msg):
    """Display an error window - non-modal"""
    d = gtk.MessageDialog(
        type = gtk.MESSAGE_INFO,
        message_format = msg,
        )
    d.set_position(gtk.WIN_POS_CENTER)
    d.set_keep_above(True)
    d.show_all()
    catchup(True)
    return d

def modal_msg(msg):
    """Display an error window - model"""
    d = gtk.MessageDialog(
        type = gtk.MESSAGE_INFO,
        buttons = gtk.BUTTONS_CLOSE,
        message_format = msg,
        )
    d.set_position(gtk.WIN_POS_CENTER)
    d.set_keep_above(True)
    d.run()
    d.destroy()

# =========================================================================
# COOT FUNCTIONS
# =========================================================================

def set_main_coot_molecule(i):
    # Have to set colour manually!
    set_molecule_bonds_colour_map_rotation(i, MOL_COLOUR);
    graphics_draw()
    # Other settings
    set_pointer_atom_molecule(i)
    set_go_to_atom_molecule(i)
    update_go_to_atom_window_on_new_mol()

def coot_customisation():
    set_nomenclature_errors_on_read("ignore")
    set_recentre_on_read_pdb(0)
    set_show_symmetry_master(1)
    set_symmetry_shift_search_size(2)
    set_colour_map_rotation_for_map(0.0)
    set_colour_map_rotation_on_read_pdb(0.0)
    set_colour_map_rotation_on_read_pdb_flag(1)
    set_colour_map_rotation_on_read_pdb_c_only_flag(1)
    #set-stop-scroll-iso-map 0

    add_key_binding("Add ligand", "a", lambda: solvent_ligands_gui())
    add_key_binding("Add water", "w", lambda: place_typed_atom_at_pointer("Water"))

def coot_setup():
    # The order these are added is important -- decides order of loading
    for k, v in pandda_column_labels.items():
        set_auto_read_column_labels(v['F'], v['PHI'], v['diff']) 

def post_coot_windows():
    
    post_display_control_window()
    # post_go_to_atom_window()
    try:
        # Future-available function
        pass
        # post_delete_item_dialog()
    except:
        pass

# =========================================================================

def format_ligand_names(ligand_name_list, new_name):
    names = [ln for ln in ligand_name_list if ln]
    if new_name not in names: 
        names.append(new_name)
    return ','.join(names)


class Notices(object):
    def __init__(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_position(gtk.WIN_POS_CENTER)
        self.window.set_border_width(10)
        self.window.set_title("Map contouring and interpretation")
        # ---------------------------
        main_vbox = gtk.VBox(spacing=5)
        # ---------------------------
        box = gtk.HBox(spacing=5)
        label = gtk.Label(
            'Maps from PanDDA are pre-scaled, so the rmsd value given by coot is not informative. The value given in e/A^3 is actually the sigma-scaled value',
            )
        label.props.width_chars = 100
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.set_line_wrap(True)
        box.pack_start(label, expand=True)
        # ---------------------------
        main_vbox.pack_start(box)
        # ---------------------------
        box = gtk.HBox(spacing=5)
        label = gtk.Label('NOTE: This only applies to the PanDDA maps (those loaded from files ending with *.ccp4)')
        label.props.width_chars = 100
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.set_line_wrap(True)
        box.pack_start(label, expand=True)
        # ---------------------------
        main_vbox.pack_start(box)
        main_vbox.pack_start(gtk.HSeparator())
        # ---------------------------
        image = gtk.Image()
        image.set_from_file(str(IMG_DIR / '1-sigma-anno.png'))
        image.show()
        label = gtk.Label('This map is contoured at 1 sigma')
        vbox_1 = gtk.VBox(spacing=5)
        vbox_1.pack_start(image)
        vbox_1.pack_start(label)
        image = gtk.Image()
        image.set_from_file(str(IMG_DIR / '0.2-sigma-anno.png'))
        image.show()
        label = gtk.Label('This map is contoured at 0.2 sigma')
        vbox_2 = gtk.VBox(spacing=5)
        vbox_2.pack_start(image, expand=False)
        vbox_2.pack_start(label, expand=False)
        box = gtk.HBox(spacing=5, homogeneous=True)
        box.pack_start(vbox_1, expand=True)
        box.pack_start(vbox_2, expand=True)
        # ---------------------------
        main_vbox.pack_start(box)
        main_vbox.pack_start(gtk.HSeparator())
        # ---------------------------
        box = gtk.HBox(spacing=5)
        label = gtk.Label(
            'However, things are further complicated when looking at EVENT MAPS. Since these are partial-difference maps, the EFFECTIVE sigma-value is calculated by dividing the contour level by the (1-BDC) value.')
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.props.width_chars = 100
        label.set_line_wrap(True)
        box.pack_start(label)
        main_vbox.pack_start(box)
        # ---------------------------
        main_vbox.pack_start(gtk.HSeparator())
        # ---------------------------
        box = gtk.HBox(spacing=5)
        label = gtk.Label('For event maps where (1-BDC) is 0.2:')
        box.pack_start(label)
        # ---------------------------
        main_vbox.pack_start(box)
        # ---------------------------
        image = gtk.Image()
        image.set_from_file(str(IMG_DIR / '1-sigma-anno.png'))
        image.show()
        label = gtk.Label(
            '(1-BDC) is 0.2, so the event map is contoured at an effective value of 5 Sigma (1.0/0.2 = 5)')
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.props.width_chars = 30
        label.set_line_wrap(True)
        vbox_1 = gtk.VBox(spacing=5)
        vbox_1.pack_start(image)
        vbox_1.pack_start(label)
        image = gtk.Image()
        image.set_from_file(str(IMG_DIR / '0.2-sigma-anno.png'))
        image.show()
        label = gtk.Label(
            '(1-BDC) is 0.2, so the event map is contoured at an effective value of 1 Sigma (0.2/0.2 = 1)')
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.props.width_chars = 30
        label.set_line_wrap(True)
        vbox_2 = gtk.VBox(spacing=5)
        vbox_2.pack_start(image)
        vbox_2.pack_start(label)
        box = gtk.HBox(spacing=5, homogeneous=True)
        box.pack_start(vbox_1, expand=True)
        box.pack_start(vbox_2, expand=True)
        # ---------------------------
        main_vbox.pack_start(box)
        # ---------------------------
        button = gtk.Button(label='Close window')
        button.child.set_padding(5, 5)
        button.child.set_justify(gtk.JUSTIFY_CENTER)
        button.child.set_line_wrap(True)
        button.connect("clicked", lambda x: self.window.destroy())
        box = gtk.HBox(spacing=5, homogeneous=False)
        box.pack_start(gtk.Label(''), expand=True)
        box.pack_start(button, expand=False)
        box.pack_start(gtk.Label(''), expand=True)
        # ---------------------------
        main_vbox.pack_end(box)
        # ---------------------------
        self.window.add(main_vbox)
        self.window.show_all()
        self.window.set_keep_above(True)


class SplashScreen(object):

    def __init__(self):
        # ---------------------------
        # Create main window
        # ---------------------------
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_position(gtk.WIN_POS_CENTER_ALWAYS)
        self.window.set_border_width(0)
        self.window.set_decorated(False)
        self.window.set_default_size(300, 300)
        # ---------------------------
        self.main_hbox = gtk.HBox(spacing=0)
        frame1 = gtk.EventBox()
        frame1.set_border_width(5)
        frame1.modify_bg(gtk.STATE_NORMAL, None)
        frame1.add(self.main_hbox)
        frame2 = gtk.EventBox()
        frame2.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("black"))
        frame2.add(frame1)
        self.window.add(frame2)
        # ---------------------------
        image = gtk.Image()
        image.set_from_file(LOGO_PATH)
        image.show()
        self.main_hbox.pack_start(image)
        # ---------------------------
        self.window.set_title("Loading PanDDA Inspect")
        self.window.show_all()
        self.window.set_keep_above(True)
        catchup(True)

    def show_menu(self):
        # ---------------------------
        box = gtk.VBox(spacing=7)
        box.set_border_width(7)
        self.main_hbox.pack_start(box)
        # ---------------------------
        box.pack_start(gtk.Label(''), expand=True)
        # ---------------------------
        label = gtk.Label('New to PanDDA Inspect?')
        label.set_justify(gtk.JUSTIFY_CENTER)
        box.pack_start(label, expand=False)
        # ---------------------------
        button = gtk.Button(label='Read this important information regarding map contouring')
        button.child.set_padding(5, 5)
        button.child.props.width_chars = 30
        button.child.set_justify(gtk.JUSTIFY_CENTER)
        button.child.set_line_wrap(True)
        button.connect("clicked", lambda x: Notices())
        box.pack_start(button, expand=False)
        # ---------------------------
        label = gtk.Label('For more information visit')
        label.set_padding(5, 5)
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.set_line_wrap(True)
        box.pack_start(label, expand=False)
        # ---------------------------
        label = gtk.Label('https://pandda.bitbucket.io')
        label.set_padding(5, 5)
        label.set_justify(gtk.JUSTIFY_CENTER)
        label.set_line_wrap(True)
        box.pack_start(label, expand=False)
        # ---------------------------
        box.pack_start(gtk.Label(''), expand=True)
        # ---------------------------
        button = gtk.Button(label='Close window')
        button.child.set_padding(5, 5)
        button.child.set_justify(gtk.JUSTIFY_CENTER)
        button.child.set_line_wrap(True)
        button.connect("clicked", lambda x: self.close())
        box.pack_end(button, expand=False)
        # ---------------------------
        self.window.set_title("Welcome to PanDDA Inspect")
        self.window.show_all()
        self.window.queue_draw()
        catchup(True)

    def close(self):
        self.window.destroy()
        post_coot_windows()


# =========================================================================


### TMP HACK ############################################
# from pandda.inspect.events import (
#     GetNextModelFile,
#     )
###
def rel_symlink(orig, link):
    """Make a relative symlink from link to orig"""
    orig = str(orig)
    link = str(link)
    assert not os.path.exists(link), 'Link already exists: {!s}'.format(link)
    orig = os.path.abspath(orig)
    link = os.path.abspath(link)
    assert not link.endswith('/'), 'LINK CANNOT END WITH /'
    os.symlink(os.path.relpath(orig, start=os.path.dirname(link)), link)
###
class Event(object):

    def __init__(self, info):

        self.index = info.name

        self.dtag = info.name[0]

        self.event_num = (
            int(info.name[1])
            if not np.isnan(info.name[1])
            else None
            )

        self.map_resolution = (
            round(info['analysed_resolution'], 2)
            )

        self.map_uncertainty = (
            round(info['map_uncertainty'], 2)
            )

        self.rwork_rfree = (
            round(info['r_work'], 3), 
            round(info['r_free'], 3),
            )

        self.site_num = (
            int(info['site_num'])
            if not np.isnan(info['site_num'])
            else None
            )

        self.est_bdc = (
            round(info['bdc'], 2)
            )

        self.est_1_bdc = (
            round(1-info['bdc'], 2)
            )

        self.z_peak = (
            round(info['z_peak'], 1)
            if not np.isnan(info['z_peak'])
            else None
            )

        self.z_mean = (
            info['z_mean']
            if not np.isnan(info['z_mean'])
            else None
            )

        self.cluster_size = (
            int(info['cluster_size'])
            if not np.isnan(info['cluster_size'])
            else None
            )

        self.xyz = tuple(
            info[c] 
            for c in ['x','y','z']
            if not np.isnan(info[c])
            )

        if len(self.xyz) != 3: 
            self.xyz = None

        self.added_ligand_names = [
            s for s in info['Ligand Names'].split(',') if s
            ]
class GetEventFiles(object):

    def __init__(self,
        pandda_directory,
        pandda_files_dict, 
        pandda_path_prefix,
        ):

        self.pandda_directory = pandda_directory
        self.pandda_files_dict = pandda_files_dict
        self.pandda_path_prefix = pandda_path_prefix

    def __call__(self,
        event,
        ):

        files_dict = self.pandda_files_dict

        # Indentify main dataset directory 

        dataset_dir = (
            self.pandda_directory / 'processed_datasets' / event.dtag
            )

        assert dataset_dir.exists(), 'Directory does not exist: {}'.format(event.dtag)

        # Identify dataset subdirectories and files

        ligand_dir = (
            dataset_dir / 'ligand_files'
            )

        if not ligand_dir.exists(): 
            ligand_dir.mkdir(parents=True)

        model_dir = (
            dataset_dir / 'modelled_structures'
            )

        if not model_dir.exists(): 
            model_dir.mkdir(parents=True)

        # The most recent model of the protein in the pandda maps

        output_model_link = (
            model_dir / '{dkey}-pandda-model.pdb'.format(dkey=event.dtag)
            )

        # Files output by pandda

        dataset_files = files_dict['dataset_files'][event.dtag]

        import json
        logger.info(json.dumps(dataset_files, indent=2))

        input_model = (
            self.pandda_directory / pl.Path(
                dataset_files['structure']
                ).relative_to(self.pandda_path_prefix)
            )

        input_data = (
            self.pandda_directory / pl.Path(
                dataset_files['data']
                ).relative_to(self.pandda_path_prefix)
            )

        output_data = (
            self.pandda_directory / pl.Path(
                dataset_files['output_data']
                ).relative_to(self.pandda_path_prefix)
            )

        # events are optional
        if ('event_data' in dataset_files) and (event.event_num is not None):
            event_data = (
                self.pandda_directory / pl.Path(
                    dataset_files['event_data'][str(event.event_num)]
                    ).relative_to(self.pandda_path_prefix)
                )
        else:
            event_data = None

        return {
            'dataset_dir' : dataset_dir,
            'ligand_dir' : ligand_dir,
            'output_model_dir' : model_dir,
            'output_model_link' : str(output_model_link),
            'input_model' : str(input_model),
            'input_data' : str(input_data),
            'output_data' : str(output_data),
            'event_data' : (
                str(event_data)
                if (event_data is not None)
                else None
                ),
            # 'ligand_pdbs' : map(str, ligand_pdbs),
            # 'ligand_cifs' : map(str, ligand_cifs),
        }
class GetNextModelFile(object):

    output_prefix = 'fitted-v'

    def __init__(self,
        model_directory,
        output_model_link, 
        ):

        self.model_directory = model_directory
        self.output_model_link = str(output_model_link)

    def __call__(self):

        return self.next_file_and_update()

    def next_file_and_update(self):

        next_path = self.next_file()

        self.update_link(next_path)

        return next_path

    def last_file(self):
        """Get the most recent saved model of this protein"""

        fitted_outputs = list(map(
            str,
            sorted(
                self.model_directory.glob(self.output_prefix+'*')
                )
            ))

        if fitted_outputs:
            logger.info('Current models: \n\t{}'.format('\n\t'.join(fitted_outputs)))
            return fitted_outputs[-1]
        else:
            logger.info('No current models')
            return None

    def next_file(self):

        current = self.last_file()

        logger.info('Most recent saved model: {!s}'.format(current))

        if current:
            last_idx = int(current[-8:-4]) # matches to {:04d} below
        else:
            last_idx = 0 # will be incremented

        new_fitted = self.output_prefix + '{:04d}.pdb'.format(last_idx + 1)

        return str(self.model_directory / new_fitted)

    def update_link(self, path=None):
        
        if path is None: 
            path = self.last_file()

        # Always want to remove the link
        if os.path.islink(self.output_model_link):
            os.remove(self.output_model_link)

        # Check if it's a file just in case...
        if os.path.exists(self.output_model_link):
            os.remove(self.output_model_link)

        # No files
        if path is None: 
            return

        logger.info('Linking {!s} -> {!s}'.format(
            os.path.basename(path), 
            os.path.basename(self.output_model_link),
            ))

        # Create new link the most recent file
        # from giant.paths import rel_symlink
        rel_symlink(
            path, 
            self.output_model_link,
            )

        return self.output_model_link
### TMP HACK ############################################

### TMP HACK ############################################
# from pandda.inspect.trackers import (
#     LigandTracker,
#     )
###
class iTracker(object):

    def __init__(self,
        n_total,
        i_start = 0,
        ):

        self.i_current = i_start
        self.n_total = n_total

    def get(self):

        return self.i_current

    def set(self, i):

        self.i_current = (
            i % self.n_total
            )

        return self.get()

    def next(self):

        self.set(self.i_current + 1)

        return self.get()

    def prev(self):

        self.set(self.i_current - 1)

        return self.get()

    def at_first(self):

        return self.get() == 0

    def at_last(self):

        return self.get() == (self.n_total - 1)
class EventTracker(iTracker):
    pass
class LigandTracker(iTracker):
    pass
class SiteTracker(object):

    def __init__(self,
        site_idxs,
        event_tracker,
        ):

        self.site_idxs = site_idxs
        self.n_total = len(set(site_idxs))
        
        if not np.isnan(self.site_idxs).all():
            assert self.site_idxs.min() == 0
            assert self.site_idxs.max() == (self.n_total - 1)

        self.event_tracker = event_tracker

    def get(self):

        i_event = self.event_tracker.get()

        return self.site_idxs[i_event]

    def next(self):

        i_event = self.event_tracker.get()

        curr_site_idx = self.site_idxs[i_event]

        next_site_idx = (curr_site_idx + 1) % self.n_total

        i_event_new = self.find_first_event(next_site_idx)

        self.event_tracker.set(i_event_new)

        return next_site_idx

    def prev(self):

        i_event = self.event_tracker.get()

        curr_site_idx = self.site_idxs[i_event]

        prev_site_idx = (curr_site_idx - 1) % self.n_total

        i_event_new = self.find_first_event(prev_site_idx)

        self.event_tracker.set(i_event_new)

        return prev_site_idx

    def find_first_event(self, site_idx):

        event_idxs = np.where(self.site_idxs == site_idx)[0]

        return event_idxs[0]
### TMP HACK ############################################


class LoadEventModelsAndMaps(object):

    default_contours = {
        'dataset' : 2.0,
        'ground' : 2.0,
        'zvalues' : 3.0,
        }

    def __init__(self, 
        event_files,
        event_map_contour = 2.0,
        ):

        self.event_files = event_files
        self.event_map_contour = event_map_contour

    def __call__(self):

        molecules = {}

        molecules.update(
            self.load_and_assign_mols_for_models(
                pdb_path = (
                    self.event_files['output_model_link']
                    if os.path.exists(
                        self.event_files['output_model_link']
                        )
                    else self.event_files['input_model']
                    ),
                )
            )

        molecules.update(
            self.load_and_assign_mols_for_maps(
                mtz_path = self.event_files['output_data'],
                )
            )

        self.dataset_map_setup(molecules.get('dataset'))
        self.ground_map_setup(molecules.get('ground'))
        self.zvalues_map_setup(molecules.get('zvalues'))

        if self.event_files['event_data'] is not None:

            molecules.update(
                self.load_event_map(
                    mtz_path = self.event_files['event_data'],
                    )
                )

            self.event_map_setup(molecules.get('event'))

        return molecules

    def load_and_assign_mols_for_models(self, 
        pdb_path,
        ):

        i = read_pdb(pdb_path)

        return {
            'model' : i,
            }

    def load_and_assign_mols_for_maps(self, 
        mtz_path,
        ):
        """

        Ordering of columns in MTZ files: 
            (not required to be present)
            -> FDATASET / PHDATASET
            -> FGROUND / PHGROUND
            -> FZVALUES / PHZVALUES
            -> FEVENT / PHEVENT

        """

        loaded_mols = list(
            auto_read_make_and_draw_maps_from_mtz(
                mtz_path,
                )
            )

        output_mols = {}

        # Order is important for these as it is the same as the order loaded!
        for m_key, m_labels in pandda_column_labels.items(): 

            is_valid = valid_labels(
                mtz_path,
                m_labels['F'],
                m_labels['PHI'],
                '', 0,
                )

            if bool(is_valid) is False:
                continue

            output_mols[m_key] = loaded_mols.pop(0)

        assert len(loaded_mols) == 0

        return output_mols

    def load_event_map(self, 
        mtz_path,
        ):

        loaded_mols = list(
            auto_read_make_and_draw_maps_from_mtz(
                mtz_path,
                )
            )
        
        if len(loaded_mols) > 1: 
            for i in loaded_mols[1:]:
                close_molecule(i)
        #assert len(loaded_mols) == 1

        return {'event' : loaded_mols[0]}

    def dataset_map_setup(self, imol):
        
        if imol is None: 
            return
        
        set_contour_level_absolute(imol, self.default_contours['dataset'])
        set_map_displayed(imol, 0)

    def ground_map_setup(self, imol):
        
        if imol is None: 
            return

        set_contour_level_absolute(imol, self.default_contours['ground'])
        set_map_displayed(imol, 0)

    def zvalues_map_setup(self, imol):
        
        if imol is None: 
            return

        set_contour_level_absolute(imol, self.default_contours['zvalues'])
        set_map_displayed(imol, 1)

    def event_map_setup(self, imol):

        if imol is None: 
            return

        set_contour_level_absolute(imol, self.event_map_contour)
        set_imol_refinement_map(imol)


class EventModelMapHandler(object):

    def __init__(self,
        event_files,
        event_map_contour,
        ): 
        
        self.event_files = event_files

        self.get_model_path = GetNextModelFile(
            model_directory = event_files['output_model_dir'],
            output_model_link = event_files['output_model_link'],
            )

        self.load_models_and_maps = LoadEventModelsAndMaps(
            event_files = self.event_files,
            event_map_contour = event_map_contour,
            )

        self.molecules = {}

    def load_all(self):

        self.close_all()

        self.molecules.update(
            self.load_models_and_maps()
            )

        if self.event_files['event_data'] is None: 
            self.make_custom_event_map(0.0)

        set_main_coot_molecule(
            self.molecules['model']
            )

        return self.molecules

    def update_model_link(self):

        self.get_model_path.update_link()

    def revert_to_last_model(self):

        close_molecule(self.molecules.pop('model'))

        last_file = self.get_model_path.last_file()

        if last_file is None: 
            last_file = self.event_files['input_model']

        imol = read_pdb(last_file)

        set_main_coot_molecule(imol)
        
        self.molecules['model'] = imol

    def revert_to_input_model(self):

        close_molecule(self.molecules.pop('model'))

        last_file = self.event_files['input_model']

        imol = read_pdb(last_file)

        set_main_coot_molecule(imol)
        
        self.molecules['model'] = imol

    def write_current_model(self):

        next_file = self.get_model_path()

        imol = self.molecules['model']

        # coot command
        write_pdb_file(imol, next_file)

    def load_full_dataset_mtz(self):

        for k, im in self.molecules.items():
            if k.startswith('input_data_'):
                close_molecule(im)

        imols = auto_read_make_and_draw_maps(
            self.event_files['input_data'],
            )

        for i, imol in enumerate(imols):
            self.molecules['input_data_'+str(i+1)] = imol

    def load_original_model_not_active(self):

        imol = self.molecules.get('original_model_reference', None)

        if imol is not None: 
            close_molecule(imol)

        imol = read_pdb(
            self.event_files['input_model']
            )

        self.molecules['original_model_reference'] = imol

    def make_custom_event_map(self, bdc_value):

        imol = self.molecules.get('custom_event_map', None)

        if imol is not None: 
            close_molecule(imol)

        imol = difference_map(
            self.molecules['dataset'], # imol1
            self.molecules['ground'], # imol2
            float(bdc_value), # map_scale
            )

        self.molecules['custom_event_map'] = imol

        set_contour_level_absolute(imol, 2.0*(1-bdc_value))
        set_imol_refinement_map(imol)
        try: 
            set_map_is_difference_map(imol, False)
        except: 
            pass
        set_map_colour(imol, 255, 0, 255)
        set_molecule_name(imol, 'custom_event_map_1_bdc_'+str(1.0-bdc_value))
        set_scrollable_map(imol)

        # Turn off event map to avoid confusion
        imol_event = self.molecules.get('event')
        if imol_event is not None: 
            set_map_displayed(imol_event, 0)

    def close_all(self):

        while self.molecules:
            k_mol, imol = self.molecules.popitem()
            try:
                close_molecule(imol)
            except Exception as e: 
                pass


class LigandHandler(object):

    def __init__(self, 
        ligand_directory,
        show_ligands = True,
        starting_occupancy = 1.0,
        ):

        self.ligand_directory = ligand_directory

        self.ligand_pdbs = None
        self.ligand_cifs = None

        self.tracker = None

        self.show_ligands = show_ligands
        self.starting_occupancy = starting_occupancy

        self.update()

        self.molecules = {}

    def update(self):

        cifs = sorted(
            self.ligand_directory.glob('*.cif')
            )

        pdbs = [
            p.with_suffix('.pdb') for p in cifs
            ]

        for p in pdbs: assert p.exists()
        for p in cifs: assert p.exists()


        # Need ordered lists for tracker
        self.ligand_pdbs = list(map(str, pdbs))
        self.ligand_cifs = list(map(str, cifs))

        # Create hashable for convenience
        self.ligand_cif_dict = {
            self.get_cif_name(p) : str(p) for p in cifs
            }

        self.tracker = LigandTracker(
            n_total = len(self.ligand_cifs),
            )

    def go_to(self, i_ligand, show=None):

        return self.load(
            self.tracker.set(i_ligand),
            show = show, # if cycling through, then show
            )

    def next(self):

        if self.tracker.n_total == 0:
            modal_msg('This dataset has no ligands')
            return

        if self.tracker.n_total == 1:
            modal_msg('This dataset has only one ligand')
            return

        return self.load(
            self.tracker.next(),
            show = True, # if cycling through, then show
            )

    def prev(self):

        if self.tracker.n_total == 0:
            modal_msg('This dataset has no ligands')
            return

        if self.tracker.n_total == 1:
            modal_msg('This dataset has only one ligand')
            return

        return self.load(
            self.tracker.prev(),
            show = True, # if cycling through, then show
            )

    def load(self, i_ligand=None, show=None):

        self.close_all()

        # Do this as a quiet return
        if self.tracker.n_total == 0: 
            return

        if i_ligand is None: 
            i_ligand = self.tracker.get()

        if show is None: 
            show = self.show_ligands

        l_dict = read_cif_dictionary(
            self.ligand_cifs[i_ligand],
            )

        imol = handle_read_draw_molecule_and_move_molecule_here(
            self.ligand_pdbs[i_ligand],
            )

        self.molecules['ligand'] = imol
        
        self.apply_settings(imol)
            
        if show is False: 
            set_mol_displayed(imol, 0)
        else: 
            set_mol_displayed(imol, 1)

        return imol

    def apply_settings(self, imol):

        set_molecule_bonds_colour_map_rotation(
            imol, LIG_COLOUR,
            )
        # set_mol_displayed(
        #     imol, 1,
        #     )
        set_b_factor_molecule(
            imol, 20,
            )

        # Set the occupancy of the ligand to 2*(1-bdc)
        all_residue_ids = all_residues(imol)

        #  TODO: What is this mystery bool from coot - I must know.
        if (all_residue_ids):

            for res_set in all_residue_ids: 

                if len(res_set) == 4:
                    mystery_bool, res_chn, res_num, res_ins = res_set
                else: 
                    res_chn, res_num, res_ins = res_set

                set_alt_conf_occ(
                    imol, res_chn, res_num, res_ins, [['', self.starting_occupancy]],
                    )

        return imol

    def get_cif_name(self, cif_file):

        return pl.Path(cif_file).stem

    def get_i_ligand_name(self, i_ligand=None):

        if i_ligand is None: 
            i_ligand = self.tracker.get()

        return self.get_cif_name(
            self.ligand_cifs[i_ligand]
            )

    def get_i_ligand_for_cif(self, cif_name=None, cif_file=None):

        assert not ((cif_name is None) and (cif_file is None))

        if cif_file is None: 
            
            cif_file = self.ligand_cif_dict[cif_name]

        i_ligand = self.ligand_cifs.index(cif_file)

        return i_ligand

    def load_cif(self, i_ligand=None, cif_name=None):

        assert not ((i_ligand is None) and (cif_name is None))

        if i_ligand is not None: 
            cif_file = self.ligand_cifs[i_ligand]
        else: 
            cif_file = self.ligand_cif_dict[cif_name]

        read_cif_dictionary(
            cif_file,
            )

    def get_imol(self):

        imol_ligand = self.molecules.get('ligand', None)

        if imol_ligand is None: 
            modal_msg(msg='No ligand has been loaded')

        return imol_ligand

    def merge_into_other(self, imol_other, imol_ligand=None):

        if imol_ligand is None:
            imol_ligand = self.get_imol()

        if imol_ligand is None: 
            return

        try:

            merge_molecules(
                [imol_ligand], 
                imol_other,
                )

        except Exception as e:

            print(e)

    def move_here(self):

        imol_ligand = self.get_imol()

        if imol_ligand is None:
            return

        try:

            move_molecule_to_screen_centre(
                imol_ligand
                )

        except Exception as e:

            print(e)

    def close_all(self):

        while self.molecules:
            k_mol, imol = self.molecules.popitem()
            close_molecule(imol)


class EventHandler(object):

    def __init__(self, 
        event,
        event_files,
        ):

        self.event = event
        self.event_files = event_files

        self.model_map_handler = EventModelMapHandler(
            event_files = self.event_files,
            event_map_contour = 2.0 * event.est_1_bdc,
            )

        self.ligand_handler = LigandHandler(
            ligand_directory = self.event_files['ligand_dir'],
            show_ligands = (
                not os.path.exists(
                    self.event_files['output_model_link']
                    )
                ),
            starting_occupancy = (2.0 * event.est_1_bdc),
            )

    def load_all(self):

        event = self.event

        if event.xyz is not None:
            set_rotation_centre(*event.xyz)

        self.model_map_handler.load_all()

        # Load cif files for any ligand that were previously added 
        for cif_name in event.added_ligand_names: 
            self.ligand_handler.load_cif(cif_name=cif_name)

        # Ensure current (first) ligand is loaded last
        self.ligand_handler.load()

    def close_all(self):

        self.model_map_handler.close_all()
        self.ligand_handler.close_all()

    def save(self):

        self.model_map_handler.write_current_model()

    def revert_to_last_model(self):

        self.model_map_handler.revert_to_last_model()

    def revert_to_input_model(self):

        self.model_map_handler.revert_to_input_model()

    def merge_ligand_with_model(self):

        imol_model = (
            self.model_map_handler.molecules.get('model')
            )

        if imol_model is None: 
            modal_msg(msg='No model has been loaded')
            return

        self.ligand_handler.merge_into_other(
            imol_other = imol_model,
            )

    def move_ligand_here(self):

        self.ligand_handler.move_here()

    def prev_ligand(self):

        self.ligand_handler.prev()

    def next_ligand(self):

        self.ligand_handler.next()

    def update_and_go_to_ligand(self, info_dict, show=None):

        self.ligand_handler.update()

        return self.go_to_ligand(
            info_dict = info_dict,
            show = show,
            )

    def go_to_ligand(self, info_dict, show=None):

        if info_dict is None: 
            return None

        i_ligand = self.ligand_handler.get_i_ligand_for_cif(
            cif_name = info_dict.get('name'),
            cif_file = info_dict.get('cif'),
            )

        if i_ligand is None: 
            return

        self.ligand_handler.go_to(
            i_ligand = i_ligand,
            show = show,
            )
        

class GetEventHandler(object):

    def __init__(self,
        pandda_directory, 
        pandda_files_dict, 
        pandda_path_prefix,
        ):

        # from pandda.inspect.events import (                
        #     GetEventFiles,
        #     )

        self.get_event_files = GetEventFiles(
            pandda_directory = pandda_directory,
            pandda_files_dict = pandda_files_dict,
            pandda_path_prefix = pandda_path_prefix,
            )

    def __call__(self, 
        event_info,
        ):

        # from pandda.inspect.events import (
        #     Event,
        #     )

        event = Event(info=event_info)

        event_files = self.get_event_files(event)

        return EventHandler(
            event = event,
            event_files = event_files,
            )


# =========================================================================

# from pandda.inspect.trackers import (
#     EventTracker,
#     SiteTracker,
#     )


class PanddaEventListController(object):
    """Responsible for loading events models and maps"""

    def __init__(self,
        event_table,
        get_event_handler,
        ):

        self.event_table = event_table

        logger.info(
            str(event_table)
            )

        self.event_counts = self.get_event_counts(
            self.event_table,
            )

        self.get_event_handler = get_event_handler
        self.event_handler = None
        
        self.event_tracker = EventTracker(
            n_total = len(self.event_table.index),
            )

        self.site_tracker = SiteTracker(
            site_idxs = (self.event_table['site_num'] - 1), # safe hack?!
            event_tracker = self.event_tracker,
            )

        self.refresh()

    def get_event_counts(self, table):

        table = table.reset_index()

        dtags = table['dtag']

        return dtags.value_counts().to_dict()

    def get_event_info(self, i_event):
        return self.event_table.iloc[i_event]

    def load_event(self, i_event):

        self.close_all()

        event_info = self.get_event_info(i_event)

        self.event_handler = self.get_event_handler(
            event_info = event_info,
            )

        self.event_handler.load_all()

    def save(self):
        # Currently here because we might want to add additional save functions? 
        if self.event_handler is not None:
            self.event_handler.save()

    def refresh(self):
        self.load_event(
            self.event_tracker.get()
            )

    def go_to(self, i_event):
        self.load_event(
            self.event_tracker.set(i_event)
            )

    def next_event(self):
        self.load_event(
            self.event_tracker.next()
            )

    def prev_event(self):
        self.load_event(
            self.event_tracker.prev()
            )

    def next_site(self):
        i_site = self.site_tracker.next()
        self.load_event(
            self.event_tracker.get()
            )

    def prev_site(self):
        i_site = self.site_tracker.prev()
        self.load_event(
            self.event_tracker.get()
            )

    def next_event_modelled(self):

        selection = (self.event_table['Ligand Placed'] == True)

        if (selection == False).all():
            modal_msg( 
                msg = 'No modelled events (none marked as "Ligand Placed")',
                )
            return

        i_next = self.get_i_next_from_selection(
            selection = selection,
            )

        if (i_next is None):
            modal_msg(
                msg = 'Current event is the only modelled event (marked as "Ligand Placed")',
                )
            return

        self.go_to(i_next)

    def next_event_unviewed(self):

        selection = (self.event_table['Viewed'] == False)

        if (selection == False).all():
            modal_msg(
                msg = 'No unviewed events',
                )
            return

        i_next = self.get_i_next_from_selection(
            selection = selection,
            )

        if (i_next is None):
            modal_msg(
                msg = 'Current event was the only unviewed event',
                )
            return

        self.go_to(i_next)

    def next_event_for_dataset(self, dataset_id):

        selection = (self.event_table.reset_index()['dtag'] == dataset_id)

        if (selection == False).all():
            modal_msg(
                msg = 'No datasets with this identifier: "{}"'.format(
                    dataset_id,
                    ),
                )
            return

        i_next = self.get_i_next_from_selection(selection) 

        if (i_next is None): 
            modal_msg(
                msg = 'Current event is the only event from this dataset ({})'.format(
                    dataset_id,
                    ),
                )
            return

        self.go_to(i_next)

    def get_i_next_from_selection(self, selection):

        i_current = self.event_tracker.get()
        n_total = self.event_tracker.n_total

        indices = np.where(selection)[0]

        if len(indices) == 0: 
            # No selected

            return None

        elif len(indices) == 1:
            # Special case -- one model

            i_next = indices[0]

            if (i_current == i_next):
                # do not return the current index
                return None

        else: 

            # Modulate the indices by n_total
            indices = sorted(
                ((indices - i_current) % n_total) + i_current
            )

            # We know that there are at least two events so this should be safe!
            i_next = (
                indices[0]
                if 
                indices[0] != i_current
                else 
                indices[1]
                )

        assert i_next != i_current

        return i_next

    def close_all(self):
        for imol in molecule_number_list():
            close_molecule(imol)


# =========================================================================

class DummyUpdateOutput(object):

    def __init__(self,
        inspector,
        mode,
        ):

        self.inspector = inspector
        self.mode = mode

    def __call__(self):

        d = nonmodal_msg(
            'Note: there is no HTML output when running in the "{mode}" mode'.format(
                mode = self.mode,
                )
            )

        catchup(True) 

        time.sleep(2)

        d.destroy()

### TMP HACK ############################################
try:
    import matplotlib as mpl
    mpl.use('agg')
    mpl.interactive(False)
    from matplotlib import pyplot as plt
    plt.switch_backend('agg') # yes I know this done twice -- for safety!
    plt.interactive(0)
except Exception as e:
    logger.info(str(e))
    import matplotlib as mpl
    from matplotlib import pyplot as plt
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
            fig = fig,
            filename = filename, 
            )

        return {self.output_key : filename}

    def setup(self, nrows=1, ncols=1, **kwargs):

        fig, axes = plt.subplots(
            nrows=nrows, ncols=ncols,
            **kwargs
            )

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
class EventSitePlotter(PanddaPlotter):

    output_key = "site_events"

    def __init__(self, output_path_template):

        self.output_path_template = str(output_path_template)

        assert ('{site_num' in output_path_template)

    def __call__(self, event_dicts, *args, **kwargs):

        event_data = self.unpack_events(event_dicts)

        output_files = collections.OrderedDict()

        for site_num, e_data in sorted(event_data.items()):

            fig = self.plot(
                site_num = site_num, 
                event_data = e_data, 
                *args, **kwargs
                )

            filename = self.get_path(
                site_num = site_num,
                )

            self.save(
                fig = fig, 
                filename = filename,
                )

            output_files[site_num] = filename

        return {self.output_key : output_files}

    def plot(self,
        site_num,
        event_data,
        *args, **kwargs
        ):

        ###

        fig, axis = self.setup()

        axis.set_title(
            'Events for Site {site_num}'.format(
                site_num = site_num,
                )
            )
        axis.set_xlabel('Event #')
        axis.set_ylabel('Z-Peak Value')

        self.make_bar_plot(
            axis = axis,
            bar_values = self.get_plot_values(event_data),
            bar_colours = self.get_colour_values(event_data),
            )

        return fig

    def make_bar_plot(self, 
        axis, 
        bar_values,
        bar_colours, 
        min_x = 5,
        ):
        """Plot set of bar graphs in one figure"""

        n = len(bar_values)

        assert n > 0

        axis.bar(
            x = np.arange(n) + 1.0,
            height = bar_values, 
            width = 0.8, 
            color = bar_colours,
            )

        axis.set_xticks(
            list(range(1, n+1))
            )

        axis.set_yticks(
            list(range(0, int(max(bar_values)+0.5)))
            )

        axis.set_xlim(
            [0.5, max(min_x, n) + 2]
            )

    def get_path(self, **kwargs):

        p = pl.Path(
            self.output_path_template.format(
                **kwargs
                )
            )
        
        if not p.parent.exists():
            p.parent.mkdir(parents=True)

        return str(p)

    def unpack_events(self, events):

        site_nums = sorted(set(
            [e['site_num'] for e in events]
            ))

        site_events = {
            i : [
                e for e in events 
                if (e['site_num'] == i)
                ]
            for i in site_nums
            }

        return site_events
        
    def get_plot_values(self, events):

        return [e['z_peak'] for e in events]

    def get_colour_values(self, events):

        return [e.get('colour','slategray') for e in events]
### TMP HACK ############################################

class UpdateOutput(object):

    def __init__(self,
        inspector,
        write_html,
        ):

        self.inspector = inspector
        self.write_html = write_html

    def __call__(self):

        d = nonmodal_msg('Updating output...')

        catchup(True)

        of = self.make_output_graphs()

        # self.write_html( 
        #     inspector = self.inspector,
        #     output_files = of,
        #     )

        # Destroy update message
        d.destroy()

    def make_output_graphs(self):

        event_dicts = self.inspector.tables.events.table.to_dict('records')

        if len(event_dicts) == 0:
            return dict()

        for e in event_dicts: 
            e['colour'] = (
                'limegreen' if e['Ligand Placed']
                else 'red' if e['Viewed']
                else 'blue'
                )

        try:
            
            # from pandda.analyse.output.graphs import EventSitePlotter
            plot = EventSitePlotter(
                'inspect/graphs/analyse_events_site_{site_num:03d}.png',
                )

            return plot(event_dicts=event_dicts)

        except Exception as e: 

            return dict()


# =========================================================================

### TMP HACK ############################################
import distutils
def is_available(program):
    """Return whether the program is available in the path"""
    from distutils.spawn import find_executable
    if find_executable(program):
        return True
    return False
import procrunner
class Dispatcher(object):

    def __init__(self, program):

        if not is_available(program):
            raise ValueError("Can't find the program '{!s}'. Is it available?".format(program))

        self.program = program
        self.command = [program]
        self.stdin = None
        self.working_directory = None
        self.timeout = None
        self.debug = False

        self.result = None

    def __str__(self):
        return (
            'Program:\n\t' + self.command[0] + '\n' +
            'Command line arguments:\n\t' + '\n\t'.join(self.command[1:] if self.command[1:] else ['None']) + '\n' +
            'Standard Input:\n\t' + (self.stdin if self.stdin is not None else 'None\n').replace('\n','\n\t')
        )

    def append_arg(self, arg):
        self.command.append(arg)

    def extend_args(self, args):
        self.command.extend(args)

    def append_stdin(self, line):
        if self.stdin is None:
            self.stdin = ''
        self.stdin += line.strip('\n') + '\n'

    def extend_stdin(self, lines):
        for line in lines:
            self.append_stdin(line)

    def as_string(self):
        return str(self)

    def as_command(self):
        out_str = ' '.join(self.command)
        if self.stdin:
            out_str += (' <<eof\n' + self.stdin + '\neof')
        return out_str

    def run(self):

        assert self.result is None, 'Program has already been run.'

        self.result = procrunner.run(
            command = self.command,
            timeout = self.timeout,
            debug = self.debug,
            stdin = (
                bytes(self.stdin, 'utf-8')
                if self.stdin is not None
                else None
                ),
            print_stdout = False,
            print_stderr = False,
            #callback_stdout = None,
            #callback_stderr = None,
            #environment = None,
            #environment_override = None,
            #win32resolve = True,
            working_directory = self.working_directory,
        )

        return self.result

    def write_output(self, log_file):

        separator = '\n----------------------\n'

        with open(log_file, 'a') as fh:

            fh.write(separator)
            fh.write('Command information:\n')
            fh.write(self.as_string()+'\n')
            fh.write(separator)
            fh.write('Command to re-run:\n')
            fh.write(self.as_command()+'\n')
            fh.write(separator)
            fh.write('Program STDOUT:\n')
            fh.write(str(self.result.stdout))
            fh.write(separator)
            fh.write('Program STDERR:\n')
            fh.write(str(self.result.stderr))
            fh.write(separator)
def generate_restraints(smiles, name='LIG', prefix='ligand', verbose=False):
    """Generate pdb and cif files from smiles string"""

    assert len(name) == 3

    # Common fixes
    smiles = smiles.replace('CL', 'Cl')
    smiles = smiles.replace('BR', 'Br')

    out_pdb = prefix+'.pdb'
    out_cif = prefix+'.cif'
    out_log = prefix+'-acedrg.log'

    if os.path.exists(out_pdb):
        raise IOError('Output PDB file already exists: {}'.format(out_pdb))
    if os.path.exists(out_cif):
        raise IOError('Output CIF file already exists: {}'.format(out_cif))

    # Run acedrg
    acedrg = Dispatcher('acedrg')

    acedrg.extend_args([
        '--smi={}'.format(smiles),
        '-r', name,
        '-o', prefix,
    ])

    if (verbose is True):
        logger.info(acedrg.as_string())

    acedrg.run()
    acedrg.write_output(out_log)

    if not (os.path.exists(out_pdb) and os.path.exists(out_cif)):
        logger.info(str(acedrg.result.stdout))
        logger.info(str(acedrg.result.stderr))
        raise Exception('acedrg failed during ligand generation')

    return out_pdb, out_cif
import string, shutil
class MakeNewLigand(object):

    def __init__(self,
        output_directory,
        ):

        self.output_directory = pl.Path(
            output_directory
            )

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        self.tmp_directory = (
            self.output_directory / 'make_ligand_tmp'
            )

        self.disallowed = set(string.punctuation)
        self.disallowed.add(' ')

    def __call__(self,
        ligand_id,
        ligand_name,
        ligand_smiles,
        ):

        self.validate(
            ligand_id = ligand_id,
            ligand_name = ligand_name,
            ligand_smiles = ligand_smiles,
            )

        real_out_pdb = (
            self.output_directory / (ligand_name + '.pdb')
            )

        real_out_cif = real_out_pdb.with_suffix('.cif')

        if real_out_pdb.exists() or real_out_cif.exists():
            raise ValueError(
                'Output ligand files already exist -- please choose a different id/name'
                )

        # Remove temporary directory if it exists
        self.cleanup()

        # This should really always evaluate to True...
        if not self.tmp_directory.exists():
            self.tmp_directory.mkdir(parents=True)

        try:

            tmp_pdb, tmp_cif = generate_restraints(
                smiles = ligand_smiles,
                name = ligand_id, 
                prefix = str(
                    self.tmp_directory / ligand_name
                    ),
                )

        except Exception as e: 
            raise Exception(
                'Error during ligand generation. See terminal output for acedrg error message.\n{}'.format(str(e))
                )

        self.move_output_files(
            ligand_pdb_tmp = tmp_pdb,
            ligand_pdb = real_out_pdb,
            ligand_cif_tmp = tmp_cif,
            ligand_cif = real_out_cif,
            )

        # Only do this upon successful completion
        self.cleanup()

        return {
            'pdb' : str(real_out_pdb),
            'cif' : str(real_out_cif),
            'id' : ligand_id,
            'name' : ligand_name,
            'smiles' : ligand_smiles,
        }

    def validate(self, 
        ligand_id,
        ligand_name,
        ligand_smiles,
        ):

        allowed_name = ['-','_']
        disallowed_name = self.disallowed.difference(allowed_name)

        allowed_id = ['_']
        disallowed_id = self.disallowed.difference(allowed_id)

        if disallowed_name.intersection(ligand_name):

            raise ValueError(
                'Ligand name cannot contain spaces or punctuation except for {}'.format(
                    ' or '.join(allowed_path),
                    )
                )

        if disallowed_id.intersection(ligand_id):
            
            raise ValueError(
                'Ligand ID cannot contain spaces or punctuation except for {}'.format(
                    ' or '.join(allowed_id),
                    )
                )

        if len(ligand_id) == 0:

            raise ValueError(
                'No ligand id provided'
                )

        if len(ligand_id) != 3:

            raise ValueError(
                'Ligand ID must be three characters'
                )

        if len(ligand_name) == 0:

            raise ValueError(
                'No ligand name provided'
                )
            
        if len(ligand_smiles) == 0:

            raise ValueError(
                'No ligand smiles provided'
                )
            
        return

    def move_output_files(self,
        ligand_pdb_tmp,
        ligand_pdb,
        ligand_cif_tmp,
        ligand_cif,
        ):

        ligand_pdb_tmp = pl.Path(ligand_pdb_tmp)
        ligand_pdb_tmp.rename(ligand_pdb)

        ligand_cif_tmp = pl.Path(ligand_cif_tmp)
        ligand_cif_tmp.rename(ligand_cif)

        assert ligand_pdb.exists()
        assert ligand_cif.exists()

    def cleanup(self):

        if self.tmp_directory.exists():
            shutil.rmtree(str(self.tmp_directory))
### TMP HACK ############################################


class MakeNewLigandModal(object):

    def __init__(self, 
        output_directory,
        ):

        # from pandda.inspect.ligands import (
        #     MakeNewLigand,
        #     )

        self.make_ligand = MakeNewLigand(
            output_directory = output_directory,
            )

    def __call__(self):

        dialog = self.make_window()

        objects = self.populate(
            vbox = dialog.vbox,
            )

        dialog.show_all()

        ligand_dict = self.handle(
            dialog = dialog, 
            objects = objects,
            )

        dialog.destroy()

        return ligand_dict

    def make_window(self):

        dialog = gtk.Dialog(
            "Create New Ligand",
            None,
            gtk.DIALOG_MODAL,
            (gtk.STOCK_CANCEL, gtk.RESPONSE_DELETE_EVENT, gtk.STOCK_OK, gtk.RESPONSE_ACCEPT),
            )

        return dialog

    def populate(self, vbox):

        # ID for the ligand
        id_hbox = gtk.HBox(homogeneous=False, spacing=5)
        vbox.pack_start(id_hbox)
        label = gtk.Label('3-letter code')
        label.props.width_chars = 20
        label.set_justify(gtk.JUSTIFY_RIGHT)
        id_hbox.pack_start(label)
        id_entry = gtk.Entry(max=3)
        id_entry.set_text('UNL')
        id_hbox.pack_start(id_entry)

        # SMILES for the ligand
        smiles_hbox = gtk.HBox(homogeneous=False, spacing=5)
        vbox.pack_start(smiles_hbox)
        label = gtk.Label('Smiles string')
        label.props.width_chars = 20
        smiles_hbox.pack_start(label)
        smiles_entry = gtk.Entry(max=300)
        smiles_entry.set_text('')
        smiles_hbox.pack_start(smiles_entry)

        # Name of the ligand
        name_hbox = gtk.HBox(homogeneous=False, spacing=5)
        vbox.pack_start(name_hbox)
        label = gtk.Label('Ligand Name (optional)')
        label.props.width_chars = 20
        label.set_justify(gtk.JUSTIFY_RIGHT)
        name_hbox.pack_start(label, expand=True)
        name_entry = gtk.Entry(max=100)
        name_entry.set_text('')
        name_hbox.pack_start(name_entry)

        status_hbox = gtk.HBox(homogeneous=False, spacing=5)
        vbox.pack_start(status_hbox)
        status_label = gtk.Label('Click OK to run acedrg')
        status_label.props.width_chars = 30
        status_label.set_line_wrap(True)
        status_hbox.pack_start(status_label)

        return {
            'id' : id_entry,
            'name' : name_entry,
            'smiles' : smiles_entry,
            'status' : status_label,
        }

    def handle(self, dialog, objects):

        success = False

        while success is False:

            ligand_pdb = ligand_cif = None

            response = dialog.run()

            # Delete window/cancel?
            if response in [int(gtk.RESPONSE_REJECT), int(gtk.RESPONSE_DELETE_EVENT)]:
                return None

            assert response is int(gtk.RESPONSE_ACCEPT), (
                'invalid response received ({} should be {})'.format(response, int(gtk.RESPONSE_ACCEPT))
                )

            self.set_status(objects, 'running... may take a moment!')

            ligand_args = self.get_args(objects)

            try: 
                ligand_info = self.make_ligand(**ligand_args)
            except Exception as e: 
                self.set_status(
                    objects, 'Error: change parameters and then click OK to re-run',
                    )
                modal_msg(str(e))
                continue

            ligand_pdb = ligand_info['pdb']
            ligand_cif = ligand_info['cif']

            if ligand_pdb.exists() and ligand_cif.exists():
                self.set_status(objects, '')
                modal_msg('Ligand generated successfully!')
                break
            
            self.set_status(
                objects, 'No files generated but no error was raised - please contact developer',
                )

        return ligand_info

    def get_args(self, objects):

        l_id = objects['id'].get_text().strip(' ')
        l_name = objects['name'].get_text().strip(' ')
        l_smiles = objects['smiles'].get_text().strip(' ')

        if len(l_name) == 0:
            l_name = l_id
            
        return {
            'ligand_id' : l_id,
            'ligand_name' : l_name,
            'ligand_smiles' : l_smiles,
        }

    def set_status(self, objects, message):

        objects['status'].set_label(message)
        catchup(True)


class GetMakeCustomEventMapWindow(object):

    def __init__(self, event_handler):

        self.event_handler = event_handler

    def __call__(self):

        window = self.make_window()

        objects = self.populate(
            window = window,
            starting_value = self.event_handler.event.est_1_bdc,
            )

        self.connect(
            objects = objects,
            )

        window.show_all()

        return window

    def make_window(self):

        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_position(gtk.WIN_POS_CENTER)
        window.set_title("Make custom event map")
        window.set_border_width(10)
        
        # dialog = gtk.Dialog(
        #     "Make custom event map",
        #     None,
        #     gtk.DIALOG_MODAL,
        #     (gtk.STOCK_OK, gtk.RESPONSE_ACCEPT),
        #     )

        return window

    def populate(self, window, starting_value):

        vbox = gtk.VBox(homogeneous=False, spacing=5)
        window.add(vbox)

        hbox = gtk.HBox(homogeneous=False, spacing=5)
        vbox.pack_start(hbox, padding=5)

        label = gtk.Label('New Occupancy (1-BDC):')
        label.props.width_chars = 20
        label.set_justify(gtk.JUSTIFY_RIGHT)
        hbox.pack_start(label)

        bdc_entry = gtk.Entry()
        bdc_entry.set_text(
            str(starting_value)
            )
        hbox.pack_start(bdc_entry)
        
        ##

        hbox = gtk.HBox(homogeneous=False, spacing=5)
        vbox.pack_start(hbox, padding=5)
        
        bdc_adjuster = gtk.Adjustment(
            lower = 0.0,
            upper = 1.0,
            value = float(starting_value),
            page_size = 0.01,
            step_incr = 0.01,
            page_incr = 0.1,
            )

        bdc_slider = gtk.HScrollbar(
            adjustment = bdc_adjuster,
            )
        bdc_slider.set_update_policy(gtk.UPDATE_DISCONTINUOUS)
        bdc_slider.set_round_digits(2)
        hbox.pack_start(bdc_slider)

        #bdc_slider.set_value_pos(gtk.POS_RIGHT)
        #bdc_slider.set_draw_value(True)

        return {
            'entry' : bdc_entry,
            'adjustment' : bdc_adjuster,
            'slider' : bdc_slider,
        }

    def connect(self, objects):
        
        objects['entry'].connect(
            "activate", 
            lambda x: [
                objects['adjustment'].set_value(
                    float(
                        objects['entry'].get_text()
                        )
                    ),
                ],
            )

        objects['adjustment'].connect(
            "value_changed",
            lambda x: [
                objects['entry'].set_text(
                    str(objects['adjustment'].get_value()),
                    ),
                self.event_handler.model_map_handler.make_custom_event_map(
                    bdc_value = (
                        1.0 - float(
                            objects['adjustment'].get_value()
                            )
                        ),
                    ),
                ],
            )


# =========================================================================


class GuiPart(object):

    def __init__(self):

        self.labels = {}
        self.buttons = {}
        self.objects = {}


class NavigationButtons1(GuiPart):

    def __call__(self):

        box = gtk.HBox(homogeneous=False, spacing=2)
        box.set_border_width(3)

        b = self.buttons.setdefault(
            'prev',
            gtk.Button(label="<<< Prev <<<\n(Don't Save Model)"),
            )
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box.pack_start(b)

        b = self.buttons.setdefault(
            'skip',
            gtk.Button(label=">>> Next >>>\n(Don't Save Model)"),
            )
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box.pack_start(b)

        b = gtk.VSeparator()
        box.pack_start(b, expand=False, padding=5)

        b = self.buttons.setdefault(
            'next',
            gtk.Button(label=">>> Next >>>\n(Save Model)")
            )
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box.pack_start(b)

        return box


class NavigationButtons2(GuiPart):

    def __call__(self):

        main_box = gtk.VBox(homogeneous=True, spacing=5)

        box1 = gtk.HBox(homogeneous=False, spacing=2)
        box1.set_border_width(3)
        frame = gtk.Frame();
        frame.add(box1)
        main_box.pack_start(frame)

        box2 = gtk.HBox(homogeneous=True, spacing=2)
        box2.set_border_width(3)
        frame = gtk.Frame();
        frame.add(box2)
        main_box.pack_start(frame)

        box3 = gtk.HBox(homogeneous=True, spacing=2)
        box3.set_border_width(3)
        frame = gtk.Frame();
        frame.add(box3)
        main_box.pack_start(frame)

        l = gtk.Label('Go to Dataset:')
        box1.pack_start(l, expand=False, fill=False, padding=5)

        e = self.objects.setdefault(
            'go-to-text',
            gtk.Entry(max=200),
            )
        box1.pack_start(e, expand=True, fill=True, padding=5)

        b = self.buttons.setdefault(
            'go-to',
            gtk.Button(label="Go"),
            )
        box1.pack_start(b, expand=False, fill=False, padding=5)

        b = self.buttons.setdefault(
            'prev-site',
            gtk.Button(label="<<< Go to Prev Site <<<"),
            )        
        box2.add(b)

        b = self.buttons.setdefault(
            'next-site',
            gtk.Button(label=">>> Go to Next Site >>>"),
            )
        box2.add(b)

        b = self.buttons.setdefault(
            'next-unviewed',
            gtk.Button(label=">>> Go to Next Unviewed >>>"),
            )
        box3.add(b)

        b = self.buttons.setdefault(
            'next-modelled',
            gtk.Button(label=">>> Go to Next Modelled >>>"),
            )
        box3.add(b)

        return main_box


class LigandButtons(GuiPart):

    def __call__(self):

        box1 = gtk.VBox(homogeneous=True, spacing=2)
        box1.set_border_width(3)

        b = self.buttons.setdefault(
            'merge-ligand',
            gtk.Button(label="Merge Ligand\nWith Model"),
            )
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box1.add(b)

        b = self.buttons.setdefault(
            'move-ligand',
            gtk.Button(label="Move New Ligand Here"),
            )
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 10
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box1.add(b)

        b = self.buttons.setdefault(
            'next-ligand',
            gtk.Button(label="Open Next Ligand"),
            )
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box1.add(b)

        box2 = gtk.VBox(homogeneous=True, spacing=2)
        box2.set_border_width(3)

        b = self.buttons.setdefault(
            'save',
            gtk.Button(label="Save Model"),
            )
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 10
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box2.add(b)

        b = self.buttons.setdefault(
            'reload',
            gtk.Button(label="Reload Last Saved Model"),
            )
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box2.add(b)

        b = self.buttons.setdefault(
            'reset',
            gtk.Button(label="Reset to Unfitted Model"),
            )
        b.child.set_line_wrap(True)
        b.child.props.width_chars = 15
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        box2.add(b)

        hbox_main = gtk.HBox(spacing=5)
        hbox_main.pack_start(box1)
        hbox_main.pack_start(gtk.VSeparator(), expand=False, padding=5)
        hbox_main.pack_start(box2)

        return hbox_main


class EventRecordButtons(GuiPart):

    def __call__(self):

        hbox_1 = gtk.HBox(homogeneous=False, spacing=5)

        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)

        l = gtk.Label('Event Comment:')
        hbox_1.pack_start(l, expand=False, fill=False, padding=5)

        e = self.objects.setdefault(
            'event comment text',
            gtk.Entry(max=200),
            )
        hbox_1.pack_start(e, expand=True, fill=True, padding=5)

        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        
        hbox_2 = gtk.HBox(homogeneous=True, spacing=5)
        hbox_2.set_border_width(3)
        vbox_1_1 = gtk.VBox(homogeneous=True, spacing=2)
        vbox_1_2 = gtk.VBox(homogeneous=True, spacing=2)
        vbox_1_3 = gtk.VBox(homogeneous=True, spacing=2)
        hbox_2.add(vbox_1_1);
        hbox_2.add(vbox_1_2);
        hbox_2.add(vbox_1_3)
        
        b = self.buttons.setdefault(
            'tp',
            gtk.RadioButton(label="Mark Event as Interesting"),
            )
        vbox_1_1.add(b)

        b = self.buttons.setdefault(
            'fp',
            gtk.RadioButton(label="Mark Event as Not Interesting", group=b),
            )
        vbox_1_1.add(b)

        b = self.buttons.setdefault(
            'placed',
            gtk.RadioButton(label="Ligand Placed"),
            )
        vbox_1_2.add(b)

        b = self.buttons.setdefault(
            'not placed',
            gtk.RadioButton(label="No Ligand Placed", group=b),
            )
        vbox_1_2.add(b)

        b = self.buttons.setdefault(
            'highconf',
            gtk.RadioButton(label="Model: High Confidence"),
            )
        vbox_1_3.add(b)

        b = self.buttons.setdefault(
            'medconf',
            gtk.RadioButton(label="Model: Medium Confidence", group=b),
            )
        vbox_1_3.add(b)

        b = self.buttons.setdefault(
            'lowconf',
            gtk.RadioButton(label="Model: Low Confidence", group=b),
            )
        vbox_1_3.add(b)

        vbox_main = gtk.VBox(spacing=0)
        vbox_main.pack_start(
            gtk.Label('Record Event Information (this event only)'), 
            expand=False, 
            fill=False, 
            padding=5,
            )
        vbox_main.pack_start(hbox_1, padding=0)
        vbox_main.pack_start(hbox_2, padding=5)

        return vbox_main


class SiteRecordButtons(GuiPart):
    
    def __call__(self):

        hbox_1 = gtk.HBox(homogeneous=False, spacing=5)
        
        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        
        l = gtk.Label('Name:')
        l.set_width_chars(10)
        hbox_1.pack_start(l, expand=False, fill=False, padding=5)
        
        e = self.objects.setdefault(
            'site name text',
            gtk.Entry(max=200),
            )
        hbox_1.pack_start(e, expand=True, fill=True, padding=5)
        
        hbox_1.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        
        hbox_2 = gtk.HBox(homogeneous=False, spacing=5)
        
        hbox_2.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        
        l = gtk.Label('Comment:')
        l.set_width_chars(10)
        hbox_2.pack_start(l, expand=False, fill=False, padding=5)
        
        e = self.objects.setdefault(
            'site comment text',
            gtk.Entry(max=200),
            )
        hbox_2.pack_start(e, expand=True, fill=True, padding=5)
        
        hbox_2.pack_start(gtk.HBox(), expand=False, fill=False, padding=10)
        
        vbox_main = gtk.VBox(spacing=0)
        vbox_main.pack_start(
            gtk.Label('Record Site Information (for all events with this site)'), 
            expand=False,
            fill=False, 
            padding=5,
            )
        vbox_main.pack_start(hbox_1, padding=5)
        vbox_main.pack_start(hbox_2, padding=5)

        return vbox_main


class MiscButtons(GuiPart):

    def __call__(self):
        
        hbox_1 = gtk.HBox(homogeneous=False, spacing=5)
        
        hbox_1.pack_start(gtk.HBox(), expand=True, fill=False, padding=10)
        
        l = gtk.Label('Miscellaneous buttons')
        hbox_1.pack_start(l, expand=False, fill=False, padding=5)
        
        b = self.buttons.setdefault(
            'load-full-dataset-mtz',
            gtk.Button('Load input mtz file'),
            )
        hbox_1.pack_start(b, expand=False, fill=False, padding=5)
        
        # b = self.buttons.setdefault(
        #     'load-ground-state-map',
        #     gtk.Button('Load average map'),
        #     )
        # hbox_1.pack_start(b, expand=False, fill=False, padding=5)
        
        b = self.buttons.setdefault(
            'load-original-model',
            gtk.Button('Load unfitted model\n(for comparison only)'),
            )
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        hbox_1.pack_start(b, expand=False, fill=False, padding=5)

        b = self.buttons.setdefault(
            'create-ligand',
            gtk.Button('Create New Ligand\n(from smiles string)'),
            )
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        hbox_1.pack_start(b, expand=False, fill=False, padding=5)

        b = self.buttons.setdefault(
            'custom-event-map',
            gtk.Button('Make custom\nevent map'),
            )
        b.child.set_justify(gtk.JUSTIFY_CENTER)
        hbox_1.pack_start(b, expand=False, fill=False, padding=5)
        
        hbox_1.pack_start(gtk.HBox(), expand=True, fill=False, padding=10)
        
        vbox_main = gtk.VBox(spacing=0)
        vbox_main.pack_start(hbox_1, padding=5)

        return vbox_main


class QuitButtons(GuiPart):
    
    def __call__(self):
        
        vbox_1 = gtk.VBox(spacing=5)
        vbox_1.set_border_width(3)
        
        b = self.buttons.setdefault(
            'quit',
            gtk.Button(label="Quit"),
            )
        vbox_1.pack_start(b)
        
        vbox_1.pack_start(gtk.HSeparator(), expand=False, padding=2)
        
        # b = self.buttons.setdefault(
        #     'summary',
        #     gtk.Button(label="Summary"),
        #     )
        # vbox_1.pack_start(b)
        
        # vbox_1.pack_start(gtk.HSeparator(), expand=False, padding=2)
        
        # b = self.buttons.setdefault(
        #     'updatehtml',
        #     gtk.Button(label="Update HTML"),
        #     )
        # vbox_1.pack_start(b)
        
        hbox_main = gtk.HBox(spacing=5)
        hbox_main.pack_start(vbox_1)

        return hbox_main


class EventInfoTable(GuiPart):

    def __call__(self):
        
        # Main box
        vbox_main = gtk.VBox()
        vbox_main.set_border_width(3)

        #  Pack sub-boxes
        hbox_sub_1 = gtk.HBox()
        hbox_sub_2 = gtk.HBox()
        vbox_main.pack_start(hbox_sub_1)
        vbox_main.pack_start(hbox_sub_2)

        ##############
        # HBOX SUB 1 #
        ##############

        # Dataset Name
        gtk_label = gtk.Label('Dataset ID')
        gtk_value = self.labels.setdefault(
            'dtag',
            gtk.Label('None'),
            )
        gtk_value.set_use_markup(gtk.TRUE)
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to first box
        hbox_sub_1.pack_start(frame)

        ##############
        # HBOX SUB 2 #
        ##############

        vbox_1 = gtk.VBox()
        vbox_2 = gtk.VBox()
        hbox_sub_2.pack_start(vbox_1)
        hbox_sub_2.pack_start(vbox_2)

        ##########
        # VBOX 1 #
        ##########

        # Create title
        title = gtk.Label('Event Information:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame();
        frame.add(title)
        # Add to first column
        vbox_1.pack_start(frame)

        # Event Number for Dataset
        gtk_label = gtk.Label('Event #')
        gtk_value = self.labels.setdefault(
            'e_num',
            gtk.Label('-'),
            )
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box);
        hbox.set_border_width(3)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)

        # Estimated Event Background Correction
        gtk_label = gtk.Label('1-BDC ("occ")')
        gtk_value = self.labels.setdefault(
            'e_1_bdc',
            gtk.Label('-'),
            )
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box);
        hbox.set_border_width(3)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)
        # Store label to allow editing

        # Z-Peak for Dataset
        gtk_label = gtk.Label('Z-blob Peak')
        gtk_value = self.labels.setdefault(
            'zpeak',
            gtk.Label('-'),
            )
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box);
        hbox.set_border_width(3)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to first column
        vbox_1.pack_start(frame)

        ##########
        # VBOX 2 #
        ##########

        # Create title
        title = gtk.Label('Dataset Information:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame();
        frame.add(title)
        # Add to second column
        vbox_2.pack_start(frame)

        # Resolution for Dataset
        gtk_label = gtk.Label('Resolution')
        gtk_value = self.labels.setdefault(
            'map_res',
            gtk.Label('-'),
            )
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box);
        hbox.set_border_width(3)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)

        # Map Uncertainty for Dataset
        gtk_label = gtk.Label('Map Uncertainty')
        gtk_value = self.labels.setdefault(
            'map_unc',
            gtk.Label('-'),
            )
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box);
        hbox.set_border_width(3)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)

        # R-Free/R-Work for Dataset
        gtk_label = gtk.Label('R-Work / R-Free')
        gtk_value = self.labels.setdefault(
            'rwork_rfree',
            gtk.Label('-'),
            )
        gtk_box = gtk.EventBox();
        gtk_box.add(gtk_value)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label);
        hbox.add(gtk_box);
        hbox.set_border_width(3)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to second column
        vbox_2.pack_start(frame)

        # # Currently Blank
        # gtk_label = gtk.Label('-')
        # gtk_value = self.labels.setdefault(
        #     'blank',
        #     gtk.Label('-'),
        #     )
        # gtk_box = gtk.EventBox();
        # gtk_box.add(gtk_value)
        # hbox = gtk.HBox(homogeneous=True);
        # hbox.add(gtk_label);
        # hbox.add(gtk_box);
        # hbox.set_border_width(3)
        # frame = gtk.Frame();
        # frame.add(hbox)
        # # Add to second column
        # vbox_2.pack_start(frame)

        return vbox_main


class ProgressTable(GuiPart):

    def __call__(self):

        # First Column
        vbox_main = gtk.VBox(spacing=5)
        vbox_main.set_border_width(3)

        # Create title
        title = gtk.Label('Overall Inspection Event/Site Progress:')
        title.set_justify(gtk.JUSTIFY_LEFT)
        frame = gtk.Frame();
        frame.add(title)
        # Add to first column
        vbox_main.pack_start(frame)

        # Event Number Name
        gtk_label_1 = gtk.Label('Event')
        gtk_value_1 = self.labels.setdefault(
            'rank_val',
            gtk.Label(0),
            )
        gtk_label_2 = gtk.Label('of')
        gtk_value_2 = self.labels.setdefault(
            'rank_tot',
            gtk.Label(0),
            )
        # Add values to boxes
        gtk_box_1 = gtk.EventBox();
        gtk_box_1.add(gtk_value_1)
        gtk_box_2 = gtk.EventBox();
        gtk_box_2.add(gtk_value_2)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label_1);
        hbox.add(gtk_box_1);
        hbox.add(gtk_label_2);
        hbox.add(gtk_box_2)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to second column
        vbox_main.pack_start(frame)

        # Event Number Name
        gtk_label_1 = gtk.Label('Site')
        gtk_value_1 = self.labels.setdefault(
            'site_val',
            gtk.Label(0),
            )
        gtk_label_2 = gtk.Label('of')
        gtk_value_2 = self.labels.setdefault(
            'site_tot',
            gtk.Label(0),
            )
        # Add values to boxes
        gtk_box_1 = gtk.EventBox();
        gtk_box_1.add(gtk_value_1)
        gtk_box_2 = gtk.EventBox();
        gtk_box_2.add(gtk_value_2)
        hbox = gtk.HBox(homogeneous=True);
        hbox.add(gtk_label_1);
        hbox.add(gtk_box_1);
        hbox.add(gtk_label_2);
        hbox.add(gtk_box_2)
        frame = gtk.Frame();
        frame.add(hbox)
        # Add to second column
        vbox_main.pack_start(frame)

        return vbox_main


class PanddaGUI(object):

    def __init__(self):

        self.navigation_buttons_1 = NavigationButtons1()
        self.navigation_buttons_2 = NavigationButtons2()
        self.ligand_buttons = LigandButtons()
        self.event_record_buttons = EventRecordButtons()
        self.site_record_buttons = SiteRecordButtons()
        self.misc_buttons = MiscButtons()
        self.quit_buttons = QuitButtons()
        self.event_info_table = EventInfoTable()
        self.progress_table = ProgressTable()

        self.widgets = [
            self.navigation_buttons_1,
            self.navigation_buttons_2,
            self.ligand_buttons,
            self.event_record_buttons,
            self.site_record_buttons,
            self.misc_buttons,
            self.quit_buttons,
            self.event_info_table,
            self.progress_table,
        ]

    def build(self):

        # Create main window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_position(gtk.WIN_POS_CENTER)
        self.window.set_border_width(10)
        self.window.set_default_size(600, 400)
        self.window.set_title("PANDDA inspect")

        # Main VBox object for the window
        main_vbox = gtk.VBox(spacing=5)
        self.window.add(main_vbox)

        # -----------------------------------------------------

        hbox = gtk.HBox(spacing=5)
        # Create buttones to allow user to quit
        quit_buttons = self.quit_buttons()
        frame = gtk.Frame();
        frame.add(quit_buttons)
        hbox.pack_start(frame)
        # Create progress summary table at the top of the window
        progress_table = self.progress_table()
        frame = gtk.Frame();
        frame.add(progress_table)
        hbox.pack_start(frame)
        # Create buttons to navigate between datasets
        nav_buttons = self.navigation_buttons_2()
        hbox.pack_start(nav_buttons)
        # Add to main vbox
        main_vbox.pack_start(hbox)

        # -----------------------------------------------------

        # Create buttons to navigate between datasets
        nav_buttons = self.navigation_buttons_1()
        frame = gtk.Frame();
        frame.add(nav_buttons)
        main_vbox.pack_start(frame)

        # -----------------------------------------------------

        hbox = gtk.HBox(homogeneous=False, spacing=5)
        # Create event summary table at the top of the window
        event_info_table = self.event_info_table()
        frame = gtk.Frame();
        frame.add(event_info_table)
        hbox.pack_start(frame)
        # Create buttons to control the ligand
        lig_buttons = self.ligand_buttons()
        frame = gtk.Frame();
        frame.add(lig_buttons)
        hbox.pack_start(frame)
        # Add to main vbox
        main_vbox.pack_start(hbox)

        # -----------------------------------------------------

        # Create buttones to record meta about the event
        rec_e_buttons = self.event_record_buttons()
        frame = gtk.Frame();
        frame.add(rec_e_buttons)
        main_vbox.pack_start(frame)

        # -----------------------------------------------------

        # Create buttones to record meta about the event
        rec_s_buttons = self.site_record_buttons()
        frame = gtk.Frame();
        frame.add(rec_s_buttons)
        main_vbox.pack_start(frame)

        # -----------------------------------------------------

        # Miscellaneous buttons
        misc_buttons = self.misc_buttons()
        frame = gtk.Frame();
        frame.add(misc_buttons)
        main_vbox.pack_start(frame)

        return self.window

    def launch(self, widget=None, *data):

        self.window.show_all()

        catchup(True)

    def quit(self):

        gtk.main_quit()

    def buttons(self):
        d = {}
        for w in self.widgets:
            d.update(
                w.buttons
                )
        return d

    def objects(self):
        d = {}
        for w in self.widgets:
            d.update(
                w.objects
                )
        return d

    def labels(self):
        d = {}
        for w in self.widgets:
            d.update(
                w.labels
                )
        return d
           

##########
def validate_path(path):

    if not os.path.exists(str(path)):
        raise MissingFile(
            'File does not exist: {}'.format(str(path))
            )
class TableHandler(object):

    require_not_empty = True

    def __init__(self, 
        input_csv_path,
        output_csv_path,
        ):

        validate_path(input_csv_path)
        
        self.input_path = pl.Path(input_csv_path)
        self.output_path = pl.Path(output_csv_path)

        table = self.read_csv(self.input_path)

        if self.require_not_empty is True:
            if len(table) == 0: 
                logger.info(str(table))
                raise SystemExit(
                    '\n\nInput csv contains no data!\n\n'
                    )

        self.initialise(table)

        if self.output_path.exists():

            self.update(
                master = table,
                other = self.read_csv(self.output_path),
                )
        
        self.table = table

        self.write()

    def get(self, index, col):

        return self.table[col].iat[index]

    def set(self, index, col, value):

        self.table[col].iat[index] = value

    def write(self):

        logger.info('Writing output csv: {}'.format(str(self.output_path)))

        self.table.to_csv(str(self.output_path))
class PanddaEventTableHandler(TableHandler):

    def read_csv(self, path):
            
        table = pd.read_csv(str(path), sep=',', dtype={'dtag': str})

        table = table.set_index(['dtag', 'event_num'])

        if self.require_not_empty is True:
            if len(table) == 0: 
                logger.info(str(table))
                raise SystemExit(
                    '\n\nNo events to inspect: Input events csv contains no data!\n\n'
                    )

        return table

    def initialise(self, table):

        table['Interesting'] = False
        table['Ligand Placed'] = False
        table['Ligand Confidence'] = 'Low'
        table['Comment'] = 'None'
        table['Viewed'] = False
        table['Ligand Names'] = ''

    def update(self, master, other):

        master.update(
            other[
                [
                    'Interesting', 
                    'Ligand Placed', 
                    'Ligand Confidence', 
                    'Comment', 
                    'Viewed',
                    'Ligand Names',
                    ]
                ]
            )
class PanddaSiteTableHandler(TableHandler):

    def read_csv(self, path):
            
        table = pd.read_csv(str(path), sep=',')

        table = table.set_index('site_num')

        if self.require_not_empty is True:
            if len(table) == 0: 
                logger.info(str(table))
                raise SystemExit(
                    '\n\nNo sites: Input sites csv contains no data!\n\n'
                    )

        return table

    def initialise(self, table):

        table['Name'] = None
        table['Comment'] = None

    def update(self, master, other):

        master.update(
            other[
                [
                    'Name', 
                    'Comment', 
                    ]
                ]
            )
class PanddaInspectTableHandler(object): 

    def __init__(self,
        files_dict,
        ):

        self.events = PanddaEventTableHandler(
            input_csv_path = files_dict['input_events'],
            output_csv_path = files_dict['output_events'],
            )

        self.sites = PanddaSiteTableHandler(
            input_csv_path = files_dict['input_sites'],
            output_csv_path = files_dict['output_sites'],
            )

    def write(self):
        self.events.write()
        self.sites.write()
###
class Div(object):

    __slots__ = (
        'id',
        'contents',
    )

    type = 'div'

    def __init__(self, **kw_args):

        # Ensure has contents attribute
        self.contents = []

        # Apply supplied args
        self.set(**kw_args)

    def append(self, other):
        self.contents.append(other)
        return other

    def extend(self, other):
        self.contents.extend(other)

    def set(self, **kw_args):
        for k, v in list(kw_args.items()):
            setattr(self, k, v)

    def get(self, key):
        assert key != 'contents'
        return getattr(self, key)

    def __getitem__(self, attr):
        return getattr(self, attr)

    def __setitem__(self, attr, value):
        setattr(self, attr, value)

    def update(self, dictionary):
        self.set(**dictionary)
class Block(Div):

    __slots__ = (
        'title',
        'text',
        'html',
        'image',
        'table',
        'footnote',
        'width',
        'colour',
        'classes',
        'styles',
        'max_width',
        'fancy_title',
        'title_size',
    )

    type = 'block'

    def __init__(
        self,
        width = 12,
        fancy_title = False,
        **kw_args
        ):

        # MUST BE FIRST
        super(Block, self).__init__(**kw_args)

        # Not captured by kw_args
        self.set(
            width = width,
            fancy_title = fancy_title,
        )
class Alert(Block):

    type = 'alert'
class ScrollX(Block):

    type = 'scroll'
class Panel(Block):

    __slots__ = Block.__slots__ + (
        'show',
    )

    type = 'panel'
class TabSet(Div):

    __slots__ = Div.__slots__ + (
        'title',
        'width',
        'title_size',
        'classes',
        'colour',
    )

    type = 'tabs'

    def set_active(
        self,
        tab = None,
        i_tab = 0,
        ):

        if not self.contents:
            return

        if tab is None:
            tab = self.contents[i_tab]

        assert tab in self.contents

        tab.active = True
class Tab(Block):

    __slots__ = Block.__slots__ + (
        'alt_title',
        'active',
    )

    type = 'tab'

    def __init__(
        self,
        alt_title = None,
        fancy_title = True,
        **kw_args
        ):

        super(Tab, self).__init__(**kw_args)

        if alt_title is None:
            alt_title = kw_args['title']

        self.set(
            alt_title = alt_title,
            fancy_title = fancy_title,
        )
class ProgressBar(Block):

    __slots__ = (
        'data',
    )

    type = 'progressbar'

    def __init__(
        self,
        data,
        width = 12,
        add_counts = True,
        add_percentages = True,
        **kw_args
        ):

        # MUST BE FIRST
        super(Block, self).__init__(**kw_args)

        self.set(
            width = width,
            )

        self.set_data(
            data = data,
            add_counts = add_counts,
            add_percentages = add_percentages,
            )

    def set_data(self, data, add_counts, add_percentages):

        total_width = float(
            np.sum(
                [d['value'] for d in data]
                )
            )

        self.data = []

        for d in data: 
            
            l = d['label']
            v = d['value']

            p = 100.0 * float(v) / total_width

            if (add_counts and add_percentages):
                l = '{}: {} ({}%)'.format(
                    str(l),
                    str(v),
                    int(p),
                    )
            elif add_counts: 
                l = '{} ({})'.format(
                    str(l), 
                    str(v),
                    )
            elif add_percentages: 
                l = '{} ({}%)'.format(
                    str(l), 
                    int(p),
                    )
            else:
                l = str(l)

            # Pass through all attributes of input dict
            d_copy = copy.deepcopy(d)
            d_copy['size'] = p
            d_copy['text'] = l

            self.data.append(
                d_copy
                )
class DivFormatter(object):

    def __init__(self,
        div_class = Block,
        *args, **kwargs
        ):

        self.div_class = div_class
        self.args = args
        self.kwargs = kwargs

    def __call__(self):

        return self.div_class()
class CitationFormatter(DivFormatter):

    def __call__(self,
        title, 
        journal,
        year,
        authors,
        link = None,
        div_class = Block,
        ):

        text_lines = [
            self.format_title(title),
            self.format_authors(authors),
            self.format_journal_year(journal, year),
            ]
            
        if link is not None:
            text_lines.append(
                self.format_link(
                    link = link,
                    link_text = link,
                    )
                )

        args, kwargs = self.get_args()

        main_div = div_class(
            text = '<br>'.join(text_lines),
            *args, **kwargs
            )

        return main_div

    def get_args(self):

        args = copy.deepcopy(self.args)
        kwargs = copy.deepcopy(self.kwargs)

        return args, kwargs

    def format_title(self, title):

        return '<strong>{}</strong>'.format(title)

    def format_journal_year(self, journal, year):

        journal_year = '{journal} ({year})'.format(
            journal = journal,
            year = year,
            )

        return '<small><strong>{}</strong></small>'.format(journal_year)

    def format_authors(self, authors):

        if isinstance(authors, list):
            if len(authors) == 1:
                authors = authors[0]
            else: 
                authors = (
                    ', '.join(authors[:-1]) + 
                    ' and ' + 
                    authors[-1]
                    )
        else: 
            authors = str(authors)

        return '<small>{}</small>'.format(authors)

    def format_link(self, link, link_text):

        link_text = (
            '<a href="{link}"><small>{link_text}</small></a>'.format(
                link = link,
                link_text = link_text,
                )
            )

        return link_text

class ImageEmbedder(object):

    def __init__(self, embed=False, relative_to=None):

        self.embed = embed
        self.relative_to = relative_to

    def __call__(self, image_path):

        # Read image and return as string
        if (self.embed is True):
            from giant.html import png2base64src_maybe
            return png2base64src_maybe(
                image_path, 
                print_on_missing = False,
                )

        # Create relative path
        if self.relative_to is not None:
            image_path = os.path.relpath(
                path = image_path, 
                start = self.relative_to,
                )

        return image_path
d_label = '<span class="label label-info">{k}</span>'
class MakePanddaInspectHtml(object):

    output_filename = 'pandda_inspect.html'
    template_name = 'pandda_page.html'

    def __init__(self, 
        output_directory,
        ):

        self.output_directory = pl.Path(output_directory)

        self.input_files = self.build_input_files_dict()

        self.input_path_prefix = self.find_input_prefix(
            filepath = self.input_files['html']['main_html'],
            )

        self.output_path = pl.Path(
            self.get_path(['html','inspect_html'])
            )

        self.image = ImageEmbedder(
            embed = False,
            relative_to = str(self.output_path.parent),
            )

    def __call__(self,
        inspector,
        output_files,
        ):

        results = self.unpack_inspector(
            inspector = inspector,
            )

        contents = self.get_contents(
            results = results,
            output_files = output_files,
            )

        self.make_output(
            contents = contents,
            )

        return self.output_path

    def build_input_files_dict(self):

        import json

        json_file = (
            self.output_directory / 'results.json'
            )

        if not json_file.exists(): 
            raise MissingFile('Json results file not found: {}'.format(str(json_file)))

        json_string = open(str(json_file), 'r').read()

        pandda_log_obj = json.loads(json_string)

        input_files = pandda_log_obj['output_files']

        return input_files

    def find_input_prefix(self, filepath):

        filepath = pl.Path(filepath)

        prefix = filepath.parent

        for _ in filepath.parts:

            test_f = filepath.relative_to(prefix)

            if (self.output_directory / test_f).exists():
                break

            prefix = prefix.parent
        
        return prefix

    def get_path(self, input_keys):

        assert isinstance(input_keys, list)

        d = self.input_files
        for k in input_keys:
            d = d[k]

        assert not hasattr(d, 'keys')

        output_path = (
            self.output_directory / pl.Path(d).relative_to(
                self.input_path_prefix
                )
            )

        return str(output_path)

    def unpack_inspector(self,
        inspector,
        ):

        results = {
            'event_table' : inspector.tables.events.table.reset_index(),
            'site_table' : inspector.tables.sites.table.reset_index(),
        }

        return results

    def get_contents(self, 
        results,
        output_files,
        ):

        contents = []

        contents.extend(
            self.get_header()
            )

        contents.extend(
            self.get_main(
                results = results,
                output_files = output_files,
                )
            )

        return contents

    def get_header(self):

        contents = [
            Block(
                title = "PanDDA inspection summary",
                fancy_title = True,
                ),
            ]

        return contents

    def get_main(self,
        results,
        output_files,
        ):

        event_table = results['event_table']

        event_table_dict = self.build_event_table_summary_dict(
            event_table = event_table,
            )

        ###

        contents = []

        contents.extend(
            self.get_inspection_progress_bars(
                summary_dict = event_table_dict,
                )
            )

        contents.extend(
            self.get_inspection_site_summary(
                summary_dict = event_table_dict,
                output_files = output_files,
                )
            )

        contents.extend(
            self.get_inspection_event_summary(
                event_table = event_table,
                )
            )

        return contents

    def get_inspection_progress_bars(self,
        summary_dict,
        ):

        block1 = Block(
            classes = ['bordered'],
            contents = [
                ProgressBar(
                    title = "Inspection Progress",
                    title_size = 4,
                    data = [
                        {
                            'label' : 'Fitted',
                            'value' : summary_dict['n_fitted'],
                            'colour' : 'success',
                            },
                        {
                            'label' : 'Unviewed',
                            'value' : summary_dict['n_unviewed'],
                            'colour' : 'info',
                            },
                        {
                            'label' : 'No Ligand Fitted',
                            'value' : summary_dict['n_empty'],
                            'colour' : 'danger',
                            },
                        ],
                    add_counts = True,
                    add_percentages = False,
                    ),
                Alert(
                    text = "Total events: {}".format(
                        summary_dict['n_blobs']
                        ),
                    colour = "info",
                    width = 12,
                    ),
                Alert(
                    text = "Fitted ligands: {}".format(
                        summary_dict['n_fitted']
                        ),
                    colour = "success",
                    width = 4,
                    ),
                Alert(
                    text = "Unviewed: {}".format(
                        summary_dict['n_unviewed']
                        ),
                    colour = "info",
                    width = 4,
                    ),
                Alert(
                    text = "No ligand fitted: {}".format(
                        summary_dict['n_empty']
                        ),
                    colour = "danger",
                    width = 4,
                    ),
                ],
            )

        block2 = Block(
            classes = ['bordered'],
            contents = [
                Alert(
                    text = "Datasets with ligands: {}".format(
                        summary_dict['n_datasets_w_hit']
                        ),
                    colour = "info",
                    width = 6,
                    ),
                Alert(
                    text = "Sites with ligands: {}".format(
                        summary_dict['n_sites_w_hit']
                        ),
                    colour = "info",
                    width = 6,
                    ),
                ],
            )

        block3 = Block(
            classes = ['bordered'],
            contents = (
                [
                    ProgressBar(
                        title = "Fitted ligands",
                        title_size = 4,
                        data = [
                            {
                                'label' : 'High confidence',
                                'value' : summary_dict['n_high_confidence'],
                                'colour' : 'success',
                                },
                            {
                                'label' : 'Medium confidence',
                                'value' : summary_dict['n_medium_confidence'],
                                'colour' : 'info',
                                },
                            {
                                'label' : 'Low confidence',
                                'value' : summary_dict['n_low_confidence'],
                                'colour' : 'danger',
                                },
                            {
                                'label' : 'Unranked',
                                'value' : summary_dict['n_unranked'],
                                'colour' : 'info',
                                },
                            ],
                        add_counts = True,
                        add_percentages = False,
                        ),
                    Alert(
                        text = "Total ligands: {}".format(
                            summary_dict['n_fitted']
                            ),
                        colour = "info",
                        width = 12,
                        ),
                    Alert(
                        text = "High Confidence: {}".format(
                            summary_dict['n_high_confidence']
                            ),
                        colour = "success",
                        width = 4,
                        ),
                    Alert(
                        text = "Medium Confidence: {}".format(
                            summary_dict['n_medium_confidence']
                            ),
                        colour = "warning",
                        width = 4,
                        ),
                    Alert(
                        text = "Low Confidence: {}".format(
                            summary_dict['n_low_confidence']
                            ),
                        colour = "danger",
                        width = 4,
                        ),
                    Alert(
                        text = "Unfitted, marked interesting: {}".format(
                            summary_dict['n_interesting_unfitted']
                            ),
                        colour = "danger",
                        width = 6,
                        ),
                    ]
                if (summary_dict['n_fitted'] > 0) else 
                [
                    Alert(
                        text = "Nothing to show here yet: no fitted ligands"
                        )
                    ]
                ),
            )

        return [block1, block2, block3]

    def get_inspection_site_summary(self,
        summary_dict,
        output_files,
        ):

        site_counts = summary_dict['n_hits_per_site']
        site_images = output_files['site_events']

        blocks = [
            Block(
                width = 6,
                text = 'Sites and Events (front)',
                image = self.get_path(['graphs','events_front']),
                ),
            Block(
                width = 6,
                text = 'Sites and Events (back)',
                image = self.get_path(['graphs','events_back']),
                ),
            (
                ProgressBar(
                    title = "Fitted ligands distribution over sites",
                    data = [
                        {
                            'label' : 'Site {}'.format(i),
                            'value' : v,
                            'colour' : ['default','info'][i%2],
                            }
                            for i,v in sorted(site_counts.items())
                        ],
                    add_counts = True,
                    add_percentages = True,
                    )
                if sum(site_counts.values()) > 0 else
                Alert(
                    text = 'Nothing to show here yet: no fitted ligands',
                    )
                ),
            ScrollX(
                contents = [
                    Block(
                        title = 'Site {}'.format(site_num),
                        text = 'Number of events: {}'.format(n_events),
                        image = self.image(site_images[site_num]),
                        width = 5,
                        )
                    for site_num, n_events
                    in sorted(site_counts.items())
                    ],
                ),
            ]

        return blocks

    def get_inspection_event_summary(self,
        event_table,
        ):

        column_labels = [] # subsample the columns

        # ['Dataset','Viewed','Interesting','Lig. Placed','Event','Site','1 - BDC','Z-Peak','Map Res.','Map Unc.','Confidence','Comment','']

        blocks = [
            self.make_table_block(
                title = 'Events/Fitted Ligands Table', 
                table = event_table, 
                width = 12,
                ),
            ]

        return blocks

    def build_event_table_summary_dict(self,
        event_table,
        ):

        # Input counts

        n_blobs = len(
            event_table.index
            )

        n_sites = len(
            set(event_table['site_num'])
            )

        n_datasets = len(
            # this assumes the index is (dtag, e_idx)
            set(event_table['dtag'])
            )

        # Blobs Inspected/Modelled/Empty

        n_fitted = sum(
            event_table['Ligand Placed']
            )

        n_viewed = sum(
            event_table['Viewed']
            )

        n_empty = (
            n_viewed - n_fitted
            )

        n_unviewed = (
            n_blobs - n_viewed
            )

        # Interesting unfitted (bookmarked)

        n_interesting_unfitted = sum(
            event_table["Interesting"][event_table['Ligand Placed']==False]
            )

        # Confidence of models

        n_high_confidence = sum(
            event_table["Ligand Placed"][event_table['Ligand Confidence']=='High']
            )

        n_medium_confidence = sum(
            event_table["Ligand Placed"][event_table['Ligand Confidence']=='Medium']
            )
        
        n_low_confidence = sum(
            event_table["Ligand Placed"][event_table['Ligand Confidence']=='Low']
            )

        n_unranked = (
            n_fitted - n_high_confidence - n_medium_confidence - n_low_confidence
            )

        # Datasets/sites with hits

        try:    
            n_datasets_w_hit = len(
                set(zip(*event_table.index[event_table['Ligand Placed'] == True])[0])
                )
        except: 
            n_datasets_w_hit = 0

        try:    
            n_sites_w_hit = len(
                set(event_table['site_num'][event_table['Ligand Placed'] == True])
                )
        except: 
            n_sites_w_hit = 0

        # Hits per site

        n_hits_per_site = {
            i_site : sum(
                event_table["Ligand Placed"][event_table['site_num']==i_site]
                )
            for i_site in range(1, n_sites+1)
        }

        ###

        s_dict = dict(
            n_blobs = n_blobs,
            n_sites = n_sites,
            n_datasets = n_datasets,
            n_fitted = n_fitted,
            n_viewed = n_viewed,
            n_empty = n_empty,
            n_unviewed = n_unviewed,
            n_interesting_unfitted = n_interesting_unfitted,
            n_high_confidence = n_high_confidence,
            n_medium_confidence = n_medium_confidence,
            n_low_confidence = n_low_confidence,
            n_unranked = n_unranked,
            n_datasets_w_hit = n_datasets_w_hit,
            n_sites_w_hit = n_sites_w_hit,
            n_hits_per_site = n_hits_per_site,
        )

        return s_dict

    def make_table_block(self, title, table, width=12):

        table_html = table.to_html(
            index = False,
            bold_rows = False,
            na_rep = '',
            classes = ['table table-striped table-hover datatable nowrap'],
            ).replace(
            '<th></th>','<th>Dataset</th>'
            ).replace(
            'border="1" ', ''
            )

        block = Alert(
            title = title,
            width = width,
            table = table_html,
            )

        return block

    def make_output(self,
        contents,
        ):

        logger.info(
            'Writing output HTML: {}'.format(self.output_filename)
            )

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

        from pandda.html import HTML_ENV
        template = HTML_ENV.get_template(self.template_name)

        header_title = 'PanDDA Inspection Summary'

        # ===========================================================>
        # Construct the data object to populate the template
        output_data = {
            'header_title' : header_title,
            #'body_header' : body_header,
            'contents' : contents,
            }
        # Jsons
        #if json_plots:
        #    output_data['json_plots'] = json_plots
        # ===========================================================>

        # Write out and format
        with open(str(self.output_path), 'w') as out_html:
            out_html.write(
                template.render(output_data)#.encode( "utf-8" )
                )

        logger.info(
            'Output HTML written to {}'.format(str(self.output_path))
            )
###

class PanddaInspect(object):
    """Main Object in pandda.inspect"""

    def __init__(self, 
        working_directory, 
        files_dict,
        mode,
        ):

        # from pandda.inspect.tables import (
        #     PanddaInspectTableHandler
        #     )

        self.tables = PanddaInspectTableHandler(
            files_dict = files_dict,
            )

        self.controller = PanddaEventListController(
            event_table = self.tables.events.table,
            get_event_handler = GetEventHandler(
                pandda_directory = working_directory, 
                pandda_files_dict = files_dict['pandda_files'], 
                pandda_path_prefix = files_dict['remove_prefix'],
                ),
            )

        if mode in ['events']:

            # from pandda.inspect.html import (
            #     MakePanddaInspectHtml,
            #     )

            write_html = MakePanddaInspectHtml(
                output_directory = (
                    working_directory
                    ),
                )

            self.update_output = UpdateOutput(
                inspector = self,
                write_html = write_html,
                )

        else:

            self.update_output = DummyUpdateOutput(
                inspector = self,
                mode = mode,
                ) 

        self._registered_windows = []

    def __call__(self):

        self.gui = PanddaGUI()
        window = self.gui.build()

        self.link_window(window)
        self.link_gui(self.gui)

        return self.gui

    def link_window(self, window):

        window.connect("delete_event", self.quit)
        window.connect("destroy_event", self.quit)

    def link_gui(self, gui):

        self.buttons = gui.buttons()
        self.objects = gui.objects()
        self.labels = gui.labels()

        self.link_navigation()
        self.link_misc_buttons()

        self.link_structure_buttons()
        self.link_ligand_buttons()
        self.link_map_buttons()
        self.link_record_meta_buttons()

        self.refresh_gui()

    def quit(self, widget=None, *data):
        self.sync_and_write_tables()
        self.update_output()
        self.gui.quit()

    def sync_and_write_tables(self):
        self.sync_tables()
        self.write_tables()

    def sync_tables(self):
        """Sync the GUI and the output tables"""

        self.tables.events.set(
            index = self.controller.event_tracker.get(),
            col ='Comment', 
            value = self.objects['event comment text'].get_text(), 
            )

        self.tables.events.set(
            index = self.controller.event_tracker.get(),
            col = 'Viewed', 
            value = True,
            )

        # Site records
        self.tables.sites.set(
            index = self.controller.site_tracker.get(),
            col = 'Name', 
            value = self.objects['site name text'].get_text(), 
            )

        self.tables.sites.set(
            index = self.controller.site_tracker.get(),
            col = 'Comment', 
            value = self.objects['site comment text'].get_text(), 
            )

    def write_tables(self):
        """Write the output tables"""
        self.tables.write()

    def refresh_gui(self):
        """Update information in the GUI from the tables"""

        self.delete_registered_windows()

        # Event/Site Tracker

        self.labels['site_val'].set_label(
            str(self.controller.site_tracker.get()+1)
            )

        self.labels['site_tot'].set_label(
            str(self.controller.site_tracker.n_total)
            ) # move elsewhere only needs doing once...

        self.labels['rank_val'].set_label(
            str(self.controller.event_tracker.get()+1)
            )

        self.labels['rank_tot'].set_label(
            str(self.controller.event_tracker.n_total)
            ) # move elsewhere only needs doing once...

        # Event and dataset information

        e = self.controller.event_handler.event

        self.labels['dtag'].set_label(
            '<b>' + str(e.dtag) + '</b>'
            )

        self.labels['e_num'].set_label(
            '{event_num}/{n_total}'.format(
                event_num = str(
                    e.event_num
                    ),
                n_total = str(
                    self.controller.event_counts.get(e.dtag, '?')
                    ),
                )
            )

        self.labels['e_1_bdc'].set_label(
            str(e.est_1_bdc)
            if e.est_1_bdc is not None
            else 'n/a'
            )

        self.labels['zpeak'].set_label(
            str(round(e.z_peak, 3))
            if e.z_peak is not None
            else 'n/a'
            )

        self.labels['map_res'].set_label(
            str(e.map_resolution)
            )

        self.labels['map_unc'].set_label(
            str(e.map_uncertainty)
            )

        self.labels['rwork_rfree'].set_label(
            '{} / {}'.format(*e.rwork_rfree)
            )

        # Reset the event comment boxes

        self.objects['event comment text'].set_text(
            str(
                self.tables.events.get(
                    index = self.controller.event_tracker.get(),
                    col = 'Comment',
                    )
                )
            )

        # Reset the site comment boxes

        self.objects['site name text'].set_text(
            str(
                self.tables.sites.get(
                    index = self.controller.site_tracker.get(), 
                    col = 'Name',
                    )
                )
            )

        self.objects['site comment text'].set_text(
            str(
                self.tables.sites.get(
                    index = self.controller.site_tracker.get(),
                    col = 'Comment',
                    )
                )
            )

        # Update the radio buttons - "Interesting"

        interesting = self.tables.events.get(
            index = self.controller.event_tracker.get(),
            col = 'Interesting',
            )

        if interesting == True:
            self.buttons['tp'].set_active(True)
        else:
            self.buttons['fp'].set_active(True)

        # Update the radio buttons - "Ligand Placed"

        ligand_placed = self.tables.events.get(
            index = self.controller.event_tracker.get(),
            col = 'Ligand Placed',
            )

        if ligand_placed == True:
            self.buttons['placed'].set_active(True)
        else:
            self.buttons['not placed'].set_active(True)

        # Update the radio buttons - "Ligand Confidence"

        ligand_confidence = self.tables.events.get(
            index = self.controller.event_tracker.get(),
            col = 'Ligand Confidence',
            )

        if ligand_confidence == 'High':
            self.buttons['highconf'].set_active(True)
        elif ligand_confidence == 'Medium':
            self.buttons['medconf'].set_active(True)
        else:
            self.buttons['lowconf'].set_active(True)

        # Reset the merge button

        self.buttons['merge-ligand'].child.set_text(
            "Merge Ligand\nWith Model"
            )

    def link_navigation(self):

        # Basic Navigation
        self.buttons['next'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.save(),
                self.controller.next_event(),
                self.refresh_gui(),
                ],
            )

        self.buttons['prev'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.prev_event(),
                self.refresh_gui(),
                ],
            )

        self.buttons['skip'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.next_event(),
                self.refresh_gui(),
                ],
            )

        self.buttons['next-unviewed'].connect(
            "clicked",
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.next_event_unviewed(),
                self.refresh_gui(),
                ],
            )

        self.buttons['next-modelled'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(),
                self.controller.next_event_modelled(),
                self.refresh_gui(),
                ],
            )

        self.buttons['next-site'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.next_site(),
                self.refresh_gui(),
                ],
            )

        self.buttons['prev-site'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.prev_site(),
                self.refresh_gui(),
                ],
            )

        self.buttons['go-to'].connect(
            "clicked", 
            lambda x: [
                self.sync_and_write_tables(), 
                self.controller.next_event_for_dataset(
                    dataset_id = self.objects['go-to-text'].get_text().strip(' '),
                    ),
                self.refresh_gui(),
                ],
            )

        self.objects['go-to-text'].connect(
            "activate", 
            lambda x: [
                self.buttons['go-to'].emit("clicked"),
                ],
            )

    def link_misc_buttons(self):

        self.buttons['quit'].connect(
            "clicked", 
            lambda x: [
                self.quit(),
                ],
            )

        # self.buttons['summary'].connect(
        #     "clicked", 
        #     lambda x: [
        #         self.sync_and_write_tables(), 
        #         self.update_output(), 
        #         os.system('ccp4-python -Qnew -m pandda.jiffies.pandda_summary &'),
        #         ],
        #     )

        # self.buttons['updatehtml'].connect(
        #     "clicked", 
        #     lambda x: [
        #         self.sync_and_write_tables(), 
        #         self.update_output(),
        #         ],
        #     )

    def link_structure_buttons(self):

        self.buttons['save'].connect(
            "clicked", 
            lambda x: [
                self.controller.save(),
                ],
            )

        self.buttons['reload'].connect(
            "clicked", 
            lambda x: [
                self.controller.event_handler.revert_to_last_model(),
                ],
            )

        self.buttons['reset'].connect(
            "clicked", 
            lambda x: [
                self.controller.event_handler.revert_to_input_model(),
                ],
            )

    def link_ligand_buttons(self):

        # Ligand Buttons
        self.buttons['merge-ligand'].connect(
            "clicked", lambda x: [
                self.buttons['merge-ligand'].child.set_text(
                    'Already Merged\n(Click to Repeat)'
                    ),
                self.controller.event_handler.merge_ligand_with_model(),
                self.buttons['placed'].clicked(), # yes, a ligand has been placed
                self.buttons['tp'].clicked(), # yes, it would seem this is interesting
                # self.buttons['highconf'].clicked(), # I don't think clicking this button automatically is justified
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Ligand Names', 
                    value = format_ligand_names(
                        ligand_name_list = self.tables.events.get(
                            index = self.controller.event_tracker.get(),
                            col = 'Ligand Names',
                            ).split(','),
                        new_name = (
                            self.controller.event_handler.ligand_handler.get_i_ligand_name()
                            ),
                        ),
                    ),
                ],
            )

        self.buttons['move-ligand'].connect(
            "clicked", lambda x: [
                self.controller.event_handler.move_ligand_here(),
                ],
            )

        self.buttons['next-ligand'].connect(
            "clicked",
            lambda x: [
                self.controller.event_handler.next_ligand(),
                ],
            )

        # TODO TODO TODO
        # buttons['prev-ligand'].connect(
        #     "clicked",
        #     lambda x: [
        #         self.controller.event_handler.prev_ligand(),
        #         ],
        #     )

        # buttons['set-ligand'].connect(
        #     "clicked",
        #     lambda x: [
        #         self.controller.event_handler.set_ligand(
        #             ligand_id = objects['...'].get_text(),
        #             ),
        #         ],
        #     )

    def link_map_buttons(self):

        # self.buttons['load-ground-state-map'].connect(
        #     "clicked", 
        #     lambda x: [
        #         modal_msg("clicking buttons is fun, but this map is already open."),
        #         #self.parent.coot.load_ground_state_map(event=self.parent.current_event),
        #         ],
        #     )

        self.buttons['load-full-dataset-mtz'].connect(
            "clicked", 
            lambda x: [
                self.controller.event_handler.model_map_handler.load_full_dataset_mtz(),
                ],
            )

        self.buttons['load-original-model'].connect(
            "clicked",
            lambda x: [
                self.controller.event_handler.model_map_handler.load_original_model_not_active(),
                ],
            )

        self.buttons['create-ligand'].connect(
            "clicked",
            lambda x: [
                self.controller.event_handler.update_and_go_to_ligand(
                    info_dict = MakeNewLigandModal(
                        output_directory = (
                            self.controller.event_handler.event_files['ligand_dir']
                            ),
                        )(),
                    show = True,
                    ),
                ],
            )

        self.buttons['custom-event-map'].connect(
            "clicked",
            lambda x: [
                self.register_window(
                    window = GetMakeCustomEventMapWindow(
                        event_handler = self.controller.event_handler,
                        )()
                    ),
                ],
            )

    def link_record_meta_buttons(self):

        self.buttons['tp'].connect(
            "clicked", 
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Interesting', 
                    value = True,
                    ),
                ],
            )

        self.buttons['fp'].connect(
            "clicked", 
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Interesting', 
                    value = False,
                    ),
                ],
            )

        self.buttons['highconf'].connect(
            "clicked", 
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Ligand Confidence',
                    value = 'High',
                    ),
                ],
            )

        self.buttons['medconf'].connect(
            "clicked", 
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Ligand Confidence',
                    value = 'Medium',
                    ),
                ],
            )

        self.buttons['lowconf'].connect(
            "clicked", 
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Ligand Confidence',
                    value = 'Low',
                    ),
                ],
            )

        self.buttons['placed'].connect(
            "clicked",
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Ligand Placed', 
                    value = True,
                    ),
                ],
            )

        self.buttons['not placed'].connect(
            "clicked",
            lambda x: [
                self.tables.events.set(
                    index = self.controller.event_tracker.get(),
                    col = 'Ligand Placed', 
                    value = False,
                    ),
                ],
            )

    def register_window(self, window):

        self._registered_windows.append(
            window
            )

        return window

    def delete_registered_windows(self):

        while self._registered_windows:

            window = self._registered_windows.pop()
            window.destroy()

###
class MakePanddaEventTable(object):

    columns = [
        "dtag",
        "event_num",
        "site_num",
        "event_fraction",
        "bdc",
        "z_peak",
        "z_mean",
        "cluster_size",
        "x",
        "y",
        "z",
        "analysed_resolution",
        "map_uncertainty",
        "global_correlation",
        "local_correlation",
        "r_work",
        "r_free",
        ]

    def __init__(self, output_path=None):
        self.output_path = output_path

    def __call__(self, event_dicts, output_path=None):

        if (output_path is None): 
            output_path = self.output_path

        event_dicts = self.populate_dicts_maybe(
            event_dicts = event_dicts,
            )

        df = pd.DataFrame(
            data = event_dicts,
            columns = self.columns,
            )

        df['interesting'] = False

        df = df.set_index(["dtag", "event_num"])

        # df = df.sort_values(
        #     by = ['site_idx','z_peak'], 
        #     ascending = [1,0],
        #     )

        if (output_path is not None): 
            df.to_csv(
                path_or_buf = str(output_path),
                )

        return df

    def populate_dicts_maybe(self, event_dicts):

        event_dicts = copy.deepcopy(event_dicts)

        for e in event_dicts:

            xyz = e.get('xyz_centroid', (None, None, None))

            if xyz is None: 
                xyz = (None, None, None)

            e.setdefault('x', xyz[0])
            e.setdefault('y', xyz[1])
            e.setdefault('z', xyz[2])

        return event_dicts
class MakePanddaSiteTable(object):

    columns = [
        "site_num",
        "n_events",
        "max_value",
        "xyz_centroid",
        "xyz_extent",
        ]

    def __init__(self, output_path=None):
        self.output_path = output_path

    def __call__(self, site_dicts, output_path=None):

        if (output_path is None):
            output_path = self.output_path

        site_dicts = self.populate_dicts_maybe(
            site_dicts = site_dicts,
            )

        df = pd.DataFrame(
            data = site_dicts,
            columns = self.columns,
            )

        df = df.set_index("site_num")

        df = df.sort_values(
            by = 'site_num',
            ascending = True,
            )

        if (output_path is not None):
            df.to_csv(
                path_or_buf = str(output_path),
                )

        return df

    def populate_dicts_maybe(self, site_dicts):

        site_dicts = copy.deepcopy(site_dicts)

        for s in site_dicts:

            xyz = s.get('xyz_centroid', (None, None, None))

            if xyz is None: 
                xyz = (None, None, None)

            s.setdefault('x', xyz[0])
            s.setdefault('y', xyz[1])
            s.setdefault('z', xyz[2])

        return site_dicts
###
class ExpandEventTable(object):

    def __init__(self,
        dataset_info_csv,
        dataset_tags,
        output_path = None,
        keep_original_events = False,
        remove_datasets_with_events = False,
        ):
        
        self.dataset_info_csv = (
            dataset_info_csv
            )

        self.dataset_tags = (
            dataset_tags
            )

        self.make_pandda_event_table = MakePanddaEventTable(
            output_path = (
                str(output_path)
                if output_path is not None
                else None
                ),
            )

        self.keep_original_events = (
            keep_original_events
            )

        self.remove_datasets_with_events = (
            remove_datasets_with_events
            )

    def __call__(self, 
        event_table = None, 
        ):

        dataset_table = pd.read_csv(
            str(self.dataset_info_csv), 
            sep = ',', 
            dtype = {'dtag': str},
            )

        dummy_event_dicts = self.make_event_dicts_from_dataset_table(
            dataset_dicts = list(
                dataset_table.to_dict('records')
                ),
            )

        #

        if self.keep_original_events is False:

            event_dicts = []

        else:

            assert event_table is not None

            event_dicts = list(
                event_table.to_records()
                )

            if self.remove_datasets_with_events is True:

                event_dtags = sorted(
                    set([
                        e['dtag']
                        for e in event_dicts
                        ])
                    )

                dummy_event_dicts = [
                    d
                    for d in dummy_event_dicts
                    if d['dtag'] not in event_dtags
                ]

        #

        df = self.make_pandda_event_table(
            event_dicts = (
                event_dicts + dummy_event_dicts
                ),
            )

        return df

    def make_event_dicts_from_dataset_table(self, dataset_dicts):

        event_dicts = []

        defaults = {
            "event_num" : 1,
            "site_num" : 1,
            "event_fraction" : 1.0,
            "bdc" : 0.0,
        }

        for d in dataset_dicts: 

            if d['dtag'] not in self.dataset_tags:
                continue

            e = {
                c : d.get(
                    c, 
                    defaults.get(c, None),
                    ) 
                for c in self.make_pandda_event_table.columns
            }

            event_dicts.append(e)

        return event_dicts
class MakeDummySiteTable(object):

    def __init__(self,
        output_path = None,
        ):

        self.make_pandda_site_table = MakePanddaSiteTable(
            output_path = str(output_path),
            )

    def __call__(self):

        site_dicts = [
            {
            "site_num" : 1,
            "n_events" : None,
            "max_value" : None,
            "xyz_centroid" : None,
            "xyz_extent" : None,
            }
        ]

        df = self.make_pandda_site_table(
            site_dicts = (
                site_dicts
                ),
            )

        return df
###
class GetPanddaInspectInputOutputFiles(object):

    def __init__(self, 
        input_directory,
        output_directory,
        mode = 'events',
        ):

        self.input_directory = pl.Path(
            input_directory
            )

        self.output_directory = pl.Path(
            output_directory
            )

        self.mode = mode

        assert self.input_directory.exists()

        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)

    def __call__(self):

        results_dict = self.get_json_data()

        output_files = results_dict['output_files']

        ret_dict = {
            'pandda_files' : output_files,
            }

        ret_dict.update(
            self.get_input_files(output_files)
            )

        ret_dict.update(
            self.get_output_files()
            )

        return ret_dict

    def get_json_data(self):

        results_json = (
            self.input_directory / "results.json"
            )

        results_dict = json.loads(
            open(str(results_json), 'r').read()
            )

        return results_dict

    def get_input_files(self, output_files):

        remove_prefix = self.find_output_prefix(
            output_files['tables']['dataset_info']
            )

        if self.mode == 'events':

            input_events_csv = (
                self.input_directory / (
                    pl.Path(
                        output_files['tables']['events_table']
                        ).relative_to(remove_prefix)
                    )
                )

            input_sites_csv = (
                self.input_directory / (
                    pl.Path(
                        output_files['tables']['sites_table']
                        ).relative_to(remove_prefix)
                    )
                )

        elif self.mode == 'datasets':

            input_events_csv = (
                self.output_directory / (
                    pl.Path(
                        'dummy_events.csv'
                        )
                    )
                )

            input_sites_csv = (
                self.output_directory / (
                    pl.Path(
                        'dummy_sites.csv'
                        )
                    )
                )

            # need to actually make these tables...

            make_event_table = ExpandEventTable(
                dataset_info_csv = str(
                    self.input_directory / (
                        pl.Path(
                            output_files['tables']['dataset_table']
                            ).relative_to(remove_prefix)
                        )
                    ),
                dataset_tags = [
                    d 
                    for (d, d_files) 
                    in output_files['dataset_files'].items()
                    if d_files.get('output_data')
                    ],
                output_path = str(input_events_csv),
                )

            make_site_table = MakeDummySiteTable(
                output_path = str(input_sites_csv),
                )

            make_event_table()
            make_site_table()
        
        assert input_events_csv.exists()
        assert input_sites_csv.exists()

        return {
            'input_events' : input_events_csv,
            'input_sites' : input_sites_csv,
            'remove_prefix' : remove_prefix,
        }

    def get_output_files(self):

        out_dir = (
            self.output_directory
            )

        if not out_dir.exists():
            out_dir.mkdir(parents=True)

        if self.mode == 'events': 

            output_events_csv = str(
                out_dir / "pandda_inspect_events.csv"
                )

            output_sites_csv = str(
                out_dir / "pandda_inspect_sites.csv"
                )

        elif self.mode == 'datasets':

            output_events_csv = str(
                out_dir / "pandda_inspect_events.datasets.csv"
                )

            output_sites_csv = str(
                out_dir / "pandda_inspect_sites.datasets.csv"
                )

        return {
            'output_events' : output_events_csv,
            'output_sites' : output_sites_csv,
        }

    def find_output_prefix(self, filepath):

        filepath = pl.Path(filepath)

        top_dir = self.input_directory

        prefix = filepath.parent

        for _ in filepath.parts:

            test_f = filepath.relative_to(prefix)

            if (top_dir / test_f).exists():
                break

            prefix = prefix.parent
        
        return prefix
###

if __name__ == '__main__':

    #############################################################################################
    #
    # CHANGE COOT SETTINGS
    #
    #############################################################################################

    coot_setup()

    try:
        coot_customisation()
    except:
        pass

    #############################################################################################
    #
    # RUN
    #
    #############################################################################################

    # args = parser.parse_args(sys.argv)
    args, unknown = parser.parse_known_args(sys.argv)

    working_directory = pl.Path(
        args.pandda_directory
        # os.getcwd()
        )
    
    try: 

        splash = SplashScreen()

        # from pandda.inspect.io import (
        #     GetPanddaInspectInputOutputFiles,
        #     )

        get_io_files = GetPanddaInspectInputOutputFiles(
            input_directory = (
                working_directory
                ),
            output_directory = (
                working_directory / 'inspect'
                ),
            mode = args.mode,
            )

        inspector = PanddaInspect(
            working_directory = working_directory,
            files_dict = get_io_files(),
            mode = args.mode,
            )

        gui = inspector()

        splash.show_menu()

        splash.window.connect("destroy", gui.launch)
        # splash.window.connect("destroy_event", gui.launch)
        
        # gui.launch()

    except MissingFile as e: 

        modal_msg(
            'Missing files! Are you in the right directory?\n\nNo {}\n\npandda.inspect will now close'.format(
                str(e)
                )
            )

        sys.exit()

