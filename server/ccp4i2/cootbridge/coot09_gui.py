# -*- coding: utf-8 -*-
"""Coot 0.9 GUI adapter for the CCP4i2 bridge (Python 2.7 / PyGTK2).

The GTK2 twin of coot1_gui: a "CCP4i2" menubar menu, insertions into
Coot's own File menu (as the legacy Qt integration did), and a
single-panel drill-down project browser - all over the same shared
data layer (api_client) and browse model as Coot 1.x.

Loaded by the coot_rebuild stub via imp.load_source (Coot 0.9's
embedded Python 2.7 cannot import the ccp4i2 package). Because an
imported module cannot see Coot's flat scripting namespace, every Coot
function it needs is handed in as the ``fns`` dict; anything missing
just disables the feature that needs it.

No f-strings, no pathlib: this runs on Python 2.7.
"""

from __future__ import absolute_import, print_function

import os
import traceback


def _log(message):
    print("[ccp4i2-cootbridge-0.9] {0}".format(message))


def install(fns, bridge_module, loader_module):
    """Entry point called from the stub: load the job's data, then build
    the menus. Returns the menu controller (kept alive by the caller)."""
    try:
        config = bridge_module.BridgeConfig()
        client = bridge_module.CootBridgeClient(config)
    except Exception:
        _log("failed to configure bridge:\n" + traceback.format_exc())
        return None

    # Initial data load (same load plan as Coot 1).
    try:
        if config.has_job_context():
            for item in bridge_module.load_plan(config, client):
                loader_module.load_item(item, fns)
    except Exception:
        _log("initial data load failed:\n" + traceback.format_exc())

    try:
        import gtk  # noqa: F401  (PyGTK2)
    except Exception:
        _log("PyGTK2 not available; menus/browser disabled (data still "
             "loaded)")
        return None

    try:
        controller = _MenuController(fns, bridge_module, loader_module,
                                     config, client)
        controller.install_menus()
        return controller
    except Exception:
        _log("menu installation failed:\n" + traceback.format_exc())
        return None


class _MenuController(object):
    def __init__(self, fns, bridge_module, loader_module, config, client):
        self.fns = fns
        self.bridge = bridge_module
        self.loader = loader_module
        self.config = config
        self.client = client
        self.browser = None

    def install_menus(self):
        menubar_menu = self.fns.get("coot_menubar_menu")
        add_item = self.fns.get("add_simple_coot_menu_menuitem")
        if menubar_menu is None or add_item is None:
            _log("coot_menubar_menu/add_simple_coot_menu_menuitem missing; "
                 "cannot build menus")
            return
        ccp4i2_menu = menubar_menu("CCP4i2")
        add_item(ccp4i2_menu, "Browse CCP4i2 projects...",
                 lambda *_a: self._open_browser())
        add_item(ccp4i2_menu, "Save model to CCP4i2 job...",
                 lambda *_a: self._save())
        # Legacy parity: also offer the save from Coot's File menu.
        try:
            file_menu = menubar_menu("File")
            add_item(file_menu, "Save to CCP4i2",
                     lambda *_a: self._save())
        except Exception:
            _log("could not extend the File menu:\n" + traceback.format_exc())
        _log("CCP4i2 menu installed")

    # -- saving -------------------------------------------------------------

    def _model_molecules(self):
        n_molecules = self.fns.get("graphics_n_molecules")
        valid = self.fns.get("is_valid_model_molecule")
        name_of = self.fns.get("molecule_name")
        molecules = []
        if n_molecules is None:
            return molecules
        for imol in range(n_molecules()):
            try:
                if valid is not None and not valid(imol):
                    continue
                label = name_of(imol) if name_of is not None else str(imol)
                molecules.append((imol, label))
            except Exception:
                continue
        return molecules

    def _save(self):
        chooser = self.fns.get("molecule_chooser_gui")
        molecules = self._model_molecules()
        if not molecules:
            _log("no model molecules to save")
            return
        if chooser is not None and len(molecules) > 1:
            chooser("Molecule to save to CCP4i2:",
                    lambda imol: self._save_molecule(imol))
        else:
            self._save_molecule(molecules[0][0])

    def _save_molecule(self, imol):
        save = self.fns.get("save_coordinates")
        name_of = self.fns.get("molecule_name")
        if save is None:
            _log("save_coordinates unavailable")
            return
        drop_dir = self.config.drop_dir or os.path.join(os.getcwd(),
                                                        "COOT_FILE_DROP")
        name = ""
        try:
            if name_of is not None:
                name = name_of(imol) or ""
        except Exception:
            pass
        extension = "cif" if name.lower().endswith((".cif", "(cif)")) \
            else "pdb"
        number = self.bridge.next_output_number(drop_dir)
        path = self.bridge.output_path(drop_dir, number, extension)
        save(imol, path)
        _log("saved molecule {0} as {1}".format(imol, path))

    # -- browser ------------------------------------------------------------

    def _open_browser(self):
        if self.browser is not None and self.browser.is_open():
            self.browser.present()
            return
        self.browser = _Browser09(self.fns, self.bridge, self.loader,
                                  self.client)


class _Browser09(object):
    """Single-panel drill-down (Projects -> Jobs -> Files) in PyGTK2,
    mirroring the Coot 1 browser. Fetches are synchronous (loopback HTTP
    is quick, and PyGTK2 threading is best avoided)."""

    def __init__(self, fns, bridge_module, loader_module, client):
        import gtk
        self.fns = fns
        self.bridge = bridge_module
        self.loader = loader_module
        self.client = client
        self.mode = "projects"
        self.projects = []
        self.job_rows = []
        self.project = None
        self.job = None
        self.row_data = []

        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_title("CCP4i2 projects")
        self.window.set_default_size(340, 480)
        self.window.connect("destroy", self._on_destroy)

        vbox = gtk.VBox(False, 4)
        vbox.set_border_width(6)

        header = gtk.HBox(False, 4)
        self.back_button = gtk.Button("<")
        self.back_button.connect("clicked", lambda *_a: self._go_back())
        self.crumb = gtk.Label("Projects")
        self.crumb.set_alignment(0.0, 0.5)
        header.pack_start(self.back_button, False, False, 0)
        header.pack_start(self.crumb, True, True, 0)
        vbox.pack_start(header, False, False, 0)

        scroller = gtk.ScrolledWindow()
        scroller.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.store = gtk.ListStore(str)
        self.tree = gtk.TreeView(self.store)
        self.tree.set_headers_visible(False)
        column = gtk.TreeViewColumn("", gtk.CellRendererText(), text=0)
        self.tree.append_column(column)
        self.tree.connect("row-activated", self._on_row_activated)
        scroller.add(self.tree)
        vbox.pack_start(scroller, True, True, 0)

        bar = gtk.HBox(True, 4)
        self.load_button = gtk.Button("Load file")
        self.load_button.connect("clicked", lambda *_a: self._load_selected())
        refresh = gtk.Button("Refresh")
        refresh.connect("clicked", lambda *_a: self._refresh())
        close = gtk.Button("Close")
        close.connect("clicked", lambda *_a: self.window.destroy())
        bar.pack_start(self.load_button, True, True, 0)
        bar.pack_start(refresh, True, True, 0)
        bar.pack_start(close, True, True, 0)
        vbox.pack_start(bar, False, False, 0)

        self.status = gtk.Label("")
        self.status.set_alignment(0.0, 0.5)
        vbox.pack_start(self.status, False, False, 0)

        self.window.add(vbox)
        self.window.show_all()
        self._show_projects()

    # -- lifecycle ----------------------------------------------------------

    def _on_destroy(self, *_args):
        self.window = None

    def is_open(self):
        return self.window is not None

    def present(self):
        if self.window is not None:
            self.window.present()

    # -- list helpers -------------------------------------------------------

    def _fill(self, labels_and_data):
        self.store.clear()
        self.row_data = []
        for text, payload in labels_and_data:
            self.store.append([text])
            self.row_data.append(payload)

    def _selected_payload(self):
        model, tree_iter = self.tree.get_selection().get_selected()
        if tree_iter is None:
            return None
        path = model.get_path(tree_iter)
        index = path[0]
        if 0 <= index < len(self.row_data):
            return self.row_data[index]
        return None

    def _say(self, message):
        self.status.set_text(message)

    def _set_mode(self, mode, crumb):
        self.mode = mode
        self.crumb.set_text(crumb)
        self.back_button.set_sensitive(mode != "projects")
        self.load_button.set_label(
            "Load job" if mode == "jobs" else "Load file")

    # -- navigation ---------------------------------------------------------

    def _show_projects(self):
        self._say("Fetching projects...")
        try:
            self.projects = sorted(
                self.client.projects() or [],
                key=lambda p: (p.get("name") or "").lower())
        except Exception:
            self._say("Could not list projects (is the CCP4i2 server "
                      "running?)")
            _log(traceback.format_exc())
            return
        self._set_mode("projects", "Projects")
        self._fill([((p.get("name") or str(p.get("id"))) + "  >", p)
                    for p in self.projects])
        self._say("{0} projects".format(len(self.projects)))

    def _show_jobs(self, project):
        self.project = project
        self._say("Fetching jobs...")
        try:
            self.job_rows = self.bridge.browse_model(
                self.client.job_tree(project["id"]))
        except Exception:
            self._say("Could not fetch the job tree")
            _log(traceback.format_exc())
            return
        self._set_mode("jobs", project.get("name") or "Jobs")
        self._fill([
            ("{0}{1}  >".format("  " * job["depth"], job["label"]), job)
            for job in self.job_rows])
        self._say("{0} jobs with loadable files".format(len(self.job_rows)))

    def _show_files(self, job):
        self.job = job
        self._set_mode("files", "{0} / {1}".format(
            self.project.get("name") or "", job["label"]))
        self._fill([("{0}  [{1}]".format(f["label"], f["kind"]), f)
                    for f in job["files"]])
        self._say("{0} files - activate one to load it".format(
            len(job["files"])))

    def _go_back(self):
        if self.mode == "files":
            self._set_mode("jobs", self.project.get("name") or "Jobs")
            self._fill([
                ("{0}{1}  >".format("  " * job["depth"], job["label"]), job)
                for job in self.job_rows])
        elif self.mode == "jobs":
            self._set_mode("projects", "Projects")
            self._fill([((p.get("name") or str(p.get("id"))) + "  >", p)
                        for p in self.projects])

    def _refresh(self):
        if self.mode == "projects":
            self._show_projects()
        elif self.mode == "jobs" and self.project is not None:
            self._show_jobs(self.project)
        elif self.mode == "files" and self.job is not None:
            self._show_files(self.job)

    def _on_row_activated(self, _tree, _path, _column):
        payload = self._selected_payload()
        if payload is None:
            return
        if self.mode == "projects":
            self._show_jobs(payload)
        elif self.mode == "jobs":
            self._show_files(payload)
        elif self.mode == "files":
            self._load_files([payload])

    def _load_selected(self):
        payload = self._selected_payload()
        if payload is None:
            self._say("Select a row first")
            return
        if self.mode == "jobs":
            self._load_files(payload["files"])
        elif self.mode == "files":
            self._load_files([payload])

    def _load_files(self, records):
        self._say("Loading...")
        own_project = self.client.config.project_name
        source_project = self.project.get("name") if self.project else None
        loaded = 0
        for record in records:
            try:
                path = self.client.file_local_path(record["file_id"])
                if path is None:
                    destination = os.path.join(
                        self.client.config.cache_dir(),
                        "file_{0}".format(record["file_id"]))
                    path = self.client.download_file(record["file_id"],
                                                     destination)
                self.loader.load_item(
                    {"kind": record["kind"], "path": path,
                     "label": self.bridge.display_label(
                         record["label"], source_project, own_project)},
                    self.fns)
                loaded += 1
            except Exception:
                _log("load failed:\n" + traceback.format_exc())
        self._say("Loaded {0} file(s)".format(loaded))
