"""Coot 1.x adapter for the CCP4i2 bridge (Python 3 / GTK4).

Loaded into Coot's embedded interpreter by the stub script the coot1
task writes (``--script``). Given only the environment contract (see
api_client), it:

* loads the launching job's input data (fetched via the load plan);
* installs a "CCP4i2" menu with save-to-job and browse actions;
* provides a project-hierarchy browser (projects -> jobs -> files)
  driven entirely by the shared browse model.

UI surface is deliberately restricted to what is verified reliable in
Coot 1.3.x: ``coot_gui.attach_module_menu_button`` + ``Gio.SimpleAction``
for the menu (with a raw-Gtk fallback), plain ``Gtk.Window``/``Gtk.ListBox``
for the browser, and ``GLib.idle_add`` to marshal worker-thread results
back to the GUI thread. All Coot calls happen on the GUI thread.

Import of this module must never crash Coot: ``start()`` wraps every
stage and prints failures to the terminal instead of raising.
"""

import os
import sys
import threading
import traceback

from ccp4i2.cootbridge import api_client


def _log(message):
    print("[ccp4i2-cootbridge] {0}".format(message))


# ---------------------------------------------------------------------------
# Entry point (called from the generated stub script at Coot startup)
# ---------------------------------------------------------------------------

_state = {"config": None, "client": None}


def start():
    try:
        config = api_client.BridgeConfig()
        client = api_client.CootBridgeClient(config)
        _state["config"] = config
        _state["client"] = client
        _log("API: {0}  job: {1}".format(config.api_url,
                                         config.job_id or config.job_uuid))
    except Exception:
        _log("failed to configure bridge:\n" + traceback.format_exc())
        return
    try:
        if config.has_job_context():
            for item in api_client.load_plan(config, client):
                _load_item(item)
    except Exception:
        _log("initial data load failed:\n" + traceback.format_exc())
    try:
        _install_menu()
    except Exception:
        _log("menu installation failed:\n" + traceback.format_exc())


# ---------------------------------------------------------------------------
# Load dispatch: browse-model kinds -> Coot 1 API calls
# ---------------------------------------------------------------------------


def _coot():
    import coot
    return coot


def _load_item(item):
    """Load one {kind, path, label} item. GUI thread only."""
    coot = _coot()
    kind, path, label = item["kind"], item["path"], item.get("label")
    _log("loading {0}: {1}".format(kind, path))
    if kind == "coordinates":
        imol = coot.read_pdb(path)
        if label and imol >= 0 and hasattr(coot, "set_molecule_name"):
            coot.set_molecule_name(imol, label)
    elif kind == "map_2fofc":
        coot.read_mtz(path, "F", "PHI", "", False, False)
    elif kind == "map_fofc":
        coot.read_mtz(path, "F", "PHI", "", False, True)
    elif kind == "map_anom":
        imap = coot.read_mtz(path, "F", "PHI", "", False, True)
        if imap >= 0 and hasattr(coot, "set_map_colour"):
            coot.set_map_colour(imap, 0.75, 0.9, 0.75)
    elif kind == "map":
        if hasattr(coot, "handle_read_ccp4_map"):
            coot.handle_read_ccp4_map(path, 0)
        elif hasattr(coot, "read_ccp4_map"):
            coot.read_ccp4_map(path, 0)
    elif kind == "dictionary":
        coot.read_cif_dictionary(path)
    else:
        _log("unknown kind {0!r}, skipped".format(kind))


# ---------------------------------------------------------------------------
# Saving back to the job (the COOT_FILE_DROP harvest contract)
# ---------------------------------------------------------------------------


def _model_molecules():
    coot = _coot()
    molecules = []
    for imol in range(coot.graphics_n_molecules()):
        try:
            if coot.is_valid_model_molecule(imol):
                molecules.append((imol, coot.molecule_name(imol)))
        except Exception:
            continue
    return molecules


def _save_molecule(imol):
    coot = _coot()
    config = _state["config"]
    drop_dir = config.drop_dir if config else None
    if not drop_dir:
        drop_dir = os.path.join(os.getcwd(), "COOT_FILE_DROP")
    name = ""
    try:
        name = coot.molecule_name(imol) or ""
    except Exception:
        pass
    extension = "cif" if name.lower().endswith((".cif", "(cif)")) else "pdb"
    number = api_client.next_output_number(drop_dir)
    path = api_client.output_path(drop_dir, number, extension)
    coot.save_coordinates(imol, path)
    _log("saved molecule {0} as {1}".format(imol, path))
    if hasattr(coot, "add_status_bar_text"):
        coot.add_status_bar_text(
            "Saved to CCP4i2 job as output{0}.{1}".format(number, extension))


def _on_save_action():
    """One model -> save it; several -> chooser dialog."""
    molecules = _model_molecules()
    if not molecules:
        _log("no model molecules to save")
        return
    if len(molecules) == 1:
        _save_molecule(molecules[0][0])
        return
    _molecule_chooser(molecules)


def _molecule_chooser(molecules):
    from gi.repository import Gtk
    window = Gtk.Window(title="Save molecule to CCP4i2")
    window.set_default_size(420, 300)
    box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
    box.set_margin_top(8)
    box.set_margin_bottom(8)
    box.set_margin_start(8)
    box.set_margin_end(8)
    listbox = Gtk.ListBox()
    for imol, name in molecules:
        row = Gtk.ListBoxRow()
        label = Gtk.Label(label="{0}: {1}".format(imol, name), xalign=0)
        row.set_child(label)
        row.imol = imol
        listbox.append(row)
    scroller = Gtk.ScrolledWindow()
    scroller.set_child(listbox)
    scroller.set_vexpand(True)
    box.append(scroller)
    button = Gtk.Button(label="Save to CCP4i2 job")

    def on_clicked(_button):
        row = listbox.get_selected_row()
        if row is not None:
            _save_molecule(row.imol)
            window.close()

    button.connect("clicked", on_clicked)
    box.append(button)
    window.set_child(box)
    window.show()


# ---------------------------------------------------------------------------
# Menu installation
# ---------------------------------------------------------------------------


def _install_menu():
    """A "CCP4i2" toolbar menu-button. Prefers coot_gui's blessed helper;
    falls back to building the same structure by hand."""
    try:
        import coot_gui
        menu = coot_gui.attach_module_menu_button("CCP4i2")
        coot_gui.add_simple_action_to_menu(
            menu, "Save model to CCP4i2 job...", "ccp4i2_save",
            lambda *_args: _on_save_action())
        coot_gui.add_simple_action_to_menu(
            menu, "Browse CCP4i2 projects...", "ccp4i2_browse",
            lambda *_args: _open_browser())
        _log("CCP4i2 menu installed (coot_gui)")
        return
    except Exception:
        _log("coot_gui menu helpers unavailable, using fallback:\n" +
             traceback.format_exc())
    import coot_gui_api
    from gi.repository import Gio, Gtk
    menu = Gio.Menu.new()
    application = coot_gui_api.application()
    for label, action_name, callback in (
            ("Save model to CCP4i2 job...", "ccp4i2_save",
             lambda *_args: _on_save_action()),
            ("Browse CCP4i2 projects...", "ccp4i2_browse",
             lambda *_args: _open_browser())):
        action = Gio.SimpleAction.new(action_name, None)
        action.connect("activate", lambda _a, _p, cb=callback: cb())
        application.add_action(action)
        menu.append(label, "app." + action_name)
    button = Gtk.MenuButton(label="CCP4i2")
    button.set_popover(Gtk.PopoverMenu.new_from_model(menu))
    coot_gui_api.main_toolbar().append(button)
    _log("CCP4i2 menu installed (fallback)")


# ---------------------------------------------------------------------------
# Worker-thread plumbing
# ---------------------------------------------------------------------------


# Background threads are unreliable in Coot's embedded interpreter (Coot 1's
# own thread-using GUI modules are broken/unported, and its native remote
# paths all run Python on the GUI thread). Loopback HTTP is milliseconds, so
# the browser fetches synchronously by default; flip this once threading
# inside Coot is verified.
USE_THREADS = False

# Dock the browser into Coot's main window rather than floating it. Off by
# default: on macOS the GTK-in-Coot window resizes unreliably, so an
# in-window dock ends up coupled to the main window's geometry and may not
# appear until the user manually grows the window. A floating window opens
# at its own size, travels with Coot (transient-for), and just works. The
# paned dock (see _try_dock) is kept for Linux/experimentation.
DOCK_BROWSER = False


def _in_background(work, on_done):
    """Run ``work()`` and deliver (result, error) to ``on_done`` on the
    GUI thread. Synchronous by default (see USE_THREADS above)."""
    if not USE_THREADS:
        result, error = None, None
        try:
            _log("fetch: {0}".format(getattr(work, "__name__", "work")))
            result = work()
        except Exception:
            error = traceback.format_exc()
        _deliver(on_done, result, error)
        return

    from gi.repository import GLib

    def runner():
        result, error = None, None
        try:
            _log("worker fetch: {0}".format(getattr(work, "__name__", "work")))
            result = work()
        except Exception:
            error = traceback.format_exc()
        GLib.idle_add(_deliver, on_done, result, error)

    threading.Thread(target=runner, daemon=True).start()


def _deliver(on_done, result, error):
    try:
        on_done(result, error)
    except Exception:
        _log("callback failed:\n" + traceback.format_exc())
    return False  # one-shot idle handler


# ---------------------------------------------------------------------------
# The project-hierarchy browser
# ---------------------------------------------------------------------------

_browser = None


def _open_browser():
    """Toggle the browser: docked panel preferred, floating fallback."""
    global _browser
    if _browser is not None:
        if _browser.is_open():
            _browser.close()
            _browser = None
            return
        _browser = None
    _browser = _Browser(_state["client"])


def _close_browser():
    global _browser
    if _browser is not None:
        _browser.close()
        _browser = None


class _Browser(object):
    """Moorhen-style single-level drill-down over the shared browse
    model: one list shows the current level (Projects, then a project's
    jobs, then a job's files), with a back arrow and breadcrumb.
    Activating a row (double-click / Enter) descends - or, at the file
    level, loads the file.

    Docks as a slim column in Coot's main_window_hbox, inserted right
    after the vertical modelling toolbar so it reads as part of it.
    (Appending *inside* the toolbar column was tried and overflows the
    screen: the buttons already consume the full window height, so any
    addition below them grows the window's minimum size. Stacking there
    needs Coot's own sidebar to become scrollable - an upstream ask.)
    Falls back to a floating window when the dock point is unavailable.
    """

    def __init__(self, client, prefer_dock=None):
        self.client = client
        self.mode = "projects"    # "projects" | "jobs" | "files"
        self.projects = []
        self.job_rows = []        # browse_model output for current project
        self.project = None       # selected project dict
        self.job = None           # selected job row (with files)
        self.window = None
        self.panel = None
        self._dock_parent = None
        self._dock_kind = None

        if prefer_dock is None:
            prefer_dock = DOCK_BROWSER
        docked = prefer_dock and self._try_dock()
        if not docked:
            self._build_floating()
        self._show_projects()

    # -- construction -------------------------------------------------------

    def _build_content(self):
        from gi.repository import Gtk
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
        for edge in ("top", "bottom", "start", "end"):
            getattr(outer, "set_margin_" + edge)(6)

        header = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=4)
        self.back_button = Gtk.Button(label="←")
        self.back_button.set_has_frame(False)
        self.back_button.set_tooltip_text("Back")
        self.back_button.connect("clicked", lambda *_: self._go_back())
        self.crumb = Gtk.Label(label="CCP4i2", xalign=0)
        self.crumb.set_hexpand(True)
        self._ellipsize(self.crumb)
        close = Gtk.Button(label="✕")
        close.set_has_frame(False)
        close.set_tooltip_text("Close browser")
        close.connect("clicked", lambda *_: _close_browser())
        header.append(self.back_button)
        header.append(self.crumb)
        header.append(close)
        outer.append(header)

        self.listbox = Gtk.ListBox()
        self.listbox.connect("row-activated",
                             lambda _l, row: self._row_activated(row))
        scroller = Gtk.ScrolledWindow()
        scroller.set_child(self.listbox)
        scroller.set_vexpand(True)
        # Keep the minimum small: when docked under the modelling
        # toolbar, every pixel of minimum height here adds to the main
        # window's minimum height. vexpand grows it into any slack.
        scroller.set_min_content_height(100)
        outer.append(scroller)

        bar = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=4)
        self.load_button = Gtk.Button(label="Load")
        self.load_button.connect("clicked", lambda *_: self._load_clicked())
        self.load_all_button = Gtk.Button(label="Load all")
        self.load_all_button.connect(
            "clicked", lambda *_: self._load_all_clicked())
        refresh = Gtk.Button(label="Refresh")
        refresh.connect("clicked", lambda *_: self._refresh())
        bar.append(self.load_button)
        bar.append(self.load_all_button)
        bar.append(refresh)
        outer.append(bar)

        self.status = Gtk.Label(label="", xalign=0)
        self.status.set_wrap(True)
        outer.append(self.status)
        return outer

    @staticmethod
    def _ellipsize(label):
        try:
            from gi.repository import Pango
            label.set_ellipsize(Pango.EllipsizeMode.END)
        except Exception:
            pass

    @staticmethod
    def _buildable_id(widget):
        try:
            return widget.get_buildable_id()
        except AttributeError:
            try:
                from gi.repository import Gtk
                return Gtk.Buildable.get_buildable_id(widget)
            except Exception:
                return None

    def _try_dock(self):
        """Dock the browser under the modelling toolbar as the lower half
        of a vertical Gtk.Paned. True on success.

        The toolbar column (main_window_vbox_inner) is not scrollable, so
        *appending* into it raises the column's minimum height and GTK
        grows the whole window to satisfy it - the "alarming jump". A
        Paned negotiates by divider position instead, and a shrinkable
        end child can be sized below its content's minimum, so the window
        never has to grow. The user also gets a draggable divider to size
        or collapse the browser. Only the small toolbar box is reparented
        (never the GL canvas), so the graphics context is untouched.
        Falls back to a plain side column, then to a floating window.
        """
        try:
            import coot_gui_api
            from gi.repository import Gtk
            hbox = coot_gui_api.main_hbox()
            panel = self._build_content()

            sidebar = None
            child = hbox.get_first_child()
            while child is not None:
                if self._buildable_id(child) == "main_window_vbox_inner":
                    sidebar = child
                    break
                child = child.get_next_sibling()

            if sidebar is not None and hasattr(Gtk, "Paned"):
                paned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
                # Swap the toolbar's slot in the hbox for the paned, then
                # move the toolbar into the paned's top half.
                sib = sidebar.get_next_sibling()
                hbox.remove(sidebar)
                paned.set_start_child(sidebar)
                paned.set_resize_start_child(False)
                paned.set_shrink_start_child(False)
                paned.set_end_child(panel)
                paned.set_resize_end_child(True)
                # Shrinkable end child = the browser can be sized below its
                # own minimum, so it never forces the window taller.
                paned.set_shrink_end_child(True)
                panel.set_size_request(210, -1)
                if sib is not None:
                    hbox.insert_child_after(paned, sib)
                else:
                    hbox.append(paned)
                # Give the toolbar its natural height, the browser the rest.
                try:
                    natural = sidebar.measure(Gtk.Orientation.VERTICAL, -1)[1]
                    if natural and natural > 0:
                        paned.set_position(natural)
                except Exception:
                    pass
                self.panel = panel
                self._dock_parent = paned
                self._dock_kind = "paned"
                self._paned_reparented = sidebar
                self._paned_hbox = hbox
                _log("browser docked under the toolbar (paned)")
                return True

            # Fallback: a plain slim side column (stable, but squeezes the
            # graphics rather than sharing the toolbar column).
            panel.set_size_request(250, -1)
            hbox.append(panel)
            self.panel = panel
            self._dock_parent = hbox
            self._dock_kind = "column"
            _log("browser docked as a side column")
            return True
        except Exception:
            _log("could not dock into main window, floating instead:\n" +
                 traceback.format_exc())
            return False

    def _build_floating(self):
        from gi.repository import Gtk
        self.window = Gtk.Window(title="CCP4i2 projects")
        self.window.set_default_size(360, 560)
        self.window.set_child(self._build_content())
        # Tie the window to Coot's main window so it stays above it and
        # travels with the app rather than becoming a stray top-level.
        try:
            import coot_gui_api
            main_window = coot_gui_api.application().get_active_window()
            if main_window is not None:
                self.window.set_transient_for(main_window)
                self.window.set_destroy_with_parent(True)
        except Exception:
            pass
        self.window.show()

    # -- lifecycle ----------------------------------------------------------

    def is_open(self):
        if self.panel is not None:
            return self.panel.get_parent() is not None
        if self.window is not None:
            return self.window.get_visible()
        return False

    def close(self):
        if getattr(self, "_dock_kind", None) == "paned":
            # Undo the reparenting: return the toolbar to the hbox slot the
            # paned occupies, then drop the paned.
            try:
                from gi.repository import Gtk  # noqa: F401
                paned = self._dock_parent
                toolbar = self._paned_reparented
                hbox = self._paned_hbox
                after = paned.get_next_sibling()
                paned.set_start_child(None)
                paned.set_end_child(None)
                hbox.remove(paned)
                if after is not None:
                    hbox.insert_child_after(toolbar, after)
                else:
                    hbox.append(toolbar)
            except Exception:
                _log(traceback.format_exc())
            self.panel = None
            self._dock_kind = None
            return
        if self.panel is not None and self._dock_parent is not None:
            try:
                self._dock_parent.remove(self.panel)
            except Exception:
                _log(traceback.format_exc())
            self.panel = None
            return
        if self.window is not None:
            self.window.close()

    # -- list helpers -------------------------------------------------------

    def _fill(self, labels_and_data):
        from gi.repository import Gtk
        while True:
            row = self.listbox.get_row_at_index(0)
            if row is None:
                break
            self.listbox.remove(row)
        for text, payload in labels_and_data:
            row = Gtk.ListBoxRow()
            label = Gtk.Label(label=text, xalign=0)
            self._ellipsize(label)
            row.set_child(label)
            row.payload = payload
            self.listbox.append(row)

    def _say(self, message):
        self.status.set_text(message)

    def _set_mode(self, mode, crumb):
        self.mode = mode
        self.crumb.set_text(crumb)
        self.back_button.set_sensitive(mode != "projects")
        self.load_button.set_visible(mode in ("jobs", "files"))
        self.load_button.set_label(
            "Load job" if mode == "jobs" else "Load file")
        self.load_all_button.set_visible(mode == "files")

    # -- navigation ---------------------------------------------------------

    def _show_projects(self, refetch=True):
        if refetch or not self.projects:
            self._say("Fetching projects...")
            _in_background(self.client.projects, self._projects_arrived)
        else:
            self._projects_arrived(self.projects, None)

    def _projects_arrived(self, result, error):
        if error:
            self._say("Could not list projects (is the CCP4i2 server "
                      "running?)")
            _log(error)
            return
        try:
            self.projects = sorted(
                result or [], key=lambda p: (p.get("name") or "").lower())
            self._set_mode("projects", "Projects")
            self._fill([((p.get("name") or str(p.get("id"))) + "  ▸", p)
                        for p in self.projects])
            self._say("{0} projects".format(len(self.projects)))
        except Exception:
            self._say("Browser error rendering projects - see terminal")
            _log(traceback.format_exc())

    def _show_jobs(self, project):
        self.project = project
        self._say("Fetching jobs...")
        _in_background(
            lambda: api_client.browse_model(
                self.client.job_tree(project["id"])),
            self._jobs_arrived)

    def _jobs_arrived(self, result, error):
        if error:
            self._say("Could not fetch the job tree")
            _log(error)
            return
        try:
            self.job_rows = result or []
            self._set_mode("jobs", self.project.get("name") or "Jobs")
            self._fill([
                ("{0}{1}  ▸".format("  " * job["depth"], job["label"]),
                 job) for job in self.job_rows])
            self._say("{0} jobs with loadable files".format(
                len(self.job_rows)))
        except Exception:
            self._say("Browser error rendering jobs - see terminal")
            _log(traceback.format_exc())

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
            self._jobs_arrived(self.job_rows, None)
        elif self.mode == "jobs":
            self._show_projects(refetch=False)

    def _refresh(self):
        if self.mode == "projects":
            self._show_projects(refetch=True)
        elif self.mode == "jobs" and self.project is not None:
            self._show_jobs(self.project)
        elif self.mode == "files" and self.job is not None:
            self._show_files(self.job)

    def _row_activated(self, row):
        if row is None:
            return
        if self.mode == "projects":
            self._show_jobs(row.payload)
        elif self.mode == "jobs":
            self._show_files(row.payload)
        elif self.mode == "files":
            self._resolve_and_load([row.payload])

    # -- loading ------------------------------------------------------------

    def _load_clicked(self):
        row = self.listbox.get_selected_row()
        if row is None:
            self._say("Select a row first")
            return
        if self.mode == "jobs":
            self._resolve_and_load(row.payload["files"])
        elif self.mode == "files":
            self._resolve_and_load([row.payload])

    def _load_all_clicked(self):
        if self.mode == "files" and self.job is not None:
            self._resolve_and_load(self.job["files"])

    def _resolve_and_load(self, files):
        """Resolve paths, then load. Runs via _in_background (synchronous
        by default; see USE_THREADS)."""
        client = self.client

        own_project = client.config.project_name
        source_project = self.project.get("name") if self.project else None

        def work():
            plan = []
            for record in files:
                path = client.file_local_path(record["file_id"])
                if path is None:
                    destination = os.path.join(
                        client.config.cache_dir(),
                        "file_{0}".format(record["file_id"]))
                    path = client.download_file(record["file_id"],
                                                destination)
                plan.append({
                    "kind": record["kind"], "path": path,
                    "label": api_client.display_label(
                        record["label"], source_project, own_project)})
            return plan

        def done(plan, error):
            if error:
                self._say("Load failed")
                _log(error)
                return
            for item in plan:
                try:
                    _load_item(item)
                except Exception:
                    _log("load failed:\n" + traceback.format_exc())
            self._say("Loaded {0} file(s)".format(len(plan)))

        self._say("Loading...")
        _in_background(work, done)
