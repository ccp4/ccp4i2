"""
Build a fragment-screening campaign from public data, for demos and teaching.

Why this exists: the interactive parts of fragment analysis -- tagging a
dataset at a site, recording a verdict, comparing siblings -- are awkward to
develop against a real campaign. Real campaigns are federated, so every
round trip is slow, and they hold live data that nobody wants to experiment
on. This command makes a throwaway campaign that behaves like a real one.

The datasets are a genuine BAZ2B bromodomain fragment series from the PDB:
same protein, same space group, cells within a couple of Angstroms, and a
DIFFERENT fragment bound in each. That is what makes it a useful fixture --
the siblings really are isomorphous and really do differ at the site, so
site-based comparison has something to compare.

Everything comes from PDBe, which serves both coordinates and structure
factors per entry. The structure factors arrive as mmCIF and are converted to
MTZ by the ordinary import path (gemmi), so no external conversion step is
needed. Note that PDB-REDO is NOT used: it was unreachable when this was
written, and PDBe carries everything required.

The jobs refine with DIMPLE: these are isomorphous crystals of a known
structure being re-refined against the campaign reference, which is what
dimple is for.

The SubstituteLigand jobs are created and fully parameterised but NOT run --
running them is slow and needs CCP4, and is left to the operator. Each job is
set to merged mode (OBSAS=MERGED), which is what lets this work without any
unmerged data: in merged mode the pipeline takes F_SIGF_IN and skips aimless
entirely, so no Diamond/Zenodo sweep is needed.

Usage:
    python manage.py make_demo_campaign
    python manage.py make_demo_campaign --name BAZ2B_teaching
    python manage.py make_demo_campaign --entries 5e9i,5dyu,5e9k
    python manage.py make_demo_campaign --dry-run
"""

import json
import logging
import shutil
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

logger = logging.getLogger(__name__)

PDBE = "https://www.ebi.ac.uk/pdbe/entry-files/download"

#: A real BAZ2B bromodomain fragment series. The first entry is the campaign
#: reference (the parent project); the rest become members. Each was checked to
#: serve both coordinates and structure factors, and to carry its own fragment.
DEFAULT_ENTRIES = ["5e9i", "5dyu", "5e9k", "5e9l", "5e9m", "5e9y"]

#: The acetyl-lysine pocket, the site these fragments bind. Filled in from the
#: reference structure at build time (see _site_origin_from_ligand), so the
#: marker lands on the real pocket rather than at a guessed coordinate.
FALLBACK_SITE_NAME = "Acetyl-lysine pocket"

#: Ligand codes that are solvent/cryoprotectant rather than the bound fragment.
NOT_A_FRAGMENT = {
    "HOH", "GOL", "EDO", "SO4", "PO4", "CL", "NA", "MG", "CA", "ZN",
    "DMS", "PEG", "PGE", "ACT", "MPD", "TRS", "IMD", "FMT", "NO3",
}


def coordinate_url(entry):
    return f"{PDBE}/pdb{entry}.ent"


def structure_factor_url(entry):
    return f"{PDBE}/r{entry}sf.ent"


def _download(url, destination):
    """Fetch one file, or raise CommandError naming the URL."""
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            with open(destination, "wb") as handle:
                shutil.copyfileobj(response, handle)
    except urllib.error.HTTPError as err:
        raise CommandError(f"{url} answered {err.code}")
    except (urllib.error.URLError, OSError) as err:
        raise CommandError(f"Could not download {url}: {err}")
    return destination


def _fragment_smiles(code):
    """The fragment's SMILES from PDBe's chemical component API, or None.

    SubstituteLigand needs to know WHAT to substitute, not just where. Without
    this the job is created with an empty SMILES field and cannot run -- the
    ligand panel shows no structure and the pipeline has nothing to build.

    Returns None rather than raising: a fragment whose SMILES cannot be found
    still makes a usable member project (its coordinates and reflections are
    fine), it just needs the ligand filling in by hand.
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{code}"
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            payload = json.load(response)
    except (urllib.error.URLError, OSError, ValueError):
        return None

    # PDBe returns each SMILES as {"program", "version", "name"} -- the string
    # itself is under "name", not "smiles", which is easy to get wrong.
    for entries in payload.values():
        for entry in entries:
            for item in entry.get("smiles", []):
                smiles = item.get("name")
                if smiles:
                    return smiles
    return None


def _fragment_code(pdb_path):
    """The bound fragment's residue code, or None.

    Picks the most common non-solvent HETATM residue, which for a fragment
    structure is the fragment. Solvent and cryoprotectants are excluded by
    name because they appear in nearly every entry and are never the point.
    """
    counts = {}
    for line in Path(pdb_path).read_text(errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        code = line[17:20].strip()
        if code and code not in NOT_A_FRAGMENT:
            counts[code] = counts.get(code, 0) + 1
    if not counts:
        return None
    return max(counts, key=counts.get)


def _site_origin_from_ligand(pdb_path, code):
    """Centroid of the named residue's atoms: where the site marker goes.

    Taken from the reference structure's own fragment, so the site sits in the
    real pocket. Returns None if the code is absent or the coordinates do not
    parse, and the caller then declines to create a site rather than placing
    one somewhere arbitrary.
    """
    xs, ys, zs = [], [], []
    for line in Path(pdb_path).read_text(errors="replace").splitlines():
        if not line.startswith("HETATM") or line[17:20].strip() != code:
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue
    if not xs:
        return None
    return [
        round(sum(xs) / len(xs), 3),
        round(sum(ys) / len(ys), 3),
        round(sum(zs) / len(zs), 3),
    ]


class Command(BaseCommand):
    help = (
        "Create a demo fragment-screening campaign from public PDBe data, with "
        "SubstituteLigand jobs configured but not run."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--name",
            default="BAZ2B_demo_campaign",
            help="Campaign (project group) name. Default: BAZ2B_demo_campaign",
        )
        parser.add_argument(
            "--entries",
            default=",".join(DEFAULT_ENTRIES),
            help=(
                "Comma-separated PDB entries. The FIRST becomes the campaign "
                "reference (parent project); the rest become members. "
                f"Default: {','.join(DEFAULT_ENTRIES)}"
            ),
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Check the data is fetchable and report the plan; create nothing.",
        )
        parser.add_argument(
            "--no-jobs",
            action="store_true",
            help="Create projects and the campaign, but no SubstituteLigand jobs.",
        )

    # -- helpers ---------------------------------------------------------

    def _import_into_job(self, job, object_path, path, description):
        """Import a local file into a job's project at ``object_path``.

        Uses the ordinary import path, which is what converts the structure
        factor mmCIF into an MTZ (gemmi) and registers the file with the
        project -- so the demo data arrives exactly as a user's would.
        """
        from ccp4i2.lib.utils.files.upload_param import (
            ImportSpec,
            _LocalPathUpload,
            import_file_for_param,
        )

        spec = ImportSpec(
            object_path=object_path,
            files=[_LocalPathUpload(path)],
            provenance_description=description,
        )
        return import_file_for_param(job, spec)

    def _set(self, job, path, value):
        """Set one parameter, or fail loudly.

        set_parameter returns a Result and does NOT raise, so an unchecked call
        leaves the job quietly holding its default -- which is how this command
        first shipped jobs whose OBSAS was still UNMERGED while the log claimed
        otherwise.
        """
        from ccp4i2.lib.utils.parameters.set_param import set_parameter

        result = set_parameter(job, path, value)
        if not getattr(result, "success", False):
            raise CommandError(
                f"Could not set {path} = {value!r} on job {job.number}: "
                f"{getattr(result, 'error', 'unknown error')}"
            )
        return result

    def _make_project(self, name, description):
        from ccp4i2.api.serializers import ProjectSerializer

        serializer = ProjectSerializer(data={"name": name, "description": description})
        serializer.is_valid(raise_exception=True)
        return serializer.save()

    # -- main ------------------------------------------------------------

    def handle(self, *args, **options):
        from ccp4i2.db import models
        from ccp4i2.lib.utils.jobs.create import create_job
        from ccp4i2.lib.utils.parameters.set_param import set_parameter

        entries = [e.strip().lower() for e in options["entries"].split(",") if e.strip()]
        if len(entries) < 2:
            raise CommandError(
                "Need at least two entries: one reference plus one member."
            )

        name = options["name"]
        reference, members = entries[0], entries[1:]

        self.stdout.write(f"Reference (parent): {reference}")
        self.stdout.write(f"Members:            {', '.join(members)}")

        # Not under --dry-run: a dry run is a check of the DATA SOURCES and
        # must work without a database (a fresh checkout, a CI box), so it
        # touches no tables.
        if not options["dry_run"] and models.ProjectGroup.objects.filter(name=name).exists():
            raise CommandError(
                f"A campaign called '{name}' already exists. "
                "Use --name, or delete it first."
            )

        work = Path(tempfile.mkdtemp(prefix="demo-campaign-"))
        try:
            # Fetch everything first: a half-built campaign is worse than none,
            # and this also makes --dry-run a genuine check of the sources.
            fetched = {}
            for entry in entries:
                coords = _download(coordinate_url(entry), work / f"{entry}.pdb")
                sfs = _download(structure_factor_url(entry), work / f"{entry}-sf.cif")
                code = _fragment_code(coords)
                smiles = _fragment_smiles(code) if code else None
                fetched[entry] = {
                    "coords": coords,
                    "sfs": sfs,
                    "fragment": code,
                    "smiles": smiles,
                }
                if code and smiles:
                    detail = f"fragment {code} ({smiles})"
                elif code:
                    detail = f"fragment {code} (no SMILES found)"
                else:
                    detail = "no fragment found"
                self.stdout.write(
                    f"  {entry}: coordinates + structure factors, {detail}"
                )

            site_origin = None
            reference_fragment = fetched[reference]["fragment"]
            if reference_fragment:
                site_origin = _site_origin_from_ligand(
                    fetched[reference]["coords"], reference_fragment
                )
            if site_origin:
                self.stdout.write(
                    f"  site '{FALLBACK_SITE_NAME}' at {site_origin} "
                    f"(centroid of {reference_fragment} in {reference})"
                )
            else:
                self.stdout.write(
                    self.style.WARNING(
                        "  no fragment found in the reference; no site will be created"
                    )
                )

            if options["dry_run"]:
                self.stdout.write(
                    self.style.SUCCESS(
                        "\nDry run: all data fetchable, nothing created."
                    )
                )
                return

            with transaction.atomic():
                group = models.ProjectGroup.objects.create(
                    name=name,
                    type=models.ProjectGroup.GroupType.FRAGMENT_SET,
                )

                parent_project = self._make_project(
                    f"{name}_{reference}",
                    f"Campaign reference: PDB {reference.upper()}",
                )
                models.ProjectGroupMembership.objects.create(
                    group=group,
                    project=parent_project,
                    type=models.ProjectGroupMembership.MembershipType.PARENT,
                )

                if site_origin:
                    models.CampaignSite.objects.create(
                        group=group,
                        name=FALLBACK_SITE_NAME,
                        origin_x=site_origin[0],
                        origin_y=site_origin[1],
                        origin_z=site_origin[2],
                        order=0,
                    )

                member_projects = []
                for entry in members:
                    project = self._make_project(
                        f"{name}_{entry}",
                        f"PDB {entry.upper()} "
                        f"(fragment {fetched[entry]['fragment'] or 'unknown'})",
                    )
                    models.ProjectGroupMembership.objects.create(
                        group=group,
                        project=project,
                        type=models.ProjectGroupMembership.MembershipType.MEMBER,
                    )
                    member_projects.append((entry, project))

            self.stdout.write(
                self.style.SUCCESS(
                    f"\nCampaign '{name}' created with "
                    f"{len(member_projects)} member(s)."
                )
            )

            if options["no_jobs"]:
                self.stdout.write("Skipping job creation (--no-jobs).")
                return

            # One SubstituteLigand per member, configured but not run.
            for entry, project in member_projects:
                job_id = create_job(
                    projectId=str(project.uuid),
                    taskName="SubstituteLigand",
                    title=f"SubstituteLigand {entry.upper()}",
                )
                job = models.Job.objects.get(uuid=job_id)

                # Merged mode: takes F_SIGF_IN directly and skips aimless, so
                # the demo needs no unmerged sweep from Diamond or Zenodo.
                self._set(job, "container.controlParameters.OBSAS", "MERGED")

                # Dimple, which is also the def.xml default now -- set
                # explicitly so the demo does not silently change route if
                # that default ever moves back.
                self._set(job, "container.inputData.PIPELINE", "DIMPLE")

                # What to substitute. Without these the job is created with an
                # empty SMILES field: the ligand panel shows no structure and
                # the pipeline has nothing to build, so the demo stops at the
                # first job. The fragment each entry actually contains is the
                # obvious thing to rebuild, and PDBe publishes its SMILES.
                smiles = fetched[entry]["smiles"]
                code = fetched[entry]["fragment"]
                if smiles:
                    self._set(job, "container.controlParameters.LIGANDAS", "SMILES")
                    # SMILESIN, not SMILES: both exist in the def.xml, but
                    # SMILESIN is what the GUI binds and what the pipeline
                    # reads (processInputFiles copies inputData.SMILESIN).
                    # SMILES is a legacy field carrying a hardcoded default,
                    # so writing to it leaves runTimeValidity reporting "a
                    # SMILES string has not been given" over a job that looks
                    # configured.
                    self._set(job, "container.inputData.SMILESIN", smiles)
                    # The three-letter code the rebuilt ligand is given. Using
                    # the entry's own code keeps the job self-describing;
                    # without one the pipeline defaults to DRG.
                    if code:
                        self._set(job, "container.controlParameters.LIGAND_CODE", code)
                else:
                    # Nothing to build from: say so rather than leaving a job
                    # that looks configured but cannot run.
                    self._set(job, "container.controlParameters.LIGANDAS", "NONE")

                self._import_into_job(
                    job,
                    "inputData.XYZIN",
                    fetched[reference]["coords"],
                    f"Campaign reference coordinates (PDB {reference.upper()})",
                )
                # Strip the reference's own fragment and waters.
                #
                # 5E9I is itself a fragment structure, so without this every
                # job starts from a model that already has a ligand sitting in
                # the pocket it is being asked to rebuild -- the density is
                # pre-explained, and what comes out says more about the start
                # point than about the dataset. XYZIN declares
                # ifAtomSelection, and its own tooltip recommends excluding
                # waters and ligands. 'protein' takes 5E9I from 1128 atoms to
                # 952, dropping F60 and the waters.
                self._set(job, "container.inputData.XYZIN.selection.text", "protein")
                self._import_into_job(
                    job,
                    "inputData.F_SIGF_IN",
                    fetched[entry]["sfs"],
                    f"Observed reflections (PDB {entry.upper()})",
                )

                # One free set for the whole campaign, taken from the
                # reference.
                #
                # Without this the pipeline warns that no free-R set was given
                # and nothing on the merged route generates one -- aimless
                # would, but aimless only runs on unmerged data. Refining every
                # member against its own deposited free set would be worse than
                # no set at all: R-free would not be comparable between
                # siblings, which is the whole point of a campaign.
                #
                # The deposited structure factors carry _refln.status, which
                # gemmi turns into a FreeR_flag column on import (833 free of
                # 16551 for 5E9I), so the reference's own file is the free set.
                self._import_into_job(
                    job,
                    "inputData.FREERFLAG_IN",
                    fetched[reference]["sfs"],
                    f"Shared campaign free-R set (PDB {reference.upper()})",
                )

                # These crystals are isomorphous but not identical -- the cells
                # differ by up to ~2 A across the series. aimless_pipe only
                # extends an input free set whose cell agrees to within
                # Clipper's 1 A default, so a shared set needs that check
                # relaxed. This is what PR #548 added for exactly this case.
                self._set(
                    job,
                    "container.controlParameters.OVERRIDE_CELL_DIFFERENCE",
                    True,
                )

                ligand_note = (
                    f"ligand {code} from SMILES"
                    if fetched[entry]["smiles"]
                    else "NO LIGAND -- set one by hand before running"
                )
                self.stdout.write(
                    f"  {entry}: SubstituteLigand job {job.number}, {ligand_note}"
                )

            self.stdout.write(
                self.style.SUCCESS(
                    "\nJobs are configured but NOT run. Run them from the UI, "
                    "or with: manage.py run_job <job-uuid>"
                )
            )

        finally:
            shutil.rmtree(work, ignore_errors=True)
