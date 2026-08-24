# Organising projects: tags and groups

CCP4i2 has two ways to say "these projects belong together", and they are not
alternatives to each other. This document records what each is for, so that the
next person adding a grouping feature does not have to guess.

Short version:

* **`ProjectTag`** is a label. It has a name and nothing else. Use it to
  organise and find projects.
* **`ProjectGroup`** is a set that *owns something* — state and behaviour that
  depend on the group's type. Use it when the collection is itself a working
  object.

## The distinction that matters

Both models are many-to-many, and neither is containment: a project can carry
several tags and belong to several groups, and neither has any bearing on where
the project sits on disk. So the difference is not shape, it is payload.

A tag has no payload. [`ProjectTag`](../server/ccp4i2/db/models.py) is
`text`, an optional self-referential `parent`, and the M2M to projects. There is
nothing else it could carry, and that is the point — a label that accrued state
would stop being a label.

A group has a payload keyed on `type`. `ProjectGroup.sites` holds stored view
states for fragment campaigns; `ProjectGroupViewSet` adds `pandda_data`,
`export_pandda`, `summary_scene` and `parent_files`, all of which only mean
anything for a `fragment_set`. A group also distinguishes one member as
`MembershipType.PARENT` — the reference structure a campaign is measured
against.

|  | `ProjectTag` | `ProjectGroup` |
|---|---|---|
| Identity | its text (plus parent) | a globally unique `name` |
| Carries state | no | yes — `sites`, and whatever a future type needs |
| Type-dependent behaviour | no | yes — PanDDA export, summary scene, parent files |
| Distinguished member | no | yes — `MembershipType.PARENT` |
| Hierarchical | yes — `parent` FK, arbitrary depth | no — groups are flat |
| Exclusive | no | no |
| Has its own page | no — it filters the project list | yes — the campaigns UI |
| REST | `/projecttags/` (+ `tree/`, `{id}/add_projects/`, `{id}/remove_projects/`), `/projects/{id}/tags/`, `/projects/?tag=` | `/projectgroups/` + type-specific actions |

## Which one do I use?

Ask whether the collection owns anything. If the answer is "it is just a way of
finding these projects again", it is a tag. If the collection has state of its
own, a view, or an operation that acts on the set as a whole, it is a group.

A useful check: could you delete the collection and lose nothing but
organisation? Then it is a tag.

The promotion path runs one way. A tag that turns out to deserve state can
become a group; a group demoted to a tag has to throw its payload away. When in
doubt, start with a tag.

## Hierarchy lives on tags, not groups

`ProjectTag.parent` makes tags a forest, so the nested project browser —
`SARS-CoV-2 → Mpro → fragment-soaks → projects` — is built on tags. It is the
pane to the left of the project list: clicking a node filters, dropping a
selection onto a node files those projects there. `ProjectGroup` has no parent
FK and is deliberately flat: a campaign is not inside another campaign, and
groups do not appear in the tree.

A tag node is not a location. Because tags are non-exclusive a project can
appear under several nodes at once, which is what makes the tree useful for
discovery and is also why it can never be the answer to "where does this
project live?".

This is also why Qt-i2's "folders" belong with tags rather than groups. A legacy
folder was a parent `Project` row that organised the project list and did
nothing else — it held no state, had no view of its own, and could not even be
opened while it had children. It is a label that was implemented as a project,
and the resulting ambiguity (should the child project's directory live inside
the parent's?) is exactly why parent-project relationships were dropped from the
Django models. A tag node has no directory, so the question cannot be asked.

So `import_sqlite.py` maps a legacy folder to a tag, preserving nesting through
`ProjectTag.parent`. Two details follow from what a folder was:

* A folder with **no jobs of its own** is not imported as a project at all. It
  was never openable in Qt-i2 while it had children, so importing it would only
  add an empty project to the list. It survives as a tag.
* A folder that **does** hold jobs is a real project as well, so it is imported
  *and* tagged with its own node — which is where a user will look for it.

Each project is tagged with its immediate folder only. Ancestors are reached by
rolling up at read time (see below), never by writing implied tags into the
database.

> Before this, folders were mapped onto `ProjectGroup`s of type `fragment_set`
> named `Legacy_<name>_group` — wrong on both counts, and invisible except as
> spurious entries in the campaigns list. Databases migrated with an older
> build still carry those groups; nothing removes them automatically.

## Consequences worth knowing

**Tags are non-exclusive, folders were not.** A legacy project sat in exactly
one folder; a tagged project can sit in five. In a nested browser it therefore
appears in several branches. This is a deliberate relaxation rather than a bug,
but it means "where does this project live?" no longer has a single answer, and
the UI has to show a project's full set of tags rather than a location.

**Roll-up is a browser concern, not a storage one.** If the tree shows a project
tagged `Mpro` under its ancestor `SARS-CoV-2`, that is a descendant query, not a
second row. The stored M2M holds only what the user actually applied. Keep it
that way — implied ancestor tags in the database are impossible to un-imply
later. `ProjectTag.tagged_projects()`, `GET projecttags/tree/` and
`GET projects/?tag=<id>` all roll up by default; the tree reports
`project_count` and `total_project_count` separately so both readings are
available, and the totals are distinct-project counts rather than a sum of the
children's.

**Uniqueness is enforced by `path`, not by `unique_together`.**
`unique_together = ["parent", "text"]` does not constrain rows with
`parent = NULL`, because SQL treats NULLs as distinct — so nothing at the
database level prevented two top-level tags both called `SARS`. Django 5.2's
`nulls_distinct=False` is PostgreSQL-only and the desktop app ships SQLite, so
the portable fix is the unique materialised `path`, added in migration
`0017_projecttag_path`. That migration merges any duplicates it finds on the way
through: two nodes with the same ancestry and the same text are
indistinguishable to a user, so folding them together loses nothing, whereas
renaming one would invent a tag nobody asked for.

## Where the code is

| | |
|---|---|
| Models | [`server/ccp4i2/db/models.py`](../server/ccp4i2/db/models.py) — `ProjectTag`, `ProjectGroup`, `ProjectGroupMembership` |
| Tag API | [`server/ccp4i2/api/ProjectTagViewSet.py`](../server/ccp4i2/api/ProjectTagViewSet.py), plus `tags` actions on `ProjectViewSet` |
| Group API | [`server/ccp4i2/api/ProjectGroupViewSet.py`](../server/ccp4i2/api/ProjectGroupViewSet.py) — documented in the [Service Contract](CCP4I2_SERVICE_CONTRACT.md) |
| Tag UI | `client/renderer/components/project-tag-tree.tsx` (the browser), `tag-selection-dialog.tsx` (bulk tagging), `edit-tags.tsx`, `tags-of-project.tsx`, `project-tag-chips.tsx` |
| Group UI | `client/renderer/app/ccp4i2/(authed)/campaigns/` |
| Legacy migration | [`server/ccp4i2/db/import_sqlite.py`](../server/ccp4i2/db/import_sqlite.py) — `_import_tags`, `_import_projecttags`, `_create_project_groups` |
