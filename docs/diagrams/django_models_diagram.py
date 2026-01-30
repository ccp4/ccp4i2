#!/usr/bin/env python3
"""
CCP4i2 Django Data Models Diagram
Generates a professional ER diagram showing all Django models and their relationships.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, ConnectionPatch
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(24, 18))
ax.set_xlim(0, 24)
ax.set_ylim(0, 18)
ax.set_aspect('equal')
ax.axis('off')

# Color palette
colors = {
    'registry': '#3B82F6',      # Blue - Compound Registry
    'assays': '#10B981',        # Green - Assays
    'constructs': '#8B5CF6',    # Purple - Constructs
    'server': '#F59E0B',        # Orange - Server/Project
    'users': '#EF4444',         # Red - Users
    'reference': '#6B7280',     # Gray - Reference data
    'text_dark': '#1F2937',
    'white': '#FFFFFF',
    'section_registry': '#DBEAFE',
    'section_assays': '#D1FAE5',
    'section_constructs': '#EDE9FE',
    'section_server': '#FEF3C7',
    'arrow': '#374151',
}

def draw_entity(ax, x, y, width, height, color, title, fields, pk_field=None):
    """Draw an entity box with fields"""
    # Main box
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.05",
                         facecolor='white', edgecolor=color,
                         linewidth=2)
    ax.add_patch(box)

    # Title bar
    title_height = 0.4
    title_box = FancyBboxPatch((x, y + height - title_height), width, title_height,
                                boxstyle="round,pad=0.01,rounding_size=0.05",
                                facecolor=color, edgecolor='none')
    ax.add_patch(title_box)

    # Title text
    ax.text(x + width/2, y + height - title_height/2, title,
            fontsize=8, ha='center', va='center', color='white', fontweight='bold')

    # Fields
    field_y = y + height - title_height - 0.25
    for field in fields:
        if field == pk_field:
            ax.text(x + 0.1, field_y, f'{field}', fontsize=6.5, va='center',
                    color=colors['text_dark'], fontweight='bold')
        elif field.startswith('FK:'):
            ax.text(x + 0.1, field_y, f'> {field[3:]}', fontsize=6.5, va='center',
                    color=color, fontweight='bold')
        else:
            ax.text(x + 0.1, field_y, field, fontsize=6.5, va='center',
                    color=colors['text_dark'])
        field_y -= 0.22

def draw_section(ax, x, y, width, height, color, title):
    """Draw a section background"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.1",
                         facecolor=color, edgecolor='#9CA3AF',
                         linewidth=1, alpha=0.4)
    ax.add_patch(box)
    ax.text(x + 0.15, y + height - 0.2, title, fontsize=11,
            ha='left', va='top', color=colors['text_dark'], fontweight='bold')

def draw_relationship(ax, start, end, label=None, style='1:N'):
    """Draw a relationship line between entities"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='->', color=colors['arrow'],
                                lw=1.2, connectionstyle='arc3,rad=0.1'))
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(mid_x, mid_y + 0.15, label, fontsize=6, ha='center',
                color=colors['arrow'], style='italic')

# Title
ax.text(12, 17.6, 'CCP4i2 Django Data Models', fontsize=20,
        ha='center', va='center', fontweight='bold', color=colors['text_dark'])
ax.text(12, 17.2, 'Entity Relationship Diagram', fontsize=12,
        ha='center', va='center', color='#6B7280')

# ==================== COMPOUND REGISTRY SECTION ====================

draw_section(ax, 0.2, 10.5, 7.5, 6.5, colors['section_registry'], 'Compound Registry')

# Supplier
draw_entity(ax, 0.4, 15.2, 2.2, 1.5, colors['registry'], 'Supplier',
            ['id (UUID)', 'name', 'initials', 'FK:user'], pk_field='id (UUID)')

# Target
draw_entity(ax, 2.9, 15.2, 2.2, 1.7, colors['registry'], 'Target',
            ['id (UUID)', 'name', 'FK:parent', 'image', 'saved_view'], pk_field='id (UUID)')

# Compound
draw_entity(ax, 0.4, 12.5, 3.5, 2.4, colors['registry'], 'Compound',
            ['id (UUID)', 'reg_number', 'FK:target', 'FK:supplier',
             'smiles', 'rdkit_smiles', 'molecular_weight',
             'stereo_comment', 'aliases[]'], pk_field='id (UUID)')

# MolecularProperties
draw_entity(ax, 4.2, 12.5, 2.8, 2.2, colors['registry'], 'MolecularProperties',
            ['FK:compound (PK)', 'molecular_weight', 'hbd, hba',
             'clogp, tpsa', 'rotatable_bonds', 'fraction_sp3'], pk_field='FK:compound (PK)')

# Batch
draw_entity(ax, 0.4, 10.7, 3.2, 1.5, colors['registry'], 'Batch',
            ['id (UUID)', 'FK:compound', 'batch_number',
             'amount', 'salt_code'], pk_field='id (UUID)')

# BatchQCFile
draw_entity(ax, 4, 10.7, 2.5, 1.2, colors['registry'], 'BatchQCFile',
            ['id (UUID)', 'FK:batch', 'file'], pk_field='id (UUID)')

# CompoundTemplate
draw_entity(ax, 5.4, 15.2, 2.2, 1.4, colors['registry'], 'CompoundTemplate',
            ['id (UUID)', 'FK:target', 'name', 'mol2d'], pk_field='id (UUID)')

# Registry Relationships
draw_relationship(ax, (2.9, 13.5), (3.9, 13.5), 'has')  # Compound -> MolecularProperties
draw_relationship(ax, (1.9, 12.5), (1.9, 12.2), 'has')  # Compound -> Batch
draw_relationship(ax, (3.6, 11.3), (4, 11.3), 'has')    # Batch -> BatchQCFile
draw_relationship(ax, (2.1, 15), (2.1, 14.9), 'belongs')  # Compound -> Target
draw_relationship(ax, (1.2, 14.9), (1.2, 15.2), 'from')   # Compound -> Supplier

# ==================== ASSAYS SECTION ====================

draw_section(ax, 7.9, 10.5, 7.8, 6.5, colors['section_assays'], 'Assays')

# Protocol
draw_entity(ax, 8.1, 15, 2.8, 1.8, colors['assays'], 'Protocol',
            ['id (UUID)', 'name', 'FK:fitting_method',
             'FK:plate_layout', 'FK:target',
             'fitting_parameters{}'], pk_field='id (UUID)')

# FittingMethod
draw_entity(ax, 11.2, 15.3, 2.3, 1.5, colors['assays'], 'FittingMethod',
            ['id (UUID)', 'name, slug', 'version', 'script'], pk_field='id (UUID)')

# PlateLayout
draw_entity(ax, 13.7, 15.3, 2, 1.2, colors['assays'], 'PlateLayout',
            ['id (UUID)', 'name', 'config{}'], pk_field='id (UUID)')

# Assay
draw_entity(ax, 8.1, 12.8, 2.5, 1.8, colors['assays'], 'Assay',
            ['id (UUID)', 'FK:protocol', 'FK:target',
             'data_file', 'labbook_number'], pk_field='id (UUID)')

# DataSeries
draw_entity(ax, 11, 12.5, 2.8, 2.1, colors['assays'], 'DataSeries',
            ['id (UUID)', 'FK:assay', 'FK:compound',
             'FK:batch', 'row, columns',
             'extracted_data{}', 'skip_points[]'], pk_field='id (UUID)')

# AnalysisResult
draw_entity(ax, 14.1, 12.8, 2.2, 1.3, colors['assays'], 'AnalysisResult',
            ['id (UUID)', 'status', 'results{}'], pk_field='id (UUID)')

# Hypothesis
draw_entity(ax, 8.1, 10.7, 2.8, 1.7, colors['assays'], 'Hypothesis',
            ['id (UUID)', 'FK:target', 'FK:parent',
             'FK:product_compound', 'smiles', 'status'], pk_field='id (UUID)')

# DilutionSeries
draw_entity(ax, 11.2, 10.7, 2.2, 1.2, colors['assays'], 'DilutionSeries',
            ['id (UUID)', 'concentrations[]', 'unit'], pk_field='id (UUID)')

# Assay Relationships
draw_relationship(ax, (10.6, 13.5), (11, 13.5), 'contains')  # Assay -> DataSeries
draw_relationship(ax, (13.8, 13.2), (14.1, 13.2), 'analysis')  # DataSeries -> AnalysisResult
draw_relationship(ax, (9.4, 14.6), (9.4, 15), 'defines')  # Assay -> Protocol
draw_relationship(ax, (10.9, 15.8), (11.2, 15.8), 'uses')  # Protocol -> FittingMethod

# Cross-domain relationships (to Registry)
ax.annotate('', xy=(7.9, 13.2), xytext=(11.8, 13.2),
            arrowprops=dict(arrowstyle='->', color='#DC2626', lw=1.2,
                           connectionstyle='arc3,rad=-0.2', linestyle='--'))
ax.text(9.8, 12.7, 'compound', fontsize=6, ha='center', color='#DC2626', style='italic')

# ==================== CONSTRUCTS SECTION ====================

draw_section(ax, 15.9, 10.5, 7.9, 6.5, colors['section_constructs'], 'Constructs')

# ConstructProject
draw_entity(ax, 16.1, 15.2, 2.4, 1.4, colors['constructs'], 'ConstructProject',
            ['id (UUID)', 'name', 'FK:parent'], pk_field='id (UUID)')

# Protein
draw_entity(ax, 18.8, 15.2, 2.2, 1.4, colors['constructs'], 'Protein',
            ['id (UUID)', 'uniprot_id'], pk_field='id (UUID)')

# ProteinSynonym
draw_entity(ax, 21.3, 15.2, 2.2, 1.2, colors['constructs'], 'ProteinSynonym',
            ['id (UUID)', 'FK:protein', 'name'], pk_field='id (UUID)')

# Plasmid
draw_entity(ax, 16.1, 12.8, 2.5, 1.8, colors['constructs'], 'Plasmid',
            ['id (UUID)', 'ncn_id', 'name', 'FK:project',
             'FK:parent', 'genbank_file'], pk_field='id (UUID)')

# Cassette
draw_entity(ax, 19, 12.8, 2.2, 1.4, colors['constructs'], 'Cassette',
            ['id (UUID)', 'FK:protein', 'start', 'end'], pk_field='id (UUID)')

# CassetteUse
draw_entity(ax, 16.1, 10.7, 2.4, 1.6, colors['constructs'], 'CassetteUse',
            ['id (UUID)', 'FK:cassette', 'FK:plasmid',
             'alignment_file'], pk_field='id (UUID)')

# ExpressionTag
draw_entity(ax, 18.8, 10.7, 2.6, 1.6, colors['constructs'], 'ExpressionTag',
            ['id (UUID)', 'FK:cassette_use', 'FK:type',
             'FK:location', 'FK:protease'], pk_field='id (UUID)')

# SequencingResult
draw_entity(ax, 21.6, 10.7, 2.2, 1.4, colors['constructs'], 'SequencingResult',
            ['id (UUID)', 'FK:cassette_use', 'FK:plasmid', 'file'], pk_field='id (UUID)')

# Construct Relationships
draw_relationship(ax, (17.3, 14.6), (17.3, 15.2), 'in')  # Plasmid -> Project
draw_relationship(ax, (19.9, 14.2), (19.9, 15.2), 'of')  # Cassette -> Protein
draw_relationship(ax, (17.3, 12.3), (17.3, 12.8), 'uses')  # CassetteUse -> Plasmid
draw_relationship(ax, (18.5, 11.5), (19, 12.8), 'cassette')  # CassetteUse -> Cassette
draw_relationship(ax, (18.8, 11.3), (18.8, 11.5), 'tags')  # ExpressionTag -> CassetteUse
draw_relationship(ax, (21, 15.7), (21.3, 15.7), 'alias')  # Protein -> ProteinSynonym

# ==================== SERVER/PROJECT SECTION ====================

draw_section(ax, 0.2, 3.5, 11.5, 6.5, colors['section_server'], 'Projects & Jobs')

# Project
draw_entity(ax, 0.4, 7.8, 2.8, 1.8, colors['server'], 'Project',
            ['id', 'uuid', 'name', 'directory',
             'creation_time', 'last_job_number'], pk_field='id')

# ProjectGroup
draw_entity(ax, 3.5, 8, 2.2, 1.4, colors['server'], 'ProjectGroup',
            ['id', 'name', 'type'], pk_field='id')

# ProjectGroupMembership
draw_entity(ax, 3.5, 6.3, 2.5, 1.4, colors['server'], 'ProjectGroupMembership',
            ['id', 'FK:group', 'FK:project', 'type'], pk_field='id')

# ProjectTag
draw_entity(ax, 6, 8, 2.2, 1.4, colors['server'], 'ProjectTag',
            ['id', 'text', 'FK:parent', 'M2M:projects'], pk_field='id')

# Job
draw_entity(ax, 0.4, 4.2, 3.2, 2.4, colors['server'], 'Job',
            ['id', 'uuid', 'number', 'FK:project',
             'FK:parent', 'task_name', 'status',
             'evaluation', 'creation_time'], pk_field='id')

# ServerJob
draw_entity(ax, 3.9, 4.2, 2.4, 1.8, colors['server'], 'ServerJob',
            ['FK:job (PK)', 'server_process_id',
             'machine', 'mechanism'], pk_field='FK:job (PK)')

# File
draw_entity(ax, 6.6, 4.2, 2.6, 2, colors['server'], 'File',
            ['id', 'uuid', 'name', 'FK:type',
             'FK:job', 'directory'], pk_field='id')

# FileType
draw_entity(ax, 9.5, 4.5, 2, 1.2, colors['server'], 'FileType',
            ['name (PK)', 'description'], pk_field='name (PK)')

# FileUse
draw_entity(ax, 6.6, 6.5, 2.4, 1.4, colors['server'], 'FileUse',
            ['FK:file', 'FK:job', 'role', 'job_param_name'], pk_field='FK:file')

# XData
draw_entity(ax, 9.5, 6.5, 2, 1.2, colors['server'], 'XData',
            ['id', 'data_class', 'xml', 'FK:job'], pk_field='id')

# Server Relationships
draw_relationship(ax, (2.0, 7.8), (2.0, 6.6), 'has')  # Project -> Job
draw_relationship(ax, (3.2, 8.5), (3.5, 8.5), 'in')  # Project -> ProjectGroup
draw_relationship(ax, (3.6, 5.3), (3.9, 5.3), 'extends')  # Job -> ServerJob
draw_relationship(ax, (3.6, 4.8), (6.6, 4.8), 'outputs')  # Job -> File
draw_relationship(ax, (9.2, 5.2), (9.5, 5.2), 'typed')  # File -> FileType
draw_relationship(ax, (6.8, 6.2), (6.8, 6.5), 'usage')  # File -> FileUse

# ==================== USERS SECTION ====================

draw_section(ax, 12, 3.5, 5.5, 3, colors['section_registry'], 'Users')

# User (Django)
draw_entity(ax, 12.2, 4, 2.4, 2, colors['users'], 'User (Django)',
            ['id', 'username', 'email', 'first_name',
             'last_name', 'is_staff'], pk_field='id')

# UserProfile
draw_entity(ax, 14.9, 4, 2.4, 2, colors['users'], 'UserProfile',
            ['FK:user (PK)', 'role', 'operating_level',
             'legacy_username', 'last_seen_at'], pk_field='FK:user (PK)')

draw_relationship(ax, (14.6, 5), (14.9, 5), 'profile')

# ==================== REFERENCE DATA ====================

draw_section(ax, 17.8, 3.5, 6, 6.5, colors['section_constructs'], 'Reference Data')

# Assay Reference
ax.text(18, 9.6, 'Assay Config', fontsize=8, fontweight='bold', color=colors['assays'])

draw_entity(ax, 18, 8.2, 1.8, 1, colors['reference'], 'PlateLayout',
            ['id', 'name', 'config{}'], pk_field='id')

draw_entity(ax, 20, 8.2, 1.8, 1, colors['reference'], 'DilutionSeries',
            ['id', 'concentrations[]'], pk_field='id')

draw_entity(ax, 22, 8.2, 1.6, 1, colors['reference'], 'FittingMethod',
            ['id', 'slug'], pk_field='id')

# Construct Reference
ax.text(18, 7.6, 'Construct Config', fontsize=8, fontweight='bold', color=colors['constructs'])

draw_entity(ax, 18, 6.2, 2, 1, colors['reference'], 'ExpressionTagType',
            ['id', 'name'], pk_field='id')

draw_entity(ax, 20.2, 6.2, 1.8, 1, colors['reference'], 'Protease',
            ['id', 'name'], pk_field='id')

draw_entity(ax, 22.2, 6.2, 1.5, 1, colors['reference'], 'ExprTagLocation',
            ['id', 'name'], pk_field='id')

# Compound Reference
ax.text(18, 5.6, 'Compound Config', fontsize=8, fontweight='bold', color=colors['registry'])

draw_entity(ax, 18, 4, 2.6, 1.2, colors['reference'], 'MolecularPropertyThreshold',
            ['id', 'property_name', 'direction',
             'amber_threshold', 'red_threshold'], pk_field='id')

draw_entity(ax, 21, 4, 2.4, 1.2, colors['reference'], 'StereoComment',
            ['(choices on Compound)', 'achiral', 'racemate',
             'single_unknown...'], pk_field=None)

# ==================== LEGEND ====================

ax.text(0.5, 2.8, 'Legend', fontsize=11, fontweight='bold', color=colors['text_dark'])

# Domain colors
legend_items = [
    (colors['registry'], 'Compound Registry'),
    (colors['assays'], 'Assays'),
    (colors['constructs'], 'Constructs'),
    (colors['server'], 'Projects/Jobs'),
    (colors['users'], 'Users'),
    (colors['reference'], 'Reference Data'),
]

for i, (color, label) in enumerate(legend_items):
    x = 0.5 + (i % 3) * 3.5
    y = 2.3 - (i // 3) * 0.45
    rect = Rectangle((x, y - 0.15), 0.35, 0.3, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(x + 0.45, y, label, fontsize=8, va='center', color=colors['text_dark'])

# Relationship types
ax.text(11, 2.8, 'Relationships', fontsize=11, fontweight='bold', color=colors['text_dark'])
ax.text(11, 2.3, '>  Foreign Key (1:N)', fontsize=8, color=colors['text_dark'])
ax.text(11, 1.9, '-->  Cross-domain FK', fontsize=8, color='#DC2626')
ax.text(11, 1.5, 'Primary Key', fontsize=8, color=colors['text_dark'])

# Key counts
ax.text(16, 2.8, 'Model Counts', fontsize=11, fontweight='bold', color=colors['text_dark'])
counts = [
    'Registry: 7 models',
    'Assays: 9 models',
    'Constructs: 11 models',
    'Server: 14 models',
    'Users: 2 models',
]
for i, count in enumerate(counts):
    ax.text(16, 2.3 - i * 0.35, f'• {count}', fontsize=8, color=colors['text_dark'])

# Audit note
ax.text(21, 2.8, 'Common Audit Fields', fontsize=9, fontweight='bold', color=colors['text_dark'])
ax.text(21, 2.3, '• created_by (FK>User)', fontsize=7, color=colors['text_dark'])
ax.text(21, 2.0, '• modified_by (FK>User)', fontsize=7, color=colors['text_dark'])
ax.text(21, 1.7, '• created_at (auto)', fontsize=7, color=colors['text_dark'])
ax.text(21, 1.4, '• modified_at (auto)', fontsize=7, color=colors['text_dark'])

plt.tight_layout()
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/django_models_diagram.png',
            dpi=200, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/django_models_diagram.svg',
            format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: django_models_diagram.png and .svg")
