#!/usr/bin/env python3
"""
CCP4i2 REST API Architecture Diagram
Shows all API endpoints, ViewSets, and request/response patterns.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(24, 18))
ax.set_xlim(0, 24)
ax.set_ylim(0, 18)
ax.set_aspect('equal')
ax.axis('off')

# Color palette
colors = {
    'get_green': '#10B981',
    'post_blue': '#3B82F6',
    'put_orange': '#F59E0B',
    'delete_red': '#EF4444',
    'viewset_purple': '#8B5CF6',
    'text_dark': '#1F2937',
    'white': '#FFFFFF',
    'section_bg': '#F3F4F6',
    'arrow': '#6B7280',
}

def draw_rounded_box(ax, x, y, width, height, color, label, sublabel=None,
                     text_color='white', fontsize=9):
    """Draw a rounded rectangle with label"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.08",
                         facecolor=color, edgecolor='none', alpha=0.95)
    ax.add_patch(box)

    label_y = y + height*0.65 if sublabel else y + height/2
    ax.text(x + width/2, label_y, label, fontsize=fontsize,
            ha='center', va='center', color=text_color, fontweight='bold')

    if sublabel:
        ax.text(x + width/2, y + height*0.3, sublabel, fontsize=fontsize-2,
                ha='center', va='center', color=text_color, alpha=0.9)

def draw_section_box(ax, x, y, width, height, color, title):
    """Draw a section box with title"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.1",
                         facecolor=color, edgecolor='#9CA3AF',
                         linewidth=1, alpha=0.4)
    ax.add_patch(box)
    ax.text(x + 0.15, y + height - 0.2, title, fontsize=11,
            ha='left', va='top', color=colors['text_dark'], fontweight='bold')

def draw_method_badge(ax, x, y, method):
    """Draw a small HTTP method badge"""
    method_colors = {
        'GET': colors['get_green'],
        'POST': colors['post_blue'],
        'PUT': colors['put_orange'],
        'DELETE': colors['delete_red'],
        'PATCH': '#EC4899',
    }
    color = method_colors.get(method, colors['arrow'])
    rect = Rectangle((x, y), 0.55, 0.25, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(x + 0.275, y + 0.125, method, fontsize=5.5, ha='center', va='center',
            color='white', fontweight='bold')

def draw_endpoint(ax, x, y, method, path, description):
    """Draw an endpoint with method badge"""
    draw_method_badge(ax, x, y, method)
    ax.text(x + 0.65, y + 0.125, path, fontsize=6.5, va='center',
            color=colors['text_dark'], family='monospace')
    ax.text(x + 0.65, y - 0.15, description, fontsize=5.5, va='center',
            color=colors['arrow'])
    return y - 0.45

# Title
ax.text(12, 17.6, 'CCP4i2 REST API Architecture', fontsize=20,
        ha='center', va='center', fontweight='bold', color=colors['text_dark'])
ax.text(12, 17.2, '/api/ccp4i2/* - ~135 Endpoints across 10 ViewSets', fontsize=12,
        ha='center', va='center', color=colors['arrow'])

# ==================== PROJECTS VIEWSET ====================

draw_section_box(ax, 0.2, 10.5, 7.8, 6.5, '#DBEAFE', 'ProjectViewSet - /projects/')

# Standard CRUD
ax.text(0.4, 16.5, 'Standard CRUD:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 16.1
y = draw_endpoint(ax, 0.4, y, 'GET', '/', 'List all projects')
y = draw_endpoint(ax, 0.4, y, 'POST', '/', 'Create new project')
y = draw_endpoint(ax, 0.4, y, 'GET', '/{id}/', 'Get project details')
y = draw_endpoint(ax, 0.4, y, 'PUT', '/{id}/', 'Update project')
y = draw_endpoint(ax, 0.4, y, 'DELETE', '/{id}/', 'Delete project')

# Custom actions
ax.text(0.4, y - 0.1, 'Custom Actions:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = y - 0.45
y = draw_endpoint(ax, 0.4, y, 'GET', '/{id}/job_tree/', 'Consolidated job hierarchy')
y = draw_endpoint(ax, 0.4, y, 'POST', '/{id}/create_task/', 'Create new job/task')
y = draw_endpoint(ax, 0.4, y, 'GET', '/{id}/directory/', 'List project directory')
y = draw_endpoint(ax, 0.4, y, 'GET', '/{id}/project_file/', 'Download project file')
y = draw_endpoint(ax, 0.4, y, 'POST', '/{id}/export/', 'Export as ZIP')
y = draw_endpoint(ax, 0.4, y, 'POST', '/import_project/', 'Import from ZIP')

# More endpoints in second column
ax.text(4.2, 16.5, 'More Actions:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 16.1
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/jobs/', 'List project jobs')
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/files/', 'List project files')
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/tags/', 'Get project tags')
y = draw_endpoint(ax, 4.2, y, 'POST', '/{id}/tags/', 'Add tag to project')
y = draw_endpoint(ax, 4.2, y, 'DELETE', '/{id}/tags/{tag_id}/', 'Remove tag')
y = draw_endpoint(ax, 4.2, y, 'POST', '/{id}/preview_file/', 'Open in viewer')
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/job_float_values/', 'Get job KPIs')
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/exports/', 'List exports')

# ==================== JOBS VIEWSET ====================

draw_section_box(ax, 8.2, 10.5, 7.8, 6.5, '#D1FAE5', 'JobViewSet - /jobs/')

ax.text(8.4, 16.5, 'Standard CRUD:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 16.1
y = draw_endpoint(ax, 8.4, y, 'GET', '/', 'List jobs (filter by project)')
y = draw_endpoint(ax, 8.4, y, 'POST', '/', 'Create new job')
y = draw_endpoint(ax, 8.4, y, 'GET', '/{id}/', 'Get job details')
y = draw_endpoint(ax, 8.4, y, 'PUT', '/{id}/', 'Update job')
y = draw_endpoint(ax, 8.4, y, 'DELETE', '/{id}/', 'Delete job + dependents')

ax.text(8.4, y - 0.1, 'Execution:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = y - 0.45
y = draw_endpoint(ax, 8.4, y, 'POST', '/{id}/run/', 'Execute job')
y = draw_endpoint(ax, 8.4, y, 'POST', '/{id}/run_local/', 'Force local execution')
y = draw_endpoint(ax, 8.4, y, 'POST', '/{id}/clone/', 'Clone job')
y = draw_endpoint(ax, 8.4, y, 'GET', '/{id}/validation/', 'Validate parameters')
y = draw_endpoint(ax, 8.4, y, 'GET', '/{id}/i2run_command/', 'Get CLI command')

ax.text(12, 16.5, 'Parameters & Data:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 16.1
y = draw_endpoint(ax, 12, y, 'GET', '/{id}/params_xml/', 'Get params as XML')
y = draw_endpoint(ax, 12, y, 'PUT', '/{id}/params_xml/', 'Update params XML')
y = draw_endpoint(ax, 12, y, 'POST', '/{id}/set_parameter/', 'Set single param')
y = draw_endpoint(ax, 12, y, 'GET', '/{id}/get_parameter/', 'Get single param')
y = draw_endpoint(ax, 12, y, 'POST', '/{id}/upload_file_param/', 'Upload file param')
y = draw_endpoint(ax, 12, y, 'GET', '/{id}/container/', 'Get container JSON')

ax.text(12, y - 0.1, 'Reports & Output:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = y - 0.45
y = draw_endpoint(ax, 12, y, 'GET', '/{id}/digest/', 'Get job digest')
y = draw_endpoint(ax, 12, y, 'GET', '/{id}/report_xml/', 'Get XML report')
y = draw_endpoint(ax, 12, y, 'POST', '/{id}/regenerate_report/', 'Regenerate report')
y = draw_endpoint(ax, 12, y, 'GET', '/{id}/export_job/', 'Export as ZIP')

# ==================== FILES VIEWSET ====================

draw_section_box(ax, 16.2, 10.5, 7.6, 6.5, '#FEF3C7', 'FileViewSet - /files/')

ax.text(16.4, 16.5, 'Standard CRUD:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 16.1
y = draw_endpoint(ax, 16.4, y, 'GET', '/', 'List all files')
y = draw_endpoint(ax, 16.4, y, 'POST', '/', 'Create file record')
y = draw_endpoint(ax, 16.4, y, 'GET', '/{id}/', 'Get file metadata')
y = draw_endpoint(ax, 16.4, y, 'PUT', '/{id}/', 'Update file metadata')
y = draw_endpoint(ax, 16.4, y, 'DELETE', '/{id}/', 'Delete file record')

ax.text(16.4, y - 0.1, 'Download & Preview:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = y - 0.45
y = draw_endpoint(ax, 16.4, y, 'GET', '/{id}/download/', 'Download file by ID')
y = draw_endpoint(ax, 16.4, y, 'GET', '/{id}/download_by_uuid/', 'Download by UUID')
y = draw_endpoint(ax, 16.4, y, 'GET', '/{id}/digest/', 'Get file digest')
y = draw_endpoint(ax, 16.4, y, 'GET', '/{id}/digest_by_uuid/', 'Digest by UUID')
y = draw_endpoint(ax, 16.4, y, 'POST', '/{id}/preview/', 'Open in viewer')

ax.text(20, 16.5, 'Related ViewSets:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 16.1
ax.text(20, y, '/filetypes/ - File type definitions', fontsize=7, color=colors['text_dark'])
y -= 0.35
ax.text(20, y, '/fileimports/ - Import records', fontsize=7, color=colors['text_dark'])
y -= 0.35
ax.text(20, y, '/fileuses/ - File usage tracking', fontsize=7, color=colors['text_dark'])
y -= 0.35
ax.text(20, y, '/projectexports/ - Export records', fontsize=7, color=colors['text_dark'])

# ==================== PROJECT GROUPS ====================

draw_section_box(ax, 0.2, 5.5, 7.8, 4.5, '#EDE9FE', 'ProjectGroupViewSet - /projectgroups/')

ax.text(0.4, 9.5, 'CRUD + Custom:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 9.1
y = draw_endpoint(ax, 0.4, y, 'GET', '/', 'List groups (filter by type)')
y = draw_endpoint(ax, 0.4, y, 'POST', '/', 'Create group')
y = draw_endpoint(ax, 0.4, y, 'GET', '/{id}/', 'Get group with memberships')
y = draw_endpoint(ax, 0.4, y, 'POST', '/create_with_parent/', 'Create campaign + parent')

ax.text(4.2, 9.5, 'Member Management:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 9.1
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/member_projects/', 'Get members + KPIs')
y = draw_endpoint(ax, 4.2, y, 'GET', '/{id}/parent_project/', 'Get parent project')
y = draw_endpoint(ax, 4.2, y, 'POST', '/{id}/add_member/', 'Add project to group')
y = draw_endpoint(ax, 4.2, y, 'DELETE', '/{id}/members/{pid}/', 'Remove member')

# ==================== PROJECT TAGS ====================

draw_section_box(ax, 8.2, 5.5, 3.8, 4.5, '#FCE7F3', 'ProjectTagViewSet')

ax.text(8.4, 9.5, '/projecttags/', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 9.0
y = draw_endpoint(ax, 8.4, y, 'GET', '/', 'List all tags')
y = draw_endpoint(ax, 8.4, y, 'POST', '/', 'Create tag')
y = draw_endpoint(ax, 8.4, y, 'GET', '/{id}/', 'Get tag')
y = draw_endpoint(ax, 8.4, y, 'PUT', '/{id}/', 'Update tag')
y = draw_endpoint(ax, 8.4, y, 'DELETE', '/{id}/', 'Delete tag')
ax.text(8.4, y - 0.2, 'Hierarchical (parent FK)', fontsize=6, color=colors['arrow'])

# ==================== UTILITY ENDPOINTS ====================

draw_section_box(ax, 12.2, 5.5, 5.4, 4.5, '#E0E7FF', 'Utility Endpoints')

ax.text(12.4, 9.5, 'System:', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = 9.0
y = draw_endpoint(ax, 12.4, y, 'GET', '/health/', 'Health check (DB status)')
y = draw_endpoint(ax, 12.4, y, 'GET', '/task_tree/', 'Available plugins (176)')
y = draw_endpoint(ax, 12.4, y, 'GET', '/active_jobs/', 'Running jobs + resources')

ax.text(12.4, y - 0.2, 'Admin (Platform Admin only):', fontsize=8, fontweight='bold', color=colors['text_dark'])
y = y - 0.6
y = draw_endpoint(ax, 12.4, y, 'POST', '/admin/import-legacy/', 'Import dumpdata')
y = draw_endpoint(ax, 12.4, y, 'GET', '/admin/import-status/', 'Database counts')

# ==================== REQUEST/RESPONSE ====================

draw_section_box(ax, 17.8, 5.5, 6, 4.5, '#F3F4F6', 'Request/Response Pattern')

ax.text(18, 9.5, 'Standard Response:', fontsize=8, fontweight='bold', color=colors['text_dark'])
response_json = '''{"success": true,
 "data": {
   "id": 123,
   "name": "...",
   ...
 }}'''
ax.text(18, 9.1, response_json, fontsize=6, color=colors['text_dark'],
        family='monospace', va='top')

ax.text(18, 7.3, 'Error Response:', fontsize=8, fontweight='bold', color=colors['text_dark'])
error_json = '''{"success": false,
 "error": "Message",
 "details": {...}}'''
ax.text(18, 6.9, error_json, fontsize=6, color=colors['text_dark'],
        family='monospace', va='top')

# ==================== AUTHENTICATION ====================

draw_section_box(ax, 0.2, 0.5, 11.6, 4.5, '#FEE2E2', 'Authentication & Middleware')

ax.text(0.4, 4.5, 'Authentication Flow:', fontsize=9, fontweight='bold', color=colors['text_dark'])

# Auth flow boxes
draw_rounded_box(ax, 0.5, 2.8, 2.2, 1.2, colors['post_blue'],
                 'Request', 'Authorization: Bearer', fontsize=7)

draw_rounded_box(ax, 3.2, 2.8, 2.4, 1.2, colors['viewset_purple'],
                 'Azure AD\nMiddleware', 'Validate JWT', fontsize=7)

draw_rounded_box(ax, 6.1, 2.8, 2, 1.2, colors['get_green'],
                 'ViewSet', 'IsAuthenticated', fontsize=7)

draw_rounded_box(ax, 8.6, 2.8, 2.4, 1.2, colors['get_green'],
                 'Response', 'JSON + Status', fontsize=7)

# Arrows
ax.annotate('', xy=(3.2, 3.4), xytext=(2.7, 3.4),
            arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=1.5))
ax.annotate('', xy=(6.1, 3.4), xytext=(5.6, 3.4),
            arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=1.5))
ax.annotate('', xy=(8.6, 3.4), xytext=(8.1, 3.4),
            arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=1.5))

# Public paths
ax.text(0.4, 2.3, 'Public Paths (no auth):', fontsize=8, fontweight='bold', color=colors['text_dark'])
ax.text(0.4, 1.9, '/health/, /.auth/*, /static/*, /media/*', fontsize=7,
        color=colors['text_dark'], family='monospace')

# Modes
ax.text(0.4, 1.4, 'Auth Modes:', fontsize=8, fontweight='bold', color=colors['text_dark'])
ax.text(0.4, 1.0, 'Production: Azure AD OAuth2 + JWT validation', fontsize=7, color=colors['text_dark'])
ax.text(0.4, 0.7, 'Development: Auto-dev-user when DEBUG=True', fontsize=7, color=colors['text_dark'])

# ==================== ENDPOINT STATISTICS ====================

draw_section_box(ax, 12, 0.5, 11.8, 4.5, '#DBEAFE', 'API Statistics & Features')

# Stats
stats = [
    ('ViewSets:', '10 (Project, Job, File, etc.)'),
    ('Total Endpoints:', '~135'),
    ('Standard CRUD:', '50 (5 per ViewSet)'),
    ('Custom Actions:', '~80'),
    ('Utility:', '5'),
]

ax.text(12.2, 4.5, 'Endpoint Count:', fontsize=9, fontweight='bold', color=colors['text_dark'])
for i, (label, value) in enumerate(stats):
    ax.text(12.2, 4.0 - i * 0.4, f'{label} {value}', fontsize=7, color=colors['text_dark'])

# Features
ax.text(16.5, 4.5, 'Features:', fontsize=9, fontweight='bold', color=colors['text_dark'])
features = [
    'Pagination (PageNumberPagination)',
    'Filtering (?project=, ?type=)',
    'Ordering (?ordering=-last_access)',
    'Multipart uploads supported',
    'XML parameter round-trips',
]
for i, feat in enumerate(features):
    ax.text(16.5, 4.0 - i * 0.4, f'- {feat}', fontsize=7, color=colors['text_dark'])

ax.text(20.5, 4.5, 'Parsers:', fontsize=9, fontweight='bold', color=colors['text_dark'])
parsers = ['MultiPartParser', 'JSONParser', 'FormParser']
for i, parser in enumerate(parsers):
    ax.text(20.5, 4.0 - i * 0.4, f'- {parser}', fontsize=7, color=colors['text_dark'])

# ==================== METHOD LEGEND ====================

ax.text(12.2, 2.0, 'HTTP Methods:', fontsize=9, fontweight='bold', color=colors['text_dark'])
methods = [
    ('GET', colors['get_green'], 'Read/List'),
    ('POST', colors['post_blue'], 'Create/Action'),
    ('PUT', colors['put_orange'], 'Update'),
    ('DELETE', colors['delete_red'], 'Delete'),
]
for i, (method, color, desc) in enumerate(methods):
    x = 12.2 + (i % 2) * 4
    y = 1.5 - (i // 2) * 0.5
    rect = Rectangle((x, y), 0.6, 0.3, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(x + 0.3, y + 0.15, method, fontsize=6, ha='center', va='center',
            color='white', fontweight='bold')
    ax.text(x + 0.75, y + 0.15, desc, fontsize=7, va='center', color=colors['text_dark'])

plt.tight_layout()
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/api_architecture.png',
            dpi=200, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/api_architecture.svg',
            format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: api_architecture.png and .svg")
