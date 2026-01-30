#!/usr/bin/env python3
"""
CCP4i2 Electron App Architecture Diagram
Shows the three-process model, IPC communication, and startup sequence.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(22, 16))
ax.set_xlim(0, 22)
ax.set_ylim(0, 16)
ax.set_aspect('equal')
ax.axis('off')

# Color palette
colors = {
    'electron_purple': '#8B5CF6',
    'nextjs_black': '#000000',
    'react_cyan': '#06B6D4',
    'django_green': '#059669',
    'ipc_orange': '#F59E0B',
    'storage_blue': '#3B82F6',
    'text_dark': '#1F2937',
    'white': '#FFFFFF',
    'section_purple': '#EDE9FE',
    'section_cyan': '#CFFAFE',
    'section_green': '#D1FAE5',
    'section_orange': '#FEF3C7',
    'section_blue': '#DBEAFE',
    'arrow': '#6B7280',
}

def draw_rounded_box(ax, x, y, width, height, color, label, sublabel=None,
                     text_color='white', fontsize=9, border_color=None):
    """Draw a rounded rectangle with label"""
    ec = border_color if border_color else 'none'
    lw = 2 if border_color else 0
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.1",
                         facecolor=color, edgecolor=ec, linewidth=lw, alpha=0.95)
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
                         boxstyle="round,pad=0.01,rounding_size=0.15",
                         facecolor=color, edgecolor='#9CA3AF',
                         linewidth=1.5, alpha=0.5)
    ax.add_patch(box)
    ax.text(x + 0.15, y + height - 0.25, title, fontsize=11,
            ha='left', va='top', color=colors['text_dark'], fontweight='bold')

def draw_arrow(ax, start, end, color='#6B7280', label=None, curved=False):
    """Draw an arrow with optional label"""
    style = 'arc3,rad=0.15' if curved else 'arc3,rad=0'
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='->', color=color, lw=2,
                               connectionstyle=style))
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(mid_x, mid_y + 0.2, label, fontsize=7, ha='center',
                color=color, fontweight='bold')

def draw_bidirectional_arrow(ax, start, end, color='#6B7280', label=None):
    """Draw a bidirectional arrow"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='<->', color=color, lw=2))
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(mid_x, mid_y + 0.25, label, fontsize=7, ha='center',
                color=color, fontweight='bold')

# Title
ax.text(11, 15.7, 'CCP4i2 Electron App Architecture', fontsize=20,
        ha='center', va='center', fontweight='bold', color=colors['text_dark'])
ax.text(11, 15.3, 'Three-Process Model with IPC Communication', fontsize=12,
        ha='center', va='center', color=colors['arrow'])

# ==================== MAIN PROCESS SECTION ====================

draw_section_box(ax, 0.3, 9.5, 7, 5.5, colors['section_purple'], 'MAIN PROCESS (Node.js)')

# Main process components
draw_rounded_box(ax, 0.6, 13.3, 3, 1.3, colors['electron_purple'],
                 'ccp4i2-master.ts', 'App Entry Point', fontsize=8)

draw_rounded_box(ax, 3.9, 13.3, 3.1, 1.3, colors['electron_purple'],
                 'electron-store', 'Persistent Config', fontsize=8)

draw_rounded_box(ax, 0.6, 11.7, 2.1, 1.3, colors['electron_purple'],
                 'IPC Handlers', '15+ channels', fontsize=8)

draw_rounded_box(ax, 2.9, 11.7, 2.1, 1.3, colors['electron_purple'],
                 'Window Mgr', 'BrowserWindow', fontsize=8)

draw_rounded_box(ax, 5.2, 11.7, 1.8, 1.3, colors['electron_purple'],
                 'Menus', 'App Menu', fontsize=8)

draw_rounded_box(ax, 0.6, 10, 3.1, 1.3, colors['electron_purple'],
                 'Django Server Mgr', 'Spawn/Kill Process', fontsize=8)

draw_rounded_box(ax, 4, 10, 3, 1.3, colors['electron_purple'],
                 'Next.js Server', 'Express + Next', fontsize=8)

# Config store contents
ax.text(5.45, 12.9, 'CCP4Dir\nprojectRoot\nzoomLevel\ndevMode\ntheme',
        fontsize=6, ha='center', va='top', color=colors['text_dark'])

# ==================== RENDERER PROCESS SECTION ====================

draw_section_box(ax, 7.6, 9.5, 6.8, 5.5, colors['section_cyan'], 'RENDERER PROCESS (Chromium)')

# Renderer components
draw_rounded_box(ax, 7.9, 13.3, 3.1, 1.3, colors['react_cyan'],
                 'Next.js 15', 'App Router', fontsize=8, text_color='#0E7490')

draw_rounded_box(ax, 11.3, 13.3, 2.8, 1.3, colors['react_cyan'],
                 'React 19', 'UI Components', fontsize=8, text_color='#0E7490')

draw_rounded_box(ax, 7.9, 11.7, 2.5, 1.3, '#000000',
                 'Preload Bridge', 'electronAPI', fontsize=8)

draw_rounded_box(ax, 10.7, 11.7, 3.4, 1.3, colors['react_cyan'],
                 'Moorhen Viewer', '3D Structure Viewer', fontsize=8, text_color='#0E7490')

draw_rounded_box(ax, 7.9, 10, 3, 1.3, colors['storage_blue'],
                 'api-fetch.ts', 'HTTP Client', fontsize=8)

draw_rounded_box(ax, 11.2, 10, 2.9, 1.3, colors['storage_blue'],
                 'SWR Hooks', 'Data Caching', fontsize=8)

# Routes
ax.text(9.45, 12.9, '/ccp4i2/*\n/registry/*\n/assays/*',
        fontsize=6, ha='center', va='top', color=colors['text_dark'])

# ==================== DJANGO SERVER SECTION ====================

draw_section_box(ax, 14.7, 9.5, 7, 5.5, colors['section_green'], 'DJANGO SERVER (Python)')

# Django components
draw_rounded_box(ax, 15, 13.3, 3.2, 1.3, colors['django_green'],
                 'Uvicorn ASGI', '2 Workers', fontsize=8)

draw_rounded_box(ax, 18.5, 13.3, 2.9, 1.3, colors['django_green'],
                 'Django 4.2', 'REST Framework', fontsize=8)

draw_rounded_box(ax, 15, 11.7, 2.5, 1.3, colors['django_green'],
                 'ViewSets', 'API Endpoints', fontsize=8)

draw_rounded_box(ax, 17.8, 11.7, 2, 1.3, colors['django_green'],
                 'ORM', 'Models', fontsize=8)

draw_rounded_box(ax, 20.1, 11.7, 1.4, 1.3, colors['django_green'],
                 'Auth', 'JWT/Dev', fontsize=8)

draw_rounded_box(ax, 15, 10, 2.8, 1.3, colors['django_green'],
                 'Task Manager', '176 Plugins', fontsize=8)

draw_rounded_box(ax, 18.1, 10, 2.4, 1.3, colors['django_green'],
                 'Process Mgr', 'Job Execution', fontsize=8)

draw_rounded_box(ax, 20.8, 10, 0.9, 1.3, colors['django_green'],
                 'DB', 'SQLite', fontsize=7)

# ==================== IPC COMMUNICATION ====================

draw_section_box(ax, 0.3, 6, 14.1, 3, colors['section_orange'], 'IPC COMMUNICATION CHANNELS')

# IPC channels listed
ipc_channels = [
    ('get-config', 'Request configuration'),
    ('start-uvicorn', 'Start Django server'),
    ('locate-ccp4', 'File dialog for CCP4'),
    ('check-file-exists', 'File system check'),
    ('toggle-dev-mode', 'Toggle development'),
]

for i, (channel, desc) in enumerate(ipc_channels):
    x = 0.6 + (i % 3) * 4.7
    y = 7.8 - (i // 3) * 0.9
    draw_rounded_box(ax, x, y, 2.2, 0.7, colors['ipc_orange'],
                     channel, fontsize=7, text_color='#1F2937')
    ax.text(x + 2.35, y + 0.35, desc, fontsize=6, va='center', color=colors['text_dark'])

# Additional channels
more_channels = ['zoom-in/out', 'set-theme-mode', 'install-requirements', 'message-from-main']
ax.text(0.6, 6.35, 'More: ' + ', '.join(more_channels), fontsize=7, color=colors['text_dark'])

# ==================== HTTP COMMUNICATION ====================

draw_section_box(ax, 14.7, 6, 7, 3, colors['section_blue'], 'HTTP COMMUNICATION')

# HTTP flow
ax.text(15, 8.2, 'API Routes:', fontsize=9, fontweight='bold', color=colors['text_dark'])
routes = [
    '/api/proxy/ccp4i2/* -> localhost:8001',
    'POST /jobs/{id}/run/',
    'GET /projects/{id}/job_tree/',
    'GET /files/{id}/download/',
]
for i, route in enumerate(routes):
    ax.text(15, 7.8 - i * 0.4, route, fontsize=7, color=colors['text_dark'],
            family='monospace')

# ==================== STARTUP SEQUENCE ====================

draw_section_box(ax, 0.3, 0.5, 21.4, 5, colors['section_purple'], 'STARTUP SEQUENCE')

# Sequence steps
steps = [
    ('1', 'npm start:electron', 'Build main + preload via Vite'),
    ('2', 'electron dist/main.js', 'Main process starts'),
    ('3', 'detectPort(3000, 8001)', 'Find available ports'),
    ('4', 'startNextServer()', 'Launch Express + Next.js'),
    ('5', 'createWindow()', 'Open BrowserWindow'),
    ('6', 'load /ccp4i2/config', 'Show configuration page'),
    ('7', 'IPC: start-uvicorn', 'User clicks Start Server'),
    ('8', 'spawn(ccp4-python)', 'Start Django/Uvicorn'),
    ('9', 'Run migrations', 'Initialize database'),
    ('10', 'navigate /ccp4i2', 'User starts working'),
]

for i, (num, action, desc) in enumerate(steps):
    x = 0.6 + (i % 5) * 4.3
    y = 4.1 - (i // 5) * 1.7

    # Number circle
    circle = plt.Circle((x + 0.2, y + 0.5), 0.25, color=colors['electron_purple'])
    ax.add_patch(circle)
    ax.text(x + 0.2, y + 0.5, num, fontsize=8, ha='center', va='center',
            color='white', fontweight='bold')

    # Action and description
    ax.text(x + 0.55, y + 0.7, action, fontsize=7, fontweight='bold',
            color=colors['text_dark'])
    ax.text(x + 0.55, y + 0.3, desc, fontsize=6, color=colors['arrow'])

# ==================== ARROWS BETWEEN SECTIONS ====================

# Main -> Renderer (IPC)
draw_bidirectional_arrow(ax, (7.0, 12.35), (7.9, 12.35), colors['ipc_orange'], 'IPC')

# Renderer -> Django (HTTP)
draw_arrow(ax, (11.2, 10.65), (15, 10.65), colors['storage_blue'], 'HTTP :8001')

# Main -> Django (spawn/kill)
draw_arrow(ax, (3.65, 10), (3.65, 9.5), colors['electron_purple'])
ax.annotate('', xy=(15.6, 9.5), xytext=(3.65, 9.5),
            arrowprops=dict(arrowstyle='->', color=colors['electron_purple'], lw=2,
                           connectionstyle='arc3,rad=-0.1'))
ax.text(9.6, 9.2, 'spawn(ccp4-python, uvicorn)', fontsize=7, ha='center',
        color=colors['electron_purple'])

# ==================== LEGEND ====================

ax.text(15, 5.3, 'Port Allocation', fontsize=10, fontweight='bold', color=colors['text_dark'])
ax.text(15, 4.9, 'Next.js: auto-detect from 3000', fontsize=8, color=colors['text_dark'])
ax.text(15, 4.5, 'Django: auto-detect from 8001', fontsize=8, color=colors['text_dark'])

ax.text(18.5, 5.3, 'Build Tools', fontsize=10, fontweight='bold', color=colors['text_dark'])
ax.text(18.5, 4.9, 'Vite: main + preload', fontsize=8, color=colors['text_dark'])
ax.text(18.5, 4.5, 'Next.js: renderer', fontsize=8, color=colors['text_dark'])

# Color legend
legend_items = [
    (colors['electron_purple'], 'Electron/Main'),
    (colors['react_cyan'], 'React/Renderer'),
    (colors['django_green'], 'Django/Python'),
    (colors['ipc_orange'], 'IPC Channel'),
    (colors['storage_blue'], 'HTTP/API'),
]

for i, (color, label) in enumerate(legend_items):
    x = 15 + (i % 3) * 2.3
    y = 4.0 - (i // 3) * 0.5
    rect = Rectangle((x, y - 0.1), 0.25, 0.25, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(x + 0.35, y + 0.02, label, fontsize=7, va='center', color=colors['text_dark'])

plt.tight_layout()
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/electron_architecture.png',
            dpi=200, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/electron_architecture.svg',
            format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: electron_architecture.png and .svg")
