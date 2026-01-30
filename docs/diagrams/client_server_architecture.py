#!/usr/bin/env python3
"""
CCP4i2 Client/Server Architecture Diagram
Generates a professional diagram showing the frontend-backend communication,
deployment modes, and key components.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(20, 14))
ax.set_xlim(0, 20)
ax.set_ylim(0, 14)
ax.set_aspect('equal')
ax.axis('off')

# Color palette
colors = {
    'frontend_blue': '#3B82F6',
    'backend_green': '#10B981',
    'electron_purple': '#8B5CF6',
    'react_cyan': '#06B6D4',
    'django_green': '#059669',
    'database_amber': '#F59E0B',
    'arrow_gray': '#6B7280',
    'text_dark': '#1F2937',
    'white': '#FFFFFF',
    'section_blue': '#DBEAFE',
    'section_green': '#D1FAE5',
    'section_purple': '#EDE9FE',
    'section_orange': '#FEF3C7',
}

def draw_rounded_box(ax, x, y, width, height, color, label, sublabel=None,
                     text_color='white', fontsize=9, border_color=None):
    """Draw a rounded rectangle with label"""
    ec = border_color if border_color else 'none'
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.1",
                         facecolor=color, edgecolor=ec, linewidth=2, alpha=0.95)
    ax.add_patch(box)

    label_y = y + height*0.6 if sublabel else y + height/2
    ax.text(x + width/2, label_y, label, fontsize=fontsize,
            ha='center', va='center', color=text_color, fontweight='bold')

    if sublabel:
        ax.text(x + width/2, y + height*0.25, sublabel, fontsize=fontsize-2,
                ha='center', va='center', color=text_color, alpha=0.9)

def draw_section_box(ax, x, y, width, height, color, title, title_pos='top'):
    """Draw a section box with title"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.15",
                         facecolor=color, edgecolor='#9CA3AF',
                         linewidth=1.5, alpha=0.5)
    ax.add_patch(box)
    if title_pos == 'top':
        ax.text(x + 0.15, y + height - 0.2, title, fontsize=10,
                ha='left', va='top', color=colors['text_dark'], fontweight='bold')
    else:
        ax.text(x + width/2, y + height + 0.15, title, fontsize=10,
                ha='center', va='bottom', color=colors['text_dark'], fontweight='bold')

def draw_arrow(ax, start, end, color='#6B7280', label=None, label_pos='above'):
    """Draw an arrow with optional label"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='->', color=color, lw=2))
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        offset = 0.15 if label_pos == 'above' else -0.15
        ax.text(mid_x, mid_y + offset, label, fontsize=7, ha='center',
                va='center', color=color, fontweight='bold')

def draw_bidirectional_arrow(ax, start, end, color='#6B7280', label=None):
    """Draw a bidirectional arrow"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='<->', color=color, lw=2))
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(mid_x, mid_y + 0.2, label, fontsize=7, ha='center',
                va='center', color=color, fontweight='bold')

# Title
ax.text(10, 13.7, 'CCP4i2 Client/Server Architecture', fontsize=18,
        ha='center', va='center', fontweight='bold', color=colors['text_dark'])
ax.text(10, 13.3, 'Next.js 15 + React 19 Frontend • Django 4.2 + DRF Backend', fontsize=11,
        ha='center', va='center', color=colors['arrow_gray'])

# ==================== DEPLOYMENT MODES (Top) ====================

# Desktop Mode
draw_section_box(ax, 0.3, 11, 6, 2, colors['section_purple'], 'Desktop Mode (Electron)')
draw_rounded_box(ax, 0.6, 11.2, 2.5, 1.4, colors['electron_purple'],
                 'Electron Main', 'Process Manager', fontsize=8)
draw_rounded_box(ax, 3.4, 11.2, 2.5, 1.4, colors['electron_purple'],
                 'IPC Bridge', 'Main ↔ Renderer', fontsize=8)

# Web/Docker Mode
draw_section_box(ax, 7, 11, 6, 2, colors['section_orange'], 'Web/Docker Mode')
draw_rounded_box(ax, 7.3, 11.2, 2.5, 1.4, colors['database_amber'],
                 'Docker Web', 'Next.js Container', fontsize=8, text_color='#1F2937')
draw_rounded_box(ax, 10.1, 11.2, 2.5, 1.4, colors['database_amber'],
                 'Docker Server', 'Django Container', fontsize=8, text_color='#1F2937')

# Azure Mode
draw_section_box(ax, 13.7, 11, 6, 2, colors['section_blue'], 'Azure Mode')
draw_rounded_box(ax, 14, 11.2, 2.5, 1.4, '#0078D4',
                 'Azure AD', 'OAuth 2.0 + JWT', fontsize=8)
draw_rounded_box(ax, 16.8, 11.2, 2.5, 1.4, '#0078D4',
                 'Container Apps', 'Auto-scaling', fontsize=8)

# ==================== FRONTEND SECTION ====================

draw_section_box(ax, 0.3, 5.5, 9.2, 5, colors['section_blue'], 'Frontend (client/renderer/)')

# React/Next.js Core
draw_rounded_box(ax, 0.6, 9, 4, 1.2, colors['react_cyan'],
                 'Next.js 15 App Router', 'React 19 • TypeScript', fontsize=9, text_color='#0E7490')

# Pages/Routes
draw_rounded_box(ax, 5, 9, 4.2, 1.2, colors['frontend_blue'],
                 'App Routes', '/ccp4i2 • /registry • /assays', fontsize=9)

# Components
draw_rounded_box(ax, 0.6, 7.3, 2.8, 1.3, colors['frontend_blue'],
                 'Components', 'Task • Job • Project', fontsize=8)

# Providers
draw_rounded_box(ax, 3.6, 7.3, 2.6, 1.3, colors['frontend_blue'],
                 'Providers', 'Auth • Toast • Dialog', fontsize=8)

# API Layer
draw_rounded_box(ax, 6.4, 7.3, 2.8, 1.3, '#1D4ED8',
                 'API Layer', 'api.ts • SWR Hooks', fontsize=8)

# API Fetch
draw_rounded_box(ax, 0.6, 5.8, 4, 1.1, '#1E40AF',
                 'api-fetch.ts', 'Centralized fetch • Token injection', fontsize=8)

# Proxy Route
draw_rounded_box(ax, 5, 5.8, 4.2, 1.1, '#1E40AF',
                 'API Proxy Routes', '/api/proxy/ccp4i2/[...path]', fontsize=8)

# Arrows within frontend
draw_arrow(ax, (2.6, 9), (2.6, 8.6), colors['arrow_gray'])
draw_arrow(ax, (4.9, 7.9), (6.4, 7.9), colors['arrow_gray'])
draw_arrow(ax, (7.8, 7.3), (7.8, 6.9), colors['arrow_gray'])

# ==================== BACKEND SECTION ====================

draw_section_box(ax, 10.5, 5.5, 9.2, 5, colors['section_green'], 'Backend (server/ccp4i2/)')

# Django Core
draw_rounded_box(ax, 10.8, 9, 4, 1.2, colors['django_green'],
                 'Django 4.2 + DRF', 'Uvicorn ASGI Server', fontsize=9)

# URL Routing
draw_rounded_box(ax, 15.2, 9, 4.2, 1.2, colors['backend_green'],
                 'URL Router', '/api/ccp4i2/* endpoints', fontsize=9)

# ViewSets
draw_rounded_box(ax, 10.8, 7.3, 2.8, 1.3, colors['backend_green'],
                 'ViewSets', 'Project • Job • File', fontsize=8)

# Middleware
draw_rounded_box(ax, 13.8, 7.3, 2.6, 1.3, colors['backend_green'],
                 'Middleware', 'Auth • CORS • CORP', fontsize=8)

# Core Modules
draw_rounded_box(ax, 16.6, 7.3, 2.8, 1.3, colors['django_green'],
                 'Core Modules', '176 Task Plugins', fontsize=8)

# Models
draw_rounded_box(ax, 10.8, 5.8, 4, 1.1, '#047857',
                 'db/models.py', 'Project • Job • File • Tags', fontsize=8)

# Task Manager
draw_rounded_box(ax, 15.2, 5.8, 4.2, 1.1, '#047857',
                 'CCP4TaskManager', 'Plugin Registry • Containers', fontsize=8)

# Arrows within backend
draw_arrow(ax, (12.8, 9), (12.8, 8.6), colors['arrow_gray'])
draw_arrow(ax, (16.2, 7.9), (16.6, 7.9), colors['arrow_gray'])
draw_arrow(ax, (15.2, 7.9), (13.8, 7.3), colors['arrow_gray'])

# ==================== COMMUNICATION FLOW ====================

# Main communication arrow
draw_bidirectional_arrow(ax, (9.5, 6.35), (10.5, 6.35), '#374151', 'HTTP/REST')

# Add communication details
ax.text(10, 5.1, '← POST/GET/PATCH/DELETE >', fontsize=8, ha='center',
        color=colors['arrow_gray'], style='italic')

# ==================== BOTTOM SECTIONS ====================

# BaseLayer
draw_section_box(ax, 0.3, 2.5, 4.5, 2.5, '#FEE2E2', 'Compatibility Layer')
draw_rounded_box(ax, 0.6, 2.8, 4, 1.8, '#DC2626',
                 'BaseLayer', 'Qt-free stubs for legacy plugins\nSignal • Slot • QObject', fontsize=8)

# Database
draw_section_box(ax, 5.2, 2.5, 4.5, 2.5, colors['section_orange'], 'Database')
draw_rounded_box(ax, 5.5, 3.8, 1.9, 0.9, colors['database_amber'],
                 'SQLite', 'Desktop', fontsize=8, text_color='#1F2937')
draw_rounded_box(ax, 7.6, 3.8, 1.9, 0.9, colors['database_amber'],
                 'PostgreSQL', 'Cloud', fontsize=8, text_color='#1F2937')
ax.text(7.45, 2.9, 'Django ORM', fontsize=8, ha='center', color=colors['text_dark'])

# Job Execution
draw_section_box(ax, 10.1, 2.5, 4.8, 2.5, colors['section_green'], 'Job Execution')
draw_rounded_box(ax, 10.4, 3.5, 2.1, 1.2, colors['backend_green'],
                 'ProcessMgr', 'Subprocess', fontsize=8)
draw_rounded_box(ax, 12.7, 3.5, 2, 1.2, colors['backend_green'],
                 'i2run CLI', 'Task Runner', fontsize=8)
draw_arrow(ax, (12.5, 4.1), (12.7, 4.1), colors['arrow_gray'])
ax.text(12.35, 2.9, 'Spawn > Execute > Report', fontsize=7, ha='center',
        color=colors['text_dark'])

# Storage
draw_section_box(ax, 15.3, 2.5, 4.4, 2.5, colors['section_blue'], 'File Storage')
draw_rounded_box(ax, 15.6, 3.5, 1.9, 1.2, colors['frontend_blue'],
                 '~/.ccp4i2', 'Local', fontsize=8)
draw_rounded_box(ax, 17.7, 3.5, 1.8, 1.2, colors['frontend_blue'],
                 'Azure Files', 'Cloud', fontsize=8)
ax.text(17.5, 2.9, 'Projects • Media • Static', fontsize=7, ha='center',
        color=colors['text_dark'])

# ==================== LEGEND ====================

ax.text(0.5, 1.8, 'Data Flow', fontsize=10, fontweight='bold', color=colors['text_dark'])

# Flow description
flows = [
    '1. User action > React component',
    '2. useApi() hook > api-fetch.ts',
    '3. fetch() > /api/proxy/ccp4i2/...',
    '4. Next.js proxy > Django :8000',
    '5. ViewSet > ORM > Response',
    '6. SWR cache > UI update',
]
for i, flow in enumerate(flows):
    col = i // 3
    row = i % 3
    ax.text(0.5 + col * 4, 1.4 - row * 0.35, flow, fontsize=7, color=colors['text_dark'])

# Key features
ax.text(9, 1.8, 'Key Features', fontsize=10, fontweight='bold', color=colors['text_dark'])
features = [
    '• SWR client-side caching + polling',
    '• Azure AD JWT authentication',
    '• Namespaced APIs (/ccp4i2, /compounds)',
    '• Async job execution with progress',
    '• 176 crystallography task plugins',
    '• Qt-free via BaseLayer stubs',
]
for i, feat in enumerate(features):
    col = i // 3
    row = i % 3
    ax.text(9 + col * 4.5, 1.4 - row * 0.35, feat, fontsize=7, color=colors['text_dark'])

# Color legend
ax.text(17.5, 1.8, 'Colors', fontsize=10, fontweight='bold', color=colors['text_dark'])
legend_items = [
    (colors['frontend_blue'], 'Frontend'),
    (colors['backend_green'], 'Backend'),
    (colors['electron_purple'], 'Electron'),
]
for i, (color, label) in enumerate(legend_items):
    rect = Rectangle((17.5, 1.3 - i * 0.4), 0.25, 0.25, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(17.85, 1.42 - i * 0.4, label, fontsize=7, va='center', color=colors['text_dark'])

plt.tight_layout()
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/client_server_architecture.png',
            dpi=200, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/client_server_architecture.svg',
            format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: client_server_architecture.png and .svg")
