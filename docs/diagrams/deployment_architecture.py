#!/usr/bin/env python3
"""
CCP4i2 Azure Deployment Architecture Diagram
Generates a professional diagram showing the 3-layer container architecture,
Azure services, and networking.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
import matplotlib.lines as mlines
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(20, 14))
ax.set_xlim(0, 20)
ax.set_ylim(0, 14)
ax.set_aspect('equal')
ax.axis('off')

# Color palette
colors = {
    'azure_blue': '#0078D4',
    'container_green': '#10B981',
    'storage_orange': '#F59E0B',
    'database_purple': '#8B5CF6',
    'security_red': '#EF4444',
    'network_gray': '#6B7280',
    'background_light': '#F3F4F6',
    'text_dark': '#1F2937',
    'white': '#FFFFFF',
    'layer1': '#DBEAFE',  # Light blue
    'layer2': '#D1FAE5',  # Light green
    'layer3': '#FEF3C7',  # Light yellow
}

def draw_rounded_box(ax, x, y, width, height, color, label, sublabel=None,
                     text_color='white', fontsize=9, icon=None):
    """Draw a rounded rectangle with label"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.1",
                         facecolor=color, edgecolor='none', alpha=0.95)
    ax.add_patch(box)

    # Add icon if provided
    if icon:
        ax.text(x + width/2, y + height*0.7, icon, fontsize=16,
                ha='center', va='center', color=text_color)
        label_y = y + height*0.35
    else:
        label_y = y + height*0.6 if sublabel else y + height/2

    ax.text(x + width/2, label_y, label, fontsize=fontsize,
            ha='center', va='center', color=text_color, fontweight='bold')

    if sublabel:
        ax.text(x + width/2, y + height*0.25, sublabel, fontsize=fontsize-2,
                ha='center', va='center', color=text_color, alpha=0.9)

def draw_section_box(ax, x, y, width, height, color, title):
    """Draw a section box with title"""
    box = FancyBboxPatch((x, y), width, height,
                         boxstyle="round,pad=0.01,rounding_size=0.15",
                         facecolor=color, edgecolor='#9CA3AF',
                         linewidth=1.5, alpha=0.4)
    ax.add_patch(box)
    ax.text(x + 0.15, y + height - 0.25, title, fontsize=10,
            ha='left', va='top', color=colors['text_dark'], fontweight='bold')

def draw_arrow(ax, start, end, color='#6B7280', style='->'):
    """Draw an arrow between two points"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle=style, color=color, lw=1.5))

# Title
ax.text(10, 13.7, 'CCP4i2 Azure Deployment Architecture', fontsize=18,
        ha='center', va='center', fontweight='bold', color=colors['text_dark'])
ax.text(10, 13.3, 'UK South Region • Container Apps Environment', fontsize=11,
        ha='center', va='center', color=colors['network_gray'])

# ==================== LEFT SIDE: Container Build Layers ====================

# Section: Container Image Layers
draw_section_box(ax, 0.3, 6.5, 5.4, 6.5, colors['layer1'], 'Container Image Layers')

# Layer 1: CCP4 Base
draw_rounded_box(ax, 0.5, 11.5, 5, 1.2, '#3B82F6', 'Layer 1: CCP4 Base',
                 'ccp4i2/base:ccp4-20251105 (~10GB)', fontsize=9)
ax.text(5.7, 12.1, '• Python 3.11 slim\n• CCP4 suite\n• System deps',
        fontsize=7, va='center', color=colors['text_dark'])

# Layer 2: ARP/wARP
draw_rounded_box(ax, 0.5, 9.8, 5, 1.2, '#10B981', 'Layer 2: ARP/wARP',
                 'ccp4i2/base-arpwarp:ccp4-20251105', fontsize=9)
ax.text(5.7, 10.4, '• ARP/wARP tools\n• Builds on Layer 1',
        fontsize=7, va='center', color=colors['text_dark'])

# Layer 3: Application
draw_rounded_box(ax, 0.5, 8.1, 5, 1.2, '#F59E0B', 'Layer 3: Application',
                 'ccp4i2/server & ccp4i2/web (~500MB)', fontsize=9, text_color='#1F2937')
ax.text(5.7, 8.7, '• Django/Next.js\n• Fast rebuild',
        fontsize=7, va='center', color=colors['text_dark'])

# Arrows between layers
draw_arrow(ax, (3, 11.5), (3, 11.0), '#6B7280')
draw_arrow(ax, (3, 9.8), (3, 9.3), '#6B7280')

# ACR Box
draw_rounded_box(ax, 0.8, 6.7, 4.4, 1.1, colors['azure_blue'],
                 'Azure Container Registry', 'ccp4acrukbwmx.azurecr.io', fontsize=8)

# ==================== CENTER: Container Apps ====================

# VNet Section
draw_section_box(ax, 6, 3.5, 8, 9.5, '#E5E7EB', 'Virtual Network (10.0.0.0/16)')

# Container Apps Environment
draw_section_box(ax, 6.3, 4, 7.4, 7.8, colors['layer2'], 'Container Apps Environment')

# Web Container
draw_rounded_box(ax, 6.6, 9.8, 3.2, 1.6, colors['container_green'],
                 'Web App', 'Next.js - Port 3000\n1-5 replicas (HTTP)', fontsize=9)
ax.text(6.6, 9.6, 'External', fontsize=7, color='#059669', fontweight='bold')

# Server Container
draw_rounded_box(ax, 10.2, 9.8, 3.2, 1.6, colors['container_green'],
                 'Server App', 'Django - Port 8000\n1-10 replicas (CPU/HTTP)', fontsize=9)
ax.text(10.2, 9.6, 'Internal', fontsize=7, color='#6B7280', fontweight='bold')

# Worker Container
draw_rounded_box(ax, 8.4, 7.5, 3.2, 1.6, colors['container_green'],
                 'Worker App', 'Background Jobs\n0-20 replicas (Queue)', fontsize=9)
ax.text(8.4, 7.3, 'No Ingress', fontsize=7, color='#6B7280', fontweight='bold')

# Arrows within containers
draw_arrow(ax, (8.2, 10.6), (10.2, 10.6), '#374151')  # Web -> Server
ax.text(9.2, 10.85, 'API Proxy', fontsize=7, ha='center', color='#374151')

draw_arrow(ax, (11.8, 9.8), (10.8, 9.1), '#374151')  # Server -> Worker
ax.text(11.5, 9.35, 'Queue', fontsize=7, ha='center', color='#374151')

# Private Endpoints Section
draw_section_box(ax, 6.5, 4.2, 6.8, 2.8, '#FEE2E2', 'Private Endpoints')

# Private DNS
ax.text(9.9, 6.7, 'Private DNS Zones', fontsize=8, ha='center',
        color=colors['security_red'], fontweight='bold')
ax.text(9.9, 6.4, 'privatelink.*.azure.net', fontsize=7, ha='center',
        color=colors['text_dark'])

# Small endpoint boxes
endpoints = [
    ('PostgreSQL', 6.7, 5.2),
    ('Storage', 8.3, 5.2),
    ('Key Vault', 9.9, 5.2),
    ('Service Bus', 11.5, 5.2),
]
for name, x, y in endpoints:
    draw_rounded_box(ax, x, y, 1.4, 0.7, '#DC2626', name, fontsize=7)

# ==================== RIGHT SIDE: Azure Services ====================

# Section: Azure Services
draw_section_box(ax, 14.5, 3.5, 5.2, 9.5, colors['layer3'], 'Azure Services')

# PostgreSQL
draw_rounded_box(ax, 14.8, 11, 4.6, 1.4, colors['database_purple'],
                 'PostgreSQL', 'Flexible Server v15\nStandard_B1ms - 32GB', fontsize=9)

# Storage Account
draw_rounded_box(ax, 14.8, 9.2, 4.6, 1.4, colors['storage_orange'],
                 'Storage Account', 'Projects - Media - Static\nPrivate Endpoints',
                 fontsize=9, text_color='#1F2937')

# Key Vault
draw_rounded_box(ax, 14.8, 7.4, 4.6, 1.4, colors['security_red'],
                 'Key Vault', 'Secrets & Credentials\nManaged Identity RBAC', fontsize=9)

# Service Bus
draw_rounded_box(ax, 14.8, 5.6, 4.6, 1.4, colors['azure_blue'],
                 'Service Bus', 'Job Queue (Premium)\nScale: 5+ msgs -> workers', fontsize=9)

# Azure AD
draw_rounded_box(ax, 14.8, 3.8, 4.6, 1.4, '#0078D4',
                 'Azure AD', 'OAuth 2.0 Authentication\nJWT Token Validation', fontsize=9)

# ==================== TOP: Internet/User ====================

# Internet cloud
ax.text(8.2, 13, 'Internet', fontsize=10, ha='center',
        fontweight='bold', color=colors['text_dark'])
draw_arrow(ax, (8.2, 12.7), (8.2, 11.4), colors['azure_blue'], '->')
ax.text(8.5, 12, 'HTTPS', fontsize=8, color=colors['azure_blue'])

# ==================== BOTTOM: Legend ====================

# Legend
legend_y = 1.8
ax.text(1, legend_y + 1, 'Legend', fontsize=10, fontweight='bold', color=colors['text_dark'])

# Legend items
legend_items = [
    (colors['layer1'], 'Image Build Layers'),
    (colors['container_green'], 'Container Apps'),
    (colors['storage_orange'], 'Storage'),
    (colors['database_purple'], 'Database'),
    (colors['security_red'], 'Security/Networking'),
    (colors['azure_blue'], 'Azure Services'),
]

for i, (color, label) in enumerate(legend_items):
    x = 1 + (i % 3) * 3.5
    y = legend_y - (i // 3) * 0.5
    rect = Rectangle((x, y - 0.15), 0.3, 0.3, facecolor=color, edgecolor='none')
    ax.add_patch(rect)
    ax.text(x + 0.4, y, label, fontsize=8, va='center', color=colors['text_dark'])

# Key info
ax.text(12, legend_y + 0.5, 'Key Configuration', fontsize=10, fontweight='bold',
        color=colors['text_dark'])
ax.text(12, legend_y, '• Region: UK South', fontsize=8, color=colors['text_dark'])
ax.text(12, legend_y - 0.4, '• Resource Group: ccp4i2-bicep-rg-uksouth',
        fontsize=8, color=colors['text_dark'])
ax.text(12, legend_y - 0.8, '• ACR: ccp4acrukbwmx.azurecr.io',
        fontsize=8, color=colors['text_dark'])

# Scaling info
ax.text(12, 0.8, 'Auto-Scaling', fontsize=10, fontweight='bold', color=colors['text_dark'])
ax.text(12, 0.4, 'Web: 1->5 (HTTP) | Server: 1->10 (CPU) | Worker: 0->20 (Queue)',
        fontsize=8, color=colors['text_dark'])

plt.tight_layout()
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/deployment_architecture.png',
            dpi=200, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/Users/nmemn/Developer/ccp4i2/docs/diagrams/deployment_architecture.svg',
            format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
print("Saved: deployment_architecture.png and .svg")
