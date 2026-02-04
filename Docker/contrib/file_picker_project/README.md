# filepicker

A tool for selecting and transferring crystallography data files from XChem-style fragment screening datasets, with built-in support for remote access to Diamond Light Source via SFTP.

## Features

- **Remote SFTP access**: Connect directly to Diamond Light Source (`nx.diamond.ac.uk`) without needing to install anything on the remote server
- **Flexible crystal selection**: by index (`x0001`), ranges (`x0001-x0010`), or wildcards (`*`)
- **Template-based file selection**: predefined templates for common use cases, or custom patterns
- **Symlink dereferencing**: automatically resolves symlinks to download actual file content
- **Sibling file discovery**: grabs related files (unmerged data, xia2.mmcif.bz2) from autoprocessing directories
- **Progress reporting**: real-time download progress with speed and ETA

## Installation

```bash
# Create a virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install paramiko

# Install filepicker
pip install -e /path/to/file_picker_project

# Or run directly without installing
cd /path/to/file_picker_project
python -m filepicker --help
```

## Quick Start (Diamond Light Source)

```bash
# Show available crystals
filepicker fedid@nx.diamond.ac.uk:/dls/i04-1/data/2024/sw12345-6/processing/analysis/model_building info

# List crystals in a range
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building list -c x0001-x0010

# Preview files that would be downloaded
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building resolve -c x0001 x0002 -v

# Download selected crystals
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building pull -c x0001-x0050 -d ./local_data

# Download with manifest and skip existing files
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building pull -c '*' -d ./data --skip-existing -m manifest.txt

# Dry run (show what would be downloaded)
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building pull -c x0001-x0100 -d ./data --dry-run

# Two-stage workflow: generate manifest first, then download from it
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building manifest -c '*' -o transfer
filepicker fedid@nx.diamond.ac.uk:/dls/.../model_building pull --from-manifest transfer_manifest.txt -d ./data
```

## Commands

### Remote Commands

#### `info`
Show information about the dataset:
```bash
filepicker fedid@nx.diamond.ac.uk:/path/to/model_building info
```

#### `list`
List crystals matching selection:
```bash
filepicker fedid@nx.diamond.ac.uk:/path/to/model_building list -c x0001-x0010
```

#### `resolve`
Preview file resolution (shows what files would be selected):
```bash
filepicker fedid@nx.diamond.ac.uk:/path/to/model_building resolve -c x0001 -v
```

#### `pull`
Download files from remote server:
```bash
filepicker fedid@nx.diamond.ac.uk:/path/to/model_building pull \
    -c x0001-x0100 \
    -d ./local_data \
    --skip-existing
```

Options:
- `-c, --crystals`: Crystal specifiers (e.g., `x0001-x0010`, `*`)
- `-d, --dest`: Local destination directory (required)
- `-t, --template`: File template (`default`, `minimal`, `full`)
- `-m, --manifest`: Write transfer manifest to file
- `--from-manifest`: Read files from existing manifest (skips SFTP resolution)
- `--skip-existing`: Skip files that already exist with same size
- `-n, --dry-run`: Show what would be downloaded without downloading
- `-y, --yes`: Skip confirmation for large downloads (>1GB)
- `-q, --quiet`: Suppress progress output

#### `manifest`
Generate a transfer manifest without downloading:
```bash
filepicker fedid@nx.diamond.ac.uk:/path/to/model_building manifest \
    -c '*' \
    -o transfer \
    --rsync ./local_dest
```

### Local Commands

#### `archive`
Create a compressed archive (only for local data):
```bash
filepicker /local/model_building archive \
    -c x0001-x0100 \
    -o dataset.tar.gz
```

## Crystal Specifiers

| Specifier | Description |
|-----------|-------------|
| `x0001` | Single crystal by index |
| `x0001-x0010` | Range (inclusive) |
| `x0001:x0010` | Range (alternative syntax) |
| `*` or `all` | All crystals |
| `!x0005` | Exclude (combine with other specifiers) |

## File Templates

### Default template
Includes key files for each crystal:
- `dimple/dimple/final.pdb` and `final.mtz` - refined model
- `compound/*.cif`, `*.smiles`, `*.png` - ligand files
- Top-level `*.mtz` (dereferenced) - selected merged data
- Sibling files from autoprocessing: `*_scaled_unmerged.mtz`, `xia2.mmcif.bz2`
- Electron density maps (`2fofc.map`, `fofc.map`)

### Minimal template
Just the essentials:
- `final.pdb` and `final.mtz`
- `compound/*.cif` and `*.smiles`

### Custom patterns
```bash
filepicker ... --files "final.pdb,compound/*.cif,*.mtz"
```

## Authentication

The tool uses standard SSH authentication:

1. **SSH Agent**: If you have an SSH agent running with your key loaded, it will be used automatically
2. **SSH Keys**: Looks for keys in `~/.ssh/`
3. **Password**: If keys aren't available, you'll be prompted for your password

For Diamond, use your FedID as the username:
```bash
filepicker fedid@nx.diamond.ac.uk:/dls/... info
```

## Data Retention Note

Diamond Light Source keeps data on disk for **40 days** after your experiment. After that, data is archived and must be retrieved via DataGateway.

## Requirements

- Python 3.6+
- paramiko (for SFTP operations)

## Troubleshooting

### "paramiko is required for remote operations"
Install paramiko:
```bash
pip install paramiko
```

### Connection timeouts
Diamond's SSH service may be slow for large operations. Consider using `--skip-existing` for resumable downloads.

### Permission denied
Make sure you're using your Diamond FedID and that you have access to the visit data.
