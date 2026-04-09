'use client';

import { useState } from 'react';
import { useRouter } from 'next/navigation';
import {
  Button,
  Menu,
  MenuItem,
  ListItemIcon,
  ListItemText,
  Tooltip,
} from '@mui/material';
import {
  Add,
  GridOn,
  TableChart,
  Science,
} from '@mui/icons-material';
import { routes } from '@/lib/compounds/routes';

interface AddAssayMenuProps {
  /** Target ID to pre-select in import forms */
  targetId?: string;
  /** Whether the menu trigger is disabled (e.g., user lacks contributor permission) */
  disabled?: boolean;
}

export function AddAssayMenu({ targetId, disabled = false }: AddAssayMenuProps) {
  const router = useRouter();
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);

  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  const navigateTo = (path: string) => {
    handleClose();
    // Append target query param if provided
    const url = targetId ? `${path}${path.includes('?') ? '&' : '?'}target=${targetId}` : path;
    router.push(url);
  };

  return (
    <>
      <Tooltip title={disabled ? 'Requires Contributor or Admin operating level' : ''} arrow>
        <span>
          <Button
            variant="outlined"
            size="small"
            startIcon={<Add />}
            onClick={handleClick}
            disabled={disabled}
            aria-controls={open ? 'add-assay-menu' : undefined}
            aria-haspopup="true"
            aria-expanded={open ? 'true' : undefined}
          >
            Add Assay
          </Button>
        </span>
      </Tooltip>
      <Menu
        id="add-assay-menu"
        anchorEl={anchorEl}
        open={open}
        onClose={handleClose}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
        transformOrigin={{ vertical: 'top', horizontal: 'right' }}
      >
        <MenuItem onClick={() => navigateTo(routes.assays.import())}>
          <ListItemIcon>
            <GridOn fontSize="small" />
          </ListItemIcon>
          <ListItemText
            primary="Upload Plate Data"
            secondary="PHERAstar Excel files"
            secondaryTypographyProps={{ variant: 'caption' }}
          />
        </MenuItem>
        <MenuItem onClick={() => navigateTo(routes.assays.importTableOfValues())}>
          <ListItemIcon>
            <TableChart fontSize="small" />
          </ListItemIcon>
          <ListItemText
            primary="Import Table of Values"
            secondary="Pre-analyzed results"
            secondaryTypographyProps={{ variant: 'caption' }}
          />
        </MenuItem>
        <MenuItem onClick={() => navigateTo(routes.assays.importAdme())}>
          <ListItemIcon>
            <Science fontSize="small" />
          </ListItemIcon>
          <ListItemText
            primary="Import Pharmaron ADME"
            secondary="ADME property data"
            secondaryTypographyProps={{ variant: 'caption' }}
          />
        </MenuItem>
      </Menu>
    </>
  );
}
