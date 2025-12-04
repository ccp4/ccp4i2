import { moorhen } from "moorhen/types/moorhen";
import {
  hideMolecule,
  showMolecule,
  hideMap,
  showMap,
  removeMap,
  removeMolecule,
  setRequestDrawScene,
} from "moorhen";
import { useDispatch, useSelector } from "react-redux";
import {
  Box,
  List,
  ListItem,
  ListItemButton,
  ListItemText,
  ListItemSecondaryAction,
  IconButton,
  Menu,
  MenuItem,
  Typography,
  Paper,
  Chip,
  Stack,
  Dialog,
} from "@mui/material";
import { MoreVert, Visibility, VisibilityOff } from "@mui/icons-material";
import { useEffect, useMemo, useState } from "react";
import { useApi } from "../../api";
import { useTheme } from "../../theme/theme-provider";
import { PushToCCP4i2Panel } from "./push-to-ccp4i2-panel";
import {
  extractFileId,
  fetchItemMetadata,
  ItemMetadata,
} from "./item-metadata-utils"; // <-- Import utilities

type ContentType = "Molecule" | "Map";
type ContentItem = moorhen.Molecule | moorhen.Map;

interface ItemMenuState {
  anchorEl: HTMLElement | null;
  item: ContentItem | null;
}

interface MoorhenLoadedContentProps {
  onFileSelect: (fileId: number) => void;
  type: ContentType;
}

// Utility to get favicon URL
function getFaviconUrl(): string | undefined {
  const link = document.querySelector("link[rel~='icon']");
  const href = link ? link.getAttribute("href") : undefined;
  return href === null ? undefined : href;
}

export const MoorhenLoadedContent: React.FC<MoorhenLoadedContentProps> = ({
  onFileSelect,
  type = "Molecule",
}) => {
  const { customColors } = useTheme();
  const [menuState, setMenuState] = useState<ItemMenuState>({
    anchorEl: null,
    item: null,
  });
  const [itemMetadata, setItemMetadata] = useState<Map<number, ItemMetadata>>(
    new Map()
  );
  const [pushDialogOpen, setPushDialogOpen] = useState(false);
  const [faviconUrl, setFaviconUrl] = useState<string | undefined>(undefined);

  const dispatch = useDispatch();
  const api = useApi();

  const cootInitialized = useSelector(
    (state: moorhen.State) => state.generalStates.cootInitialized
  );

  // Molecule selectors
  const molecules = useSelector(
    (state: moorhen.State) => state.molecules.moleculeList
  );
  const visibleMolecules = useSelector(
    (state: moorhen.State) => state.molecules.visibleMolecules
  );

  // Map selectors
  const maps = useSelector((state: moorhen.State) => state.maps);
  const visibleMaps = useSelector(
    (state: moorhen.State) => state.mapContourSettings.visibleMaps
  );

  // Get the appropriate data based on type
  const items = type === "Molecule" ? molecules : maps;
  const visibleItems = type === "Molecule" ? visibleMolecules : visibleMaps;

  // Check if item was loaded from database
  const isFromDatabase = (item: ContentItem): boolean => {
    return !!(item.uniqueId && extractFileId(item.uniqueId));
  };

  // Fetch metadata for a specific item using the utility
  const fetchAndStoreItemMetadata = async (item: ContentItem) => {
    const fileId = extractFileId(item.uniqueId || "");
    if (!fileId) return;

    const molNo = item.molNo;

    // Set loading state
    setItemMetadata(
      (prev) =>
        new Map(
          prev.set(molNo, {
            fileId,
            isLoading: true,
          })
        )
    );

    try {
      const metadata = await fetchItemMetadata(item);
      if (metadata) {
        setItemMetadata(
          (prev) => new Map(prev.set(molNo, { ...metadata, isLoading: false }))
        );
      }
    } catch (error) {
      setItemMetadata(
        (prev) =>
          new Map(
            prev.set(molNo, {
              fileId,
              isLoading: false,
              error: error instanceof Error ? error.message : "Unknown error",
            })
          )
      );
    }
  };

  // Fetch metadata for all items when they change
  useEffect(() => {
    if (!items || items.length === 0) return;

    items.forEach((item) => {
      if (item.uniqueId) {
        const fileId = extractFileId(item.uniqueId);
        if (fileId && !itemMetadata.has(item.molNo)) {
          fetchAndStoreItemMetadata(item);
        }
      }
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [items]);

  useEffect(() => {
    setFaviconUrl(getFaviconUrl());
  }, []);

  const handleMenuOpen = (
    event: React.MouseEvent<HTMLElement>,
    item: ContentItem
  ) => {
    event.stopPropagation();
    setMenuState({
      anchorEl: event.currentTarget,
      item: item,
    });
  };

  const handleContextMenu = (
    event: React.MouseEvent<HTMLElement>,
    item: ContentItem
  ) => {
    event.preventDefault();
    setMenuState({
      anchorEl: event.currentTarget,
      item: item,
    });
  };

  const handleMenuClose = () => {
    setMenuState({
      anchorEl: null,
      item: null,
    });
  };

  const handleHideItem = () => {
    if (menuState.item) {
      if (type === "Molecule") {
        const molecule = menuState.item as moorhen.Molecule;
        dispatch(hideMolecule({ molNo: molecule.molNo }));
        dispatch(setRequestDrawScene(true));
      } else {
        const map = menuState.item as moorhen.Map;
        map.hideMapContour();
        dispatch(hideMap({ molNo: map.molNo }));
        dispatch(setRequestDrawScene(true));
      }
    }
    handleMenuClose();
  };

  const handleShowItem = () => {
    if (menuState.item) {
      if (type === "Molecule") {
        const molecule = menuState.item as moorhen.Molecule;
        dispatch(showMolecule({ molNo: molecule.molNo, show: true }));
        dispatch(setRequestDrawScene(true));
      } else {
        const map = menuState.item as moorhen.Map;
        dispatch(showMap({ molNo: map.molNo, show: true }));
        dispatch(setRequestDrawScene(true));
      }
    }
    handleMenuClose();
  };

  const handleDeleteItem = () => {
    if (menuState.item) {
      if (type === "Map") {
        dispatch(removeMap(menuState.item as moorhen.Map));
        menuState.item.delete();
        dispatch(setRequestDrawScene(true));
        console.log(`Delete map:`, menuState.item.name);
      } else {
        dispatch(removeMolecule(menuState.item as moorhen.Molecule));
        menuState.item.delete();
        dispatch(setRequestDrawScene(true));
        console.log(`Delete ${type.toLowerCase()}:`, menuState.item.name);
      }
    }
    handleMenuClose();
  };

  const handleCenterOnItem = () => {
    if (menuState.item) {
      if (type === "Molecule") {
        const molecule = menuState.item as moorhen.Molecule;
        molecule.centreOn("/*/*/*/*", false, true);
      } else {
        const map = menuState.item as moorhen.Map;
        map.centreOnMap();
        //console.log("Center on map:", map.name);
      }
    }
    handleMenuClose();
  };

  const handleItemClick = (item: ContentItem) => {
    if (isVisible(item)) {
      if (type === "Molecule") {
        const molecule = item as moorhen.Molecule;
        dispatch(hideMolecule({ molNo: molecule.molNo }));
      } else {
        const map = item as moorhen.Map;
        dispatch(hideMap({ molNo: map.molNo }));
      }
    } else {
      if (type === "Molecule") {
        const molecule = item as moorhen.Molecule;
        dispatch(showMolecule({ molNo: molecule.molNo, show: true }));
      } else {
        const map = item as moorhen.Map;
        dispatch(showMap({ molNo: map.molNo, show: true }));
      }
    }
  };

  const isVisible = (item: ContentItem) => {
    if (type === "Molecule") {
      const molecule = item as moorhen.Molecule;
      return visibleMolecules.includes(molecule.molNo);
    } else {
      const map = item as moorhen.Map;
      return visibleMaps.includes(map.molNo);
    }
  };

  const getItemId = (item: ContentItem) => {
    return item.molNo;
  };

  const getItemMetadata = (item: ContentItem): ItemMetadata | undefined => {
    return itemMetadata.get(item.molNo);
  };

  // Get the primary display text for an item
  const getPrimaryDisplayText = (item: ContentItem): React.ReactNode => {
    const metadata = getItemMetadata(item);
    const fromDatabase = isFromDatabase(item);

    if (fromDatabase && metadata) {
      if (metadata.isLoading) {
        return (
          <Typography variant="body2" color="text.secondary">
            Loading...
          </Typography>
        );
      }

      if (metadata.error) {
        return (
          <Typography variant="body2" color="error">
            Error loading metadata
          </Typography>
        );
      }

      if (metadata.projectName && metadata.jobNumber) {
        return (
          <Stack
            direction="row"
            alignItems="center"
            spacing={0.5}
            sx={{ flexWrap: "wrap" }}
          >
            <Typography
              variant="body2"
              color="primary.main"
              sx={{ fontWeight: "medium" }}
            >
              📁 {metadata.projectName}
            </Typography>
            <Typography variant="body2" color="text.secondary">
              •
            </Typography>
            <Typography variant="body2" color="text.secondary">
              Job {metadata.jobNumber}
            </Typography>
            {metadata.fileAnnotation && (
              <>
                <Typography variant="body2" color="text.secondary">
                  •
                </Typography>
                <Typography
                  variant="body2"
                  color="text.secondary"
                  sx={{ fontStyle: "italic" }}
                >
                  {metadata.fileAnnotation}
                </Typography>
              </>
            )}
          </Stack>
        );
      }
    }

    // Fallback to molecule/map name for non-database items or when metadata is unavailable
    return (
      <Typography
        variant="body2"
        sx={{
          fontStyle: fromDatabase ? "normal" : "italic",
          color: fromDatabase ? "text.primary" : "text.secondary",
        }}
      >
        {item.name}
        {!fromDatabase && (
          <Typography
            component="span"
            variant="caption"
            sx={{ ml: 1, color: "text.disabled" }}
          >
            (external file)
          </Typography>
        )}
      </Typography>
    );
  };

  const itemToPush = useMemo(() => {
    return menuState.item;
  }, [menuState]);
  const molNo = itemToPush?.molNo;
  const handlePushToCCP4i2 = () => {
    setPushDialogOpen(true);
    //handleMenuClose();
  };

  const handlePushDialogClose = () => {
    setPushDialogOpen(false);
  };

  if (!cootInitialized) {
    return (
      <Box
        sx={{
          height: "100%",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
        }}
      >
        <Typography variant="body2" color="text.secondary">
          Loading...
        </Typography>
      </Box>
    );
  }

  if (!items || items.length === 0) {
    return (
      <Box
        sx={{
          height: "100%",
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
        }}
      >
        <Typography variant="body2" color="text.secondary">
          No {type.toLowerCase()}s loaded.
        </Typography>
      </Box>
    );
  }

  return (
    <Box
      sx={{
        height: "100%",
        width: "100%",
        display: "flex",
        flexDirection: "column",
      }}
    >
      <Paper
        sx={{
          flex: 1,
          overflow: "auto",
          boxShadow: "none",
          border: `1px solid ${customColors.ui.mediumGray}`,
        }}
      >
        <List
          sx={{
            width: "100%",
            bgcolor: "background.paper",
            padding: 0,
          }}
        >
          {items.map((item, index) => {
            const itemIsVisible = isVisible(item);
            return (
              <ListItem
                key={getItemId(item)}
                disablePadding
                sx={{
                  borderBottom:
                    index < items.length - 1 ? "1px solid #f0f0f0" : "none",
                }}
              >
                <ListItemButton
                  onClick={() => handleItemClick(item)}
                  onContextMenu={(event) => handleContextMenu(event, item)}
                  sx={{
                    paddingY: 1,
                    paddingX: 2,
                    "&:hover": {
                      backgroundColor: customColors.ui.lightGray,
                    },
                    opacity: itemIsVisible ? 1 : 0.7,
                  }}
                >
                  <ListItemText
                    primary={
                      <Stack direction="row" alignItems="center" spacing={1}>
                        {/* Visibility Eye Icon */}
                        <Box
                          sx={{
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                            minWidth: "24px",
                            height: "24px",
                            flexShrink: 0,
                          }}
                        >
                          {itemIsVisible ? (
                            <Visibility
                              fontSize="small"
                              sx={{
                                color: "primary.main",
                                opacity: 0.8,
                              }}
                            />
                          ) : (
                            <VisibilityOff
                              fontSize="small"
                              sx={{
                                color: "text.disabled",
                                opacity: 0.6,
                              }}
                            />
                          )}
                        </Box>

                        {/* Molecule/Map ID Chip */}
                        <Chip
                          label={getItemId(item)}
                          size="small"
                          variant="outlined"
                          sx={{
                            fontSize: "0.75rem",
                            fontFamily: "monospace",
                            minWidth: "40px",
                            height: "20px",
                            flexShrink: 0,
                            opacity: itemIsVisible ? 1 : 0.6,
                          }}
                        />

                        {/* Content Text */}
                        <Box sx={{ flex: 1, minWidth: 0 }}>
                          {getPrimaryDisplayText(item)}
                        </Box>
                      </Stack>
                    }
                  />
                  <ListItemSecondaryAction>
                    <IconButton
                      edge="end"
                      size="small"
                      onClick={(event) => handleMenuOpen(event, item)}
                      sx={{
                        marginRight: 1,
                        opacity: 0.7,
                        "&:hover": {
                          opacity: 1,
                        },
                      }}
                    >
                      <MoreVert fontSize="small" />
                    </IconButton>
                  </ListItemSecondaryAction>
                </ListItemButton>
              </ListItem>
            );
          })}
        </List>
      </Paper>

      {/* Contextual Menu */}
      <Menu
        anchorEl={menuState.anchorEl}
        open={Boolean(menuState.anchorEl)}
        onClose={handleMenuClose}
        slotProps={{
          paper: {
            sx: {
              minWidth: "180px",
              boxShadow: "0 4px 12px rgba(0, 0, 0, 0.15)",
            },
          },
        }}
      >
        {menuState.item && isVisible(menuState.item) ? (
          <MenuItem onClick={handleHideItem}>
            <VisibilityOff sx={{ mr: 1 }} fontSize="small" />
            Hide {type}
          </MenuItem>
        ) : (
          <MenuItem onClick={handleShowItem}>
            <Visibility sx={{ mr: 1 }} fontSize="small" />
            Show {type}
          </MenuItem>
        )}
        <MenuItem onClick={handleCenterOnItem}>
          <Typography sx={{ mr: 1 }}>🎯</Typography>
          Centre on {type}
        </MenuItem>
        <MenuItem onClick={handleDeleteItem} sx={{ color: "error.main" }}>
          <Typography sx={{ mr: 1 }}>🗑️</Typography>
          Delete {type}
        </MenuItem>
        <MenuItem onClick={handlePushToCCP4i2}>
          {faviconUrl ? (
            <img
              src={faviconUrl}
              alt="CCP4i2"
              style={{
                width: 20,
                height: 20,
                marginRight: 8,
                verticalAlign: "middle",
              }}
            />
          ) : (
            <Typography sx={{ mr: 1 }}>🚀</Typography>
          )}
          Push to CCP4i2
        </MenuItem>
      </Menu>

      {/* Push to CCP4i2 Dialog */}
      <Dialog
        open={pushDialogOpen}
        onClose={handlePushDialogClose}
        maxWidth="md"
        fullWidth
      >
        <PushToCCP4i2Panel
          molNo={molNo}
          item={itemToPush}
          itemMetadata={itemMetadata.get(itemToPush?.molNo || 0)}
          onClose={handlePushDialogClose}
        />
      </Dialog>
    </Box>
  );
};
