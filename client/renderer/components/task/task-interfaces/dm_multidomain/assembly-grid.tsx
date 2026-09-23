import React, { useCallback, useEffect, useMemo, useState } from "react";
import {
  Box,
  Button,
  Chip,
  IconButton,
  MenuItem,
  Select,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TextField,
  Tooltip,
  Typography,
} from "@mui/material";
import { Add, AutoFixHigh, Delete, Edit } from "@mui/icons-material";

import {
  Assembly,
  ChainInfo,
  IMPLICIT_ROLE,
  Model,
  addInstance,
  addRole,
  removeInstance,
  removeRole,
  renameRole,
  roleLabel,
  setCell,
} from "./dm-spec";

/**
 * The NCS assembly as a grid: one column per entity, one row per copy.
 *
 * The parameter underneath is a list of "CDK=A cyclin=B" strings, and asking
 * for it as text is what made this task need a paragraph of explanation. As a
 * grid the three load-bearing facts stop needing saying at all: the first row
 * is labelled "reference", a partial copy is a visible hole rather than an
 * omitted token, and a role is a column heading rather than a word that has to
 * match another word somewhere else. Every cell is a chain the model actually
 * has.
 */

interface AssemblyGridProps {
  model: Model;
  chains: ChainInfo[];
  entities: string[][];
  disabled?: boolean;
  onChange: (next: Model) => void;
  onDetect?: () => void;
  detectLabel?: string;
}

/** Column heading: the role name, renameable in place. */
const RoleHeader: React.FC<{
  role: string;
  assembly: Assembly;
  disabled?: boolean;
  canRemove: boolean;
  onRename: (name: string) => void;
  onRemove: () => void;
}> = ({ role, assembly, disabled, canRemove, onRename, onRemove }) => {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(role === IMPLICIT_ROLE ? "" : role);

  useEffect(() => {
    setDraft(role === IMPLICIT_ROLE ? "" : role);
  }, [role]);

  const commit = useCallback(() => {
    setEditing(false);
    const name = draft.trim();
    if (name && name !== role) onRename(name);
  }, [draft, role, onRename]);

  if (editing) {
    return (
      <TextField
        size="small"
        variant="standard"
        autoFocus
        value={draft}
        placeholder="CDK"
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => {
          if (e.key === "Enter") (e.target as HTMLInputElement).blur();
          if (e.key === "Escape") setEditing(false);
        }}
        sx={{ width: "8rem" }}
      />
    );
  }

  return (
    <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
      <Typography
        variant="body2"
        sx={{
          fontWeight: 600,
          fontStyle: role === IMPLICIT_ROLE ? "italic" : "normal",
        }}
      >
        {roleLabel(role, assembly)}
      </Typography>
      <Tooltip title="Name this entity (e.g. CDK) — the name is what rigid bodies refer to">
        <span>
          <IconButton size="small" disabled={disabled} onClick={() => setEditing(true)}>
            <Edit fontSize="inherit" />
          </IconButton>
        </span>
      </Tooltip>
      {canRemove && (
        <Tooltip title="Remove this entity from the assembly">
          <span>
            <IconButton size="small" color="error" disabled={disabled} onClick={onRemove}>
              <Delete fontSize="inherit" />
            </IconButton>
          </span>
        </Tooltip>
      )}
    </Box>
  );
};

export const AssemblyGrid: React.FC<AssemblyGridProps> = ({
  model,
  chains,
  entities,
  disabled,
  onChange,
  onDetect,
  detectLabel = "Detect from model",
}) => {
  const { assembly } = model;
  const roles = assembly.roles.length ? assembly.roles : [IMPLICIT_ROLE];

  /** Where else a chain is already used, so the menu can say why it is out. */
  const claimedBy = useCallback(
    (chain: string, instanceIndex: number, role: string): string | null => {
      for (let i = 0; i < assembly.instances.length; i += 1) {
        for (const [r, c] of Object.entries(assembly.instances[i])) {
          if (c !== chain) continue;
          if (i === instanceIndex && r === role) continue;
          return i === 0 ? "the reference copy" : `copy ${i + 1}`;
        }
      }
      return null;
    },
    [assembly.instances]
  );

  const entityOf = useCallback(
    (chain: string) => entities.findIndex((group) => group.includes(chain)),
    [entities]
  );

  // Names to offer the new column, in order: the model's own chain ids, so an
  // entity is called after something the user can see in the grid.
  const roleCandidates = useMemo(() => chains.map((c) => c.id), [chains]);

  return (
    <Box>
      <Box sx={{ display: "flex", alignItems: "center", gap: 1, mb: 1, flexWrap: "wrap" }}>
        <Typography variant="body2" color="text.secondary">
          The first row is the reference: masks are cut from it and every other
          copy is superposed onto it. Leave a cell empty where a copy is missing
          that entity.
        </Typography>
        {onDetect && (
          <Button
            size="small"
            startIcon={<AutoFixHigh />}
            onClick={onDetect}
            disabled={disabled}
          >
            {detectLabel}
          </Button>
        )}
      </Box>

      <Table size="small" sx={{ width: "auto" }}>
        <TableHead>
          <TableRow>
            <TableCell sx={{ borderBottom: "none" }} />
            {roles.map((role) => (
              <TableCell key={role} sx={{ borderBottom: "none" }}>
                <RoleHeader
                  role={role}
                  assembly={assembly}
                  disabled={disabled}
                  canRemove={roles.length > 1}
                  onRename={(name) => onChange(renameRole(model, role, name))}
                  onRemove={() => onChange(removeRole(model, role))}
                />
              </TableCell>
            ))}
            <TableCell sx={{ borderBottom: "none" }}>
              <Tooltip title="Add another entity — a second, different protein in each copy">
                <span>
                  <IconButton
                    size="small"
                    color="primary"
                    disabled={disabled}
                    onClick={() => onChange(addRole(model, roleCandidates))}
                    aria-label="Add entity"
                  >
                    <Add fontSize="small" />
                  </IconButton>
                </span>
              </Tooltip>
            </TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {assembly.instances.map((instance, index) => (
            <TableRow key={index}>
              <TableCell sx={{ borderBottom: "none", pr: 2 }}>
                {index === 0 ? (
                  <Chip size="small" color="primary" variant="outlined" label="reference" />
                ) : (
                  <Typography variant="body2" color="text.secondary">
                    copy {index + 1}
                  </Typography>
                )}
              </TableCell>
              {roles.map((role) => {
                const value = instance[role] ?? "";
                return (
                  <TableCell key={role} sx={{ borderBottom: "none" }}>
                    <Select
                      size="small"
                      variant="standard"
                      displayEmpty
                      value={value}
                      disabled={disabled}
                      onChange={(e) =>
                        onChange(setCell(model, index, role, String(e.target.value)))
                      }
                      sx={{ minWidth: "6rem" }}
                      renderValue={(chain) =>
                        chain ? (
                          String(chain)
                        ) : (
                          <Typography component="span" variant="body2" color="text.disabled">
                            —
                          </Typography>
                        )
                      }
                    >
                      <MenuItem value="">
                        <em>— not in this copy</em>
                      </MenuItem>
                      {chains.map((chain) => {
                        const claimed = claimedBy(chain.id, index, role);
                        const entity = entityOf(chain.id);
                        return (
                          <MenuItem
                            key={chain.id}
                            value={chain.id}
                            disabled={Boolean(claimed)}
                          >
                            chain {chain.id}
                            {chain.first !== null && chain.last !== null && (
                              <Typography
                                component="span"
                                variant="caption"
                                color="text.secondary"
                                sx={{ ml: 1 }}
                              >
                                {chain.first}–{chain.last}
                                {entity >= 0 ? ` · entity ${entity + 1}` : ""}
                              </Typography>
                            )}
                            {claimed && (
                              <Typography
                                component="span"
                                variant="caption"
                                color="text.secondary"
                                sx={{ ml: 1 }}
                              >
                                (already in {claimed})
                              </Typography>
                            )}
                          </MenuItem>
                        );
                      })}
                    </Select>
                  </TableCell>
                );
              })}
              <TableCell sx={{ borderBottom: "none" }}>
                <Tooltip title={index === 0 ? "Remove the reference copy" : `Remove copy ${index + 1}`}>
                  <span>
                    <IconButton
                      size="small"
                      color="error"
                      disabled={disabled}
                      onClick={() => onChange(removeInstance(model, index))}
                      aria-label={`Remove copy ${index + 1}`}
                    >
                      <Delete fontSize="small" />
                    </IconButton>
                  </span>
                </Tooltip>
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>

      <Button
        size="small"
        startIcon={<Add />}
        disabled={disabled}
        onClick={() => onChange(addInstance(model))}
        sx={{ mt: 0.5 }}
      >
        Add copy
      </Button>
    </Box>
  );
};

export default AssemblyGrid;
