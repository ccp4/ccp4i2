import {
  Stack,
  Table,
  TableBody,
  TableCell,
  TableRow,
} from "@mui/material";

interface CCell {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

interface CResolutionRange {
  low: number;
  high: number;
}

/**
 * Digest data from MTZ/reflection files.
 * Supports both formats:
 * - Legacy: resolutionRange.low / resolutionRange.high
 * - Current: lowRes / highRes (from digest API)
 */
interface CObsData {
  cell?: CCell;
  spaceGroup?: string;
  resolutionRange?: CResolutionRange;
  // Alternative resolution format from digest API
  lowRes?: number;
  highRes?: number;
}

interface BaseSpacegroupCellElementProps {
  data?: CObsData;
}

export const BaseSpacegroupCellElement: React.FC<
  BaseSpacegroupCellElementProps
> = (props) => {
  return props.data?.cell ? (
    <Stack direction="column">
      <Table>
        <TableBody>
          <TableRow>
            <TableCell variant="head">Spacegroup</TableCell>
            <TableCell
              variant="body"
              key="Spacegroup"
              colSpan={6}
              sx={{ textAlign: "justify" }}
            >
              {typeof props.data?.spaceGroup === "string"
                ? JSON.stringify(props.data?.spaceGroup)
                : "?"}
            </TableCell>
          </TableRow>
          {props.data.cell && (
            <TableRow>
              <TableCell variant="head">Cell</TableCell>
              {["a", "b", "c", "alpha", "beta", "gamma"].map((key: string) => {
                if (props?.data?.cell) {
                  const { a, b, c, alpha, beta, gamma } = props.data.cell;
                  const cell: CCell = { a, b, c, alpha, beta, gamma };
                  const constrainedKey = key as
                    | "a"
                    | "b"
                    | "c"
                    | "alpha"
                    | "beta"
                    | "gamma";
                  return (
                    <TableCell variant="body" key={key}>
                      {key}={cell[constrainedKey].toPrecision(4)}
                    </TableCell>
                  );
                }
              })}
            </TableRow>
          )}
          <TableRow>
            <TableCell variant="head">Resolution</TableCell>
            <TableCell variant="body" key="low">
              {(() => {
                // Support both formats: resolutionRange.low (legacy) and lowRes (digest API)
                const low = props.data?.resolutionRange?.low ?? props.data?.lowRes;
                return typeof low === "number" ? low.toPrecision(4) : "?";
              })()}
            </TableCell>
            <TableCell variant="body" key="high">
              {(() => {
                // Support both formats: resolutionRange.high (legacy) and highRes (digest API)
                const high = props.data?.resolutionRange?.high ?? props.data?.highRes;
                return typeof high === "number" ? high.toPrecision(4) : "?";
              })()}
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </Stack>
  ) : (
    <div></div>
  );
};
