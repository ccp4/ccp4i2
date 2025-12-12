import { Editor } from "@monaco-editor/react";
import { useEffect, useMemo, useRef, useState } from "react";
import { CCP4Table, parseXML, Plot } from "./ChartLib";

import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  BarElement,
  BarController,
  ScatterController,
  LineElement,
  Title,
  Tooltip,
  Legend,
  ChartOptions,
  ChartData,
} from "chart.js";
import { Scatter } from "react-chartjs-2";

import annotationPlugin from "chartjs-plugin-annotation";
import {
  Autocomplete,
  Button,
  Card,
  CardContent,
  CardHeader,
  Dialog,
  DialogContent,
  TextField,
  Typography,
} from "@mui/material";
import { OpenInNew, Preview } from "@mui/icons-material";
import { useTheme } from "../../theme/theme-provider";

// Register components for both Scatter and Bar charts
ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  BarElement,
  BarController,
  ScatterController,
  LineElement,
  Title,
  Tooltip,
  Legend,
  annotationPlugin
);

interface CCP4i2ApplicationOutputViewProps {
  output: Element;
  height?: string | number;
  jobId?: number;
  graphId?: string;
}

const colours = [
  "rgba(255, 99, 132, 1)",
  "rgba(54, 162, 235, 1)",
  "rgba(255, 206, 86, 1)",
  "rgba(75, 192, 192, 1)",
  "rgba(153, 102, 255, 1)",
  "rgba(255, 159, 64, 1)",
];
interface ChartArgs {
  options: ChartOptions;
  data: any;
}
export const CCP4i2ApplicationOutputView: React.FC<
  CCP4i2ApplicationOutputViewProps
> = ({ output, height = "400px", jobId, graphId }) => {
  const [table, setTable] = useState<CCP4Table | null>(null);
  const [selectedPlot, setSelectedPlot] = useState<Plot | null>(null);
  const [showJson, setShowJson] = useState(false);
  const { mode } = useTheme();

  const chartRef = useRef(null);

  useEffect(() => {
    const asyncEffect = async () => {
      const serializedOutput = new XMLSerializer().serializeToString(output);
      try {
        const parsedOutput = await parseXML(serializedOutput);
        if (
          parsedOutput &&
          parsedOutput.CCP4Table &&
          !Array.isArray(parsedOutput.CCP4Table)
        ) {
          const graph = new CCP4Table(parsedOutput.CCP4Table);
          if (graph) {
            setTable(graph);
            setSelectedPlot(
              Array.isArray(graph.plot)
                ? graph.plot[0]
                : graph.plot
                  ? graph.plot
                  : null
            );
          }
        }
      } catch (err) {
        console.log(err);
      }
    };
    asyncEffect();
  }, [output]);

  const allPlots = useMemo(() => {
    return table ? table.allPlots : [];
  }, [table]);

  const chartArgs: ChartArgs | null = useMemo(() => {
    if (selectedPlot && table) {
      const options: ChartOptions | null = table.plotOptions(selectedPlot);
      const data: ChartData<"bar" | "scatter"> | null =
        table.plotData(selectedPlot);
      if (!options || !data) {
        return null;
      }
      return { data, options };
    }
    return null;
  }, [selectedPlot, table]);

  return (
    <>
      <Card variant="outlined">
        <CardHeader
          title={
            allPlots &&
            selectedPlot &&
            allPlots.length > 0 && (
              <Autocomplete
                sx={{
                  px: 0,
                  py: 0,
                }}
                options={allPlots}
                getOptionKey={(option) => allPlots.indexOf(option)}
                getOptionLabel={(option) => option.title || "Unnamed plot"}
                onChange={(event, newValue) => {
                  setSelectedPlot(newValue);
                }}
                disabled={allPlots.length === 1}
                value={selectedPlot}
                renderInput={(params) => (
                  <TextField {...params} size="small" label="Plot" />
                )}
              />
            )
          }
          action={
            <>
              {jobId && graphId && (
                <Button
                  onClick={() => {
                    window.open(
                      `/graph-viewer/${jobId}/${graphId}`,
                      "_blank",
                      "width=1000,height=700"
                    );
                  }}
                  title="Open in standalone viewer"
                >
                  <OpenInNew />
                </Button>
              )}
              <Button
                onClick={() => {
                  setShowJson(true);
                }}
                title="Show JSON data"
              >
                <Preview />
              </Button>
            </>
          }
        />
        <CardContent sx={{ height }}>
          {chartArgs && (
            <Scatter
              key={
                allPlots && selectedPlot
                  ? `plot-${allPlots.indexOf(selectedPlot)}`
                  : "scatter-plot"
              }
              options={chartArgs.options as any}
              data={chartArgs.data}
            />
          )}
        </CardContent>
      </Card>
      <Dialog
        fullWidth
        maxWidth="xl"
        open={showJson}
        onClose={() => {
          setShowJson(false);
        }}
      >
        <DialogContent>
          <Typography variant="h6">JSON</Typography>
          <Editor
            height="200px"
            defaultLanguage="json"
            defaultValue={JSON.stringify(selectedPlot, null, 2)}
            theme={mode === "dark" ? "vs-dark" : "light"}
          />
        </DialogContent>
      </Dialog>
    </>
  );
};
