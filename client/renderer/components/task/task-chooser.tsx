// filepath: client/renderer/components/task/task-chooser.tsx
import { useCallback, useContext, useEffect, useMemo, useState } from "react";
import {
  Avatar,
  Card,
  CardContent,
  CardHeader,
  Chip,
  Collapse,
  Grid2,
  Paper,
  Skeleton,
  Toolbar,
  Box,
} from "@mui/material";
import { ElaborateSearch } from "../General/SearchObjects";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import { useApi } from "../../api";
import { MyExpandMore } from "../expand-more";

interface TaskTree {
  lookup: any;
  tree: any[];
  iconLookup: any;
}
interface CCP4i2TaskTreeProps {
  onTaskSelect?: (taskName: string) => void;
}
interface CCP4i2TaskTreeFolderProps {
  category: [name: string, title: string, taskNames: string[]];
  searchText: string | null;
  taskTree: TaskTree;
  iconLookup: any;
  onTaskSelect?: (taskName: string) => void;
}
interface CCP4i2TaskCardProps {
  taskName: string;
  task: any;
  onTaskSelect?: (taskName: string) => void;
}

export const CCP4i2TaskTree: React.FC<CCP4i2TaskTreeProps> = ({
  onTaskSelect,
}) => {
  const api = useApi();
  const [searchText, setSearchText] = useState<string | null>(null);
  const { data: taskTreeResult } = api.get<any>(`task_tree/`);

  // Handle new API response format: {success: true, data: {task_tree: {...}}}
  const taskTree = taskTreeResult?.success
    ? taskTreeResult?.data?.task_tree
    : taskTreeResult?.task_tree; // Fallback for legacy format
  const iconLookup = taskTree?.iconLookup;

  useEffect(() => {
    console.log("taskTreeResult", taskTreeResult);
  }, [taskTreeResult]);

  return (
    <Paper
      sx={{
        maxHeight: "calc(100vh - 10rem)",
        overflowY: "auto",
      }}
    >
      <Toolbar>
        {taskTree?.lookup &&
          `Tasks (Currently numbering ${Object.keys(taskTree?.lookup).length})`}
        <ElaborateSearch
          searchValue={searchText}
          setSearchValue={setSearchText}
        />
      </Toolbar>
      <Paper>
        {taskTree?.tree?.map(
          (
            category: [name: string, title: string, taskNames: string[]],
            iCategory: number
          ) => (
            <CCP4i2TaskTreeFolder
              key={`${iCategory}`}
              {...{ category, taskTree, searchText, onTaskSelect, iconLookup }}
            />
          )
        )}
      </Paper>
    </Paper>
  );
};

const CCP4i2TaskTreeFolder: React.FC<CCP4i2TaskTreeFolderProps> = ({
  category,
  searchText,
  taskTree,
  iconLookup,
  onTaskSelect,
}) => {
  const [tasksExpanded, setTasksExpanded] = useState<boolean>(false);

  const filterFunc = useCallback(
    (taskName: string) => {
      if (searchText == null || searchText.trim().length == 0) return true;
      if (!taskTree?.lookup) return false;
      if (!Object.keys(taskTree?.lookup).includes(taskName)) return false;
      const task = taskTree.lookup[taskName];
      if (!task) return false;
      return (
        task?.TASKTITLE?.toUpperCase().includes(
          searchText.toUpperCase()
        ) ||
        taskName
          .toUpperCase()
          .includes(searchText.toUpperCase()) ||
        task?.DESCRIPTION?.toUpperCase().includes(
          searchText.toUpperCase()
        )
      );
    },
    [taskTree, searchText]
  );

  const filteredTasks = useMemo(() => {
    if (!category || !taskTree) return [];
    if (!searchText || searchText.trim().length == 0) return category[2];
    return category[2].filter(filterFunc);
  }, [searchText, taskTree, category]);

  const searchActive = useMemo(() => {
    return searchText != null && searchText.trim().length > 0;
  }, [searchText]);

  return filteredTasks && filteredTasks.length > 0 ? (
    <Card
      key={category[0]}
      onClick={(ev) => {
        ev.stopPropagation();
        setTasksExpanded(!tasksExpanded);
      }}
      sx={{
        mb: 1,
        ":hover": { boxShadow: 4 },
      }}
    >
      <CardHeader
        slotProps={{
          title: { variant: "body2", my: 0, py: 0 },
          subheader: { variant: "caption", my: 0, py: 0 },
        }}
        sx={{
          py: 1,
          "& .MuiCardHeader-action": { alignSelf: "center", margin: 0 },
        }}
        avatar={
          <Avatar
            src={
              iconLookup
                ? `/${iconLookup[category[0]]}`
                : `/qticons/ccp4i2.png`
            }
            alt="/qticons/ccp4i2.png"
          />
        }
        title={category[1]}
        subheader={category[0]}
        action={
          <MyExpandMore
            expand={tasksExpanded}
            onClick={(ev) => {
              ev.stopPropagation();
              setTasksExpanded(!tasksExpanded);
            }}
            aria-expanded={tasksExpanded}
            aria-label="Show subjobs"
          >
            <ExpandMoreIcon />
          </MyExpandMore>
        }
      />
      {(tasksExpanded || searchActive) && (
        <CardContent sx={{ py: 1 }}>
          <Collapse
            in={tasksExpanded || searchActive}
            timeout="auto"
            unmountOnExit
          >
            <Box
              sx={{
                display: "grid",
                // 4 cards per row on wide screens, falling back to 3, 2, 1
                gridTemplateColumns:
                  "repeat(auto-fill, minmax(min(280px, 100%), 1fr))",
                gap: 2,
                p: 1,
              }}
            >
              {filteredTasks.map(
                (taskName: string) =>
                  Object.keys(taskTree.lookup).includes(taskName) && (
                    <CCP4i2TaskCard
                      key={taskName}
                      taskName={taskName}
                      task={taskTree.lookup[taskName]}
                      onTaskSelect={onTaskSelect}
                    />
                  )
              )}
            </Box>
          </Collapse>
        </CardContent>
      )}
    </Card>
  ) : null;  // Don't render empty folders when filtering
};

const CCP4i2TaskCard: React.FC<CCP4i2TaskCardProps> = ({
  taskName,
  task,
  onTaskSelect,
}) => {
  const handleTaskSelect = useCallback(() => {
    if (task && onTaskSelect) {
      onTaskSelect(taskName);
    }
  }, [task, onTaskSelect]);

  //useEffect(() => { console.log(task) }, [])
  return (
    <Card
      sx={(theme) => ({
        minHeight: "14rem",
        maxHeight: "24rem",
        overflowY: "auto",
        ":hover": { boxShadow: 24 },
        minWidth: "280px",
        // Theme-aware scrollbar styling
        scrollbarColor: `${theme.palette.action.disabled} transparent`,
        scrollbarWidth: "thin",
        "&::-webkit-scrollbar": {
          width: "8px",
        },
        "&::-webkit-scrollbar-track": {
          background: "transparent",
        },
        "&::-webkit-scrollbar-thumb": {
          backgroundColor: theme.palette.action.disabled,
          borderRadius: "4px",
          "&:hover": {
            backgroundColor: theme.palette.action.active,
          },
        },
      })}
      onClick={handleTaskSelect}
    >
      <CardHeader
        slotProps={{ title: { variant: "button", my: 0, py: 0 } }}
        title={
          <>
            <Avatar
              src={`/svgicons/${taskName}.svg`}
              alt={`/qticons/${taskName}.png`}
            />
            {task.shortTitle || task.TASKTITLE || taskName}
          </>
        }
        subheader={task.TASKTITLE || taskName}
      />
      <CardContent>
        {task.isAutogenerated && (
          <Chip
            label="Auto-generated Interface"
            size="small"
            color="info"
            variant="outlined"
            sx={{ mb: 1, fontSize: "0.7rem" }}
          />
        )}
        <p>{`${taskName}`}</p>
        <p>{`${task.DESCRIPTION || ""}`}</p>
      </CardContent>
    </Card>
  );
};
