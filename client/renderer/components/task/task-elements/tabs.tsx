import { Box, Tab, Tabs } from "@mui/material";
import { PropsWithChildren, useState, Children } from "react";

export function a11yProps(index: Number) {
  return {
    id: `simple-tab-${index}`,
    "aria-controls": `simple-tabpanel-${index}`,
  };
}

export interface CCP4i2TabsProps {
  visibility?: boolean | (() => boolean);
}

export const CCP4i2Tabs: React.FC<PropsWithChildren<CCP4i2TabsProps>> = ({
  children,
  visibility,
}) => {
  const [value, setValue] = useState<number>(0);
  return (
    <>
      <Tabs
        value={value}
        onChange={(event, newValue) => {
          setValue(newValue);
        }}
        aria-label="basic tabs example"
      >
        {Children.map(children, (child: any, iChild) => (
          <Tab
            key={child.props.tab}
            label={child.props.label}
            {...a11yProps(iChild)}
          />
        ))}
      </Tabs>
      {Children.map(children, (child: any, iChild) => (
        <CCP4i2TabPanel key={iChild} value={value} index={iChild}>
          {child.props.children}
        </CCP4i2TabPanel>
      ))}
    </>
  );
};

export interface CCP4i2TabProps {
  label: string;
}

export const CCP4i2Tab: React.FC<PropsWithChildren<CCP4i2TabProps>> = ({
  children,
}) => {
  return <>{children}</>;
};

interface CCP4i2TabPanelProps {
  value: Number;
  index: Number;
  other?: any;
}
export const CCP4i2TabPanel: React.FC<
  PropsWithChildren<CCP4i2TabPanelProps>
> = ({ children, value, index, other }) => {
  return (
    <div
      role="tabpanel"
      hidden={value !== index}
      id={`simple-tabpanel-${index}`}
      aria-labelledby={`simple-tab-${index}`}
      sx={{ p: 0 }}
      {...other}
    >
      {value === index && <Box sx={{ p: 1 }}>{children}</Box>}
    </div>
  );
};
