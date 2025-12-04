import { Components, Theme } from "@mui/material/styles";

// Module declarations for custom variants
declare module "@mui/material/Paper" {
  interface PaperPropsVariantOverrides {
    validation: true;
  }
}

declare module "@mui/material/CardHeader" {
  interface CardHeaderPropsVariantOverrides {
    lightGrey: true;
  }
}

declare module "@mui/material/Box" {
  interface BoxPropsVariantOverrides {
    validationContainer: true;
  }
}

// Add this module declaration for Toolbar
declare module "@mui/material/Toolbar" {
  interface ToolbarPropsVariantOverrides {
    lightGrey: true;
  }
}

export const createComponentVariants = (theme: Theme): Components => ({
  MuiTooltip: {
    defaultProps: {
      disableInteractive: true,
    },
  },
  MuiTableCell: {
    styleOverrides: {
      root: {
        variants: [
          {
            props: { variant: "head" },
            style: {
              backgroundColor: "rgb(92, 149, 165)", // Use the primary color for the header
              marginTop: "0.5rem",
              paddingBottom: 0,
            },
          },
          {
            props: { variant: "body" },
            style: {
              paddingTop: 0,
              paddingBottom: 0,
            },
          },
        ],
      },
    },
    defaultProps: {
      variant: "body", // Default variant if none is specified
    },
  },
  MuiTextField: {
    styleOverrides: {
      root: {
        marginTop: "0.5rem", // Add a gap of 0.5rem at the top
        marginBottom: "0.5rem",
      },
    },
  },
  MuiAutocomplete: {
    styleOverrides: {
      root: {
        marginTop: "0.5rem", // Add a gap of 0.5rem at the top
        marginBottom: "0.5rem",
      },
    },
  },
  MuiIcon: {
    styleOverrides: {
      root: { mx: 0, my: 0 },
    },
  },
  MuiCard: {
    variants: [
      {
        props: { variant: "validation" },
        style: {
          border: 2,
          borderColor: "divider", // Default, will be overridden by sx prop
          "&:hover": {
            borderColor: theme.palette.primary.light,
          },
        },
      },
    ],
  },
  // MuiBox is not a valid key in MUI Components, so it has been removed.
  MuiCardHeader: {
    styleOverrides: {
      root: {
        variants: [
          {
            props: { variant: "primary" },
            defaultProps: {
              disableTypography: false,
              titleTypographyProps: {
                variant: "body1",
              },
              subheaderTypographyProps: {
                variant: "subtitle2",
              },
            },
            style: {
              paddingTop: theme.spacing(0.5),
              paddingBottom: theme.spacing(0.5),
              fontWeight: "bold",
              backgroundColor: theme.palette.primary.main, // Use the primary color
              color: theme.palette.primary.contrastText,
            },
          },
          {
            props: { variant: "secondary" },
            defaultProps: {
              disableTypography: false,
              titleTypographyProps: {
                variant: "body1",
              },
              subheaderTypographyProps: {
                variant: "subtitle2",
              },
            },
            style: {
              paddingTop: theme.spacing(0.5),
              paddingBottom: theme.spacing(0.5),
              fontWeight: "bolder",
              color: theme.palette.primary.dark,
            },
          },
          {
            props: { variant: "lightGrey" },
            style: {
              backgroundColor: theme.palette.action.hover,
              color: theme.palette.text.primary,
              cursor: "pointer",
              "& .MuiCardHeader-title": {
                color: theme.palette.text.primary,
              },
              "&:hover": {
                backgroundColor: theme.palette.action.selected,
              },
            },
          },
        ],
      },
    },
    defaultProps: {
      disableTypography: false,
      titleTypographyProps: {
        variant: "body1",
      },
      subheaderTypographyProps: {
        variant: "subtitle2",
      },
    },
  },
  MuiToolbar: {
    variants: [
      {
        props: { variant: "lightGrey" },
        style: {
          backgroundColor: theme.palette.action.hover,
          color: theme.palette.text.primary,
          "& .MuiTypography-root": {
            color: theme.palette.text.primary,
          },
          "& .MuiIconButton-root": {
            color: theme.palette.text.primary,
          },
          "&:hover": {
            backgroundColor: theme.palette.action.selected,
          },
        },
      },
    ],
  },
});
