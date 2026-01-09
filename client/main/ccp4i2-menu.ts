import { Menu, MenuItem, BrowserWindow } from "electron";
import { createWindow } from "./ccp4i2-create-window";
import { store } from "./ccp4i2-master";

// Function to toggle theme and notify all renderer windows
function toggleTheme() {
  const currentTheme = store.get("theme") as string;
  const newTheme = currentTheme === "light" ? "dark" : "light";
  store.set("theme", newTheme);
  console.log("Theme toggled to", newTheme);

  // Notify all windows of the theme change
  BrowserWindow.getAllWindows().forEach((win) => {
    win.webContents.send("message-from-main", {
      message: "theme-changed",
      theme: newTheme,
    });
  });

  // Rebuild menu to update the checkmark
  return newTheme;
}

// Function to add "New Window" item to the default menu
export function addNewWindowMenuItem(NEXT_PORT: number, DJANGO_PORT: number) {
  // Modify the default menu by adding a "New Window" option
  const menu = Menu.getApplicationMenu();
  if (!menu) {
    console.error("Menu not found");
    return;
  }
  // Find the File menu (usually at index 0)
  const fileMenu = menu.items[0];
  // If fileMenu is found, insert a "New Window" item right after the "New Tab" or similar item
  if (fileMenu) {
    fileMenu.submenu?.append(
      new MenuItem({
        label: "New Window",
        accelerator: "CmdOrCtrl+N", // Optional: add a keyboard shortcut (Cmd+N / Ctrl+N)
        click: () => {
          // Use /ccp4i2 base path for multi-app integration
          createWindow(`http://localhost:${NEXT_PORT}/ccp4i2`, store); // Create new window when this option is clicked
        },
      })
    );
    fileMenu.submenu?.append(
      new MenuItem({
        label: "Django restful API",
        accelerator: "CmdOrCtrl+D", // Optional: add a keyboard shortcut (Cmd+D / Ctrl+D)
        click: () => {
          createWindow(`http://localhost:${DJANGO_PORT}`, store); // Create new window when this option is clicked
        },
      })
    );
  }

  // Find the View menu and add theme toggle
  const viewMenu = menu.items.find((item) => item.label === "View");
  if (viewMenu && viewMenu.submenu) {
    viewMenu.submenu.append(new MenuItem({ type: "separator" }));
    viewMenu.submenu.append(
      new MenuItem({
        label: "Toggle Dark Mode",
        accelerator: "CmdOrCtrl+Shift+D",
        click: () => {
          toggleTheme();
        },
      })
    );
  }

  Menu.setApplicationMenu(menu);
}
