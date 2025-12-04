import React, {
  createContext,
  useContext,
  useState,
  useCallback,
  ReactNode,
} from "react";
import Popcorn from "@mui/material/Snackbar"; // Replace with your actual popcorn/snackbar component import

type PopcornContextType = {
  setMessage: (msg: string) => void;
};

const PopcornContext = createContext<PopcornContextType | undefined>(undefined);

export const usePopcorn = () => {
  const ctx = useContext(PopcornContext);
  if (!ctx) throw new Error("usePopcorn must be used within a PopcornProvider");
  return ctx;
};

export const PopcornProvider: React.FC<{ children: ReactNode }> = ({
  children,
}) => {
  const [message, setMessage] = useState<string | null>(null);

  // Show message for 3 seconds
  const showMessage = useCallback((msg: string) => {
    setMessage(msg);
    setTimeout(() => setMessage(null), 3000);
  }, []);

  return (
    <PopcornContext.Provider value={{ setMessage: showMessage }}>
      {/* Popcorn/Snackbar at the top of the page */}
      <Popcorn
        open={!!message}
        message={message}
        anchorOrigin={{ vertical: "top", horizontal: "center" }}
        onClose={() => setMessage(null)}
        autoHideDuration={3000}
      />
      {children}
    </PopcornContext.Provider>
  );
};

// Usage in a child component:
// const { setMessage } = usePopcorn();
// setMessage("Hello, world!");
