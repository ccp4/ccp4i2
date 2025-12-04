import { useEffect, useRef, useState } from "react";
import { useCCP4i2Window } from "../app-context";
import { Typography } from "@mui/material";

interface RDKitViewProps {
  smiles: string;
}
export const RDKitView: React.FC<RDKitViewProps> = ({ smiles }) => {
  const { rdkitModule } = useCCP4i2Window();
  const [dataURI, setDataURI] = useState<string | null>(null);
  const lineSkipper: RegExp = /\>([\s\S]*)$/;
  const theBlob = useRef<Blob | undefined>(undefined);

  useEffect(() => {
    if (!rdkitModule) return;
    try {
      const mol = rdkitModule.get_mol(smiles);
      console.log("RDKit module is ready and SMILES processed:", mol);
      let longSvg = mol.get_svg();
      mol.delete();
      //let trimmedSvg = lineSkipper.exec(longSvg)[1]
      theBlob.current = new Blob([longSvg], { type: "image/svg+xml" });
      setDataURI(URL.createObjectURL(theBlob.current));
    } catch (error) {
      console.error("RDKit error:", error);
      setDataURI(null);
    }
  }, [smiles, rdkitModule]);

  return dataURI ? (
    <img src={dataURI} style={{ height: 250, width: 250 }} />
  ) : (
    <Typography sx={{ color: "text.error" }}>No valid SMILES </Typography>
  );
};
